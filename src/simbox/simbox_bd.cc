/**
 * @file   simbox_bd.cc
 * @brief  BDSimBox class implementation
 */

#include <limits.h>
#include <float.h>
#include <string.h>
#include <sys/time.h>

#include "log.h"
#include "profile.h"
#include "simbox_bd.h"
#include "molecule_io.h"
#include "mob_spme.h"
#include "mob_ewald.h"
#include "brwn_lanczos.h"
#include "brwn_chol.h"
#include "force_bonded.h"
#include "force_steric.h"
#include "rnd_stream.h"


namespace stokesdt {

using namespace std;


BDSimBox::BDSimBox(char *config_file,
                   char *model_file,
                   char *xyz_file)
    : config_file_(config_file),
      model_file_(model_file),
      xyz_file_(xyz_file)
{

}


BDSimBox::~BDSimBox()
{
    detail::PrintProfiler();
    delete mol_io_;
    delete mob_;
    delete brwn_;
    delete rnd_stream_;
    int num_forces = force_.size();
    for (int i = 0; i < num_forces; i++) {
        delete force_[i];
    }
    if (traj_fp_ != NULL) {
        fclose(traj_fp_);
    }
}


bool BDSimBox::ImportMolecule()
{
    mol_io_ = new MoleculeIO();
    // parse config files
    if (!mol_io_->ParseConfig(config_file_.c_str(), " \t", "#")) {
        return false;
    }
    if (!mol_io_->ParseModel(model_file_.c_str(), " \t", "#")) {
        return false;
    }
    if (!mol_io_->ParseXYZ(xyz_file_.c_str(), &npos_, &Lx_, &Ly_, &Lz_,
                     &start_step_, &start_time_)) {
        return false;
    }
    
    // import particles
    pos_ = (double *)detail::AlignMalloc(sizeof(double) * npos_ * 3);
    rdi_ = (double *)detail::AlignMalloc(sizeof(double) * npos_);
    if (NULL == pos_ || NULL == rdi_) {
        LOG_ERROR("Failed to allocate memory:\n");
        return false;
    }
    mol_io_->GetParticles(pos_, rdi_);

    return true;
}


bool BDSimBox::InitMob()
{
    // create mobility matrix
    const char* str_mob = mol_io_->GetStringKey("mob-method");
    if (strcmp(str_mob, "ewald") == 0) {
        mob_method_ = MOB_EWALD;
    } else if (strcmp(str_mob, "spme") == 0) {
        mob_method_ = MOB_SPME;
    } else {
        mob_method_ = MOB_EWALD;
    }    
    if (MOB_EWALD == mob_method_) {
        double ewald_tol = mol_io_->GetFloatKey("ewald-tol", 0.0, DBL_MAX, 1e-4);
        mob_ = new MobEwald(npos_, rdi_, Lx_, ewald_tol);
    } else if (MOB_SPME == mob_method_) {
        double spme_xi = mol_io_->GetFloatKey("spme-xi", 0.0, DBL_MAX, 0.5);
        double spme_rmax = mol_io_->GetFloatKey("spme-rmax", 0.0, DBL_MAX, 4.0);
        int spme_dim = mol_io_->GetIntKey("spme-nfft-grids", 1, 1024, 128);
        int spme_porder = mol_io_->GetIntKey("spme-porder", 4, 6, 4);
        mob_ = new MobSpme(npos_, rdi_, Lx_,
                           spme_xi, spme_rmax,
                           spme_dim, spme_porder);
    }

    // Init mobility matrix
    if (!mob_->Init()) {
        LOG_ERROR("Initialization of mobility matrix failed.\n");
        return false;
    }

    mob_interval_ = mol_io_->GetIntKey("mob-interval", 1, INT_MAX, 100);

    dim_mob_ = mob_->dim();
    ldm_ = detail::PadLen(dim_mob_, sizeof(double));

    
    return true;
}


bool BDSimBox::InitBrwn()
{
    // create Brownian
    const char* str_brwn = mol_io_->GetStringKey("brownian-method");
    if (strcmp(str_brwn, "chol") == 0) {
        brwn_method_ = BRWN_CHOL;
    } else if (strcmp(str_brwn, "lanczos") == 0) {
        brwn_method_ = BRWN_LANCZOS;
    } else {
        brwn_method_ = BRWN_CHOL;
    }
    if (brwn_method_ == BRWN_CHOL && mob_method_ == MOB_SPME) {
        LOG_ERROR("Cholesky requires using the explicit matrix approaches.\n");
        return false;
    }
    if (BRWN_CHOL == brwn_method_) {
        brwn_ = new BrwnChol(dim_mob_);
    } else if (BRWN_LANCZOS == brwn_method_) {
        int max_iters = mol_io_->GetIntKey("lanczos-maxiters", 1, INT_MAX, 50);
        int max_nrhs = mol_io_->GetIntKey("lanczos-maxnrhs", 1, INT_MAX, 20);
        double tol = mol_io_->GetFloatKey("lanczos-tol", 0.0, DBL_MAX, 1.0e-3);
        brwn_ = new BrwnLanczos(dim_mob_, max_iters, max_nrhs, tol);
    }

    if (!brwn_->Init()) {
        LOG_ERROR("Initialization of Brownian displacements failed.\n");
        return false;
    }

    // random stream
    rnd_stream_ = new RndStream(1234);
    if (!rnd_stream_->Init()) {
        LOG_ERROR("Initialization of random stream failed.\n");
        return false;
    }

    // create buffers
    brwn_vec_ = (double *)detail::AlignMalloc(sizeof(double) * mob_interval_ *
                                              dim_mob_ * ldm_);
    disp_vec_ = (double *)detail::AlignMalloc(sizeof(double) * mob_interval_ *
                                              dim_mob_ * ldm_);
    if (NULL == brwn_vec_ ||
        NULL == disp_vec_) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * dim_mob_ * ldm_ * mob_interval_);
        return false;
    }

    return true;
}


bool BDSimBox::InitForces()
{
    ForceBase *f;

    // bond force
    int num_bonds = mol_io_->LineCount("bond"); 
    std::vector<int> bond_id0(num_bonds);
    std::vector<int> bond_id1(num_bonds);
    std::vector<double> bond_cutoff(num_bonds);
    std::vector<double> bond_k0(num_bonds);
    std::string str_line;
    for (int i = 0; i < num_bonds; i++) {
        const char *line = mol_io_->GetLine("bond", i);
        sscanf(line, "%*s %d %d %le %le", &(bond_id0[i]), &(bond_id1[i]),
            &(bond_cutoff[i]), &(bond_k0[i]));
    }
    if (num_bonds > 0) {
      f = new BondedForce(npos_, num_bonds,
                        &bond_id0[0], &bond_id1[0], 
                        &bond_cutoff[0], &bond_k0[0]);
      force_.push_back(f);
    }
    
    // steric force
    if (1 == mol_io_->LineCount("steric")) {
        double steric_cutoff;
        double steric_k0;
        const char *line = mol_io_->GetLine("steric", 0);
        sscanf(line, "%le %le", &steric_cutoff, &steric_k0);
        f = new StericForce(npos_, rdi_, Lx_, steric_cutoff, steric_k0);
        force_.push_back(f);
    } else if (mol_io_->LineCount("steric") > 1) {
        LOG_ERROR("Duplicate \"steric\" lines found\n");
        return false;
    }

    // initialize forces
    int num_forces = force_.size();
    for (int i = 0; i < num_forces; i++) {
        if (!force_[i]->Init()) {
            LOG_ERROR("Initialization of forces failed.\n");
            return false;         
        }
    }

    return true;
}


bool BDSimBox::Init()
{
    if (!ImportMolecule()) {
        return false;
    }

    LOG(3, "\n        Initializes BDSimBox\n");
    LOG(3, "        --------------------\n");
    LOG(3, "Number-particles = %d\n", npos_);
    LOG(3, "Lx               = %g\n", Lx_);
    LOG(3, "Ly               = %g\n", Ly_);
    LOG(3, "Lz               = %g\n", Lz_);
    LOG(3, "Starting step    = %d\n", start_step_);
    LOG(3, "Starting time    = %g\n", start_time_);
    
    if (!InitMob()) {
        return false;
    }

    if (!InitBrwn()) {
        return false;
    }

    if (!InitForces()) {
        return false;
    }

    delta_t_ = mol_io_->GetFloatKey("delta-t", 0.0, DBL_MAX, 0.001);
    
    traj_interval_ = mol_io_->GetIntKey("traj-interval", INT_MIN, INT_MAX, -1);
    traj_start_ = mol_io_->GetIntKey("traj-startstep", 0, INT_MAX, 0);
    if (traj_interval_ > 0) {
        const char* traj_file = mol_io_->GetStringKey("traj-output");
        if (NULL != traj_file) {
            traj_fp_ = fopen(traj_file, "w+");
        } else {
            traj_fp_ = fopen("./trajectory.out", "w+");
        }
        if (traj_fp_ == NULL) {
            LOG_ERROR("Open file %s failed.\n", "trajectory.out");
            return false;
        }
    }
    
    // resume simulations
    istep_ = start_step_;
    cur_time_ = start_time_;

    detail::InitProfiler();
        
    return true;
}


void BDSimBox::Advance(int nsteps)
{
    struct timeval tv1;
    struct timeval tv2;
    
    LOG(3, "\n        Advance BDSimBox\n");
    LOG(3, "        ----------------\n");
    LOG(3, "Num-steps = %d\n\n", nsteps);
    
    int irhs = 0;
    int num_forces = force_.size();   
    __declspec(align(detail::kAlignLen)) double f[ldm_];
    for (int i = 0; i < nsteps; i++)
    {
        LOG(3, "Start step %10d\n", istep_);
        gettimeofday (&tv1, NULL);
        // each mob interval
        if (istep_ % mob_interval_ == 0)
        {
            irhs = 0;
            // generate random vectors
            rnd_stream_->Gaussian(0.0, sqrt(2.0 * delta_t_),
                                 dim_mob_, mob_interval_, ldm_, brwn_vec_);
            // build mobility matrix
            mob_->Update(pos_, rdi_);
            // compute Brownian displacements
            brwn_->Compute(mob_, mob_interval_,
                           ldm_, brwn_vec_, ldm_, disp_vec_);
        }
        double *cur_disp = &(disp_vec_[irhs * ldm_]);

        // compute external forces
        memset(f, 0, sizeof(double) * ldm_);
        for (int i = 0; i < num_forces; i++) {
            force_[i]->Accumulate(pos_, rdi_, f);
        }
        mob_->MulVector(1, delta_t_, ldm_, f, 1.0, ldm_, cur_disp);
        
        START_TIMER(detail::PARTICLE_TICKS);

        // update pos
        cblas_daxpy(npos_ * 3, 1.0, cur_disp, 1, pos_, 1);

        STOP_TIMER(detail::PARTICLE_TICKS);

        // advance time and steps        
        irhs++;
        istep_++;
        cur_time_ += delta_t_;
        
        gettimeofday (&tv2, NULL);
        timersub (&tv2, &tv1, &tv1);
        double timepass = tv1.tv_sec + tv1.tv_usec/1e6;
        LOG(3, "    Elapsed time: %.3le secs\n", timepass);

        if (traj_interval_ > 0) {
            mol_io_->WriteXYZ(npos_, Lx_, Ly_, Lz_,
                              istep_, cur_time_, NULL,
                              pos_, rdi_, traj_fp_);
        }
    }        
}

} // namespace stokesdt
