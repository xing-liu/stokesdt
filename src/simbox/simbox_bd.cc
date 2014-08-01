/**
 * @file   simbox_bd.cc
 * @brief  BDSimBox class implementation
 */

#include <limits.h>
#include <float.h>
#include <string.h>

#include "log.h"
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
    delete mob_;
    delete brwn_;
    delete rnd_stream_;
    int num_forces = force_.size();
    for (int i = 0; i < num_forces; i++) {
        delete force_[i];
    }
}


bool BDSimBox::Init()
{  
    MoleculeIO *parser = new MoleculeIO();
    // parse config files
    if (!parser->ParseConfig(config_file_.c_str(), " \t", "#")) {
        return false;
    }
        printf("before\n");
    if (!parser->ParseModel(model_file_.c_str(), " \t", "#")) {
        return false;
    }
        printf("before\n");
    if (!parser->ParseXYZ(xyz_file_.c_str(), &npos_, &Lx_, &Ly_, &Lz_,
                     &start_step_, &start_time_)) {
        return false;
    }
    printf("before\n");
    // import particles
    pos_ = (double *)detail::AlignMalloc(sizeof(double) * npos_ * 3);
    rdi_ = (double *)detail::AlignMalloc(sizeof(double) * npos_);
    if (NULL == pos_ || NULL == rdi_) {
        LOG_ERROR("Failed to allocate memory:\n");
        return false;
    }
    parser->GetParticles(pos_, rdi_);

        printf("mobility\n");
    // create mobility matrix
    const char* str_mob = parser->GetStringKey("mob-method");
    if (str_mob == "ewald") {
        mob_method_ = MOB_EWALD;
    } else if (str_mob == "spme") {
        mob_method_ = MOB_SPME;
    } else {
        mob_method_ = MOB_EWALD;
    }    
    if (MOB_EWALD == mob_method_) {
        double ewald_tol = parser->GetFloatKey("ewald-tol", 0.0, DBL_MAX, 1e-4);
        mob_ = new MobEwald(npos_, rdi_, Lx_, ewald_tol);
    } else if (MOB_SPME == mob_method_) {
        double spme_xi = parser->GetFloatKey("spme-xi", 0.0, DBL_MAX, 0.5);
        double spme_rmax = parser->GetFloatKey("spme-rmax", 0.0, DBL_MAX, 4.0);
        int spme_dim = parser->GetIntKey("spme-nfft-grids", 1, 1024, 128);
        int spme_porder = parser->GetIntKey("spme-porder", 4, 6, 4);
        mob_ = new MobSpme(npos_, rdi_, Lx_,
                           spme_xi, spme_rmax,
                           spme_dim, spme_porder);
    }
    dim_mob_ = mob_->dim();
    printf("ok\n");
    // create Brownian
    const char* str_brwn = parser->GetStringKey("brownian-method");
    if (str_brwn == "chol" == 0) {
        brwn_method_ = BRWN_CHOL;
    } else if (str_brwn == "lanczos" == 0) {
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
        int max_iters = parser->GetIntKey("lanczos-maxiters", 1, INT_MAX, 50);
        int max_nrhs = parser->GetIntKey("lanczos-maxnrhs", 1, INT_MAX, 20);
        double tol = parser->GetFloatKey("lanczos-tol", 0.0, DBL_MAX, 1.0e-3);
        brwn_ = new BrwnLanczos(dim_mob_, max_iters, max_nrhs, tol);
    }
    printf("haha\n");
    ForceBase *f;
    // bond force
    int num_bonds = parser->LineCount("bond"); 
    std::vector<int> bond_id0(num_bonds);
    std::vector<int> bond_id1(num_bonds);
    std::vector<double> bond_cutoff(num_bonds);
    std::vector<double> bond_k0(num_bonds);
    std::string str_line;
    for (int i = 0; i < num_bonds; i++) {
        const char *line = parser->GetLine("bond", i);
        sscanf(line, "%*s %d %d %le %le", &(bond_id0[i]), &(bond_id1[i]),
            &(bond_cutoff[i]), &(bond_k0[i]));
    }
    if (num_bonds > 0) {
      f = new BondedForce(npos_, num_bonds,
                        &bond_id0[0], &bond_id1[0], 
                        &bond_cutoff[0], &bond_k0[0]);
      force_.push_back(f);
    }
    printf("cool\n");
    
    // steric force
    if (1 == parser->LineCount("steric")) {
        double steric_cutoff;
        double steric_k0;
        const char *line = parser->GetLine("steric", 0);
        sscanf(line, "%le %le", &steric_cutoff, &steric_k0);
        printf("push %d %lf %lf %lf\n", npos_, Lx_, steric_cutoff, steric_k0);
        f = new StericForce(npos_, rdi_, Lx_, steric_cutoff, steric_k0);
        force_.push_back(f);
    } else if (parser->LineCount("steric") > 1) {
        LOG_ERROR("Duplicate \"steric\" lines found\n");
        return false;
    }
    printf("force done\n");
    delta_t_ = parser->GetFloatKey("delta-t", 0.0, DBL_MAX, 0.001);
    // random stream
    rnd_stream_ = new RndStream(1234);    
    printf("force done\n");
    // Initilization
    if (!mob_->Init()) {
        LOG_ERROR("Initialization of mobility matrix failed.\n");
        return false;
    }
    if (!brwn_->Init()) {
        LOG_ERROR("Initialization of Brownian displacements failed.\n");
        return false;
    }
    if (!rnd_stream_->Init()) {
        LOG_ERROR("Initialization of random stream failed.\n");
        return false;
    }
    printf("force done\n");
    int num_forces = force_.size();
    for (int i = 0; i < num_forces; i++) {
        printf("init who %d\n", i);
        if (!force_[i]->Init()) {
            LOG_ERROR("Initialization of forces failed.\n");
            return false;         
        }
    }
   printf("force init\n");
    traj_interval_ = parser->GetIntKey("traj-interval", INT_MIN, INT_MAX, -1);
    traj_start_ = parser->GetIntKey("traj-startstep", 0, INT_MAX, 0);
    if (traj_interval_ > 0) {
        const char* traj_file = parser->GetStringKey("traj-output");
        if (0) {
            traj_fs_.open(traj_file);
        } else {
            traj_fs_.open("./trajectory.out");
        }
    }
     printf("done init\n");
    
    // resume simulations
    istep_ = start_step_;
    cur_time_ = start_time_;

    delete parser;
    
    return true;
}


void BDSimBox::Advance(int nsteps)
{
    int irhs = 0;
    int ldm = detail::PadLen(dim_mob_, sizeof(double));
    int num_forces = force_.size();
    __declspec(align(detail::kAlignLen)) double f[ldm];
    for (int i = 0; i < nsteps; i++)
    {   
        // each mob interval
        if (istep_ % mob_interval_ == 0)
        {
            irhs = 0;
            // generate random vectors
            rnd_stream_->Gaussian(0.0, sqrt(2.0/delta_t_),
                                 dim_mob_, mob_interval_, ldm, brwn_vec_);
            // build mobility matrix
            mob_->Update(pos_, rdi_);
            // compute Brownian displacements
            brwn_->Compute(mob_, mob_interval_, ldm, brwn_vec_, ldm, disp_vec_);            
        }
        double *cur_disp = &(disp_vec_[irhs * ldm]);

        // compute external forces
        memset(f, 0, sizeof(double) * ldm);
        for (int i = 0; i < num_forces; i++) {
            force_[i]->Accumulate(pos_, rdi_, f);
        }
        mob_->MulVector(1, 1.0, ldm, f, 1.0, ldm, cur_disp);
            
        // update pos
        cblas_daxpy(npos_ * 3, 1.0, cur_disp, 1, pos_, 1);

        // advance indices and time
        irhs++;
        istep_++;
        cur_time_ += delta_t_;
    }        
}

} // namespace stokesdt