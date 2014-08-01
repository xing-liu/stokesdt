/**
 * @file   simbox_bd.h
 * @brief  BDSimBox class definition
 */

#ifndef BD_SIMBOX_H
#define BD_SIMBOX_H


#include <vector>
#include <string>
#include <fstream>
#include "simbox_base.h"


namespace stokesdt {

class ForceBase;
class RndStream;
class MobBase;
class BrwnBase;

enum MobMethod {
    MOB_EWALD    = 0,
    MOB_SPME     = 1
};

enum BrwnbMethod {
    BRWN_CHOL    = 0,
    BRWN_LANCZOS = 1
};

/** @class  BDSimBox
 *  @brief  Container of Brownian dynamics simulations
 */
class BDSimBox : public SimBoxBase {
  public:
    BDSimBox(char *config_file,
             char *model_file,
             char *xyz_file);

    /// @copydoc SimBoxBase::~SimBoxBase()
    virtual ~BDSimBox();

    bool Init();

    void Advance(int nsteps);
        
  private:
    /// Disable copy constructor and assignement operator
    DISALLOW_COPY_AND_ASSIGN(BDSimBox);

  private:
    std::string config_file_;
    std::string model_file_;
    std::string xyz_file_;
    double delta_t_;    
    int npos_;
    double Lx_;
    double Ly_;
    double Lz_;
    double *pos_;
    double *rdi_;
    MobMethod mob_method_;
    int mob_interval_;     
    BrwnbMethod brwn_method_;
    
    int start_step_;
    double start_time_;
    int istep_;
    double cur_time_;

    // stats
    int traj_interval_;
    int traj_start_;
    std::ofstream traj_fs_;
    
    MobBase *mob_;
    int dim_mob_;
    std::vector<ForceBase *> force_;
    BrwnBase *brwn_;
    RndStream *rnd_stream_;
    double *brwn_vec_;
    double *disp_vec_;
};

} // namespace stokesdt


#endif // BD_SIMBOX_H