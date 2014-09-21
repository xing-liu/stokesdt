/**
 * @file   brwn_random.h
 * @brief  BrwnRandom class definition
 */
 
#ifndef BRWN_RANDOM_H_
#define BRWN_RANDOM_H_


#include "brwn_base.h"
#include "rnd_stream.h"
#include "mob_base.h"

namespace stokesdt {

/** 
 *  @class  BrwnRandom
 *  @brief  Computes Brownian displacement vectors using randomized sampling.
 */
class BrwnRandom : public BrwnBase {
  public:
    /** 
     * @brief  Class constructor
     *
     * Constructs a new BrwnRandom instance that uses randomized sampling 
     * to compute Brownian displacement vectors for a given MobBase instance.
     * <p>
     * The <code>dim</code> argument specifies the dimension of
     * the mobility matrix represented by the MobBase instance.
     *
     * @param[in] dim       the dimesnion of the MobBase instance to be computed
     * @param[in] num_vecs  the number of random vectors to be used
     * @param[in] max_nrhs  the maximum number of right-hand sides
     */
    BrwnRandom(const int dim, const int num_vecs, const int max_nrhs);

    /// \copydoc  BrwnBase::~BrwnBase()
    virtual ~BrwnRandom();

    /// \copydoc  BrwnBase::Init()
    virtual bool Init();

    /// \copydoc  BrwnBase::Compute()
    virtual void Compute(MobBase *mob, const int num_rhs,
                         const int ldz, const double *z,
                         const int ldy, double *y);

  private:
    DISALLOW_COPY_AND_ASSIGN(BrwnRandom);
    
  private:
    /// the seed
    enum { SEED = 5678 };
    /// the dimension of the mobility matrix
    int dim_;
    /// the leading dimension of the mobility matrix
    int ldm_;
    /// the number of random vectors
    int num_vecs_;
    /// the leading dimension of b
    int ldb_;
    /// the maximum number of right-hand sides
    int max_nrhs_;
    /// the buffer for random vectors
    double *u_ = NULL;
    /// the buffer for s
    double *s_ = NULL;
    /// the buffer for V
    double *v_ = NULL;
    /// the buffer for the first term
    double *first_term_ = NULL;
    /// the buffer for the second term
    double *second_term_ = NULL;
    /// the random number generator
    RndStream *rnd_stream_ = NULL;
    /// the mobility matrix
    MobBase *mob_ = NULL;

    double *tau_;
    double *b_;
    double *temp_;
};

} // namespace stokesdt


#endif // BRWN_RANDOM_H_