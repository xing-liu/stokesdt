/**
 *  @file   rnd_stream.h
 *  @brief  RndStream class definition
 */

#ifndef RND_STREAM_H_
#define RND_STREAM_H_


#include "common.h"


namespace stokesdt {

/** @class  RndStream
 *  @brief  Random number generator
 */
class RndStream {
  public:
    /**
     * @brief  Class constructor
     *
     * Construct a new RndStream instance.
     *
     * @param[in] seed  the seed used for generating random numbers
     */
    RndStream(const unsigned int seed);

    /// Class deconstructor
    virtual ~RndStream();

    /**
     * @brief  Initializes the instance.
     * @return <code>true</code> if the instance is initialized successfully;
     *         <code>false</code> otherwise. 
     */
    bool Init();

    /**
     * @brief  Generates a block of random vectors with normal distribution 
     *
     * @param[in]  mean     the mean value of the normal distribution
     * @param[in]  sigma    the standard deviation the normal distribution
     * @param[in]  num_rhs  the number of vectors to be generate
     * @param[in]  len_v    the length of each vector
     * @param[in]  ldv      the leading dimension of each vector
     * @param[out] v        the pointer to the random vectors 
     */
    void Gaussian(const double mean,
                  const double sigma,
                  const int num_rhs,
                  const int len_v,
                  const int ldv,
                  double *v);

    /**
     * @brief  Generates a block of random vectors with uniform distribution 
     *
     * @param[in]  low      the lower bound of the uniform distribution
     * @param[in]  high     the upper bound of the uniform distribution
     * @param[in]  num_rhs  the number of vectors to be generate
     * @param[in]  len_v    the length of each vector
     * @param[in]  ldv      the leading dimension of each vector
     * @param[out] v        the pointer to the random vectors 
     */
    void Uniform(const double low,
                 const double high,
                 const int num_rhs,
                 const int len_v,
                 const int ldv,
                 double *v);

  private:
    DISALLOW_COPY_AND_ASSIGN(RndStream);
        
  private:
    /// the seed
    unsigned int seed_;
    /// the pointer to random stream engine
    void *stream_;    
};


} // namespace stokesdt


#endif // RND_STREAM_H_ 