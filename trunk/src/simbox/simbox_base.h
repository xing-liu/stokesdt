/**
 * @file   simbox_base.h
 * @brief  SimBoxBase class definition.
 */


#ifndef SIMBOX_BASE_H_
#define SIMBOX_BASE_H_


#include "common.h"


namespace stokesdt {

/** @class  SimBoxBase
 *  @brief  Abstract base class for performing simulations
 */
class SimBoxBase {
  public:
    /// Class deconstructor
    virtual ~SimBoxBase();

  protected:
    /** 
     * @brief  Sole constructor. (For invocation by
     * derived class constructors, typically implicit.)
     */
    SimBoxBase();
    
  private:
    DISALLOW_COPY_AND_ASSIGN(SimBoxBase);
};

} // namespace stokesdt


#endif // STOKES_SIMBOX_H_
