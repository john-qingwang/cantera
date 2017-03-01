#ifndef CT_MODEL_GENERIC_H
#define CT_MODEL_GENERIC_H

#include "cantera/oneD/StFlow.h"

namespace Cantera
{

class model_generic
{

public:
    model_generic(StFlow* const sf_) : sf(sf_)
    {
        m_wdot.resize(sf_->phase().nSpecies(),sf_->nPoints(),0.0);
    }

    ~model_generic() {};

    virtual void model_diff() {};

    virtual void model_src(const doublereal* x) {};

protected:
    // a pointer to which the rdm operator is belonging to
    StFlow* const sf;
    // species production term
    Array2D m_wdot;

};

}
#endif
