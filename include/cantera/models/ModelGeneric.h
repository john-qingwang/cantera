#ifndef CT_MODEL_GENERIC_H
#define CT_MODEL_GENERIC_H

#include <iomanip>
#include "cantera/oneD/StFlow.h"

namespace Cantera
{

class ModelGeneric
{

public:
    ModelGeneric() : sf(nullptr) {};

    ModelGeneric(StFlow* const sf_) : sf(sf_)
    {
        if (sf_!=nullptr) {
            m_wdot.resize(sf_->phase().nSpecies(),sf_->nPoints(),0.0);
        }
    }

    virtual ~ModelGeneric() {};

    virtual void model_diff() {};

    virtual void getWdot(const doublereal* x) {};

    virtual void update(StFlow* sf_) {
        std::cout << "CAUSION: calling the update routine from generic class" << std::endl;
        sf = sf_;
        m_wdot.resize(sf_->phase().nSpecies(),sf_->nPoints());
    }

    virtual void model_summary() {};

    doublereal wdot_orig(size_t k, size_t j)
    {
        return sf->wdot(k,j);
    }

    doublereal wdot(size_t k, size_t j)
    {
        return m_wdot(k,j);
    }

protected:
    // a pointer to which the rdm operator is belonging to
    StFlow* sf;
    // species production term
    Array2D m_wdot;

};

}
#endif
