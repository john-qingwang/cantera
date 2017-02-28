#ifndef CT_MODEL_GENERIC_H
#define CT_MODEL_GENERIC_H

#include "cantera/oneD/StFlow.h"

namespace Cantera
{

class model_generic
{

public:
    model_generic(StFlow* const sf_) : sf(sf_) {};

    ~model_generic() {};

    virtual void model_diff() {};

    virtual void model_src(const doublereal* x) {};

protected:
    // a pointer to which the rdm operator is belonging to
    StFlow* const sf;

};

}
#endif
