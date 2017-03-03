#ifndef MODEL_FACTORY_H
#define MODEL_FACTORY_H

#include "ModelGeneric.h"
#include "RDM.h"

namespace Cantera
{

class StFlow;

class ModelFactory
{

public:
    static ModelGeneric* make_model(StFlow* sf, std::string model_name)
    {
        if (model_name.compare("RDM")==0) {
            return new RDM(sf,5);
        }
        return new ModelGeneric(sf);
    }

};

} // namespace

#endif
