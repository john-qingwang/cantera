#include "cantera/numerics/DenseMatrix.h"
#include "cantera/oneD/StFlow.h"
#include "model_generic.h"

#ifndef CT_RDM_H
#define CT_RDM_H

namespace Cantera
{

class rdm: public model_generic
{
public:
    rdm(StFlow* const sf_) : model_generic(sf_) {};

    rdm(StFlow* const sf_, size_t delta, vector_fp& z);

    ~rdm() {};

    void update_grid(const vector_fp& z);

    size_t nPoints() {
        return m_nz;
    }

    vector_fp grid() {
        return m_z;
    }

    doublereal grid(size_t j) {
        return m_z[j];
    }

    doublereal getWdot(size_t k, size_t j) {
        return m_wdot(k,j);
    }

    void interp_factory(const vector_fp& z);

    void operator_init();

    vector_fp interpolation(const vector_fp& x) {
        assert(m_interp.nColumns()==x.size());
        return mat_vec_multiplication(m_interp,x);
    }

    vector_fp deconvolution(const vector_fp& x) {
        assert(x.size()==m_nz);
        return mat_vec_multiplication(m_Q,x);
    }

    vector_fp filtering(const vector_fp& x) {
        assert(x.size()==m_nz);
        return mat_vec_multiplication(m_G,x);
    }

    void subgrid_reconstruction(const vector_fp& xbar, vector_fp& xdcv);

    vector_fp down_sampling(const vector_fp& x) {
        vector_fp xnew;
        size_t half_delta = (m_delta+1)/2;
        for (size_t i=0; i<x.size(); i+=half_delta) {
            xnew.push_back(x[i]);
        }
        return xnew;
    }

    void model_diff() {};

    void model_src(const doublereal* x);

    // for verification
    void wdot_orig(const doublereal* x,size_t j,doublereal* wdot_) {
        sf->setGas(x,j);
        sf->kinetics().getNetProductionRates(wdot_);
    }

    doublereal getWdot_fine(size_t k, size_t j) {
        return wdot_fine(k,j);
    }

private:
    vector_fp mat_vec_multiplication(const DenseMatrix& A,const vector_fp& x);

    // filter width
    size_t m_delta;
    // number of grid points in the original grid
    size_t m_npts;
    // number of grid points in the interpolation grid
    size_t m_nz;
    // interpolation grid
    vector_fp m_z;
    // interpolation option
    int opt_interp;
    std::string opt_interp_s;
    // Spline interpolation
    int spline_bc;
    std::string spline_bc_s;

    // interpolation matrix
    DenseMatrix m_interp;
    // Filter matrix
    DenseMatrix m_G;
    // deconvolution matrix
    DenseMatrix m_Q;
    // base subgrid contribution matrix (I - Q*G)
    DenseMatrix m_SG;

    // species production term
    Array2D wdot_fine;

};

}
#endif
