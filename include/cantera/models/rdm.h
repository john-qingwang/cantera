#ifndef CT_RDM_H
#define CT_RDM_H

#include "cantera/numerics/DenseMatrix.h"
#include "ModelGeneric.h"

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// choose exact integral type
// #ifdef CGAL_USE_GMP
// #include <CGAL/Gmpz.h>
// typedef CGAL::Gmpz ET;
// #else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
// #endif

namespace Cantera
{

class RDM: public ModelGeneric
{

public:
    // for optimization
    typedef CGAL::Quadratic_program_from_iterators
    <double**,                                                // for A
     double*,                                                 // for b
     CGAL::Const_oneset_iterator<CGAL::Comparison_result>,    // for r
     bool*,                                                   // for fl
     double*,                                                 // for l
     bool*,                                                   // for fu
     double*,                                                 // for u
     double**,                                                // for D
     double*>                                                 // for c 
     Program;
    typedef CGAL::Quadratic_program_solution<ET> Solution;

    RDM() : ModelGeneric(nullptr), opt_interp(1), spline_bc(1) {};

    RDM(StFlow* const sf_) : ModelGeneric(sf_), opt_interp(1), spline_bc(1)  {};

    RDM(StFlow* const sf_, size_t delta);

    ~RDM() 
    {
        delete[] qp_l;
        delete[] qp_u;
        delete[] qp_fl;
        delete[] qp_fu;
        delete[] qp_A[0];
        delete[] qp_A;
        for (size_t j = 0; j < m_nz; j++) {
            delete[] qp_D[j];
            delete[] qp_C[j];
        }
        delete[] qp_D;
        delete[] qp_C;
    };

    // derived function from generic class
    void model_diff() {};

    void getWdot(const doublereal* x);

    // helper functions
    void update(StFlow* const sf_)
    {
        sf = sf_;
        update_grid();
    }

    void update_grid();

    void update_filter_type(std::string filter_type_) {
        filter_type = filter_type_;
        operator_init();
    }

    void update_filter_width(size_t delta_new)
    {
        m_delta = delta_new;
        update_grid();
    }

    void update_interp_method(std::string interp_method)
    {
        opt_interp_s = interp_method;
        interp_factory(sf->grid());
    }

    void update_spline_bc(std::string spline_bc_new)
    {
        spline_bc_s = spline_bc_new;
        interp_factory(sf->grid());
    }

    void update_regularization_amplification(doublereal alpha_amp_)
    {
        alpha_amp = alpha_amp_;
        operator_init();
    }

    void update_deconvolution_method(std::string deconv_opt)
    {
        if (deconv_opt.compare("constrained")==0) {
            use_constrained_qp = true;
            use_Landweber = false;
        } else if (deconv_opt.compare("Landweber")==0) {
            use_constrained_qp = false;
            use_Landweber = true;
        } else {
            use_constrained_qp = false;
            use_Landweber = false;
        }
    }

    size_t nPoints() {
        return m_nz;
    }

    vector_fp grid()
    {
        return m_z;
    }

    doublereal grid(size_t j)
    {
        return m_z[j];
    }

    doublereal wdot_fine(size_t k, size_t j)
    {
        return m_wdot_fine(k,j);
    }

    void model_summary()
    {
        std::cout << " RDM operator summary " << std::endl;
        std::cout << " ==================== " << std::endl;
        std::cout << " Filter type  : " << filter_type << std::endl;
        std::cout << " Filter width : " << std::setw(5) << m_delta << std::endl;
        std::cout << " Reg. Amp.    : " << std::setw(5) << alpha_amp << std::endl;
        std::cout << " Interp size  : " << std::setw(5) << m_nz << std::endl;
        std::cout << " Interp method: " << opt_interp_s << std::endl;
        std::cout << " Deconvolution: ";
        if (use_constrained_qp) {
            std::cout << "constrained" << std::endl;
        } else if (use_Landweber) {
            std::cout << "Landweber" << std::endl;
        } else {
            std::cout << "unconstrained" << std::endl;
        }
    }

    void interp_factory(const vector_fp& z);

    void operator_init();

    vector_fp interpolation(const vector_fp& x)
    {
        assert(m_interp.nColumns()==x.size());
        if (opt_interp_s.compare("ENO")==0) {
            update_ENO_interp(x);
        }
        return mat_vec_multiplication(m_interp,x);
    }

    vector_fp deconvolution(const vector_fp& x, size_t k)
    {
        assert(x.size()==m_nz);
        if (use_constrained_qp) {
            return constrained_deconvolution(x,k);
        } else if (use_Landweber) {
            return Landweber_deconvolution(x,k);
        } else {
            vector_fp sol = mat_vec_multiplication(m_Q,x);
            if (k >= c_offset_Y) {
                for (auto it = sol.begin(); it != sol.end(); ++it) {
                    *it = std::min(std::max(*it, 0.0),1.0);
                }
            } else if (k == c_offset_T) {
                for (auto it = sol.begin(); it != sol.end(); ++it) {
                    *it = std::min(std::max(*it, 300.0),2500.0);
                }
            }
            return sol;
        }
    }

    vector_fp filtering(const vector_fp& x)
    {
        assert(x.size()==m_nz);
        return mat_vec_multiplication(m_G,x);
    }

    void subgrid_reconstruction(const vector_fp& xbar, vector_fp& xdcv);

    vector_fp constrained_deconvolution(const vector_fp& x, size_t k); 

    vector_fp Landweber_deconvolution(const vector_fp& x, size_t k);

    vector_fp down_sampling(const vector_fp& x)
    {
        vector_fp xnew;
        size_t half_delta = (m_delta+1)/2;
        for (size_t i=0; i<x.size(); i+=half_delta) {
            xnew.push_back(x[i]);
        }
        return xnew;
    }

private:
    vector_fp mat_vec_multiplication(const DenseMatrix& A,const vector_fp& x);

    vector_fp Lagrange_polynomial(const vector_fp& z, const doublereal& m_z_i, const size_t idx_l, const size_t idx_r)
    {
        vector_fp coeff;
        for (size_t j = idx_l; j <= idx_r; j++) {
            doublereal num = 1.0;
            doublereal den = 1.0;
            for (size_t k = idx_l; k <= idx_r; k++) {
                if (k == j) continue;
                num *= (m_z_i -z[k]);
                den *= (  z[j]-z[k]);
            }
            coeff.push_back(num/den);
        }

        return coeff;
    }

    void update_ENO_interp(const vector_fp& x);

    // filter width
    size_t m_delta;
    // regularization factor
    doublereal alpha,alpha_amp;
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
    // ENO interpolation
    size_t order;
    std::vector<DenseMatrix> interp_pool;
    DenseMatrix stencil_idx;

    // interpolation matrix
    DenseMatrix m_interp;
    // Filter matrix
    std::string filter_type;
    DenseMatrix m_G;
    // deconvolution matrix
    DenseMatrix m_Q;
    // base subgrid contribution matrix (I - Q*G)
    DenseMatrix m_SG;
    // constrained optimization
    doublereal **qp_D, **qp_C, **qp_A;
    doublereal *qp_l, *qp_u;
    bool *qp_fl, *qp_fu;
    // Landweber deconvolution
    DenseMatrix lw_Greg;
    size_t lw_niter;
    // deconvolution option
    bool use_constrained_qp;
    bool use_Landweber;

    // species production term
    Array2D m_wdot_fine;

};

}
#endif
