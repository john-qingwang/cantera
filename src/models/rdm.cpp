#include "cantera/models/rdm.h"

namespace Cantera
{

rdm::rdm(StFlow* const sf_, size_t delta, vector_fp& z) : 
     model_generic(sf_),m_delta(delta)
{
    // fine mesh for interpolation
    size_t refine_factor = (m_delta+1)/2;
    m_npts = z.size();
    m_nz = (m_npts-1)*refine_factor+1;
    m_z.resize(m_nz);
    for (size_t i=0; i<m_npts-1; i++) {
        doublereal dz = (z[i+1]-z[i])/(doublereal)refine_factor;
        for (size_t j=0; j<refine_factor; j++) {
            m_z[j+i*refine_factor] = z[i] + (doublereal)j*dz; 
        }
    }
    m_z[m_nz-1] = z[m_npts-1];

    // interpolation matrix
    opt_interp_s = "spline";
    spline_bc_s = "parabolic-run-out";
    if (opt_interp_s.compare("spline")==0) {
        opt_interp = 1;
        if (spline_bc_s.compare("free-run-out")==0) {
            spline_bc = 1;
        } else if(spline_bc_s.compare("parabolic-run-out")==0) {
            spline_bc = 2;
        } else if(spline_bc_s.compare("periodic")==0) {
            spline_bc = 3;
        } else {
            spline_bc = 1;
        }
    }
    m_interp.resize(m_nz,m_npts,0.0);
    interp_factory(z);
    
    // filter and deconvolution operators
    m_G.resize(m_nz,m_nz,0.0);
    m_Q.resize(m_nz,m_nz,0.0);
    m_SG.resize(m_nz,m_nz,0.0);
    operator_init();

    // solution arrays
    m_wdot.resize(sf->phase().nSpecies(),m_npts,0.0);
    wdot_fine.resize(sf->phase().nSpecies(),m_nz,0.0);

}

void rdm::update_grid(const vector_fp& z)
{
    size_t refine_factor = (m_delta+1)/2;

    m_npts = z.size();
    m_nz = (m_npts-1)*refine_factor+1;
    m_z.resize(m_nz);
    for (size_t i=0; i<m_npts-1; i++) {
        doublereal dz = (z[i+1]-z[i])/(doublereal)refine_factor;
        for (size_t j=0; j<refine_factor; j++) {
            m_z[j+i*refine_factor] = z[i] + (doublereal)j*dz; 
        }
    }
    m_z[m_nz-1] = z[m_npts-1];

    m_interp.resize(m_nz,m_npts,0.0);
    m_G.resize(m_nz,m_nz,0.0);
    m_Q.resize(m_nz,m_nz,0.0);

    interp_factory(z);
    operator_init();

    m_wdot.resize(sf->phase().nSpecies(),m_npts,0.0);
    wdot_fine.resize(sf->phase().nSpecies(),m_nz,0.0);

}

void rdm::interp_factory(const vector_fp& z)
{
    switch (opt_interp) {
    case 1: {
                // Spline interpolation
                vector_fp dz;
                dz.resize(m_npts-1);
                for (size_t i=0; i<m_npts-1; i++) {
                    dz[i] = z[i+1]-z[i];
                }
                
                DenseMatrix G(m_npts,m_npts,0.0);
                DenseMatrix rhs(m_npts,m_npts,0.0);
                G(0,0) = 1.0;
                G(m_npts-1,m_npts-1) = 1.0;
                for (size_t i = 1; i<m_npts-1; i++) {
                    G(i,i-1) = dz[i-1]/6.0;
                    G(i,i  ) = (dz[i-1]+dz[i])/3.0;
                    G(i,i+1) = dz[i]/6.0;
                    rhs(i,i-1) = 1.0/dz[i-1];
                    rhs(i,i  ) = -(1.0/dz[i-1]+1.0/dz[i]);
                    rhs(i,i+1) = 1.0/dz[i];
                }

                switch (spline_bc) {
                case 1: {
                            // free-run-out
                            // bc remain at 0
                            break;
                        }
                case 2: {
                            // parabolic-run-out
                            G(0,0) = G(1,0)+G(1,1);
                            G(0,2) = G(1,2);
                            doublereal* buf = new doublereal[m_npts];
                            rhs.getRow(1,buf);
                            rhs.setRow(0,buf);
                            G(m_npts-1,m_npts-1) = G(m_npts-2,m_npts-1)+G(m_npts-2,m_npts-2);
                            G(m_npts-1,m_npts-3) = G(m_npts-2,m_npts-3);
                            rhs.getRow(m_npts-2,buf);
                            rhs.setRow(m_npts-1,buf);
                            delete[] buf;
                            break;
                        }
                case 3: {
                            // periodic
                            break;
                        }
                default: break;
                }

                DenseMatrix C1(m_nz,m_npts,0.0);
                DenseMatrix C2(m_nz,m_npts,0.0);
                size_t i_src = 0;
                for (size_t i=0; i<m_nz; i++) {
                    doublereal z_ = m_z[i];
                    if (z_>z[i_src+1]) {
                        i_src++;
                    }
                    doublereal z_i  = z[i_src];
                    doublereal z_ip = z[i_src+1];
                    doublereal dz_i = dz[i_src];
                    C1(i,i_src  ) = (std::pow(z_ip-z_,3.0)/dz_i-dz_i*(z_ip-z_))/6.0;
                    C1(i,i_src+1) = (std::pow(z_-z_i ,3.0)/dz_i-dz_i*(z_-z_i ))/6.0;
                    C2(i,i_src  ) = (z_ip-z_)/dz_i;
                    C2(i,i_src+1) = (z_-z_i )/dz_i;
                }

                solve(G,rhs);
                C1.mult(rhs,m_interp);
                m_interp+=C2;

                break;
            }
    default: {
                 break;
             }
    }
    /* for debug
    std::ofstream debug_file;
    debug_file.open("mat_debug_interp.dat");
    for (size_t i=0; i<m_interp.nRows(); i++) {
        for (size_t j=0; j<m_interp.nColumns(); j++) {
            debug_file << m_interp(i,j) << ",";
        }
        debug_file << std::endl;
    }
    debug_file.close();
    */

}

void rdm::operator_init()
{
    size_t half_delta = (m_delta-1)/2;
    // filter operator
    for (size_t i=0; i<m_nz; i++) {
        int l_ = i-half_delta;
        int r_ = i+half_delta;
        for (size_t j=std::max(0,l_); j<=std::min((int)(m_nz)-1,r_); j++) {
            m_G(i,j) = 1.0/(doublereal)m_delta;
        }
        if (l_<0) {
            m_G(i,0) = (1.0-(doublereal)l_)*m_G(i,0);
        } else if (r_>m_nz-1) {
            m_G(i,m_nz-1) = (1.0+(doublereal)(r_-m_nz+1))*m_G(i,m_nz-1);
        }
    }
    // deconvolution operator
    DenseMatrix I(m_nz,m_nz,0.0);
    DenseMatrix GT(m_nz,m_nz,0.0);
    for (size_t i=0; i<m_nz; i++) {
        I(i,i) = 1.0;
        for (size_t j=0; j<m_nz; j++) {
            GT(j,i) = m_G(i,j);
        }
    }
    DenseMatrix G2(m_nz,m_nz,0.0);
    GT.mult(m_G,G2);
    doublereal alpha = 0.0;
    for (size_t i=0; i<m_nz; i++) {
        alpha += G2(i,i);
    }
    alpha /= (doublereal)m_nz;
    m_Q = m_G;
    m_Q *= -1.0;
    m_Q += I;
    DenseMatrix rhs(m_nz,m_nz,0.0);
    GT.mult(m_Q,rhs);
    DenseMatrix buf(m_nz,m_nz,0.0);
    buf = I;
    buf *= alpha;
    buf += G2;
    solve(buf,rhs);
    rhs += I;
    m_Q = rhs;

    // subgrid reconstruction operator
    m_Q.mult(m_G,m_SG);
    m_SG *= -1.0;
    m_SG += I;

    /* for debug
    std::ofstream debug_file;
    debug_file.open("mat_debug_Q.dat");
    for (size_t i=0; i<m_Q.nRows(); i++) {
        for (size_t j=0; j<m_Q.nColumns(); j++) {
            debug_file << m_Q(i,j) << ",";
        }
        debug_file << std::endl;
    }
    debug_file.close();
    */

}

vector_fp rdm::mat_vec_multiplication(const DenseMatrix& A,const vector_fp& x)
{
    assert(A.nColumns()==x.size());
    double* x_ = new double[x.size()];
    for (size_t i=0; i<x.size(); i++) {
        x_[i] = x[i];
    }
    double* sol = new double[A.nRows()];
    multiply(A,x_,sol);
    vector_fp sol_;
    for (size_t i=0; i<A.nRows(); i++) {
        sol_.push_back(sol[i]);
    }
    delete[] x_;
    delete[] sol;
    return sol_;
}

void rdm::subgrid_reconstruction(const vector_fp& xbar, vector_fp& xdcv)
{
    assert(xbar.size()==xdcv.size());

    vector_fp xvar1,xvar2;

    // variance 1: compute from algebric model
    doublereal Cz = 0.8;
    for (vector_fp::const_iterator it = xbar.begin(); it != xbar.end(); ++it) {
        xvar1.push_back(std::pow(*it,2.0));
    }
    vector_fp buf1 = filtering(xvar1);
    vector_fp buf2 = filtering(xbar);
    vector_fp::iterator it1 = buf1.begin();
    vector_fp::iterator it2 = buf2.begin();
    for (vector_fp::iterator it = xvar1.begin(); it != xvar1.end(); ++it) {
        *it = Cz*std::max(*it1 - std::pow(*it2,2.0), 0.0);
        ++it1;
        ++it2;
    }
    
    // variance 2: compute from deconvolution results
    vector_fp::const_iterator it0 = xbar.begin();
    it2 = xdcv.begin();
    for (vector_fp::iterator it = xvar2.begin(); it != xvar2.end(); ++it) {
        *it = std::pow(*it0-*it2,2.0);
        ++it0;
        ++it2;
    }

    // compute the subgrid contribution
    buf1 = mat_vec_multiplication(m_SG,xbar); 
    it1 = xvar1.begin();
    it2 = xvar2.begin();
    vector_fp::iterator it3 = buf1.begin();
    for (vector_fp::iterator it = xdcv.begin(); it != xdcv.end(); ++it) {
        *it += std::sqrt((*it1)/(*it2))*(*it3);
        ++it1;
        ++it2;
        ++it3;
    }
       
}

void rdm::model_src(const doublereal* x)
{
    size_t nsp = sf->phase().nSpecies();
    size_t nvar = c_offset_Y + nsp;

    doublereal* xnew = new doublereal[nvar*m_nz];

    // Deconvolving the variables
    for (size_t k=0; k<nvar; k++) {
        // Assemble the scalar k into one vector
        vector_fp sc;
        for (size_t i=0; i<m_npts; i++) {
            sc.push_back(*(x+k+i*nvar));
        }
        // Interpolation
        vector_fp sc_bar = interpolation(sc);
        // Deconvolution
        vector_fp sc_dcv = deconvolution(sc_bar);
        // Subgrid reconstruction
        subgrid_reconstruction(sc_bar,sc_dcv);
        // Assemble the scalar k back into the mixed vector
        for (size_t i=0; i<m_nz; i++) {
            xnew[k+i*nvar] = sc_dcv[i];
        }
    }
    // Compute the source term
    for (size_t j=0; j<m_nz; j++) {
        sf->setGas(xnew,j);
        sf->kinetics().getNetProductionRates(&wdot_fine(0,j));
    }
    // Filter the source term and down-scale it to the original size
    for (size_t k=0; k<nsp; k++) {
        vector_fp src;
        for (size_t j=0; j<m_nz; j++) {
            src.push_back(wdot_fine(k,j));
        }
        vector_fp src_buf = filtering(src);
        // /* for development purpose only
        for (size_t j=0; j<m_nz; j++) {
            wdot_fine(k,j) = src_buf[j];
        }
        // */
        vector_fp src_buf2 = down_sampling(src_buf);
        for (size_t j=0; j<m_npts; j++) {
            m_wdot(k,j) = src_buf2[j];
        }
    }

    delete[] xnew;

}

}

