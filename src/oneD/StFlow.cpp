//! @file StFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/StFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"
#include "cantera/models/ModelFactory.h"

using namespace std;

namespace Cantera
{

StFlow::StFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    Domain1D(nsp+c_offset_Y, points),
    m_press(-1.0),
    m_nsp(nsp),
    m_thermo(0),
    m_kin(0),
    m_trans(0),
    m_epsilon_left(0.0),
    m_epsilon_right(0.0),
    m_do_soret(false),
    m_do_multicomponent(false),
    m_do_radiation(false),
    m_kExcessLeft(0),
    m_kExcessRight(0),
    model(nullptr)
{
    m_type = cFlowType;
    m_points = points;
    m_thermo = ph;

    if (ph == 0) {
        return; // used to create a dummy object
    }

    size_t nsp2 = m_thermo->nSpecies();
    if (nsp2 != m_nsp) {
        m_nsp = nsp2;
        Domain1D::resize(m_nsp+c_offset_Y, points);
    }

    // make a local copy of the species molecular weight vector
    m_wt = m_thermo->molecularWeights();

    // the species mass fractions are the last components in the solution
    // vector, so the total number of components is the number of species
    // plus the offset of the first mass fraction.
    m_nv = c_offset_Y + m_nsp;

    // enable all species equations by default
    m_do_species.resize(m_nsp, true);

    // but turn off the energy equation at all points
    m_do_energy.resize(m_points,false);

    m_diff.resize(m_nsp*m_points);
    m_multidiff.resize(m_nsp*m_nsp*m_points);
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_ybar.resize(m_nsp);
    m_qdotRadiation.resize(m_points, 0.0);

    //-------------- default solution bounds --------------------
    setBounds(0, -1e20, 1e20); // no bounds on u
    setBounds(1, -1e20, 1e20); // V
    setBounds(2, 200.0, 1e9); // temperature bounds
    setBounds(3, -1e20, 1e20); // lambda should be negative

    // mass fraction bounds
    for (size_t k = 0; k < m_nsp; k++) {
        setBounds(c_offset_Y+k, -1.0e-7, 1.0e5);
    }

    //-------------------- grid refinement -------------------------
    m_refiner->setActive(0, false);
    m_refiner->setActive(1, false);
    m_refiner->setActive(2, false);
    m_refiner->setActive(3, false);

    vector_fp gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, gr.data());
    setID("stagnation flow");

    // Find indices for radiating species
    m_kRadiating.resize(2, npos);
    m_kRadiating[0] = m_thermo->speciesIndex("CO2");
    m_kRadiating[1] = m_thermo->speciesIndex("H2O");

    // For model development
    enable_model = false;
    if (nsp>1 && points>1) {
        std::cout << "updating model from StFlow::StFlow" << std::endl;
        set_model("RDM");
    }
}

void StFlow::resize(size_t ncomponents, size_t points)
{
    Domain1D::resize(ncomponents, points);
    m_rho.resize(m_points, 0.0);
    m_wtm.resize(m_points, 0.0);
    m_cp.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_do_energy.resize(m_points,false);
    m_qdotRadiation.resize(m_points, 0.0);
    m_fixedtemp.resize(m_points);

    m_dz.resize(m_points-1);
    m_z.resize(m_points);

    // For model development
    if (model!=nullptr) {
        std::cout << "updating model from StFlow::resize" << std::endl;
        model->update(this);
    }
}

void StFlow::setupGrid(size_t n, const doublereal* z)
{
    resize(m_nv, n);

    m_z[0] = z[0];
    for (size_t j = 1; j < m_points; j++) {
        if (z[j] <= z[j-1]) {
            throw CanteraError("StFlow::setupGrid",
                               "grid points must be monotonically increasing");
        }
        m_z[j] = z[j];
        m_dz[j-1] = m_z[j] - m_z[j-1];
    }

    // For model development: update the model whenever there's a change of the grid
    if (model!=nullptr) {
        std::cout << "updating model from StFlow::setupGrid" << std::endl;
        model->update(this);
    }
}

void StFlow::resetBadValues(double* xg) {
    double* x = xg + loc();
    for (size_t j = 0; j < m_points; j++) {
        double* Y = x + m_nv*j + c_offset_Y;
        m_thermo->setMassFractions(Y);
        m_thermo->getMassFractions(Y);
    }
}

void StFlow::setTransport(Transport& trans)
{
    m_trans = &trans;
    m_do_multicomponent = (m_trans->transportType() == "Multi");

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
}

void StFlow::enableSoret(bool withSoret)
{
    if (m_do_multicomponent) {
        m_do_soret = withSoret;
    } else {
        throw CanteraError("setTransport",
                           "Thermal diffusion (the Soret effect) "
                           "requires using a multicomponent transport model.");
    }
}

void StFlow::_getInitialSoln(double* x)
{
    for (size_t j = 0; j < m_points; j++) {
        T(x,j) = m_thermo->temperature();
        m_thermo->getMassFractions(&Y(x, 0, j));
    }
}

void StFlow::setGas(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(T(x,j));
    const doublereal* yy = x + m_nv*j + c_offset_Y;
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press);
}

void StFlow::setGasAtMidpoint(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
    const doublereal* yyj = x + m_nv*j + c_offset_Y;
    const doublereal* yyjp = x + m_nv*(j+1) + c_offset_Y;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo->setMassFractions_NoNorm(m_ybar.data());
    m_thermo->setPressure(m_press);
}

void StFlow::_finalize(const doublereal* x)
{
    size_t nz = m_zfix.size();
    bool e = m_do_energy[0];
    for (size_t j = 0; j < m_points; j++) {
        if (e || nz == 0) {
            m_fixedtemp[j] = T(x, j);
        } else {
            double zz = (z(j) - z(0))/(z(m_points - 1) - z(0));
            double tt = linearInterp(zz, m_zfix, m_tfix);
            m_fixedtemp[j] = tt;
        }
    }
    if (e) {
        solveEnergyEqn();
    }
}

void StFlow::eval(size_t jg, doublereal* xg,
                  doublereal* rg, integer* diagg, doublereal rdt)
{
    // if evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // if evaluating a Jacobian, compute the steady-state residual
    if (jg != npos) {
        rdt = 0.0;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    // properties are computed for grid points from j0 to j1
    size_t j0 = std::max<size_t>(jmin, 1) - 1;
    size_t j1 = std::min(jmax+1,m_points-1);

    // ------------ update properties ------------

    updateThermo(x, j0, j1);
    if (jg == npos || m_force_full_update) {
        // update transport properties only if a Jacobian is not being
        // evaluated, or if specifically requested
        updateTransport(x, j0, j1);
    }
    if (jg == npos) {
        double* Yleft = x + index(c_offset_Y, jmin);
        m_kExcessLeft = distance(Yleft, max_element(Yleft, Yleft + m_nsp));
        double* Yright = x + index(c_offset_Y, jmax);
        m_kExcessRight = distance(Yright, max_element(Yright, Yright + m_nsp));
    }

    // update the species diffusive mass fluxes whether or not a
    // Jacobian is being evaluated
    updateDiffFluxes(x, j0, j1);

    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    // calculation of qdotRadiation

    // The simple radiation model used was established by Y. Liu and B. Rogg [Y.
    // Liu and B. Rogg, Modelling of thermally radiating diffusion flames with
    // detailed chemistry and transport, EUROTHERM Seminars, 17:114-127, 1991].
    // This model uses the optically thin limit and the gray-gas approximation
    // to simply calculate a volume specified heat flux out of the Planck
    // absorption coefficients, the boundary emissivities and the temperature.
    // The model considers only CO2 and H2O as radiating species. Polynomial
    // lines calculate the species Planck coefficients for H2O and CO2. The data
    // for the lines is taken from the RADCAL program [Grosshandler, W. L.,
    // RADCAL: A Narrow-Band Model for Radiation Calculations in a Combustion
    // Environment, NIST technical note 1402, 1993]. The coefficients for the
    // polynomials are taken from [http://www.sandia.gov/TNF/radiation.html].

    if (m_do_radiation) {
        // variable definitions for the Planck absorption coefficient and the
        // radiation calculation:
        doublereal k_P_ref = 1.0*OneAtm;

        // polynomial coefficients:
        const doublereal c_H2O[6] = {-0.23093, -1.12390, 9.41530, -2.99880,
                                     0.51382, -1.86840e-5};
        const doublereal c_CO2[6] = {18.741, -121.310, 273.500, -194.050,
                                     56.310, -5.8169};

        // calculation of the two boundary values
        double boundary_Rad_left = m_epsilon_left * StefanBoltz * pow(T(x, 0), 4);
        double boundary_Rad_right = m_epsilon_right * StefanBoltz * pow(T(x, m_points - 1), 4);

        // loop over all grid points
        for (size_t j = jmin; j < jmax; j++) {
            // helping variable for the calculation
            double radiative_heat_loss = 0;

            // calculation of the mean Planck absorption coefficient
            double k_P = 0;
            // absorption coefficient for H2O
            if (m_kRadiating[1] != npos) {
                double k_P_H2O = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_H2O += c_H2O[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_H2O /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[1], j) * k_P_H2O;
            }
            // absorption coefficient for CO2
            if (m_kRadiating[0] != npos) {
                double k_P_CO2 = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_CO2 += c_CO2[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_CO2 /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[0], j) * k_P_CO2;
            }

            // calculation of the radiative heat loss term
            radiative_heat_loss = 2 * k_P *(2 * StefanBoltz * pow(T(x, j), 4)
            - boundary_Rad_left - boundary_Rad_right);

            // set the radiative heat loss vector
            m_qdotRadiation[j] = radiative_heat_loss;
        }
    }

    // Compute the modeled terms if necessary
    if (enable_model && jg==npos) {
        model->getWdot(x);
    }

    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left, since
            // rho_u at point 0 is dependent on rho_u at point 1, but not on
            // mdot from the inlet.
            rsd[index(c_offset_U,0)] =
                -(rho_u(x,1) - rho_u(x,0))/m_dz[0]
                -(density(1)*V(x,1) + density(0)*V(x,0));

            // the inlet (or other) object connected to this one will modify
            // these equations by subtracting its values for V, T, and mdot. As
            // a result, these residual equations will force the solution
            // variables to the values for the boundary object
            rsd[index(c_offset_V,0)] = V(x,0);
            rsd[index(c_offset_T,0)] = T(x,0);
            rsd[index(c_offset_L,0)] = -rho_u(x,0);

            // The default boundary condition for species is zero flux. However,
            // the boundary object may modify this.
            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                sum += Y(x,k,0);
                rsd[index(c_offset_Y + k, 0)] =
                    -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
            }
            rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;
        } else if (j == m_points - 1) {
            evalRightBoundary(x, rsd, diag, rdt);
        } else { // interior points
            evalContinuity(j, x, rsd, diag, rdt);

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho dV/dt + \rho u dV/dz + \rho V^2
            //       = d(\mu dV/dz)/dz - lambda
            //-------------------------------------------------
            rsd[index(c_offset_V,j)]
            = (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j)
               - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
              - rdt*(V(x,j) - V_prev(j));
            diag[index(c_offset_V, j)] = 1;

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //-------------------------------------------------
            getWdot(x,j);
            for (size_t k = 0; k < m_nsp; k++) {
                double convec = rho_u(x,j)*dYdz(x,k,j);
                double diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                                / (z(j+1) - z(j-1));
                double wdot_;
                if (enable_model && jg==npos) {
                    wdot_ = model->wdot(k,j);
                } else {
                    wdot_ = this->wdot(k,j);
                }
                rsd[index(c_offset_Y + k, j)]
                = (m_wt[k]*(wdot_)
                   - convec - diffus)/m_rho[j]
                  - rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Y + k, j)] = 1;
            }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //-----------------------------------------------
            if (m_do_energy[j]) {
                setGas(x,j);

                // heat release term
                const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
                const vector_fp& cp_R = m_thermo->cp_R_ref();
                double sum = 0.0;
                double sum2 = 0.0;
                for (size_t k = 0; k < m_nsp; k++) {
                    double flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    sum += wdot(k,j)*h_RT[k];
                    sum2 += flxk*cp_R[k]/m_wt[k];
                }
                sum *= GasConstant * T(x,j);
                double dtdzj = dTdz(x,j);
                sum2 *= GasConstant * dtdzj;

                rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x,j)*dtdzj
                                            - divHeatFlux(x,j) - sum - sum2;
                rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
                rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
                diag[index(c_offset_T, j)] = 1;
            } else {
                // residual equations if the energy equation is disabled
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
            diag[index(c_offset_L, j)] = 0;
        }
    }
}

void StFlow::updateTransport(doublereal* x, size_t j0, size_t j1)
{
     if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            doublereal wtm = m_thermo->meanMolecularWeight();
            doublereal rho = m_thermo->density();
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,j)]);

            // Use m_diff as storage for the factor outside the summation
            for (size_t k = 0; k < m_nsp; k++) {
                m_diff[k+j*m_nsp] = m_wt[k] * rho / (wtm*wtm);
            }

            m_tcon[j] = m_trans->thermalConductivity();
            if (m_do_soret) {
                m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + j*m_nsp);
            }
        }
    } else { // mixture averaged transport
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMixDiffCoeffs(&m_diff[j*m_nsp]);
            m_tcon[j] = m_trans->thermalConductivity();
        }
    }
}

void StFlow::showSolution(const doublereal* x)
{
    writelog("    Pressure:  {:10.4g} Pa\n", m_press);

    Domain1D::showSolution(x);

    if (m_do_radiation) {
        writeline('-', 79, false, true);
        writelog("\n          z      radiative heat loss");
        writeline('-', 79, false, true);
        for (size_t j = 0; j < m_points; j++) {
            writelog("\n {:10.4g}        {:10.4g}", m_z[j], m_qdotRadiation[j]);
        }
        writelog("\n");
    }
}

void StFlow::updateDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                doublereal sum = 0.0;
                for (size_t m = 0; m < m_nsp; m++) {
                    sum += m_wt[m] * m_multidiff[mindex(k,m,j)] * (X(x,m,j+1)-X(x,m,j));
                }
                m_flux(k,j) = sum * m_diff[k+j*m_nsp] / dz;
            }
        }
    } else {
        for (size_t j = j0; j < j1; j++) {
            double sum = 0.0;
            double wtm = m_wtm[j];
            double rho = density(j);
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
                m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
                sum -= m_flux(k,j);
            }
            // correction flux to insure that \sum_k Y_k V_k = 0.
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) += sum*Y(x,k,j);
            }
        }
    }

    if (m_do_soret) {
        for (size_t m = j0; m < j1; m++) {
            double gradlogT = 2.0 * (T(x,m+1) - T(x,m)) /
                              ((T(x,m+1) + T(x,m)) * (z(m+1) - z(m)));
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,m) -= m_dthermal(k,m)*gradlogT;
            }
        }
    }
}

string StFlow::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "u";
    case 1:
        return "V";
    case 2:
        return "T";
    case 3:
        return "lambda";
    default:
        if (n >= c_offset_Y && n < (c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else {
            return "<unknown>";
        }
    }
}

size_t StFlow::componentIndex(const std::string& name) const
{
    if (name=="u") {
        return 0;
    } else if (name=="V") {
        return 1;
    } else if (name=="T") {
        return 2;
    } else if (name=="lambda") {
        return 3;
    } else {
        for (size_t n=c_offset_Y; n<m_nsp+c_offset_Y; n++) {
            if (componentName(n)==name) {
                return n;
            }
        }
    }
    return npos;
}

void StFlow::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    vector<string> ignored;
    size_t nsp = m_thermo->nSpecies();
    vector_int did_species(nsp, 0);

    vector<XML_Node*> str = dom.getChildren("string");
    for (size_t istr = 0; istr < str.size(); istr++) {
        const XML_Node& nd = *str[istr];
        writelog(nd["title"]+": "+nd.value()+"\n");
    }

    double pp = getFloat(dom, "pressure", "pressure");
    setPressure(pp);
    vector<XML_Node*> d = dom.child("grid_data").getChildren("floatArray");
    vector_fp x;
    size_t np = 0;
    bool readgrid = false, wrote_header = false;
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        if (nm == "z") {
            getFloatArray(fa,x,false);
            np = x.size();
            if (loglevel >= 2) {
                writelog("Grid contains {} points.\n", np);
            }
            readgrid = true;
            setupGrid(np, x.data());
        }
    }
    if (!readgrid) {
        throw CanteraError("StFlow::restore",
                           "domain contains no grid points.");
    }

    debuglog("Importing datasets:\n", loglevel >= 2);
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        getFloatArray(fa,x,false);
        if (nm == "u") {
            debuglog("axial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "axial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(0,j)] = x[j];
            }
        } else if (nm == "z") {
            ; // already read grid
        } else if (nm == "V") {
            debuglog("radial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "radial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(1,j)] = x[j];
            }
        } else if (nm == "T") {
            debuglog("temperature   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "temperature array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(2,j)] = x[j];
            }

            // For fixed-temperature simulations, use the imported temperature
            // profile by default.  If this is not desired, call
            // setFixedTempProfile *after* restoring the solution.
            vector_fp zz(np);
            for (size_t jj = 0; jj < np; jj++) {
                zz[jj] = (grid(jj) - zmin())/(zmax() - zmin());
            }
            setFixedTempProfile(zz, x);
        } else if (nm == "L") {
            debuglog("lambda   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "lambda arary size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(3,j)] = x[j];
            }
        } else if (m_thermo->speciesIndex(nm) != npos) {
            debuglog(nm+"   ", loglevel >= 2);
            if (x.size() == np) {
                size_t k = m_thermo->speciesIndex(nm);
                did_species[k] = 1;
                for (size_t j = 0; j < np; j++) {
                    soln[index(k+c_offset_Y,j)] = x[j];
                }
            }
        } else {
            ignored.push_back(nm);
        }
    }

    if (loglevel >=2 && !ignored.empty()) {
        writelog("\n\n");
        writelog("Ignoring datasets:\n");
        size_t nn = ignored.size();
        for (size_t n = 0; n < nn; n++) {
            writelog(ignored[n]+"   ");
        }
    }

    if (loglevel >= 1) {
        for (size_t ks = 0; ks < nsp; ks++) {
            if (did_species[ks] == 0) {
                if (!wrote_header) {
                    writelog("Missing data for species:\n");
                    wrote_header = true;
                }
                writelog(m_thermo->speciesName(ks)+" ");
            }
        }
    }

    if (dom.hasChild("energy_enabled")) {
        getFloatArray(dom, x, false, "", "energy_enabled");
        if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_energy[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            throw CanteraError("StFlow::restore", "energy_enabled is length {}"
                               "but should be length {}", x.size(), nPoints());
        }
    }

    if (dom.hasChild("species_enabled")) {
        getFloatArray(dom, x, false, "", "species_enabled");
        if (x.size() == m_nsp) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_species[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            // This may occur when restoring from a mechanism with a different
            // number of species.
            if (loglevel > 0) {
                writelog("\nWarning: StFlow::restore: species_enabled is "
                    "length {} but should be length {}. Enabling all species "
                    "equations by default.", x.size(), m_nsp);
            }
            m_do_species.assign(m_nsp, true);
        }
    }

    if (dom.hasChild("refine_criteria")) {
        XML_Node& ref = dom.child("refine_criteria");
        refiner().setCriteria(getFloat(ref, "ratio"), getFloat(ref, "slope"),
                              getFloat(ref, "curve"), getFloat(ref, "prune"));
        refiner().setGridMin(getFloat(ref, "grid_min"));
    }
}

XML_Node& StFlow::save(XML_Node& o, const doublereal* const sol)
{
    Array2D soln(m_nv, m_points, sol + loc());
    XML_Node& flow = Domain1D::save(o, sol);
    flow.addAttribute("type",flowType());

    if (m_desc != "") {
        addString(flow,"description",m_desc);
    }
    XML_Node& gv = flow.addChild("grid_data");
    addFloat(flow, "pressure", m_press, "Pa", "pressure");

    addFloatArray(gv,"z",m_z.size(), m_z.data(),
                  "m","length");
    vector_fp x(soln.nColumns());

    soln.getRow(0, x.data());
    addFloatArray(gv,"u",x.size(),x.data(),"m/s","velocity");

    soln.getRow(1, x.data());
    addFloatArray(gv,"V",x.size(),x.data(),"1/s","rate");

    soln.getRow(2, x.data());
    addFloatArray(gv,"T",x.size(),x.data(),"K","temperature");

    soln.getRow(3, x.data());
    addFloatArray(gv,"L",x.size(),x.data(),"N/m^4");

    for (size_t k = 0; k < m_nsp; k++) {
        soln.getRow(c_offset_Y+k, x.data());
        addFloatArray(gv,m_thermo->speciesName(k),
                      x.size(),x.data(),"","massFraction");
    }
    if (m_do_radiation) {
        addFloatArray(gv, "radiative_heat_loss", m_z.size(),
            m_qdotRadiation.data(), "W/m^3", "specificPower");
    }
    vector_fp values(nPoints());
    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = m_do_energy[i];
    }
    addNamedFloatArray(flow, "energy_enabled", nPoints(), &values[0]);

    values.resize(m_nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        values[i] = m_do_species[i];
    }
    addNamedFloatArray(flow, "species_enabled", m_nsp, &values[0]);

    XML_Node& ref = flow.addChild("refine_criteria");
    addFloat(ref, "ratio", refiner().maxRatio());
    addFloat(ref, "slope", refiner().maxDelta());
    addFloat(ref, "curve", refiner().maxSlope());
    addFloat(ref, "prune", refiner().prune());
    addFloat(ref, "grid_min", refiner().gridMin());
    return flow;
}

void StFlow::solveEnergyEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = true;
        }
    } else {
        if (!m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = true;
    }
    m_refiner->setActive(0, true);
    m_refiner->setActive(1, true);
    m_refiner->setActive(2, true);
    if (changed) {
        needJacUpdate();
    }
}

void StFlow::setBoundaryEmissivities(doublereal e_left, doublereal e_right)
{
    if (e_left < 0 || e_left > 1) {
        throw CanteraError("setBoundaryEmissivities",
            "The left boundary emissivity must be between 0.0 and 1.0!");
    } else if (e_right < 0 || e_right > 1) {
        throw CanteraError("setBoundaryEmissivities",
            "The right boundary emissivity must be between 0.0 and 1.0!");
    } else {
        m_epsilon_left = e_left;
        m_epsilon_right = e_right;
    }
}

void StFlow::fixTemperature(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = false;
        }
    } else {
        if (m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = false;
    }
    m_refiner->setActive(0, false);
    m_refiner->setActive(1, false);
    m_refiner->setActive(2, false);
    if (changed) {
        needJacUpdate();
    }
}

void AxiStagnFlow::evalRightBoundary(doublereal* x, doublereal* rsd,
                                     integer* diag, doublereal rdt)
{
    size_t j = m_points - 1;
    // the boundary object connected to the right of this one may modify or
    // replace these equations. The default boundary conditions are zero u, V,
    // and T, and zero diffusive flux for all species.
    rsd[index(0,j)] = rho_u(x,j);
    rsd[index(1,j)] = V(x,j);
    rsd[index(2,j)] = T(x,j);
    rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
    diag[index(c_offset_L, j)] = 0;
    doublereal sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[index(c_offset_Y + rightExcessSpecies(), j)] = 1.0 - sum;
    diag[index(c_offset_Y + rightExcessSpecies(), j)] = 0;
}

void AxiStagnFlow::evalContinuity(size_t j, doublereal* x, doublereal* rsd,
                                  integer* diag, doublereal rdt)
{
    //----------------------------------------------
    //    Continuity equation
    //
    //    Note that this propagates the mass flow rate information to the left
    //    (j+1 -> j) from the value specified at the right boundary. The
    //    lambda information propagates in the opposite direction.
    //
    //    d(\rho u)/dz + 2\rho V = 0
    //------------------------------------------------
    rsd[index(c_offset_U,j)] =
        -(rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
        -(density(j+1)*V(x,j+1) + density(j)*V(x,j));

    //algebraic constraint
    diag[index(c_offset_U, j)] = 0;
}

FreeFlame::FreeFlame(IdealGasPhase* ph, size_t nsp, size_t points) :
    StFlow(ph, nsp, points),
    m_zfixed(Undef),
    m_tfixed(Undef)
{
    m_dovisc = false;
    setID("flame");
}

void FreeFlame::evalRightBoundary(doublereal* x, doublereal* rsd,
                                  integer* diag, doublereal rdt)
{
    size_t j = m_points - 1;

    // the boundary object connected to the right of this one may modify or
    // replace these equations. The default boundary conditions are zero u, V,
    // and T, and zero diffusive flux for all species.

    // zero gradient
    rsd[index(0,j)] = rho_u(x,j) - rho_u(x,j-1);
    rsd[index(1,j)] = V(x,j);
    rsd[index(2,j)] = T(x,j) - T(x,j-1);
    doublereal sum = 0.0;
    rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
    diag[index(c_offset_L, j)] = 0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[index(c_offset_Y + rightExcessSpecies(), j)] = 1.0 - sum;
    diag[index(c_offset_Y + rightExcessSpecies(), j)] = 0;
}

void FreeFlame::evalContinuity(size_t j, doublereal* x, doublereal* rsd,
                               integer* diag, doublereal rdt)
{
    //----------------------------------------------
    //    Continuity equation
    //
    //    d(\rho u)/dz + 2\rho V = 0
    //----------------------------------------------
    if (grid(j) > m_zfixed) {
        rsd[index(c_offset_U,j)] =
            - (rho_u(x,j) - rho_u(x,j-1))/m_dz[j-1]
            - (density(j-1)*V(x,j-1) + density(j)*V(x,j));
    } else if (grid(j) == m_zfixed) {
        if (m_do_energy[j]) {
            rsd[index(c_offset_U,j)] = (T(x,j) - m_tfixed);
        } else {
            rsd[index(c_offset_U,j)] = (rho_u(x,j)
                                        - m_rho[0]*0.3);
        }
    } else if (grid(j) < m_zfixed) {
        rsd[index(c_offset_U,j)] =
            - (rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
            - (density(j+1)*V(x,j+1) + density(j)*V(x,j));
    }
    //algebraic constraint
    diag[index(c_offset_U, j)] = 0;
}

void FreeFlame::_finalize(const doublereal* x)
{
    StFlow::_finalize(x);
    // If the domain contains the temperature fixed point, make sure that it
    // is correctly set. This may be necessary when the grid has been modified
    // externally.
    if (m_tfixed != Undef) {
        for (size_t j = 0; j < m_points; j++) {
            if (z(j) == m_zfixed) {
                return; // fixed point is already set correctly
            }
        }

        for (size_t j = 0; j < m_points - 1; j++) {
            // Find where the temperature profile crosses the current
            // fixed temperature.
            if ((T(x, j) - m_tfixed) * (T(x, j+1) - m_tfixed) <= 0.0) {
                m_tfixed = T(x, j+1);
                m_zfixed = z(j+1);
                return;
            }
        }
    }
}

void FreeFlame::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    StFlow::restore(dom, soln, loglevel);
    getOptionalFloat(dom, "t_fixed", m_tfixed);
    getOptionalFloat(dom, "z_fixed", m_zfixed);
}

XML_Node& FreeFlame::save(XML_Node& o, const doublereal* const sol)
{
    XML_Node& flow = StFlow::save(o, sol);
    if (m_zfixed != Undef) {
        addFloat(flow, "z_fixed", m_zfixed, "m");
        addFloat(flow, "t_fixed", m_tfixed, "K");
    }
    return flow;
}

void StFlow::set_model(std::string model_name)
{
    std::cout << "updating model from StFlow::set_model" << std::endl;
    model = ModelFactory::make_model(this,model_name);
}

SprayFlame::SprayFlame(IdealGasPhase* ph, size_t nsp, size_t points) :
    AxiStagnFlow(ph, nsp, points)
{
    m_nv = c_offset_Y+m_nsp+c_offset_nl+1;
    Domain1D::resize(m_nv,points);

    setBounds(c_offset_Y+m_nsp+c_offset_Ul, -1e20, 1e20); // no bounds on Ul
    setBounds(c_offset_Y+m_nsp+c_offset_vl, -1e20, 1e20); // no bounds on vl
    setBounds(c_offset_Y+m_nsp+c_offset_Tl, 200.0, 5000.0); // bounds on Tl
    setBounds(c_offset_Y+m_nsp+c_offset_ml, -1e-7, 1e20); // bounds on ml
    setBounds(c_offset_Y+m_nsp+c_offset_nl, -1e-7, 1e20); // positivity for nl

    setID("spray flame");
    // CHANGE HERE 
    // (default is water)
    updateFuelSpecies("H2O");
    // vapor pressure coefficients
    m_prs_A = 8.07131;
    m_prs_B = 1730.63;
    m_prs_C = 233.426-273.15;
    // liquid density
    m_rhol_A = 0.14395;
    m_rhol_B = 0.0112;
    m_rhol_C = 649.727;
    m_rhol_D = 0.05107;
    // heat capacity [J/kmol/K]
    m_cpl = 76.0e+03;
    
}

void SprayFlame::eval(size_t jg, doublereal* xg,
                  doublereal* rg, integer* diagg, doublereal rdt)
{
    // Get the residual from the gaseous phase equations
    AxiStagnFlow::eval(jg,xg,rg,diagg,rdt);

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    // Gaseous phase
    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left, since
            // rho_u at point 0 is dependent on rho_u at point 1, but not on
            // mdot from the inlet.
            continue;
            // rsd[index(c_offset_U,0)] += (nl(x,0)*mdot(x,0) + nl(x,1)*mdot(x,1))/2.0;

        } else if (j == m_points - 1) {
            continue;

        } else { // interior points

            //------------------------------------------------
            //    Coninuity equation
            //------------------------------------------------
            // rsd[index(c_offset_U,j)] += (nl(x,j)*mdot(x,j) + nl(x,j+1)*mdot(x,j+1))/2.0;

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho dV/dt + \rho u dV/dz + \rho V^2
            //       = d(\mu dV/dz)/dz - lambda
            //         + nl mdot (Ul - Ug) - nl Fr
            //-------------------------------------------------
            rsd[index(c_offset_V,j)] -= ( nl(x,j) * Fr(x,j) / m_rho[j] );
            // rsd[index(c_offset_V,j)] += 
            //     (nl(x,j) * mdot(x,j) * (Ul(x,j)-V(x,j)) - nl(x,j) * Fr(x,j)) / m_rho[j];

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //     + (\delta_kf - Y_k) nl mdot
            //-------------------------------------------------
            // for (size_t k = 0; k < m_nsp; k++) {
            //     doublereal delta_kf;
            //     if (k == c_offset_fuel) {
            //         delta_kf = 1.0;
            //     } else {
            //         delta_kf = 0.0;
            //     }
            //     rsd[index(c_offset_Y + k, j)] += 
            //         (delta_kf - Y(x,k,j)) * nl(x,j) * mdot(x,j) / m_rho[j];
            // }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //      + nl mdot cp (Tl - Tg) - nl mdot q
            //-----------------------------------------------
            // rsd[index(c_offset_T, j)] += 
            //     (nl(x,j) * mdot(x,j) * cpgf(x,j) * (Tl(x,j) - T(x,j)) - 
            //      nl(x,j) * mdot(x,j) * q(x,j)) / (m_rho[j]*m_cp[j]);

        }
    }

    // Liquid phase
    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object
            
            // Number density. This propagates information right-to-left, since
            // ml_vl at point 0 is dependent on ml_vl at point 1, but not on
            // mdot from the inlet.
            // rsd[index(c_offset_Y+m_nsp+c_offset_nl,0)] = 
            //     -(nl_vl(x,1) - nl_vl(x,0))/m_dz[0]
            //     -(nl_Ul(x,1) + nl_Ul(x,0));
            rsd[index(c_offset_Y+m_nsp+c_offset_nl,0)] = nl(x,0); 
            diag[index(c_offset_Y+m_nsp+c_offset_nl, 0)] = 0;

            // the inlet (or other) object connected to this one will modify
            // these equations by subtracting its values for V, T, and mdot. As
            // a result, these residual equations will force the solution
            // variables to the values for the boundary object
            rsd[index(c_offset_Y+m_nsp+c_offset_vl,0)] = vl(x,0);
            rsd[index(c_offset_Y+m_nsp+c_offset_Ul,0)] = Ul(x,0);
            rsd[index(c_offset_Y+m_nsp+c_offset_Tl,0)] = Tl(x,0);
            rsd[index(c_offset_Y+m_nsp+c_offset_ml,0)] = ml(x,0);

            diag[index(c_offset_Y+m_nsp+c_offset_Ul, 0)] = 0;
            diag[index(c_offset_Y+m_nsp+c_offset_vl, 0)] = 0;
            diag[index(c_offset_Y+m_nsp+c_offset_Tl, 0)] = 0;
            diag[index(c_offset_Y+m_nsp+c_offset_ml, 0)] = 0;

        } else if (j == m_points - 1) {
            evalRightBoundaryLiquid(x, rsd, diag, rdt);
        } else { // interior points

            //------------------------------------------------
            //    Number density equation
            //------------------------------------------------
            evalNumberDensity(j, x, rsd, diag, rdt);

            //------------------------------------------------
            //    Mass equation
            //
            //    dm_l/dt + v_l dm_l/dz = -mdot
            //-------------------------------------------------
            rsd[index(c_offset_Y+m_nsp+c_offset_ml,j)] = -vl(x,j) * dmldz(x,j) - 
                rdt * (ml(x,j) - ml_prev(j));
            diag[index(c_offset_Y+m_nsp+c_offset_ml, j)] = 1;

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    m_l dU_l/dt + m_l v_l dU_l/dz + m_l U_l^2
            //       = f_r/r
            //-------------------------------------------------
            rsd[index(c_offset_Y+m_nsp+c_offset_Ul,j)] = -vl(x,j) * dUldz(x,j) -
                Ul(x,j)*Ul(x,j) - rdt*(Ul(x,j)-Ul_prev(j));
            if (ml(x,j)>sqrt(numeric_limits<double>::min())) {
                rsd[index(c_offset_Y+m_nsp+c_offset_Ul,j)] += Fr(x,j)/ml(x,j);
            }
            diag[index(c_offset_Y+m_nsp+c_offset_Ul, j)] = 1;

            //------------------------------------------------
            //    Axial momentum equation
            //
            //    m_l dv_l/dt + m_l v_l dv_l/dz = f_z
            //-------------------------------------------------
            rsd[index(c_offset_Y+m_nsp+c_offset_vl,j)] = -vl(x,j)*dvldz(x,j) - 
                    rdt*(vl(x,j)-vl_prev(j));
            if (ml(x,j)>sqrt(numeric_limits<double>::min())) {
                rsd[index(c_offset_Y+m_nsp+c_offset_vl,j)] += fz(x,j)/ml(x,j);
            }
            diag[index(c_offset_Y+m_nsp+c_offset_vl, j)] = 1;

            //-----------------------------------------------
            //    energy equation
            //
            //    m_l c_p_l dT_l/dt + m_l c_p_l v_l dT_l/dz
            //    = mdot_l (q - L) 
            //-----------------------------------------------
            rsd[index(c_offset_Y+m_nsp+c_offset_Tl,j)] = -vl(x,j)*dTldz(x,j) - 
                    rdt * (Tl(x,j) - Tl_prev(j));
            // if (ml(x,j)>sqrt(numeric_limits<double>::min())) {
            //     rsd[index(c_offset_Y+m_nsp+c_offset_Tl,j)] += 
            //         comp_mdot * mdot(x,j)*(q(x,j)-Lv()) / ml(x,j) / cpl(x,j);
            // }
            diag[index(c_offset_Y+m_nsp+c_offset_Tl, j)] = 1;

        }

    }

}

void SprayFlame::evalNumberDensity(size_t j, doublereal* x, doublereal* rsd,
                                  integer* diag, doublereal rdt)
{
     //----------------------------------------------
     //    Number density equation
     //
     //    Note that this propagates the liquid mass flow rate information to the right
     //    (j+1 -> j) from the value specified at the left boundary.
     //
     //    d(n_l v_l)/dz + 2n_l U_l = 0
     //------------------------------------------------
     rsd[index(c_offset_Y+m_nsp+c_offset_nl,j)] = -vl(x,j) * dnldz(x,j) -
         nl(x,j) * (vl(x,j) - vl(x,j-1))/m_dz[j-1] -
         2.0 * nl_Ul(x,j) - rdt * (nl(x,j) - nl_prev(j));
     // rsd[index(c_offset_Y+m_nsp+c_offset_nl,j)] =
     //     -(nl_vl(x,j) - nl_vl(x,j-1))/m_dz[j-1]
     //     -(nl_Ul(x,j) + nl_Ul(x,j-1))
     //     -rdt * (nl(x,j) - nl_prev(j));

     diag[index(c_offset_Y+m_nsp+c_offset_nl, j)] = 1;
}

void SprayFlame::evalRightBoundaryLiquid(doublereal* x, doublereal* rsd,
                                         integer* diag, doublereal rdt)
{
    size_t j = m_points - 1;

    rsd[index(c_offset_Y+m_nsp+c_offset_nl,j)] = nl(x,j) - nl(x,j-1);
    // rsd[index(c_offset_Y+m_nsp+c_offset_nl,j)] =
    //      -(nl_vl(x,j) - nl_vl(x,j-1))/m_dz[j-1]
    //      -(nl_Ul(x,j) + nl_Ul(x,j-1));

    diag[index(c_offset_Y+m_nsp+c_offset_nl, j)] = 0;

    // Neumann boundary condition for Ul, vl, Tl and ml
    rsd[index(c_offset_Y+m_nsp+c_offset_Ul,j)] = Ul(x,j) - Ul(x,j-1);
    rsd[index(c_offset_Y+m_nsp+c_offset_vl,j)] = vl(x,j) - vl(x,j-1);
    rsd[index(c_offset_Y+m_nsp+c_offset_Tl,j)] = Tl(x,j) - Tl(x,j-1);
    rsd[index(c_offset_Y+m_nsp+c_offset_ml,j)] = ml(x,j) - ml(x,j-1);
    diag[index(c_offset_Y+m_nsp+c_offset_Ul, j)] = 0;
    diag[index(c_offset_Y+m_nsp+c_offset_vl, j)] = 0;
    diag[index(c_offset_Y+m_nsp+c_offset_Tl, j)] = 0;
    diag[index(c_offset_Y+m_nsp+c_offset_ml, j)] = 0;
}

string SprayFlame::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "u";
    case 1:
        return "V";
    case 2:
        return "T";
    case 3:
        return "lambda";
    default:
        if (n >= c_offset_Y && n < (c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else if (n >= (c_offset_Y + m_nsp) && n <= (c_offset_Y + m_nsp + c_offset_nl)) {
            switch (n - c_offset_Y - m_nsp) {
                case 0:
                    return "Ul";
                case 1:
                    return "vl";
                case 2:
                    return "Tl";
                case 3:
                    return "ml";
                case 4:
                    return "nl";
            }
        } else {
            return "<unknown>";
        }
    }
}

size_t SprayFlame::componentIndex(const std::string& name) const
{
    if (name=="u") {
        return 0;
    } else if (name=="V") {
        return 1;
    } else if (name=="T") {
        return 2;
    } else if (name=="lambda") {
        return 3;
    } else if (name=="Ul") {
        return c_offset_Y + m_nsp + c_offset_Ul;
    } else if (name=="vl") {
        return c_offset_Y + m_nsp + c_offset_vl;
    } else if (name=="Tl") {
        return c_offset_Y + m_nsp + c_offset_Tl;
    } else if (name=="ml") {
        return c_offset_Y + m_nsp + c_offset_ml;
    } else if (name=="nl") {
        return c_offset_Y + m_nsp + c_offset_nl;
    } else {
        for (size_t n=c_offset_Y; n<m_nsp+c_offset_Y; n++) {
            if (componentName(n)==name) {
                return n;
            }
        }
    }
    return npos;
}

XML_Node& SprayFlame::save(XML_Node& o, const doublereal* const sol)
{
    Array2D soln(m_nv, m_points, sol + loc());
    XML_Node& flow = StFlow::save(o, sol);

    XML_Node& gv = flow.addChild("liq_grid_data");
    vector_fp x(soln.nColumns());

    soln.getRow(componentIndex("Ul"), x.data());
    addFloatArray(gv,"Ul",x.size(),x.data(),"m/s","liq. r velocity");

    soln.getRow(componentIndex("vl"), x.data());
    addFloatArray(gv,"vl",x.size(),x.data(),"m/s","liq. x velocity");

    soln.getRow(componentIndex("Tl"), x.data());
    addFloatArray(gv,"Tl",x.size(),x.data(),"K","liq. temperature");

    soln.getRow(componentIndex("ml"), x.data());
    addFloatArray(gv,"ml",x.size(),x.data(),"kg","droplet mass");

    soln.getRow(componentIndex("nl"), x.data());
    addFloatArray(gv,"nl",x.size(),x.data(),"/m^3","number density");

    return flow;
}

} // namespace
