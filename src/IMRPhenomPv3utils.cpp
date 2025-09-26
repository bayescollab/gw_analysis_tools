#include "IMRPhenomPv3utils.h"
#include <math.h>
#include <string>
#include <sstream>
#include <gsl/gsl_math.h>

/*! \file
 * Utilitites for IMRPhenomPv3
 * Much taken from LALSim
*/


/**
 * Initialize spin -2 l=2 spherical harmonics for spherical angles theta and phi
*/
template <class T>
void IMRPhenomPv3InitY2m(sph_harm<T> *Ylm, T theta, T phi)
{
    const std::complex<T> zero (0.,0.);
    *Ylm = {zero, zero, zero, zero, zero};

    Ylm->Y22  = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, 2, 2);
    Ylm->Y21  = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, 2, 1);
	Ylm->Y20  = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, 2, 0);
	Ylm->Y2m1 = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, 2, -1);
	Ylm->Y2m2 = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, 2, -2);
}

/** Compute the Wigner D^{2}_{+/- 2, m}(beta) elements 
 *  Stolen from XLALSimIMRPhenomPv3HMComputeWignerdElements
*/
template <class T>
void IMRPhenomPv3ComputeWignerD(
    T (*WignerD22)[2][5], /** [out] D^{2}_{2, m} elements */
    const T b /**< beta angle */
)
{
    T b2 = 2. * b;
    T cosb = cos(b);
    T cos2b = cos(b2);
    T sinb = sin(b);

    T cosb_over_two = cosb * 0.5;
    T cos2b_fac_1 = cos2b * ONE_OVER_EIGHT + THREE_OVER_EIGHT;
    T cosb_minus_1 = cosb - 1.0;
    T cosb_plus_1 = cosb + 1.0;

    T sinb_over_two = sinb * 0.5;

    //mprime == 2
    (*WignerD22)[0][0] = -cosb_over_two + cos2b_fac_1;           //m=-2
    (*WignerD22)[0][1] = cosb_minus_1 * sinb_over_two;           //m=-1
    (*WignerD22)[0][2] = SQRT_6 * (1. - cos2b) * ONE_OVER_EIGHT; //m=0
    (*WignerD22)[0][3] = -cosb_plus_1 * sinb_over_two;           //m=1
    (*WignerD22)[0][4] = cosb_over_two + cos2b_fac_1;            //m=2

    //mprime == -2
    (*WignerD22)[1][0] = (*WignerD22)[0][4];
    (*WignerD22)[1][1] = -(*WignerD22)[0][3];
    (*WignerD22)[1][2] = (*WignerD22)[0][2];
    (*WignerD22)[1][3] = -(*WignerD22)[0][1];
    (*WignerD22)[1][4] = (*WignerD22)[0][0];
}

/**
 * Given m1 with spins (chi1x, chi1y, chi1z) and m2 with spins (chi2x,chi2y,chi2z).
 * Enforce that m1 >= m2 and swap spins accordingly.
 * Enforce that the primary object (heavier) is indexed by 1.
 * To be used with precessing-spin waveform models.
 */
template <class T> void PhenomPrecessingSpinEnforcePrimary(
    T *m1,    /**< [out] mass of body 1 */
    T *m2,    /**< [out] mass of body 2 */
    T *chi1x, /**< [out] x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    T *chi1y, /**< [out] y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    T *chi1z, /**< [out] z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    T *chi2x, /**< [out] x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    T *chi2y, /**< [out] y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    T *chi2z  /**< [out] z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
)
{
    T m1_tmp, m2_tmp;
    T chi1x_tmp, chi1y_tmp, chi1z_tmp;
    T chi2x_tmp, chi2y_tmp, chi2z_tmp;
    if (*m1 > *m2)
    {
        chi1x_tmp = *chi1x;
        chi1y_tmp = *chi1y;
        chi1z_tmp = *chi1z;

        chi2x_tmp = *chi2x;
        chi2y_tmp = *chi2y;
        chi2z_tmp = *chi2z;

        m1_tmp = *m1;
        m2_tmp = *m2;
    }
    else
    { /* swap spins and masses */
        chi1x_tmp = *chi2x;
        chi1y_tmp = *chi2y;
        chi1z_tmp = *chi2z;

        chi2x_tmp = *chi1x;
        chi2y_tmp = *chi1y;
        chi2z_tmp = *chi1z;

        m1_tmp = *m2;
        m2_tmp = *m1;
    }
    *m1 = m1_tmp;
    *m2 = m2_tmp;
    *chi1x = chi1x_tmp;
    *chi1y = chi1y_tmp;
    *chi1z = chi1z_tmp;

    *chi2x = chi2x_tmp;
    *chi2y = chi2y_tmp;
    *chi2z = chi2z_tmp;

    if (*m1 < *m2)
        throw std::runtime_error("Unable to enforce m1 as the larger mass.");
}

template <class T> int approx_equal(T x, T X, T epsilon)
{
    return !(gsl_fcmp(x, X, epsilon));
}

template<> int approx_equal<adouble>(adouble x, adouble X, adouble epsilon)
{
    double x_val = x.value();
    double X_val = X.value();
    double eps = epsilon.value();

    return approx_equal(x_val, X_val, eps);
}

/**
 * If x and X are relatively equal within epsilon, set x = X.
 * If X = 0, compare absolute values.
*/
template <class T> void nudge(T *x, T X, T epsilon)
{
    if (X != 0)
    {
        if (approx_equal(*x, X, epsilon))
        {
            *x = X;
        }
    }
    else
    {
        if (fabs(*x-X) < epsilon)
        {
            *x = X;
        }
    }
}

/*
 * Adapted from LALSimInspiralFDPrecAngles_internals
*/
template <class T> int checkOmegaz5(const T Omegaz5)
{
    if (Omegaz5 >= 1000.0)
    {
        std::ostringstream str;
        str << "Omegaz5 = " << Omegaz5 << ", which is larger than expected.";

        throw std::domain_error(str.str());
    }

    return 0;
}

/*
 * Adapted from LALSimInspiralFDPrecAngles_internals
*/
void invalidExpansionOrder(const int ExpansionOrder)
{
    std::cerr << "ExpansionOrder = " << ExpansionOrder << " not recognised. \
    Defaulting to keeping all orders in expansion.\n";
}

/**
 * Create a vector in spherical coordinates using its magnitude r, polar angle th, and azimuthal angle ph 
*/
template <class T> vector3D<T> CreateSphVector(const T r, const T th, const T ph)
{
    const T fact = r*sin(th);
    vector3D<T> out {fact*cos(ph), fact*sin(ph), r*cos(th)};

    return out;
}

/**
 * Scale a vector vec by a scalar c.
*/
template <class T> vector3D<T> ScaleVector(T c, vector3D<T> vec)
{
    vector3D<T> out {c*vec.x, c*vec.y, c*vec.z};

    return out;
}

/**
 * Sum of two vectors
*/
template <class T> vector3D<T> VectorSum(vector3D<T> vec1, vector3D<T> vec2)
{
    vector3D<T> out;
    out.x = vec1.x+vec2.x;
    out.y = vec1.y+vec2.y;
    out.z = vec1.z+vec2.z;

    return out;
}

/**
 * Dot product of two vectors.
*/
template <class T> T DotProduct(const vector3D<T> vec1, const vector3D<T> vec2)
{
    return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
}

/**
 * Magnitude of a vector.
*/
template <class T> T VectorNorm(const vector3D<T> vec)
{
    return sqrt(DotProduct(vec, vec));
}

/**
 * Cross product of two vectors.
*/
template <class T> vector3D<T> CrossProduct(const vector3D<T> vec1, const vector3D<T> vec2)
{
    vector3D<T> out;
    out.x = vec1.y*vec2.z-vec1.z*vec2.y;
    out.y = vec1.z*vec2.x-vec1.x*vec2.z;
    out.z = vec1.x*vec2.y-vec1.y*vec2.x;

    return out;
}

/**
 * The roots of (dS^2/dt)^2 in a vector
 * See Eq. 22 in arxiv:1703.03967
 * 
 * out.x = A1 = S_{3}^2
 * out.y = A2 = S_{-}^2
 * out.z = A3 = S_{+}^2
 */
template <class T> vector3D<T> Roots(const T L_norm, const T J_norm, const sysprecquant<T> *system)
{
    vector3D<T> out;
    vector3D<T> coeffs = BCDcoeff(L_norm, J_norm, system);
    T A1, A2, A3;
    const T B_2 = coeffs.x*coeffs.x;
    const T B_3 = B_2*coeffs.x;
    const T B_C =  coeffs.x*coeffs.y;

    const T p = coeffs.y - ((*system).onethird)*B_2;
    const T qc = 2./27.*B_3 - ((*system).onethird)*B_C + coeffs.z;
    const T sqrtarg = sqrt(-p/3.);
    T acosarg = 1.5*qc/p/sqrtarg;
    if(acosarg < -1) acosarg = -1;
    if(acosarg > 1) acosarg = 1;
    const T theta = ((*system).onethird)*acos(acosarg);


    /* Case when the total spin is constant */
    /* when spin1 or spin2 is zero then the norm and the total spin is constant. But there is still precession. */
    if(theta!=theta || sqrtarg!=sqrtarg || (*system).dot1n == 1 || (*system).dot2n == 1 || (*system).dot1n == -1 || (*system).dot2n == -1|| (*system).S1_norm_2 == 0 || (*system).S2_norm_2 == 0) {
        out.x = 0;
        out.y = ((*system).S0_norm)*((*system).S0_norm);
        out.z = out.y + 1e-9;
        /* The 1e-9 is a nudge because sometimes the azimuthal
         * precession angle has a term that makes it -inf.
         * This small perturbation was to remedy these cases. */
    }
    else{
        out.z = 2.*sqrtarg*cos(theta) - ((*system).onethird)*coeffs.x;
        out.y = 2.*sqrtarg*cos(theta - GWAT_TWOPI*((*system).onethird)) - ((*system).onethird)*coeffs.x;
        out.x = 2.*sqrtarg*cos(theta - 2.*GWAT_TWOPI*((*system).onethird)) - ((*system).onethird)*coeffs.x;

        A3 = fmax(fmax(out.x,out.y),out.z);
        A1 = fmin(fmin(out.x,out.y),out.z);

        if ((A3 - out.z) > 0 && (A1 - out.z) < 0)
            A2 = out.z;
        else if ((A3 - out.x) > 0 && (A1 - out.x) < 0)
            A2 = out.x;
        else
            A2 = out.y;

        out.x = A1;
        out.y = A2;
        out.z = A3;

    }
    return out;
}

/**
 * The B,C,D coefficients of Eq. 21 in arxiv:1703.03967 in a vector.
 * The expressions are given in Appendix B.
 * 
 * B coefficient = eq. B2
 * C coefficient = eq. B3
 * D coefficient = eq. B4
 */
template <class T> vector3D<T> BCDcoeff(const T L_norm, const T J_norm, const sysprecquant<T> *system)
{
    vector3D<T> out;
    const T J_norm_2 = J_norm*J_norm;
    const T L_norm_2 = L_norm*L_norm;

    out.x = (L_norm_2+((*system).S1_norm_2))*((*system).q) + (L_norm_2+((*system).S2_norm_2))/((*system).q) + 2.*L_norm*((*system).Seff) - 2.*J_norm_2 - ((*system).S1_norm_2) - ((*system).S2_norm_2);

    out.y = (L_norm_2 - J_norm_2)*(L_norm_2 - J_norm_2) + 2.*((*system).Seff)*L_norm*(L_norm_2 - J_norm_2) - 2.*L_norm_2*(((*system).S1_norm_2) - ((*system).S2_norm_2)*((*system).q))*(1./((*system).q)-1) + 4.*((*system).nu)*((*system).Seff)*((*system).Seff)*L_norm_2 - 2.*((*system).Seff)*L_norm*(((*system).S1_norm_2)-((*system).S2_norm_2))*((*system).deltam_over_M) + 2.*J_norm_2*(((*system).S1_norm_2)*((*system).q) - ((*system).S2_norm_2))*(1./((*system).q)-1);

    out.z = -(L_norm_2 - J_norm_2)*(L_norm_2 - J_norm_2)*(((*system).S1_norm_2)*((*system).q) - ((*system).S2_norm_2))*(1./((*system).q)-1) - 2.*((*system).Seff)*((*system).deltam_over_M)*(((*system).S1_norm_2) - ((*system).S2_norm_2))*L_norm*(L_norm_2 - J_norm_2) + (((*system).S1_norm_2) - ((*system).S2_norm_2))*(((*system).S1_norm_2) - ((*system).S2_norm_2))*L_norm_2*((*system).deltam_over_M)*((*system).deltam_over_M)/((*system).nu);

    return out;
}

/**
 * Compute the spin-orbit couplings
 */
template <class T> T beta(const T a, const T b, const sysprecquant<T> *system)
{
    return ((((*system).dot1)*(a + b*((*system).q))) + (((*system).dot2)*(a + b/((*system).q))));
}

/**
 * Compute the spin-spin couplings
 */
template <class T> T sigma(const T a, const T b, const sysprecquant<T> *system)
{
    return (a*((*system).dot12) - b*((*system).dot1)*((*system).dot2))/((*system).nu);
}

/**
 * Compute the spin-spin couplings
 */
template <class T> T tau(const T a, const T b, const sysprecquant<T> *system)
{
    return (((*system).q)*((*system).S1_norm_2*a - b*((*system).dot1)*((*system).dot1)) + ((*system).S2_norm_2*a - b*((*system).dot2)*((*system).dot2))/((*system).q))/((*system).nu);
}

/**
 * Magnitude of L divided by GMsquare_over_c to
 * 3PN order
 * from 0605140 and Blanchet LRR and 1212.5520 Eq. 4.7
 */
template <class T> T L_norm_3PN_of_xi(const T xi, const T xi_2, const T L_norm, const sysprecquant<T> *system)
{
    return L_norm*(1. + xi_2*((*system).constants_L[0] + xi*(*system).constants_L[1] + xi_2*((*system).constants_L[2] + xi*(*system).constants_L[3] + xi_2*((*system).constants_L[4]))));
}

/**
 * Magnitude of J to Newtonian order
 * Equation 41 (1703.03967)
 */
template <class T> T J_norm_of_xi(const T L_norm, const sysprecquant<T> *system)
{
    return sqrt(L_norm*L_norm + 2.*L_norm*((*system).c_1_over_nu) + (*system).Ssqave);
}

template <class T> int elljacobi(T u, T m, T *sn, T *cn, T *dn)
{
    return gsl_sf_elljac_e(u, m, sn, cn, dn);
}

template <> int elljacobi<adouble>(adouble u, adouble m, adouble *sn, adouble *cn, adouble *dn)
{
    int gsl_ret;
    double u_val = u.value();
    double m_val = m.value();

    double sn_val, cn_val, dn_val;

    gsl_ret = gsl_sf_elljac_e(u_val, m_val, &sn_val, &cn_val, &dn_val);

    *sn = (adouble)sn_val;
    *cn = (adouble)cn_val;
    *dn = (adouble)dn_val;

    return gsl_ret;
}

/**
 * Magnitude of S divided by GMsquare_over_c
 * Equation 23 (1703.03967)
 */
template <class T> T S_norm_of_xi(const T xi, const T xi_2, const vector3D<T> roots, const sysprecquant<T> *system)
{
    T sn, cn, dn, m, u;

    if(fabs(roots.y-roots.z)<1.e-5) sn = 0.;
    else {
        m = (roots.y - roots.z)/(roots.x - roots.z);
        u = u_of_xi(xi,xi_2,system)+(*system).constant_of_S;
        elljacobi(u, m, &sn, &cn, &dn);
    }
    const T S_norm_square_bar = roots.z + (roots.y - roots.z)*sn*sn;
    return sqrt(S_norm_square_bar);
}

/**
 * Cosine of the angle between L and J
 * equation 8 1703.03967
 */
template <class T> T costhetaL(const T J_norm, const T L_norm, const T S_norm)
{
    T out = 0.5*(J_norm*J_norm + L_norm*L_norm - S_norm*S_norm)/L_norm/J_norm;

    if (out>1) out = 1;
    if (out<-1) out = -1;

    return out;
}

/**
 * Second term in phase of the magnitude of S
 * Eqn. 51 in arxiv:1703.03967
*/
template <class T> T u_of_xi(const T xi, const T xi_2, const sysprecquant<T> *system)
{

    /*
    dictionary from code variables to paper symbols arxiv:1703.03967
    (*system).deltam_over_M = (m1-m2)/(m1+m2)
    (*system).constants_u[0] = $-g_0$ in eq. 51. Where g_0 is given in eq. A1.
    xi corresponds to v so
    1./xi_2/xi corresponds to the $v^{-3}$ factor
    (*system).constants_u[1] = $\psi_1$ in eq. 51, given in eq. C1
    (*system).constants_u[2] = $\psi_2$ in eq. 51, given in eq. C2
    */

    return 0.75*((*system).deltam_over_M)*(*system).constants_u[0]/xi_2/xi*(1. + xi*((*system).constants_u[1] + xi*((*system).constants_u[2])));
}

/**
 * "c_0", "c_2" and "c_4" from 1703.03967
 * corresponding to equations
 * B6, B7 and B8 in 1703.03967 respectively
 * out.x = c_0
 * out.y = c_2
 * out.z = c_4
 */
template <class T> vector3D<T> c_coeffs(const T xi, const T xi_2, const T J_norm, const vector3D<T> roots, const sysprecquant<T> *system)
{
    const T xi_3 = xi_2*xi;
    const T xi_4 = xi_3*xi;
    const T xi_6 = xi_3*xi_3;

    const T J_norm_2 = J_norm*J_norm;

    vector3D<T> out;

    out.x = -0.75*((J_norm_2-roots.z)*(J_norm_2-roots.z)*xi_4/((*system).nu) - 4.*((*system).nu)*((*system).Seff)*(J_norm_2-roots.z)*xi_3-2.*(J_norm_2-roots.z+2*(((*system).S1_norm_2)-((*system).S2_norm_2))*((*system).deltam_over_M))*((*system).nu)*xi_2+(4.*((*system).Seff)*xi+1)*((*system).nu)*((*system).nu_2))*J_norm*xi_2*(((*system).Seff)*xi-1.);
    out.y = 1.5*(roots.y-roots.z)*J_norm*((J_norm_2-roots.z)/((*system).nu)*xi_2-2.*((*system).nu)*((*system).Seff)*xi-((*system).nu))*(((*system).Seff)*xi-1.)*xi_4;
    out.z = -0.75*J_norm*(((*system).Seff)*xi-1.)*(roots.z-roots.y)*(roots.z-roots.y)*xi_6/((*system).nu);

    return out;
}

/**
 * "d_0", "d_2" and "d_4" from 1703.03967
 * corresponding to equations
 * B9, B10 and B11 in 1703.03967 respectively
 * out.x = d_0
 * out.y = d_2
 * out.z = d_4
 */
template <class T> vector3D<T> d_coeffs(const T L_norm, const T J_norm, const vector3D<T> roots)
{
    vector3D<T> out;

    out.x = -((L_norm+J_norm)*(L_norm+J_norm)-roots.z)*((L_norm-J_norm)*(L_norm-J_norm)-roots.z);
    out.y = 2.*(roots.y-roots.z)*(J_norm*J_norm+L_norm*L_norm-roots.z);
    out.z = -(roots.z-roots.y)*(roots.z-roots.y);

    return out;
}

/**
 * The derivative of the phase of the magnitude of S
 * Eqn 24 in arxiv:1703.03967
 */
template <class T> T psidot(const T xi, const T xi_2, const vector3D<T> roots, const sysprecquant<T> *system)
{
    const T xi_3 = xi_2*xi;
    const T xi_6 = xi_3*xi_3;

    /*
    dictionary from code variables to paper symbols arxiv:1703.03967
    roots.z = S_+^2 as computed in the 'Roots' function
    roots.x = S_3^2 as computed in the 'Roots' function
    (*system).nu = symmetric mass-ratio
    (*system).Seff = "$\xi$" (eq. 7) i.e. the effective spin,
    which is NOT the same as the first argument to this function which is
    confusingly also called 'xi'
    xi and xi_6 are the v and v^6 variables in eq. B1 (in appendix B)

    the numerical factor is 0.75 which comes from a factor of 3/2 in eq.B1
    and a factor of 1/2 in eq. 24
    */

    return -0.75/sqrt((*system).nu)*(1.-(*system).Seff*xi)*xi_6*sqrt(roots.z-roots.x);
}

/**
 * Computes the MS corrections for phiz and zeta
 * eq. 67 (for phiz) and F19 (for zeta) in arxiv:1703.03967
 */
template <class T> vector3D<T> computeMScorrections(const T xi, const T xi_2,
    const T L_norm, const T J_norm,
    const vector3D<T> roots, const sysprecquant<T> *system)
{
    vector3D<T> out;
    const vector3D<T> c_return = c_coeffs(xi,xi_2,J_norm,roots,system);
    const vector3D<T> d_return = d_coeffs(L_norm,J_norm,roots);
    const T c0 = c_return.x; //eq. B6 in 1703.03967
    const T c2 = c_return.y; //eq. B7 in 1703.03967
    const T c4 = c_return.z; //eq. B8 in 1703.03967
    const T d0 = d_return.x; //eq. B9 in 1703.03967
    const T d2 = d_return.y; //eq. B10 in 1703.03967
    const T d4 = d_return.z; //eq. B11 in 1703.03967

    //"s_d" in the paper: eq. B20 in 1703.03967
    const T sqt = sqrt(fabs(d2 * d2 - 4. * d0 * d4));

    //related to eq. F20 in 1703.03967
    const T Aa = 0.5*(J_norm/L_norm + L_norm/J_norm - roots.z/J_norm/L_norm);
    //eq. F21 in 1703.03967
    const T Bb = 0.5*(roots.z - roots.y)/L_norm/J_norm;

    const T twod0_d2 = 2*d0+d2;
    const T twod0d4 = 2*d0*d4;
    const T d2_twod4 = d2+2.*d4;
    const T d2_d4 =d2+d4;
    const T d0_d2_d4 = d0+d2+d4;

    // "n3" in the code is "n_c" in the paper eq. B16 in 1703.03967
    const T n3 = (2.*d0_d2_d4/(twod0_d2+sqt));
    // "n4" in the code is "n_d" in the paper eq. B17 in 1703.03967
    const T n4 = ((twod0_d2+sqt)/(2.*d0));
    const T sqtn3 = sqrt(fabs(n3));
    const T sqtn4 = sqrt(fabs(n4));
    const T den = 2.*d0*d0_d2_d4*sqt;

    const T u = u_of_xi(xi,xi_2,system)+(*system).constant_of_S;
    const T tanu = tan(u);
    const T atanu = atan(tan(u));
    const T psdot = psidot(xi,xi_2,roots,system);

    T L3, L4;
    //L3 is XX (1703.03967)
    if (n3 == 1){
        L3 = 0;
    } else {
        L3 = fabs((c4 * d0 * (twod0_d2 + sqt) - c2 * d0 * (d2_twod4 - sqt) - c0 * (twod0d4 - d2_d4 * (d2 - sqt))) / (den)) * (sqtn3 / (n3 - 1.) * (atanu - atan(sqtn3 * tanu))) / psdot;
    }

    //L4 is XX (1703.03967)
    if (n4 == 1){
        L4 = 0;
    } else {
        L4 = fabs((-c4 * d0 * (twod0_d2 - sqt) + c2 * d0 * (d2_twod4 + sqt) - c0 * (-twod0d4 + d2_d4 * (d2 + sqt)))) / (den) * (sqtn4 / (n4 - 1.) * (atanu - atan(sqtn4 * tanu))) / psdot;
    }

    out.x = L3+L4;
    if (out.x != out.x) out.x=0;

    // eq. F19 in 1703.03967
    out.y = Aa*out.x + 2.*Bb*d0*(L3/(sqt-d2)-L4/(sqt+d2));
    if (out.y != out.y) out.y=0;

    out.z = 0;

    return out;
}

/**
 * phiz
 * equation 66 in 1703.03967
 */
template <class T> T phiz_of_xi(const T xi, const T xi_2, const T J_norm, const sysprecquant<T> *system)
{
    const T inside_log1 = ((*system).c_1) + J_norm*((*system).nu)+((*system).nu_2)/xi;
    // eq. D28 in 1703.03967
    const T log1 = log(fabs(inside_log1));
    const T inside_log2 = ((*system).c_1) + J_norm*((*system).sqrtSsqave)*xi + ((*system).Ssqave)*xi;
    // eq. D29 in 1703.03967
    const T log2 = log(fabs(inside_log2));

    // eq. D22 in 1703.03967
    const T phiz0 = J_norm*(0.5*((*system).c1_2)/((*system).nu_4) - 0.5*((*system).onethird)*((*system).c_1)/xi/((*system).nu_2) - ((*system).onethird)*((*system).Ssqave)/((*system).nu_2) - ((*system).onethird)/xi_2) - 0.5*((*system).c_1)*(((*system).c1_2)/((*system).nu_4) - ((*system).Ssqave)/((*system).nu_2))/((*system).nu)*log1;
    // eq. D23 in 1703.03967
    const T phiz1 = -J_norm*0.5*(((*system).c_1)/((*system).nu_2) + 1./xi) + 0.5*(((*system).c1_2)/((*system).nu_2) - ((*system).Ssqave))/((*system).nu)*log1 ;
    // eq. D24 in 1703.03967
    const T phiz2 = -J_norm + ((*system).sqrtSsqave)*log2 - ((*system).c_1)/((*system).nu)*log1;
    // eq. D25 in 1703.03967
    const T phiz3 = J_norm*xi -((*system).nu)*log1 + ((*system).c_1)/((*system).sqrtSsqave)*log2;
    // eq. D26 in 1703.03967
    const T phiz4 = J_norm*0.5*xi*(((*system).c_1)/((*system).Ssqave) + xi) - 0.5*(((*system).c1_2)/((*system).Ssqave) - ((*system).nu_2))/((*system).sqrtSsqave)*log2;
    // eq. D27 in 1703.03967
    const T phiz5 = J_norm*xi*(-0.5*((*system).c1_2)/((*system).Ssqave)/((*system).Ssqave) + 0.5*((*system).onethird)*((*system).c_1)*xi/((*system).Ssqave) + ((*system).onethird)*(xi_2 + ((*system).nu_2)/((*system).Ssqave))) + 0.5*((*system).c_1)*(((*system).c1_2)/((*system).Ssqave) - ((*system).nu_2))/((*system).Ssqave)/((*system).sqrtSsqave)*log2;

    // perform the summation in eq. 66 in 1703.03967
    T phizout = phiz0*(*system).constants_phiz[0] + phiz1*(*system).constants_phiz[1] + phiz2*(*system).constants_phiz[2] + phiz3*(*system).constants_phiz[3] + phiz4*(*system).constants_phiz[4] + phiz5*(*system).constants_phiz[5] + (*system).phiz_0;
    if (phizout!=phizout) phizout = 0;

    return phizout;
}

/**
 * zeta
 * eq. F5 in 1703.03967
 */
template <class T> T zeta_of_xi(const T xi, const T xi_2, const sysprecquant<T> *system)
{
    const T logxi = log(xi);
    const T xi_3 = xi_2*xi;

    // summation in eq. F5 in 1703.03967
    T zetaout = ((*system).nu)*(-(*system).onethird*(*system).constants_zeta[0]/xi_3 - 0.5*(*system).constants_zeta[1]/xi_2 - (*system).constants_zeta[2]/xi + (*system).constants_zeta[3]*logxi + (*system).constants_zeta[4]*xi  + 0.5*(*system).constants_zeta[5]*xi_2) + (*system).zeta_0;
    if (zetaout!=zetaout) zetaout = 0;

    return zetaout;
}

/**
 * Computes phiz, zeta, and costhetaL at 3PN
 */
template <class T> vector3D<T> compute_phiz_zeta_costhetaL3PN(const T xi, const sysprecquant<T> *system)
{
    vector3D<T> angles;
    const T xi_2 = xi*xi;
    const T L_norm = ((*system).nu)/xi;
    const T L_norm3PN = L_norm_3PN_of_xi(xi,xi_2,L_norm,system);
    const T J_norm3PN = J_norm_of_xi(L_norm3PN,system);
    const T J_norm = J_norm_of_xi(L_norm,system);
    const vector3D<T> roots = Roots(L_norm,J_norm,system);
    const T S_norm = S_norm_of_xi(xi,xi_2,roots,system);

    vector3D<T> MScorrections = {0.,0.,0.};
    if(fabs(roots.y-roots.z)>1.e-5){
        MScorrections = computeMScorrections(xi,xi_2,L_norm,J_norm,roots,system);
    }

    angles.x = phiz_of_xi(xi,xi_2,J_norm,system) + MScorrections.x;
    angles.y = zeta_of_xi(xi,xi_2,system) + MScorrections.y;
    angles.z = costhetaL(J_norm3PN,L_norm3PN,S_norm);

    return angles;
}

/**
 * Simple 2PN version of the orbital angular momentum L,
 * without any spin terms expressed as a function of v.
 *
 *  Reference:
 *  - Bohe et al, 1212.5520v2 Eq 4.7 first line
 */
template <class T> T L2PNR(
  const T v,   /**< Cubic root of (Pi * Frequency (geometric)) */
  const T eta) /**< Symmetric mass-ratio */
{
  const T eta2 = eta*eta;
  const T x = v*v;
  const T x2 = x*x;
  return (eta*(1.0 + (1.5 + eta/6.0)*x + (3.375 - (19.0*eta)/8. - eta2/24.0)*x2)) / sqrt(x);
}

/**
 * Compute the orbital angular momentum at 3PN order at a single frequency.
 * We assume that Lhat = (0,0,1)
 */
template <class T> T L3PN(
    const T f_orb_hz,   /**< Orbital frequency (Hz)  */
    const T m1,         /**< mass of primary in SI (kg) */
    const T m2,         /**< mass of secondary in SI (kg) */
    const T s1x,        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const T s1y,        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const T s1z,        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const T s2x,        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const T s2y,        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const T s2z,        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const T f_0,        /**< Reference Gravitational Wave frequency (Hz) */
    const int ExpansionOrder /**< Keep terms upto ExpansionOrder in precession angles phi_z and zeta (1,2,3,4,5 or -1 for all orders) */
)
{
    const T mul = 1.0; /* Cosine of Polar angle of orbital angular momentum */
    const T phl = 0.0; /* Azimuthal angle of orbital angular momentum */
    T mu1 = 0.0;
    T ph1 = 0.0;
    T ch1 = 0.0;
    T mu2 = 0.0;
    T ph2 = 0.0;
    T ch2 = 0.0;

    CartesianToPolar(&mu1, &ph1, &ch1, s1x, s1y, s1z);
    CartesianToPolar(&mu2, &ph2, &ch2, s2x, s2y, s2z);

    T cosmu1 = cos(mu1);
    T cosmu2 = cos(mu2);
    return OrbitalAngMom3PNSpinning(
            f_orb_hz,
            m1, m2,
            mul, phl,
            cosmu1, ph1, ch1,
            cosmu2, ph2, ch2,
            f_0, ExpansionOrder);
}

/**
 * atan2 wrapper that returns 0 when both magnitudes of x and y are
 * below tol, otherwise it returns atan2(x, y)
 */
template <class T> T atan2tol(T a, T b, T tol)
{
    T out;
    if (fabs(a) < tol && fabs(b) < tol)
    {
        out = 0.;
    }
    else
    {
        out = atan2(a, b);
    }

    return out;
}

/**
 * Compute the polar coordinates of a cartesian point (x,y,z)
*/
template <class T> void CartesianToPolar(
    T *polar,
    T *azimuthal,
    T *magnitude,
    T x,
    T y,
    T z
)
{
    T zero = 0.;
    T eps = 1e-9;
    T tol_atan = MAX_TOL_ATAN;

    *magnitude = sqrt( x*x + y*y + z*z );
    if (approx_equal(*magnitude, zero, eps)){
        *polar = zero;
        *azimuthal = zero;
    } else {
        T along_z = acos(z / *magnitude);
        *polar = along_z;
        *azimuthal = atan2tol(y, x, tol_atan);
    }
}

/**
 * \brief Calculate the Incomplete Elliptic Integral
*/
template <class T> T ellint_F(T phi, T k)
{
    return gsl_sf_ellint_F(phi, k, GSL_PREC_DOUBLE);
}

template <> adouble ellint_F<adouble>(adouble phi, adouble k)
{
    double phi_val = phi.value();
    double k_val = k.value();

    double out_val = ellint_F(phi_val, k_val);

    return (adouble)out_val;
}

/*! \brief Initialize quantities needed for Pv3 precession
 *
 * Adapted from LALSimIMRPhenomPv3HM
 */
template <class T> void InitializePrecession(sysprecquant<T>* system, /** [out] Pointer to sysprecquant struct */
    const T m1,  /**< Primary mass in SI (kg) */
    const T m2,  /**< Secondary mass in SI (kg) */
    const T mul, /**< Cosine of Polar angle of orbital angular momentum */
    const T phl, /**< Azimuthal angle of orbital angular momentum  */
    const T mu1, /**< Cosine of Polar angle of primary spin w.r.t. orbital angular momentum */
    const T ph1, /**< Azimuthal angle of primary spin  */
    const T ch1, /**< Dimensionless spin magnitude of primary spin */
    const T mu2, /**< Cosine of Polar angle of secondary spin w.r.t. orbital angular momentum */
    const T ph2, /**< Azimuthal angle of secondary spin  */
    const T ch2, /**< Dimensionless spin magnitude of secondary spin */
    const T f_0, /**< Reference Gravitational Wave frequency (Hz) */
    const int ExpansionOrder /**< Keep terms upto ExpansionOrder in precession angles phi_z and zeta */
)
{
    system->onethird = 1./3.;

    const T domegadt_csts_nonspin[17] = {96./5.,-1486./35.,-264./5.,384.*M_PI/5.,34103./945.,13661./105.,944./15.,M_PI*(-4159./35.),M_PI*(-2268./5.),(16447322263./7276500. + M_PI*M_PI*512./5. - GWAT_LN2*109568./175. -gamma_E*54784./175.),(-56198689./11340. + M_PI*M_PI*902./5.),1623./140.,-1121./27.,-54784./525.,-M_PI*883./42.,M_PI*71735./63.,M_PI*73196./63.};
    const T domegadt_csts_spinorbit[18] = {-904./5.,-120.,-62638./105.,4636./5.,-6472./35.,3372./5.,-M_PI*720.,-M_PI*2416./5.,-208520./63.,796069./105.,-100019./45.,-1195759./945.,514046./105.,-8709./5.,-M_PI*307708./105.,M_PI*44011./7.,-M_PI*7992./7.,M_PI*151449./35.};
    const T domegadt_csts_spinspin[4] = {-494./5.,-1442./5.,-233./5.,-719./5.};
    const T L_csts_nonspin[9] = {3./2.,1./6.,27./8.,-19./8.,1./24.,135./16.,-6889/144.+ 41./24.*M_PI*M_PI,31./24.,7./1296.};
    const T L_csts_spinorbit[6] = {-14./6.,-3./2.,-11./2.,133./72.,-33./8.,7./4.};

    const T G_over_c = GWAT_G_SI/c;
    const T M = m1 + m2;
    system->q = m2/m1;
    const T q = system->q;
    system->nu = q/(1+q)/(1+q);
    const T nu = system->nu;
    system->nu_2 = nu*nu;
    const T nu_2 = system->nu_2;
    const T piGM_over_cthree = M_PI*G_over_c/c/c * (m1+m2);
    system->deltam_over_M = (1-q)/(1+q);
    const T deltam_over_M = system->deltam_over_M;

    /* Set up initial vectors and values */
    const T S1_norm = fabs(ch1)*nu/q;//in geometric units and with M=1
    const T S2_norm = fabs(ch2)*nu*q;

    const T xi_0 = pow(piGM_over_cthree*f_0, system->onethird);
    const T xi0_2 = xi_0*xi_0;
    T acosmul = acos(mul);
    T acosmu1 = acos(mu1);
    T acosmu2 = acos(mu2);
    const T one = 1.;
    const vector3D<T> Lhat_0 = CreateSphVector(one,acosmul,phl);
    const vector3D<T> S1_0 = CreateSphVector(S1_norm,acosmu1,ph1);//in geometric units and with M=1
    const vector3D<T> S2_0 = CreateSphVector(S2_norm,acosmu2,ph2);
    T Lnorm_0 = nu/xi_0;
    const vector3D<T> L_0 = ScaleVector(Lnorm_0, Lhat_0);

    system->dot1 = DotProduct(S1_0,Lhat_0);
    system->dot2 = DotProduct(S2_0,Lhat_0);
    system->dot1n = system->dot1/S1_norm;
    system->dot2n = system->dot2/S2_norm;
    system->dot12 = DotProduct(S1_0,S2_0);

    system->Seff = ((system->dot1)/m1 +(system->dot2)/m2)*M;

    const vector3D<T> S_0 = VectorSum(S1_0,S2_0);
    const vector3D<T> J_0 = VectorSum(L_0,S_0);

    const T S_0_norm = VectorNorm(S_0);
    const T J_0_norm = VectorNorm(J_0);
    const T L_0_norm = VectorNorm(L_0);
    system->S1_norm_2 = S1_norm*S1_norm;
    system->S2_norm_2 = S2_norm*S2_norm;
    system->S0_norm = S_0_norm;

    const vector3D<T> roots = Roots(L_0_norm,J_0_norm,system);

    const T q_2 = q*q;
    const T one_m_q_sq = (1.-q)*(1.-q);
    const T one_m_q_4 = one_m_q_sq*one_m_q_sq;
    const T one_p_q_sq = (1.+q)*(1.+q);

    const T Save_square = 0.5*(roots.z+roots.y);
    const T S1_norm_2 = system->S1_norm_2;
    const T S2_norm_2 = system->S2_norm_2;
    const T Seff = system->Seff;
    const T Seff_2 = Seff*Seff;

    //computed with initial spin couplings, they're not exactly accurate for generic precession, but the correction should be 4PN
    //these constants are used in TaylorT1 where domega/dt is expressed as a polynomial
    T beta_a = domegadt_csts_spinorbit[2] + nu*(domegadt_csts_spinorbit[3]);
    T beta_b = domegadt_csts_spinorbit[4] + nu*(domegadt_csts_spinorbit[5]);
    const T a0 = nu*domegadt_csts_nonspin[0];
    const T a2 = nu*(domegadt_csts_nonspin[1] + nu*(domegadt_csts_nonspin[2]));
    const T a3 = nu*(domegadt_csts_nonspin[3] + beta(domegadt_csts_spinorbit[0], domegadt_csts_spinorbit[1], system));
    const T a4 = nu*(domegadt_csts_nonspin[4] + nu*(domegadt_csts_nonspin[5] + nu*(domegadt_csts_nonspin[6])) + sigma(domegadt_csts_spinspin[0], domegadt_csts_spinspin[1], system) + tau(domegadt_csts_spinspin[2], domegadt_csts_spinspin[3], system));
    const T a5 = nu*(domegadt_csts_nonspin[7] + nu*(domegadt_csts_nonspin[8]) + beta(beta_a, beta_b, system));

    const T a0_2 = a0*a0;
    const T a0_3 = a0_2*a0;
    const T a2_2 = a2*a2;

    //these constants are used in TaylorT2 where domega/dt is expressed as an inverse polynomial
    //eq.A1 (1703.03967)
    const T c0 = 1./a0;
    //eq.A2 (1703.03967)
    const T c2 = -a2/a0_2;
    //eq.A3 (1703.03967)
    const T c3 = -a3/a0_2;
    //eq.A4 (1703.03967)
    const T c4 = (a2_2 - a0*a4)/a0_3;
    //eq.A5 (1703.03967)
    const T c5 = (2.*a2*a3 - a0*a5)/a0_3;

    system->nu_4 = nu_2*nu_2;
    const T nu_4 = system->nu_4;

    //eq.41 (1703.03967)
    system->c_1 =0.5*(J_0_norm*J_0_norm - L_0_norm*L_0_norm - Save_square)/L_0_norm*nu;
    T c_1 = system->c_1;
    system->c_1_over_nu = system->c_1/nu;
    T c_1_over_nu = system->c_1_over_nu;
    system->c1_2 = c_1*c_1;
    T c1_2 = system->c1_2;
    T c1_over_nu_2 = c_1_over_nu*c_1_over_nu;

    const T one_m_q2_2 = (1. - q_2) * (1. - q_2);

    // Calculate the Jacobi sine phase
    //eq.C3 (1703.03967) - note this is actually 2*\Delta
    const T Del1 = 4. * c1_over_nu_2 * one_p_q_sq;
    const T Del2 = 8. * c_1_over_nu * q * (1. + q) * Seff;
    const T Del3 = 4. * (one_m_q2_2 * S1_norm_2 - q_2 * Seff_2);
    const T Del4 = 4. * c1_over_nu_2 * q_2 * one_p_q_sq;
    const T Del5 = 8. * c_1_over_nu * q_2 * (1. + q) * Seff;
    const T Del6 = 4. * (one_m_q2_2 * S2_norm_2 - q_2 * Seff_2);
    const T Delta = sqrt((Del1 - Del2 - Del3) * (Del4 - Del5 - Del6));

    //this is g_0 in eq.51 (1703.03967)
    system->constants_u[0] = -c0;
    //eq.C1 (1703.03967)
    system->constants_u[1] = (6.*Seff*nu - 3.*c_1_over_nu)/deltam_over_M/deltam_over_M;

    const T u1 = 3. * c2 / c0;
    const T u2 = 0.75 * one_p_q_sq / one_m_q_4;
    const T u3 = -20. * c1_over_nu_2 * q_2 * one_p_q_sq;
    const T u4 = 2. * one_m_q2_2 * (q * (2. + q) * S1_norm_2 + (1. + 2. * q) * S2_norm_2 - 2. * q * Save_square);
    const T u5 = 2. * q_2 * (7. + 6. * q + 7. * q_2) * 2. * c_1_over_nu * Seff;
    const T u6 = 2. * q_2 * (3. + 4. * q + 3. * q_2) * Seff_2;
    const T u7 = q * Delta;
    //eq.C2 (1703.03967)
    system->constants_u[2] = u1 + u2*(u3 + u4 + u5 - u6 + u7);

    //eq.45 (1703.03967)
    system->Ssqave = 0.5*(roots.z+roots.y);
    system->sqrtSsqave = sqrt(system->Ssqave);

    // Calculate phi_z
    //eq.D1 (1703.03967)
    const T Rm = roots.z - roots.y;
    const T Rm_2 = Rm*Rm;
    //eq.D2 (1703.03967)
    const T cp = roots.z*nu_2 - c1_2;
    //eq.D3 (1703.03967)
    const T cm = cp-Rm*nu_2;
    //difference of spin norm squared used in eq.D6 (1703.03967)
    const T S0m = S1_norm_2 - S2_norm_2;
    const T cpcm = fabs(cp*cm);
    const T sqrt_cpcm = sqrt(cpcm);

    //eq.D4 (1703.03967)
    const T A1t = 0.5+0.75/nu;//xi^6
    //eq.D5 (1703.03967)
    const T A2t = -0.75*Seff/nu;//xi^7
    //eq.E3, called $D_2$ in paper (1703.03967)
    const T A1ave = (cp-sqrt_cpcm)/nu_2 ;
    //eq.E4, called $D_4$ in paper (1703.03967)
    const T Bave = -0.5*Rm*sqrt_cpcm/nu_2 - cp/nu_4*(sqrt_cpcm-cp);

    const T aw = (-3.*(1. + q)/q*(2.*(1. + q)*nu_2*Seff*c_1 - (1. + q)*c1_2 + (1. - q)*nu_2*S0m));
    const T cw = 3./32./nu*Rm_2;
    const T dw = 4.*cp - 4.*A1ave*nu_2;
    const T hw = -2*(2*A1ave - Rm)*c_1;
    const T fw = Rm*A1ave-Bave-0.25*Rm_2;

    //eq.D6 (1703.03967)
    const T ad = aw/dw;
    //eq.XX (1703.03967)
    const T hd = hw/dw;
    //eq.D7 (1703.03967)
    const T cd = cw/dw;
    //eq.XX (1703.03967)
    const T fd = fw/dw;

    const T hd_2 = hd*hd;
    const T hd_4 = hd_2*hd_2;
    const T adfd = ad*fd;
    const T adfdhd = adfd*hd;
    const T adfdhd_2 = adfd*hd_2;
    const T adhd = ad*hd;
    const T adhd_2 = ad*hd_2;
    const T adhd_3 = ad*hd_2*hd;
    const T adhd_4 = ad*hd_4;
    const T cdfd = cd*fd;
    const T cdhd = cd*hd;
    const T cdhd_2 = cd*hd_2;

    //eq.D10 (1703.03967)
    T Omegaz0 = A1t + ad;
    //eq.D11 (1703.03967)
    T Omegaz1 = A2t - ad*Seff - adhd;
    //eq.D12 (1703.03967)
    T Omegaz2 = cd - adfd + adhd_2 + adhd*Seff;
    //eq.D13 (1703.03967)
    T Omegaz3 = (adfd - cd - adhd_2)*Seff + 2*adfdhd - adhd_3 - cdhd;
    //eq.D14 (1703.03967)
    T Omegaz4 = -(2*adfdhd - adhd_3 - cdhd)*Seff + adfd*fd - cdfd + cdhd_2 - 3*adfdhd_2 + adhd_4;
    //eq.D15 (1703.03967)
    T Omegaz5 = -(adfd*fd - cdfd + cdhd_2 - 3*adfdhd_2 + adhd_4)*Seff + hd*(2*cdfd - 3*adfd*fd - cdhd_2 + 4*adfdhd_2 - adhd_4);

    /* We check the size of the Omegaz5 coefficient to try and catch
    cases where we think the precession model is breaking down */
    // Commented out as it isn't really significant.
    // TODO: Activate the check with a flag.
    //try
    //{
    //    checkOmegaz5(Omegaz5);
    //}
    //catch(const std::exception& e)
    //{
    //    std::cerr << e.what() << '\n';
    //    std::cerr << "f0 = " << f_0 << std::endl;
    //    std::cerr << "(m1, cos(s1z), chi1) = ";
    //    std::cerr << "(" << (m1/GWAT_MSUN_SI) << ", " << mu1 << ", " << ch1 << ")\n";
    //    std::cerr << "(m2, cos(s2z), chi2) = ";
    //    std::cerr << "(" << (m2/GWAT_MSUN_SI) << ", " << mu2 << ", " << ch2 << ")\n";
    //}
    
    //eq.D16 (1703.03967), note that "c0" in the code is "$g_0$" in the paper and likewise for the other "c's" and "g's"
    system->constants_phiz[0] = 3.*Omegaz0*c0;
    //eq.D17 (1703.03967)
    system->constants_phiz[1] = 3.*Omegaz1*c0;
    //eq.D18 (1703.03967)
    system->constants_phiz[2] = 3.*(Omegaz2*c0 + Omegaz0*c2);
    //eq.D19 (1703.03967)
    system->constants_phiz[3] = 3.*(Omegaz0*c3 + Omegaz1*c2 + Omegaz3*c0);
    //eq.D20 (1703.03967)
    system->constants_phiz[4] = 3.*(Omegaz0*c4 + Omegaz1*c3 + Omegaz2*c2 + Omegaz4*c0);
    //eq.D21 (1703.03967)
    system->constants_phiz[5] = 3.*(c5*Omegaz0 + c4*Omegaz1 + c3*Omegaz2 + c2*Omegaz3 + c0*Omegaz5);

    const T gw = 3./16./nu_2/nu*Rm_2*(c_1 - nu_2*Seff);
    //eq.F18 (1703.03967)
    const T gd = gw/dw;

    //eq.F17 (1703.03967)
    Omegaz5 += Omegaz4*c_1/nu_2 - fd*gd + gd*hd_2 + gd*hd*Seff;
    //eq.F16 (1703.03967)
    Omegaz4 += Omegaz3*c_1/nu_2 - gd*hd - gd*Seff;
    //eq.F15 (1703.03967)
    Omegaz3 += Omegaz2*c_1/nu_2 + gd;
    //eq.F14 (1703.03967)
    Omegaz2 += Omegaz1*c_1/nu_2;
    //eq.F13 (1703.03967)
    Omegaz1 += Omegaz0*c_1/nu_2;

    //eq.F12 (1703.03967)
    system->constants_zeta[0] = 3.*Omegaz0*c0;
    //eq.F13 (1703.03967)
    system->constants_zeta[1] = 3.*Omegaz1*c0;
    //eq.F14 (1703.03967)
    system->constants_zeta[2] = 3.*(Omegaz2*c0 + Omegaz0*c2);
    //eq.F15 (1703.03967)
    system->constants_zeta[3] = 3.*(Omegaz0*c3 + Omegaz1*c2 + Omegaz3*c0);
    //eq.F16 (1703.03967)
    system->constants_zeta[4] = 3.*(Omegaz0*c4 + Omegaz1*c3 + Omegaz2*c2 + Omegaz4*c0);
    //eq.F17 (1703.03967)
    system->constants_zeta[5] = 3.*(Omegaz0*c5 + Omegaz1*c4 + Omegaz2*c3 + Omegaz3*c2 + Omegaz5*c0);

    switch (ExpansionOrder)
    {
    case -1:
        /* default case, keep all orders. */
        break;
    case 1:
        /* Only keep 1st order term, i.e. only keep system->constants_phiz[0] and system->constants_zeta[0] */
        system->constants_phiz[1] *= 0.;
        system->constants_zeta[1] *= 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
    case 2:
        system->constants_phiz[2] *= 0.;
        system->constants_zeta[2] *= 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
    case 3:
        system->constants_phiz[3] *= 0.;
        system->constants_zeta[3] *= 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
    case 4:
        system->constants_phiz[4] *= 0.;
        system->constants_zeta[4] *= 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
    case 5:
        system->constants_phiz[5] *= 0.;
        system->constants_zeta[5] *= 0.;
        break;
    default:
        invalidExpansionOrder(ExpansionOrder);
    }

    T m, B, signelement;
    int sign_num;

    // I think 'constant_of_S' is the "$\psi_0$" term in equation 51 (1703.03967)
    if(fabs(roots.y-roots.z)<1.e-5)
    {
        system->constant_of_S = 0;
    }
    else
    {
        T asin_sign;
        m = sqrt((roots.y - roots.z)/(roots.x - roots.z));
        B = (S_0_norm*S_0_norm-roots.z)/(roots.y-roots.z);
        signelement = DotProduct(CrossProduct(L_0,S1_0),S2_0);
        sign_num = (signelement > 0) - (signelement < 0);

        if(B < 0. || B > 1.)
        {
            if(B > 1 && B-1. < 0.00001)
            {
                asin_sign = asin(sign_num*sqrt(1.));
                system->constant_of_S = ellint_F(asin_sign, m) - u_of_xi(xi_0,xi0_2,system);
            }
            if(B < 0 && B > -0.00001)
            {
                asin_sign = asin(sign_num*sqrt(0.));
                system->constant_of_S = ellint_F(asin_sign, m) - u_of_xi(xi_0,xi0_2,system);
            }
        }
        else
        {
            asin_sign = asin(sign_num*sqrt(B));
            system->constant_of_S = ellint_F(asin_sign, m) - u_of_xi(xi_0,xi0_2,system);
        }
    }

    beta_a = L_csts_spinorbit[2]+L_csts_spinorbit[3]*nu;
    beta_b = L_csts_spinorbit[4]+L_csts_spinorbit[5]*nu;
    system->constants_L[0] = (L_csts_nonspin[0] + nu*L_csts_nonspin[1]);
    system->constants_L[1] = beta(L_csts_spinorbit[0], L_csts_spinorbit[1], system);
    system->constants_L[2] = (L_csts_nonspin[2] + nu*L_csts_nonspin[3] + nu*nu*L_csts_nonspin[4]);
    system->constants_L[3] = beta(beta_a, beta_b, system);
    system->constants_L[4] = (L_csts_nonspin[5]+L_csts_nonspin[6]*nu +L_csts_nonspin[7]*nu*nu+L_csts_nonspin[8]*nu*nu*nu);

    vector3D<T> MScorrections = {0.,0.,0.};
    if(fabs(roots.y-roots.z)>1.e-5){
        MScorrections = computeMScorrections(xi_0,xi0_2,L_0_norm,J_0_norm,roots,system);
    }

    system->phiz_0 = 0.;
    system->phiz_0 = - phiz_of_xi(xi_0,xi0_2,J_0_norm,system) - MScorrections.x;
    system->zeta_0 = 0.;
    system->zeta_0 = - zeta_of_xi(xi_0,xi0_2,system) - MScorrections.y;
}

/**
 * Set the parameters for IMRPhenomPv3 source_parameters
*/
template <class T> void PhenomPv3_Param_Transform(source_parameters<T> *out, gen_params_base<T> *in)
{
    // chi1L, chi2L for intermediate calculations
    T chi1L, chi2L;

    // Define parameters used by Pv3
    PhenomP_ParametersFromSourceFrame(&chi1L, &chi2L,
        &(out->chip), &(out->thetaJN), &(out->alpha0),
        &(out->phi_aligned), &(out->zeta_polariz),
        in->mass1, in->mass2, in->f_ref, in->phiRef, in->incl_angle,
        in->spin1[0], in->spin1[1], in->spin1[2],
        in->spin2[0], in->spin2[1], in->spin2[2],
        IMRPhenomPv3_V);

    // Set parameters used by source_parameters
    // to avoid issues (if any) down the line
    T q = in->mass1/in->mass2;
    T m1_M_2 = q/(1.+q);
    T m2_M_2 = m1_M_2 / q;
    // squared masses
    m1_M_2 *= m1_M_2; 
    m2_M_2 *= m2_M_2;

    out->chil = chi1L + chi2L/q;
    out->SP = out->chip * m1_M_2;
    out->SL = chi1L * m1_M_2 + chi2L * m2_M_2;
}

/**
 * Function to map parameters
 * (masses, 6 spin components, phiRef and inclination at f_ref)
 * (assumed to be in the source frame
 *  where LN points in the z direction
 *  i.e. lnhat = (0,0,1)
 *  and the separation vector n is in the x direction
 *  and the spherical angles of the line of sight N are (incl,Pi/2-phiRef))
 * into IMRPhenomP intrinsic parameters
 * (chi1_l, chi2_l, chip, thetaJN, alpha0 and phi_aligned).
 *
 * See Fig. 1. in arxiv:1408.1810 for a diagram of the angles.
 */
template <class T> void PhenomP_ParametersFromSourceFrame(
    T *chi1_l,                  /**< [out] Dimensionless aligned spin on companion 1 */
    T *chi2_l,                  /**< [out] Dimensionless aligned spin on companion 2 */
    T *chip,                    /**< [out] Effective spin in the orbital plane */
    T *thetaJN,                  /**< [out] Angle between J0 and line of sight (z-direction) */
    T *alpha0,                  /**< [out] Initial value of alpha angle (azimuthal precession angle) */
    T *phi_aligned,                  /**< [out] Initial phase to feed the underlying aligned-spin model */
    T *zeta_polariz,                  /**< [out] Angle to rotate the polarizations */
    const T m1,              /**< Mass of companion 1 (solar masses) */
    const T m2,              /**< Mass of companion 2 (solar masses) */
    const T f_ref,              /**< Reference GW frequency (Hz) */
    const T phiRef,              /**< Reference phase */
    const T incl,              /**< Inclination : angle between LN and the line of sight */
    const T s1x,                /**< Initial value of s1x: dimensionless spin of BH 1 */
    const T s1y,                /**< Initial value of s1y: dimensionless spin of BH 1 */
    const T s1z,                /**< Initial value of s1z: dimensionless spin of BH 1 */
    const T s2x,                /**< Initial value of s2x: dimensionless spin of BH 2 */
    const T s2y,                /**< Initial value of s2y: dimensionless spin of BH 2 */
    const T s2z,                /**< Initial value of s2z: dimensionless spin of BH 2 */
    IMRPhenomP_version_type IMRPhenomP_version /**< IMRPhenomP(v1) uses IMRPhenomC, IMRPhenomPv2 uses IMRPhenomD, IMRPhenomPv2_NRTidal uses NRTidal framework with IMRPhenomPv2 */
)
{
    // Note that the angle phiJ defined below and alpha0 are degenerate. Therefore we do not output phiJ.
    
    const T M = m1+m2; /* Masses in solar masses */
    const T Msq = M*M;
    const T m1_2 = m1*m1;
    const T m2_2 = m2*m2;
    const T eta = m1 * m2 / (M*M);    /* Symmetric mass-ratio */

    /* From the components in the source frame, we can easily determine
    chi1_l, chi2_l, chip and phi_aligned, which we need to return.
    We also compute the spherical angles of J,
    which we need to transform to the J frame*/

    /* Aligned spins */
    *chi1_l = s1z; /* Dimensionless aligned spin on BH 1 */
    *chi2_l = s2z; /* Dimensionless aligned spin on BH 2 */

    /* Magnitude of the spin projections in the orbital plane */
    const T S1_perp = m1_2*sqrt(s1x*s1x + s1y*s1y);
    const T S2_perp = m2_2*sqrt(s2x*s2x + s2y*s2y);
    /* From this we can compute chip*/
    const T A1 = 2 + (3*m2) / (2*m1);
    const T A2 = 2 + (3*m1) / (2*m2);
    const T ASp1 = A1*S1_perp;
    const T ASp2 = A2*S2_perp;
    /* chip = max(A1 Sp1, A2 Sp2) / (A_i m_i^2) for i index of larger BH */
    const T num = (ASp2 > ASp1) ? ASp2 : ASp1;
    const T den = (m2 > m1) ? A2*m2_2 : A1*m1_2;
    *chip = num / den;

    /* Compute L, J0 and orientation angles */
    const T m_sec = M * GWAT_MTSUN_SI;   /* Total mass in seconds */
    const T piM = M_PI * m_sec;
    const T v_ref = cbrt(piM * f_ref);

    const int ExpansionOrder = 5; // Used in PhenomPv3 only

    T L0 = 0.0;
    switch (IMRPhenomP_version)
    {
    case IMRPhenomPv2_V:
    case IMRPhenomPv2NRTidal_V:
        L0 = Msq * L2PNR(v_ref, eta); /* 2PN */
        break;
    case IMRPhenomPv3_V:
      if ((s1x == 0. && s1y == 0. && s2x == 0. && s2y == 0.))
      { // non-precessing case
        L0 = Msq * L2PNR(v_ref, eta);
      }
      else
      { // precessing case, use 3PN spinning approx
        T f_ref_orb = 0.5*f_ref;
        T m1_SI = m1*GWAT_MSUN_SI;
        T m2_SI = m2*GWAT_MSUN_SI;
        L0 = Msq * L3PN(f_ref_orb,
            m1_SI, m2_SI,
            s1x, s1y, s1z, s2x, s2y, s2z,
            f_ref, ExpansionOrder);
      }
      break;
    default:
        throw std::runtime_error("Unknown PhenomP version.\n");
        break;
    }

    // Below, _sf indicates source frame components. We will also use _Jf for J frame components
    const T J0x_sf = m1_2*s1x + m2_2*s2x;
    const T J0y_sf = m1_2*s1y + m2_2*s2y;
    const T J0z_sf = L0 + m1_2*s1z + m2_2*s2z;
    const T J0 = sqrt(J0x_sf*J0x_sf + J0y_sf*J0y_sf + J0z_sf*J0z_sf);

    /* Compute thetaJ, the angle between J0 and LN (z-direction) */
    T thetaJ_sf;
    if (J0 < 1e-10) {
        thetaJ_sf = 0;
    } else {
        thetaJ_sf = acos(J0z_sf / J0);
    }

    T phiJ_sf;
    if (fabs(J0x_sf) < MAX_TOL_ATAN && fabs(J0y_sf) < MAX_TOL_ATAN)
        phiJ_sf = M_PI/2. - phiRef; // aligned spin limit
    else
        phiJ_sf = atan2(J0y_sf, J0x_sf); /* azimuthal angle of J0 in the source frame */

    *phi_aligned = - phiJ_sf;

    /* We now have to rotate to the "J frame" where we can easily
   compute alpha0, the azimuthal angle of LN,
   as well as thetaJ, the angle between J and N.
   The J frame is defined imposing that J points in the z direction
   and the line of sight N is in the xz plane (with positive projection along x).
   The components of any vector in the (new) J frame are obtained from those
   in the (old) source frame by multiplying by RZ[kappa].RY[-thetaJ].RZ[-phiJ]
   where kappa will be determined by rotating N with RY[-thetaJ].RZ[-phiJ]
   (which brings J to the z axis) and taking the opposite of azimuthal angle of the rotated N.
   */
    T tmp1,tmp2;
    // First we determine kappa
    // in the source frame, the components of N are given in Eq (35c) of T1500606-v6
    T Nx_sf = sin(incl)*cos(M_PI_2 - phiRef);
    T Ny_sf = sin(incl)*sin(M_PI_2 - phiRef);
    T Nz_sf = cos(incl);
    T tmp_x = Nx_sf;
    T tmp_y = Ny_sf;
    T tmp_z = Nz_sf;
    ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
    ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
    T kappa;
    T atan_tol = MAX_TOL_ATAN;
    kappa = - atan2tol(tmp_y,tmp_x, atan_tol);

    // Then we determine alpha0, by rotating LN
    tmp_x = 0.;
    tmp_y = 0.;
    tmp_z = 1.; // in the source frame, LN=(0,0,1)
    ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
    ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
    ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
    if (fabs(tmp_x) < MAX_TOL_ATAN && fabs(tmp_y) < MAX_TOL_ATAN)
        *alpha0 = M_PI; //this is the aligned spin case
    else
        *alpha0 = atan2(tmp_y,tmp_x);

    // Finally we determine thetaJ, by rotating N
    tmp_x = Nx_sf;
    tmp_y = Ny_sf;
    tmp_z = Nz_sf;
    ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
    ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
    ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
    T Nx_Jf = tmp_x; // let's store those two since we will reuse them later (we don't need the y component)
    T Nz_Jf = tmp_z;
    *thetaJN = acos(Nz_Jf); // No normalization needed, we are dealing with a unit vector

    /* Finally, we need to redefine the polarizations :
    PhenomP's polarizations are defined following Arun et al (arXiv:0810.5336)
    i.e. projecting the metric onto the P,Q,N triad defined with P=NxJ/|NxJ| (see (2.6) in there).
    By contrast, the triad X,Y,N used in LAL
    ("waveframe" in the nomenclature of T1500606-v6)
    is defined in e.g. eq (35) of this document
    (via its components in the source frame; note we use the defautl Omega=Pi/2).
    Both triads differ from each other by a rotation around N by an angle \zeta
    and we need to rotate the polarizations accordingly by 2\zeta
    */
    T Xx_sf = -cos(incl)*sin(phiRef);
    T Xy_sf = -cos(incl)*cos(phiRef);
    T Xz_sf = sin(incl);
    tmp_x = Xx_sf;
    tmp_y = Xy_sf;
    tmp_z = Xz_sf;
    ROTATEZ(-phiJ_sf, tmp_x, tmp_y, tmp_z);
    ROTATEY(-thetaJ_sf, tmp_x, tmp_y, tmp_z);
    ROTATEZ(kappa, tmp_x, tmp_y, tmp_z);
    //now the tmp_a are the components of X in the J frame
    //we need the polar angle of that vector in the P,Q basis of Arun et al
    // P=NxJ/|NxJ| and since we put N in the (pos x)z half plane of the J frame
    T PArunx_Jf = 0.;
    T PAruny_Jf = -1.;
    T PArunz_Jf = 0.;
    // Q=NxP
    T QArunx_Jf = Nz_Jf;
    T QAruny_Jf = 0.;
    T QArunz_Jf = -Nx_Jf;
    T XdotPArun = tmp_x*PArunx_Jf+tmp_y*PAruny_Jf+tmp_z*PArunz_Jf;
    T XdotQArun = tmp_x*QArunx_Jf+tmp_y*QAruny_Jf+tmp_z*QArunz_Jf;
    *zeta_polariz = atan2(XdotQArun , XdotPArun);
}

/**
 * Precomputes useful quantities and populates the
 * PhenomPv3Storage and sysprecquant (for precession angles) structs.
 * Converts from GWAT seconds to SI.
 */
template <class T> void init_PhenomPv3_Storage(
    PhenomPv3Storage<T> *p,   /**< [out] PhenomPv3Storage struct */
    sysprecquant<T> *pAngles,           /**< [out] precession angle pre-computations struct */
    T m1,             /**< mass of primary in solar masses */
    T m2,             /**< mass of secondary in solar masses */
    T S1x,               /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    T S1y,               /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    T S1z,               /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    T S2x,               /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    T S2y,               /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    T S2z,               /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const T distance,    /**< distance of source (sec) */
    const T inclination, /**< inclination of source (rad) */
    const T phiRef,      /**< reference orbital phase (rad) */
    const T deltaF,      /**< Sampling frequency (Hz) */
    const T f_min,       /**< Starting GW frequency (Hz) */
    const T f_max,       /**< End frequency; 0 defaults to ringdown cutoff freq */
    const T f_ref        /**< Reference GW frequency (Hz) */
)
{
    p->PRECESSING = 0;
    if (S1x == 0. && S1y == 0. && S2x == 0. && S2y == 0.)
    {
        p->PRECESSING = 1; // This means the system is not precessing
    }
    
    /* input parameters */
    p->m1_SI = m1*GWAT_MSUN_SI;
    p->m2_SI = m2*GWAT_MSUN_SI;
    p->chi1x = S1x;
    p->chi1y = S1y;
    p->chi1z = S1z;
    p->chi2x = S2x;
    p->chi2y = S2y;
    p->chi2z = S2z;
    p->distance_SI = distance*c;
    p->phiRef = phiRef;
    p->deltaF = deltaF;
    p->f_min = f_min;
    p->f_max = f_max;
    p->f_ref = f_ref;

    PhenomPrecessingSpinEnforcePrimary(
        &(p->m1_SI),
        &(p->m2_SI),
        &(p->chi1x),
        &(p->chi1y),
        &(p->chi1z),
        &(p->chi2x),
        &(p->chi2y),
        &(p->chi2z));

    p->m1_Msun = m1;
    p->m2_Msun = m2;
    p->Mtot_SI = p->m1_SI + p->m2_SI;
    p->Mtot_Msun = p->m1_Msun + p->m2_Msun;

    p->eta = p->m1_Msun * p->m2_Msun / (p->Mtot_Msun * p->Mtot_Msun);
    p->q = p->m1_Msun / p->m2_Msun; /* with m1>=m2 so q>=1 */

    /* check for rounding errors */
    if (p->eta > 0.25 || p->q < 1.0)
    {
        T comparison;
        T eps = 1e-6;

        comparison = 0.25;
        nudge(&(p->eta), comparison, eps);
        comparison = 1.0;
        nudge(&(p->q), comparison, eps);
    }

    p->Msec = p->Mtot_Msun * GWAT_MTSUN_SI; /* Total mass in seconds */

    p->amp0 = (p->Mtot_Msun) * GWAT_MRSUN_SI * (p->Mtot_Msun) * GWAT_MTSUN_SI / (p->distance_SI);

    /* Rotate to PhenomP frame */
    /* chi1_l == chi1z, chi2_l == chi2z for intermediate calculations*/
    T chi1_l, chi2_l;
    PhenomP_ParametersFromSourceFrame(
        &chi1_l, &chi2_l, &(p->chip), &(p->thetaJN), &(p->alpha0), &(p->phi_aligned), &(p->zeta_polariz),
        m1, m2, p->f_ref, p->phiRef, inclination,
        p->chi1x, p->chi1y, p->chi1z,
        p->chi2x, p->chi2y, p->chi2z, IMRPhenomPv3_V);

    p->inclination = p->thetaJN;

    /* compute spins in polar coordinates */
    CartesianToPolar(&(p->chi1_theta), &(p->chi1_phi), &(p->chi1_mag), p->chi1x, p->chi1y, p->chi1z);
    CartesianToPolar(&(p->chi2_theta), &(p->chi2_phi), &(p->chi2_mag), p->chi2x, p->chi2y, p->chi2z);

    if (p->PRECESSING != 1) // precessing case. compute angles
    {
        T one = 1.0;
        T zero = 0.;
        /* Initialize precession angles */
        /* evaluating the angles at the reference frequency */
        p->f_ref_Orb_Hz = 0.5 * p->f_ref; /* factor of 0.5 to go from GW to Orbital frequency */

        /* precompute everything needed to compute precession angles */
        /* ExpansionOrder specifies how many terms in the PN expansion of the precession angles to use.
        * In PhenomP3 we set this to 5, i.e. all but the highest order terms.
        * */
        int ExpansionOrder = 5;
        T coschi1 = cos(p->chi1_theta);
        T coschi2 = cos(p->chi2_theta);
        InitializePrecession(
            pAngles,
            p->m1_SI, p->m2_SI,
            one, zero,
            coschi1, p->chi1_phi, p->chi1_mag,
            coschi2, p->chi2_phi, p->chi2_mag,
            p->f_ref, ExpansionOrder);
    }
}

/**
 * Compute the magnitude of L divided by GMsquare_over_c to
 * 3PN order with spin terms as a function of the orbital frequency in Hz
*/
template <class T> T OrbitalAngMom3PNSpinning(
    const T f_orb_hz,    /**< input Orbital frequencies (Hz) */
    const T m1_SI,           /**< Primary mass in SI (kg) */
    const T m2_SI,           /**< Secondary mass in SI (kg) */
    const T mul,          /**< Cosine of Polar angle of orbital angular momentum */
    const T phl,          /**< Azimuthal angle of orbital angular momentum  */
    T mu1,          /**< Cosine of Polar angle of primary spin w.r.t. orbital angular momentum */
    T ph1,          /**< Azimuthal angle of primary spin  */
    T ch1,          /**< Dimensionless spin magnitude of primary spin */
    T mu2,          /**< Cosine of Polar angle of secondary spin w.r.t. orbital angular momentum */
    T ph2,          /**< Azimuthal angle of secondary spin  */
    T ch2,          /**< Dimensionless spin magnitude of secondary spin */
    const T f_0,          /**< Reference Gravitational Wave frequency (Hz) */
    const int ExpansionOrder   /**< Keep terms upto ExpansionOrder in precession angles phi_z and zeta (1,2,3,4,5 or -1 for all orders) */
)
{
    T L_norm_3PN = 0;
    sysprecquant<T> *system = (sysprecquant<T> *)malloc(sizeof(sysprecquant<T>));

    InitializePrecession(system, m1_SI, m2_SI, mul, phl, mu1, ph1, ch1, mu2, ph2, ch2, f_0, ExpansionOrder);

    T xi, xi_2, L_norm;
    const T twopiGM_over_cthree = GWAT_TWOPI * GWAT_G_SI * (m1_SI + m2_SI) / (c*c*c);

    xi = pow(f_orb_hz * twopiGM_over_cthree, system->onethird);
    xi_2 = xi*xi;
    L_norm = system->nu / xi;
    L_norm_3PN = L_norm_3PN_of_xi(xi, xi_2, L_norm, system);

    free(system);
    return L_norm_3PN;
}


// Template instantiations


template void IMRPhenomPv3InitY2m<double>(sph_harm<double> *Ylm, double theta, double phi);
template void IMRPhenomPv3InitY2m<adouble>(sph_harm<adouble> *Ylm, adouble theta, adouble phi);

template void PhenomPrecessingSpinEnforcePrimary<double>(double *m1, double *m2,
    double *chi1x, double *chi1y, double *chi1z,
    double *chi2x, double *chi2y, double *chi2z);
template void PhenomPrecessingSpinEnforcePrimary<adouble>(adouble *m1, adouble *m2,
    adouble *chi1x, adouble *chi1y, adouble *chi1z,
    adouble *chi2x, adouble *chi2y, adouble *chi2z);

template void IMRPhenomPv3ComputeWignerD<double>(double (*WignerD22)[2][5], const double beta);
template void IMRPhenomPv3ComputeWignerD<adouble>(adouble (*WignerD22)[2][5], const adouble beta);

template void nudge<double>(double *x, double X, double epsilon);

template vector3D<double> CreateSphVector<double>(const double r, const double th, const double ph);
template vector3D<adouble> CreateSphVector<adouble>(const adouble r, const adouble th, const adouble ph);

template vector3D<double> ScaleVector<double>(double c, vector3D<double> vec);
template vector3D<adouble> ScaleVector<adouble>(adouble c, vector3D<adouble> vec);

template vector3D<double> VectorSum<double>(vector3D<double> vec1, vector3D<double> vec2);
template vector3D<adouble> VectorSum<adouble>(vector3D<adouble> vec1, vector3D<adouble> vec2);

template double DotProduct<double>(const vector3D<double> vec1, const vector3D<double> vec2);
template adouble DotProduct<adouble>(const vector3D<adouble> vec1, const vector3D<adouble> vec2);

template double VectorNorm<double>(const vector3D<double> vec);
template adouble VectorNorm<adouble>(const vector3D<adouble> vec);

template vector3D<double> CrossProduct<double>(const vector3D<double> vec1, const vector3D<double> vec2);
template vector3D<adouble> CrossProduct<adouble>(const vector3D<adouble> vec1, const vector3D<adouble> vec2);

template vector3D<double> Roots<double>(const double L_norm, const double J_norm, const sysprecquant<double> *system);
template vector3D<adouble> Roots<adouble>(const adouble L_norm, const adouble J_norm, const sysprecquant<adouble> *system);

template vector3D<double> BCDcoeff<double>(const double L_norm, const double J_norm, const sysprecquant<double> *system);
template vector3D<adouble> BCDcoeff<adouble>(const adouble L_norm, const adouble J_norm, const sysprecquant<adouble> *system);

template double beta<double>(const double a, const double b, const sysprecquant<double> *system);
template adouble beta<adouble>(const adouble a, const adouble b, const sysprecquant<adouble> *system);

template double sigma<double>(const double a, const double b, const sysprecquant<double> *system);
template adouble sigma<adouble>(const adouble a, const adouble b, const sysprecquant<adouble> *system);

template double tau<double>(const double a, const double b, const sysprecquant<double> *system);
template adouble tau<adouble>(const adouble a, const adouble b, const sysprecquant<adouble> *system);

template double L_norm_3PN_of_xi<double>(const double xi, const double xi_2, const double L_norm, const sysprecquant<double> *system);
template adouble L_norm_3PN_of_xi<adouble>(const adouble xi, const adouble xi_2, const adouble L_norm, const sysprecquant<adouble> *system);

template double J_norm_of_xi<double>(const double L_norm, const sysprecquant<double> *system);
template adouble J_norm_of_xi<adouble>(const adouble L_norm, const sysprecquant<adouble> *system);

template double S_norm_of_xi<double>(const double xi, const double xi_2, const vector3D<double> roots, const sysprecquant<double> *system);
template adouble S_norm_of_xi<adouble>(const adouble xi, const adouble xi_2, const vector3D<adouble> roots, const sysprecquant<adouble> *system);

template double costhetaL<double>(const double J_norm, const double L_norm, const double S_norm);
template adouble costhetaL<adouble>(const adouble J_norm, const adouble L_norm, const adouble S_norm);

template double u_of_xi<double>(const double xi, const double xi_2, const sysprecquant<double> *system);
template adouble u_of_xi<adouble>(const adouble xi, const adouble xi_2, const sysprecquant<adouble> *system);

template double phiz_of_xi<double>(const double xi, const double xi_2, const double J_norm, const sysprecquant<double> *system);
template adouble phiz_of_xi<adouble>(const adouble xi, const adouble xi_2, const adouble J_norm, const sysprecquant<adouble> *system);

template double zeta_of_xi<double>(const double xi, const double xi_2, const sysprecquant<double> *system);
template adouble zeta_of_xi<adouble>(const adouble xi, const adouble xi_2, const sysprecquant<adouble> *system);

template vector3D<double> computeMScorrections<double>(const double xi, const double xi_2, const double L_norm, const double J_norm, const vector3D<double> roots, const sysprecquant<double> *system);
template vector3D<adouble> computeMScorrections<adouble>(const adouble xi, const adouble xi_2, const adouble L_norm, const adouble J_norm, const vector3D<adouble> roots, const sysprecquant<adouble> *system);

template vector3D<double> c_coeffs<double>(const double xi, const double xi_2, const double J_norm, const vector3D<double> roots, const sysprecquant<double> *system);
template vector3D<adouble> c_coeffs<adouble>(const adouble xi, const adouble xi_2, const adouble J_norm, const vector3D<adouble> roots, const sysprecquant<adouble> *system);

template vector3D<adouble> d_coeffs<adouble>(const adouble L_norm, const adouble J_norm, const vector3D<adouble> roots);
template vector3D<double> d_coeffs<double>(const double L_norm, const double J_norm, const vector3D<double> roots);

template int checkOmegaz5<double>(const double Omegaz5);
template int checkOmegaz5<adouble>(const adouble Omegaz5);

template vector3D<double> compute_phiz_zeta_costhetaL3PN<double>(const double xi, const sysprecquant<double> *system);
template vector3D<adouble> compute_phiz_zeta_costhetaL3PN<adouble>(const adouble xi, const sysprecquant<adouble> *system);

template double L2PNR<double>(const double v, const double eta);
template adouble L2PNR<adouble>(const adouble v, const adouble eta);

template double L3PN<double>(
    const double f_orb_hz,
    const double m1, const double m2,
    const double s1x, const double s1y, const double s1z,
    const double s2x, const double s2y, const double s2z,
    const double f_0, const int ExpansionOrder);
template adouble L3PN<adouble>(
    const adouble f_orb_hz,
    const adouble m1, const adouble m2,
    const adouble s1x, const adouble s1y, const adouble s1z,
    const adouble s2x, const adouble s2y, const adouble s2z,
    const adouble f_0, const int ExpansionOrder);

template void CartesianToPolar<double>(double *polar, double *azimuthal, double *magnitude, double x, double y, double z);
template void CartesianToPolar<adouble>(adouble *polar, adouble *azimuthal, adouble *magnitude, adouble x, adouble y, adouble z);

template void PhenomP_ParametersFromSourceFrame<double>(
    double *chi1_l, double *chi2_l, double *chip, double *thetaJN, double *alpha0, double *phi_aligned, double *zeta_polariz,
    const double m1, const double m2, const double f_ref, const double phiRef, const double incl,
    const double s1x, const double s1y, const double s1z, const double s2x, const double s2y, const double s2z,
    IMRPhenomP_version_type IMRPhenomP_version);
template void PhenomP_ParametersFromSourceFrame<adouble>(
    adouble *chi1_l, adouble *chi2_l, adouble *chip, adouble *thetaJN, adouble *alpha0, adouble *phi_aligned, adouble *zeta_polariz,
    const adouble m1, const adouble m2, const adouble f_ref, const adouble phiRef, const adouble incl,
    const adouble s1x, const adouble s1y, const adouble s1z, const adouble s2x, const adouble s2y, const adouble s2z,
    IMRPhenomP_version_type IMRPhenomP_version);

template void PhenomPv3_Param_Transform<double>(source_parameters<double> *out, gen_params_base<double> *in);
template void PhenomPv3_Param_Transform<adouble>(source_parameters<adouble> *out, gen_params_base<adouble> *in);

template void InitializePrecession<double>(
    sysprecquant<double>* system,
    const double m1, const double m2, const double mul, const double phl, const double mu1, const double ph1,
    const double ch1, const double mu2, const double ph2, const double ch2, const double f_0,
    const int ExpansionOrder);
template void InitializePrecession<adouble>(
    sysprecquant<adouble>* system,
    const adouble m1, const adouble m2, const adouble mul, const adouble phl, const adouble mu1, const adouble ph1,
    const adouble ch1, const adouble mu2, const adouble ph2, const adouble ch2, const adouble f_0,
    const int ExpansionOrder);

template void init_PhenomPv3_Storage<double>(PhenomPv3Storage<double> *p,  sysprecquant<double> *pAngles,
    double m1, double m2,
    double S1x, double S1y, double S1z,
    double S2x, double S2y, double S2z,
    const double distance, const double inclination,
    const double phiRef, const double deltaF, const double f_min, const double f_max, const double f_ref
);
template void init_PhenomPv3_Storage<adouble>(PhenomPv3Storage<adouble> *p,  sysprecquant<adouble> *pAngles,
    adouble m1, adouble m2,
    adouble S1x, adouble S1y, adouble S1z,
    adouble S2x, adouble S2y, adouble S2z,
    const adouble distance, const adouble inclination,
    const adouble phiRef, const adouble deltaF, const adouble f_min, const adouble f_max, const adouble f_ref
);

template double OrbitalAngMom3PNSpinning<double>(
    const double f_orb_hz,
    const double m1_SI, const double m2_SI,
    const double mul, const double phl, double mu1, double ph1, double ch1, double mu2, double ph2, double ch2,
    const double f_0, const int ExpansionOrder
);
template adouble OrbitalAngMom3PNSpinning<adouble>(
    const adouble f_orb_hz,
    const adouble m1_SI, const adouble m2_SI,
    const adouble mul, const adouble phl, adouble mu1, adouble ph1, adouble ch1, adouble mu2, adouble ph2, adouble ch2,
    const adouble f_0, const int ExpansionOrder
);
