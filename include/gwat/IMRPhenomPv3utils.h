#ifndef IMRPHENOMPV3_UTILS_H
#define IMRPHENOMPV3_UTILS_H

#include "util.h"
#include <stdexcept>
#include <gsl/gsl_sf.h>
#include <vector>
#include "IMRPhenomP.h"


#define ONE_OVER_EIGHT 0.125
#define THREE_OVER_EIGHT 0.375
#define SQRT_6 2.44948974278317788

/*! \struct
 * 3D vector
 */
template <class T> struct vector3D
{
    T x;
    T y;
    T z;
};

/*! \struct
 * Holds quantities needed for precession calculations.
 */
template <class T> struct sysprecquant
{
    T onethird;
    T constants_u[6];
    T constants_phiz[6];
    T constants_zeta[6];
    T constants_L[6];
    T phiz_0, zeta_0, constant_of_S;
    T c_1, Ssqave, sqrtSsqave, Seff, c1_2, nu_2, nu_4, c_1_over_nu, S0_norm;
    T S1_norm_2, S2_norm_2;
    T dot1, dot2, dot12, dot1n, dot2n;
    T deltam_over_M, nu, q;
};

/**
 * Structure storing initial and derived variables for IMRPhenomPv3
 */
template <class T> struct PhenomPv3Storage
{
    int PRECESSING;   /**< integer to signify if system is precessing, 1 for false (not precessing), 0 for true (precessing) */
    /* input parameters */
    T m1_SI;       /**< mass of primary in SI (kg) */
    T m2_SI;       /**< mass of secondary in SI (kg) */
    T chi1x;       /**< x-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    T chi1y;       /**< y-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    T chi1z;       /**< z-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    T chi2x;       /**< x-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    T chi2y;       /**< y-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    T chi2z;       /**< z-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    T distance_SI; /**< distance to source in SI (m) */
    T inclination; /**< inclination - used to compute the angle thetaJN (rad) */
    T phiRef;      /**< */
    T deltaF;      /**< frequency spacing (Hz) */
    T f_min;       /**< starting GW frequency (Hz) */
    T f_max;       /**< ending GW frequency (Hz) */
    T f_ref;       /**< reference GW frequency (Hz) */
    /* derived parameters */
    T m1_Msun;      /**< mass of primary in solar masses */
    T m2_Msun;      /**< mass of secondary in solar masses */
    T Mtot_SI;      /**< total mass in SI (kg) */
    T Mtot_Msun;    /**< total mass in solar masses */
    T eta;          /**< Symmetric mass ratio*/
    T q;            /* with m1>=m2 so q>=1 */
    T Msec;         /**< Total mass in seconds */
    T f_ref_Orb_Hz; /**< Reference orbital frequency (Hz) [It's the reference GW frequency converted to orbital frequency] */
    T twopi_Msec;   /**< LAL_TWOPI * Msec */
    T amp0;         /**< frequency domain physical scaling */
    /* variables used when rotating input parameters (LAL frame) into PhenomP intrinsic parameters  */
    T chip;         /**< effective precessing parameter */
    T thetaJN;      /**< Angle between J0 and line of sight (z-direction) */
    T alpha0;       /**< Initial value of alpha angle (azimuthal precession angle) */
    T phi_aligned;  /**< Initial phase to feed the underlying aligned-spin model */
    T zeta_polariz; /**< Angle to rotate the polarizations */
    /* compute spins in polar coordinates */
    T chi1_mag;   /**< dimensionless spin magnitude on primary */
    T chi1_theta; /**< polar angle w.r.t. Lhat = (0,0,1) on primary */
    T chi1_phi;   /**< azimuthal angle w.r.t. Lhat = (0,0,1) on primary */
    T chi2_mag;   /**< dimensionless spin magnitude on secondary */
    T chi2_theta; /**< polar angle w.r.t. Lhat = (0,0,1) on secondary */
    T chi2_phi;   /**< azimuthal angle w.r.t. Lhat = (0,0,1) on secondary */
    /* Precession angles at reference frequency */
    T alphaRef;   /**< azimuthal precession angle at f_ref */
    T epsilonRef; /**< epsilon precession angle at f_ref */
    T betaRef;    /**< beta (opening angle) precession angle at f_ref */
};

template <class T>
void IMRPhenomPv3InitY2m(sph_harm<T> *Ylm, T theta, T phi);

template <class T> void PhenomPrecessingSpinEnforcePrimary(T *m1, T *m2,
    T *chi1x, T *chi1y, T *chi1z,
    T *chi2x, T *chi2y, T *chi2z);

template <class T> void IMRPhenomPv3ComputeWignerD(T (*WignerD22)[2][5], const T beta);

template <class T> void nudge(T *x, T X, T epsilon);

template <class T> vector3D<T> CreateSphVector(const T r, const T th, const T ph);
template <class T> vector3D<T> ScaleVector(T c, vector3D<T> vec);
template <class T> vector3D<T> VectorSum(vector3D<T> vec1, vector3D<T> vec2);
template <class T> T DotProduct(const vector3D<T> vec1, const vector3D<T> vec2);
template <class T> T VectorNorm(const vector3D<T> vec);
template <class T> vector3D<T> CrossProduct(const vector3D<T> vec1, const vector3D<T> vec2);

template <class T> vector3D<T> Roots(const T L_norm, const T J_norm, const sysprecquant<T> *system);
template <class T> vector3D<T> BCDcoeff(const T L_norm, const T J_norm, const sysprecquant<T> *system);
template <class T> T beta(const T a, const T b, const sysprecquant<T> *system);
template <class T> T sigma(const T a, const T b, const sysprecquant<T> *system);
template <class T> T tau(const T a, const T b, const sysprecquant<T> *system);

template <class T> T L_norm_3PN_of_xi(const T xi, const T xi_2, const T L_norm, const sysprecquant<T> *system);
template <class T> T J_norm_of_xi(const T L_norm, const sysprecquant<T> *system);
template <class T> T S_norm_of_xi(const T xi, const T xi_2, const vector3D<T> roots, const sysprecquant<T> *system);
template <class T> T costhetaL(const T J_norm, const T L_norm, const T S_norm);

template <class T> T u_of_xi(const T xi, const T xi_2, const sysprecquant<T> *system);
template <class T> T phiz_of_xi(const T xi, const T xi_2, const T J_norm, const sysprecquant<T> *system);
template <class T> T zeta_of_xi(const T xi, const T xi_2, const sysprecquant<T> *system);

template <class T> vector3D<T> computeMScorrections(const T xi, const T xi_2, const T L_norm, const T J_norm, const vector3D<T> roots, const sysprecquant<T> *system);

template <class T> int checkOmegaz5(const T Omegaz5);
void invalidExpansionOrder(const int ExpansionOrder);

template <class T> vector3D<T> compute_phiz_zeta_costhetaL3PN(const T xi, const sysprecquant<T> *system);

template <class T> T L2PNR(const T v, const T eta);
template <class T> T L3PN(
    const T f_orb_hz,
    const T m1, const T m2,
    const T s1x, const T s1y, const T s1z,
    const T s2x, const T s2y, const T s2z,
    const T f_0, const int ExpansionOrder);

template <class T> void CartesianToPolar(T *polar, T *azimuthal, T *magnitude, T x, T y, T z);

template <class T> void PhenomP_ParametersFromSourceFrame(
    T *chi1_l, T *chi2_l, T *chip, T *thetaJN, T *alpha0, T *phi_aligned, T *zeta_polariz,
    const T m1, const T m2, const T f_ref, const T phiRef, const T incl,
    const T s1x, const T s1y, const T s1z, const T s2x, const T s2y, const T s2z,
    IMRPhenomP_version_type IMRPhenomP_version);

template <class T> void PhenomPv3_Param_Transform(source_parameters<T> *out, gen_params_base<T> *in);

template <class T> void InitializePrecession(
    sysprecquant<T>* system,
    const T m1, const T m2, const T mul, const T phl, const T mu1, const T ph1,
    const T ch1, const T mu2, const T ph2, const T ch2, const T f_0,
    const int ExpansionOrder);

template <class T> void init_PhenomPv3_Storage(PhenomPv3Storage<T> *p,  sysprecquant<T> *pAngles,
    T m1, T m2,
    T S1x, T S1y, T S1z,
    T S2x, T S2y, T S2z,
    const T distance, const T inclination,
    const T phiRef, const T deltaF, const T f_min, const T f_max, const T f_ref
);

template <class T> T OrbitalAngMom3PNSpinning(
    const T f_orb_hz,
    const T m1_SI, const T m2_SI,
    const T mul, const T phl, T mu1, T ph1, T ch1, T mu2, T ph2, T ch2,
    const T f_0, const int ExpansionOrder
);


#endif /* IMRPHENOMPV3_UTILS_H */