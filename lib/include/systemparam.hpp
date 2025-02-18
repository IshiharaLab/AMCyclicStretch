#ifndef SYSTEMPARAM_HPP
#define SYSTEMPARAM_HPP

#include <cmath>
#include <random>

template <typename T>
constexpr T Sqrt(T s)
{
    T x = s / 2.0;
    T prev = 0.0;

    while (x != prev)
    {
        prev = x;
        x = (x + s / x) / 2.0;
    }
    return x;
}

// コンパイル時にパラメータを変更するつもり
#ifndef N_
#define N_ 10000
#endif
#ifndef L_
#define L_ 10.
#endif
#ifndef dt_
#define dt_ 0.001
#endif
#ifndef T_MAX_
#define T_MAX_ 1000.
#endif
#ifndef E_LJ_
#define E_LJ_ 0.00000717
#endif
#ifndef Ri_
#define Ri_ 1.0
#endif
#ifndef V0_
#define V0_ 0.1
#endif
#ifndef E_v_
#define E_v_ 0.953
#endif
#ifndef epsilon_
#define epsilon_ 0.1
#endif
#ifndef nu_
#define nu_ 0.5
#endif
#ifndef sigma_
#define sigma_ 0.06
#endif
#ifndef E_s_
#define E_s_ 1.7
#endif
#ifndef VERBOSE_
#define VERBOSE_ 1000
#endif
#ifndef MARGIN_
#define MARGIN_ 0.5
#endif
#ifndef DATA_PATH
#define DATA_PATH "./"
#endif

// math constant
constexpr double PI2 = M_PI * 2.;

// system
constexpr int N = N_;
constexpr double L = L_;
constexpr double dt = dt_;
constexpr double T_MAX = T_MAX_;
constexpr int TC_MAX = (int)(T_MAX / dt) + 1;
constexpr double dt2 = dt * dt;
constexpr double Lh = L * 0.5;
constexpr int VERBOSE = VERBOSE_;
constexpr double Delta_t = 10. / dt;
constexpr int VERBOSE_10sec = VERBOSE - (int)(Delta_t);

// excluded volume
constexpr double PS = 1.0;                           // particle size
constexpr double E_LJ = E_LJ_;                       // energy coefficient
constexpr double PS2 = PS * PS;                      //
constexpr double CUTOFF = 0.3;                       // cutoff to avoid too large exclusive force
constexpr double CL2 = CUTOFF * CUTOFF;              //
constexpr double RC2 = PS2 / CL2;                    //
constexpr double RC6 = RC2 * RC2 * RC2;              //
constexpr double C0 = E_LJ * RC6 * (1. - RC6) / CL2; //

// vicsek
constexpr double Ri = Ri_;      //
constexpr double V0 = V0_;      //
constexpr double E_v = E_v_;    // energy coefficient
constexpr double Ri2 = Ri * Ri; //

// mesh / pair
constexpr double Rm = Ri > PS ? Ri : PS;
constexpr double MARGIN = MARGIN_;                   //
constexpr double MARGIN_h = MARGIN * 0.5;            //
constexpr double P2 = (Rm + MARGIN) * (Rm + MARGIN); // pair list length ^2
constexpr int M = (int)(L / (Rm + MARGIN)) - 1;      // 一辺のmesh
constexpr double MS = L / (double)M;                 // mesh size
constexpr int MN = M * M;                            // meshの総数
constexpr double IMS = 1.0 / MS;                     //
constexpr double cMS = 4.;
constexpr int cM = (int)(L / cMS);

// stretch
constexpr double E_s = E_s_;                                              //
constexpr double epsilon = epsilon_;                                      //
constexpr double nu = nu_;                                                //
constexpr double beta = (1. - nu) / (1. + nu);                            //
constexpr double alpha = epsilon * epsilon * E_s * (1. + nu) * (1. + nu); //

// random
static std::uniform_real_distribution ud(0.0, 1.0);

// periodic
inline double adjust_periodic_x(double x)
{
    if (x < 0)
        x += L;
    if (x >= L)
        x -= L;
    return x;
}

inline double adjust_periodic_theta(double theta)
{
    if (theta < -M_PI)
        theta += PI2;
    if (theta > M_PI)
        theta -= PI2;
    return theta;
}

inline double adjust_periodic_d(double d)
{
    if (d < -Lh)
        d += L;
    if (d > Lh)
        d -= L;
    return d;
}

#endif // SYSTEMPARAM_HPP
