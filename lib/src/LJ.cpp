#include "../include/systemparam.hpp"
#include "../include/LJ.hpp"
#include <algorithm>
#include <vector>
#include <array>
#include <cassert>
#include <iostream>

void calculate_LJ(double (&dU_LJ_dx)[N], double (&dU_LJ_dy)[N],
                  const double (&x)[N], const double (&y)[N])
{
    double dx, dy, d_sq;
    double rc2, rc6, c0;
    double sum_x, sum_y;

    for (int i = 0; i < N; i++)
    {
        dU_LJ_dx[i] = 0.;
        dU_LJ_dy[i] = 0.;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            dx = x[i] - x[j];
            dy = y[i] - y[j];
            dx = adjust_periodic_d(dx);
            dy = adjust_periodic_d(dy);
            d_sq = dx * dx + dy * dy;

            if (d_sq < CL2)
            {
                dU_LJ_dx[i] += C0 * dx;
                dU_LJ_dy[i] += C0 * dy;
            }
            else if (d_sq < PS2)
            {
                rc2 = PS2 / d_sq;
                rc6 = rc2 * rc2 * rc2;
                c0 = E_LJ * rc6 * (1. - rc6) / d_sq;
                dU_LJ_dx[i] += c0 * dx;
                dU_LJ_dy[i] += c0 * dy;
            }
        }
    }
};

void calculate_LJ_pair(double (&dU_LJ_dx)[N], double (&dU_LJ_dy)[N],
                       const std::vector<std::array<int, 2>>(&pairs), const int pair_num,
                       const double (&x)[N], const double (&y)[N])
{
    double dx, dy, d_sq;
    double rc2, rc6, c0;
    int p_i, p_j;

    for (int i = 0; i < N; i++)
    {
        dU_LJ_dx[i] = 0.;
        dU_LJ_dy[i] = 0.;
    }

    for (int p = 0; p < pair_num; p++)
    {
        p_i = pairs[(size_t)p][0];
        p_j = pairs[(size_t)p][1];

        dx = x[p_i] - x[p_j];
        dy = y[p_i] - y[p_j];
        dx = adjust_periodic_d(dx);
        dy = adjust_periodic_d(dy);
        d_sq = dx * dx + dy * dy;

        if (d_sq < CL2)
        {
            dU_LJ_dx[p_i] += C0 * dx;
            dU_LJ_dx[p_j] -= C0 * dx;
            dU_LJ_dy[p_i] += C0 * dy;
            dU_LJ_dy[p_j] -= C0 * dy;
        }
        else if (d_sq < PS2)
        {
            rc2 = PS2 / d_sq;
            rc6 = rc2 * rc2 * rc2;
            c0 = E_LJ * rc6 * (1. - rc6) / d_sq;
            dU_LJ_dx[p_i] += c0 * dx;
            dU_LJ_dx[p_j] -= c0 * dx;
            dU_LJ_dy[p_i] += c0 * dy;
            dU_LJ_dy[p_j] -= c0 * dy;
        }
    }
};
