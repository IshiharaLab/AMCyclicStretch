#include "../include/systemparam.hpp"
#include "../include/nematic.hpp"
#include <algorithm>
#include <vector>
#include <array>
#include <cassert>
#include <iostream>

void calculate_nematic(double (&dU_c_dtheta)[N],
                       const double (&x)[N], const double (&y)[N], const double (&theta)[N])
{
    double dx, dy, d_sq;
    double sum;
    int cnt;

    for (int i = 0; i < N; i++)
    {
        cnt = 0;
        sum = 0.0;
        for (int j = 0; j < N; j++)
        {
            dx = x[j] - x[i];
            dy = y[j] - y[i];
            dx = adjust_periodic_d(dx);
            dy = adjust_periodic_d(dy);
            d_sq = dx * dx + dy * dy;
            if (d_sq < Ri2)
            {
                cnt += 1;
                sum += std::sin(2. * (theta[i] - theta[j]));
            }
        }
        dU_c_dtheta[i] = E_v / (double)cnt * sum;
    }
};

void calculate_nematic_pair(double (&dU_c_dtheta)[N],
                            const std::vector<std::array<int, 2>>(&pairs), const int pair_num,
                            const double (&x)[N], const double (&y)[N], const double (&theta)[N])
{
    double dx, dy, d_sq;
    double sum[N] = {};
    int cnt[N] = {};
    int p_i, p_j;

    for (int i = 0; i < N; i++)
        cnt[i] = 1;

    for (int p = 0; p < pair_num; p++)
    {
        p_i = pairs[(size_t)p][0];
        p_j = pairs[(size_t)p][1];

        dx = x[p_i] - x[p_j];
        dy = y[p_i] - y[p_j];
        dx = adjust_periodic_d(dx);
        dy = adjust_periodic_d(dy);
        d_sq = dx * dx + dy * dy;

        if (d_sq < Ri2)
        {
            cnt[p_i] += 1;
            cnt[p_j] += 1;
            sum[p_i] += std::sin(2. * (theta[p_i] - theta[p_j]));
            sum[p_j] += std::sin(2. * (theta[p_j] - theta[p_i]));
        }
    }
    for (int i = 0; i < N; i++)
    {
        dU_c_dtheta[i] = E_v / (double)cnt[i] * sum[i];
    }
};

void calculate_nematic_pair2(double (&dU_c_dtheta)[N],
                             const std::vector<std::array<int, 2>>(&pairs), const int pair_num,
                             const double (&x)[N], const double (&y)[N], const double (&theta)[N])
{
    double dx, dy, d_sq;
    double sum[N] = {};
    int p_i, p_j;

    for (int p = 0; p < pair_num; p++)
    {
        p_i = pairs[(size_t)p][0];
        p_j = pairs[(size_t)p][1];

        dx = x[p_i] - x[p_j];
        dy = y[p_i] - y[p_j];
        dx = adjust_periodic_d(dx);
        dy = adjust_periodic_d(dy);
        d_sq = dx * dx + dy * dy;

        if (d_sq < Ri2)
        {
            sum[p_i] += std::sin(2. * (theta[p_i] - theta[p_j]));
            sum[p_j] += std::sin(2. * (theta[p_j] - theta[p_i]));
        }
    }
    for (int i = 0; i < N; i++)
    {
        dU_c_dtheta[i] = E_v * sum[i];
    }
};
