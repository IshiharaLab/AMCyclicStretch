#include "../lib/include/LJ.hpp"
#include "../lib/include/observe.hpp"
#include "../lib/include/pair.hpp"
#include "../lib/include/systemparam.hpp"
#include "../lib/include/nematic.hpp"
// #include "../lib/include/vicsek.hpp"
#include <iostream>
#include <fstream>
// #include <sstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>

static const std::string file_path = DATA_PATH;
static std::ofstream ofs, ofs_history;

int main()
{
    std::random_device rd{};
    long int seed = (long int)(rd());
    std::cout << DATA_PATH << "\t" << seed << "\n";
    static std::mt19937 mt(seed);
    static std::normal_distribution<double> noise(0.0, Sqrt(dt) * sigma_);

    double x[N], y[N], dx[N], dy[N];
    double theta[N];
    double dU_c_dtheta[N];
    double dU_LJ_dx[N], dU_LJ_dy[N];
    int new_indices[N] = {};
    int original_indices[N] = {};
    int count[MN] = {};
    int heads[MN] = {};
    std::vector<std::array<int, 2>> pairs;
    int pair_num = 0;
    double pair_life = -1.0;
    double d_sq, d_sq_max;
    double order1, order2, sum_cos1, sum_sin1, sum_cos2, sum_sin2;

    const std::string csv_file_name = "";
    std::string file_name;

    // initialize
    for (int i = 0; i < N; i++)
    {
        x[i] = ud(mt) * L;
        y[i] = ud(mt) * L;
        theta[i] = ud(mt) * PI2 - M_PI;
        original_indices[i] = i;
    }

    for (int tc = 0; tc < TC_MAX; tc++)
    {
        if (tc % VERBOSE == 0)
        {
            file_name = file_path + "data_" + std::to_string(tc / VERBOSE) + ".csv";
            save_data2(ofs, file_name, original_indices, x, y, theta);

            sum_cos1 = 0.0;
            sum_sin1 = 0.0;
            sum_cos2 = 0.0;
            sum_sin2 = 0.0;
            for (int i = 0; i < N; i++)
            {
                sum_cos1 += std::cos(theta[i]);
                sum_sin1 += std::sin(theta[i]);
                sum_cos2 += std::cos(2 * theta[i]);
                sum_sin2 += std::sin(2 * theta[i]);
            }
            order1 = std::sqrt(sum_cos1 * sum_cos1 + sum_sin1 * sum_sin1) / (double)N;
            order2 = std::sqrt(sum_cos2 * sum_cos2 + sum_sin2 * sum_sin2) / (double)N;

            ofs_history.open(file_path + "history.csv", std::ios::app);
            ofs_history << tc * dt << "," << order1 << "," << order2 << "\n";
            ofs_history.close();
        }

        if (pair_life < 0)
        {
            mesh_new_indices(new_indices, heads, count, x, y);
            mesh_sort(x, new_indices);
            mesh_sort(y, new_indices);
            mesh_sort(theta, new_indices);
            mesh_sort(original_indices, new_indices);
            pair_list(pairs, heads, count, x, y);
            pair_num = (int)pairs.size();
            pair_life = MARGIN_h;
        }
        calculate_nematic_pair(dU_c_dtheta, pairs, pair_num, x, y, theta);
        calculate_LJ_pair(dU_LJ_dx, dU_LJ_dy, pairs, pair_num, x, y);

        // update
        for (int i = 0; i < N; i++)
        {
            dx[i] = V0 * std::cos(theta[i]) - dU_LJ_dx[i];
            dy[i] = V0 * std::sin(theta[i]) - dU_LJ_dy[i];
        }
        for (int i = 0; i < N; i++)
        {
            x[i] += dt * dx[i];
            y[i] += dt * dy[i];
            theta[i] += noise(mt) +
                        dt * (alpha * (beta + std::cos(2 * theta[i])) * std::sin(2. * theta[i]) -
                              dU_c_dtheta[i]);
        }
        for (int i = 0; i < N; i++)
        {
            x[i] = adjust_periodic_x(x[i]);
            y[i] = adjust_periodic_x(y[i]);
            theta[i] = adjust_periodic_theta(theta[i]);
        }

        d_sq_max = 0.0;
        for (int i = 0; i < N; i++)
        {
            d_sq = dt2 * (dx[i] * dx[i] + dy[i] * dy[i]);
            d_sq_max = std::max(d_sq, d_sq_max);
        }
        pair_life -= std::sqrt(d_sq_max);
    }
}
