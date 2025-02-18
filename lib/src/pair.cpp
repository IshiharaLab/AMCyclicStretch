#include "../include/systemparam.hpp"
#include "../include/pair.hpp"
#include <algorithm>
#include <vector>
#include <array>
#include <iostream>

void pair_list(std::vector<std::array<int, 2>> &pairs,
               const int (&heads)[MN], const int (&count)[MN],
               const double (&x)[N], const double (&y)[N])
{
    int i_x, i_y;
    pairs.clear();
    for (int i = 0; i < MN; i++)
    {
        i_y = i / M;
        i_x = i - i_y * M;
        pair_search_same_mesh(pairs, i, heads, count, x, y);
        pair_search_next_mesh(pairs, i, i_x - 1, i_y - 1, heads, count, x, y);
        pair_search_next_mesh(pairs, i, i_x - 1, i_y + 1, heads, count, x, y);
        pair_search_next_mesh(pairs, i, i_x + 1, i_y - 1, heads, count, x, y);
        pair_search_next_mesh(pairs, i, i_x + 1, i_y + 1, heads, count, x, y);
        pair_search_next_mesh(pairs, i, i_x, i_y - 1, heads, count, x, y);
        pair_search_next_mesh(pairs, i, i_x, i_y + 1, heads, count, x, y);
        pair_search_next_mesh(pairs, i, i_x - 1, i_y, heads, count, x, y);
        pair_search_next_mesh(pairs, i, i_x + 1, i_y, heads, count, x, y);
    }
};

void mesh_new_indices(int (&new_indices)[N], int (&heads)[MN],
                      int (&count)[MN],
                      const double (&x)[N], const double (&y)[N])
{
    int ix, iy, index;
    int pos[N] = {};
    int registered_particle[MN] = {};
    int sum = 0;

    for (int i = 0; i < MN; i++)
        count[i] = 0;

    for (int i = 0; i < N; i++)
    {
        ix = (int)(x[i] * IMS);
        iy = (int)(y[i] * IMS);
        index = calculate_mesh_id(ix, iy);
        count[index] += 1;
        pos[i] = index;
    }

    sum = 0;
    for (int i = 0; i < MN - 1; i++)
    {
        sum += count[i];
        heads[i + 1] = sum;
    }

    for (int i = 0; i < N; i++)
    {
        index = pos[i];
        new_indices[heads[index] + registered_particle[index]] = i;
        registered_particle[index] += 1;
    }
};

void pair_search_same_mesh(std::vector<std::array<int, 2>> &pairs,
                           const int mesh_id,
                           const int (&heads)[MN], const int (&count)[MN],
                           const double (&x)[N], const double (&y)[N])
{
    int i_min, i_max, j_max;
    double dx, dy, d_sq;
    std::array<int, 2> pair;

    i_min = heads[mesh_id];
    i_max = heads[mesh_id] + count[mesh_id] - 1;
    // j_min = i + 1;
    j_max = heads[mesh_id] + count[mesh_id];

    for (int i = i_min; i < i_max; i++)
    {
        for (int j = i + 1; j < j_max; j++)
        {

            dx = x[j] - x[i];
            dy = y[j] - y[i];
            dx = adjust_periodic_d(dx);
            dy = adjust_periodic_d(dy);
            d_sq = dx * dx + dy * dy;

            if (d_sq < P2)
            {
                pair[0] = i;
                pair[1] = j;
                pairs.push_back(pair);
            }
        }
    }
};

void pair_search_next_mesh(std::vector<std::array<int, 2>> &pairs,
                           const int mesh_id, const int ix, const int iy,
                           const int (&heads)[MN], const int (&count)[MN],
                           const double (&x)[N], const double (&y)[N])
{
    int next_mesh_id;
    int i_min, i_max, j_min, j_max;
    double dx, dy, d_sq;
    std::array<int, 2> pair;

    next_mesh_id = calculate_mesh_id(ix, iy);

    i_min = heads[mesh_id];
    i_max = heads[mesh_id] + count[mesh_id];
    j_min = heads[next_mesh_id];
    j_max = heads[next_mesh_id] + count[next_mesh_id];

    for (int i = i_min; i < i_max; i++)
    {
        for (int j = j_min; j < j_max; j++)
        {
            if (i < j)
            {
                dx = x[j] - x[i];
                dy = y[j] - y[i];
                dx = adjust_periodic_d(dx);
                dy = adjust_periodic_d(dy);
                d_sq = dx * dx + dy * dy;
                if (d_sq < P2)
                {
                    pair[0] = i;
                    pair[1] = j;
                    pairs.push_back(pair);
                }
            }
        }
    }
};
