#ifndef PAIR_HPP_
#define PAIR_HPP_

#include "systemparam.hpp"
#include <vector>
#include <array>

inline int calculate_mesh_id(int i, int j)
{
    if (i < 0)
        i += M;
    if (i >= M)
        i -= M;
    if (j < 0)
        j += M;
    if (j >= M)
        j -= M;
    return i + j * M;
}

template <typename T>
inline void mesh_sort(T (&X)[N], const int (&new_indices)[N])
{
    T temp[N] = {};
    for (int i = 0; i < N; i++)
        temp[i] = X[new_indices[i]];
    for (int i = 0; i < N; i++)
        X[i] = temp[i];
};

void pair_list(std::vector<std::array<int, 2>> &pairs,
               const int (&heads)[MN], const int (&count)[MN],
               const double (&x)[N], const double (&y)[N]);

void mesh_new_indices(int (&new_indices)[N], int (&heads)[MN],
                      int (&count)[MN],
                      const double (&x)[N], const double (&y)[N]);

void pair_search_same_mesh(std::vector<std::array<int, 2>> &pairs,
                           const int mesh_id,
                           const int (&heads)[MN], const int (&count)[MN],
                           const double (&x)[N], const double (&y)[N]);

void pair_search_next_mesh(std::vector<std::array<int, 2>> &pairs,
                           const int mesh_id, const int ix, const int iy,
                           const int (&heads)[MN], const int (&count)[MN],
                           const double (&x)[N], const double (&y)[N]);

#endif // PAIR_HPP_
