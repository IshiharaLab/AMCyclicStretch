#ifndef LJ_HPP
#define LJ_HPP

#include "systemparam.hpp"
#include <vector>
#include <array>

void calculate_LJ(double (&dU_LJ_dx)[N], double (&dU_LJ_dy)[N],
                  const double (&x)[N], const double (&y)[N]);

void calculate_LJ_pair(double (&dU_LJ_dx)[N], double (&dU_LJ_dy)[N],
                       const std::vector<std::array<int, 2>>(&pairs), const int pair_num,
                       const double (&x)[N], const double (&y)[N]);

#endif // LJ_HPP
