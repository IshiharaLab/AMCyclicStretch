#ifndef VICSEK_HPP
#define VICSEK_HPP

#include "systemparam.hpp"
#include <vector>
#include <array>

void calculate_vicsek(double (&dU_c_dtheta)[N],
                      const double (&x)[N], const double (&y)[N], const double (&theta)[N]);

void calculate_vicsek_pair(double (&dU_c_dtheta)[N],
                           const std::vector<std::array<int, 2>>(&pairs), const int pair_num,
                           const double (&x)[N], const double (&y)[N], const double (&theta)[N]);

#endif // VICSEK_HPP
