#include "systemparam.hpp"
#include <string>
#include <fstream>
#include <cmath>

const double density_norm = 1. / (cMS * cMS * Delta_t);

inline double calculate_order1(const double (&theta)[N])
{
    double c = 0.;
    double s = 0.;
    for (int i = 0; i < N; i++)
    {
        c += std::cos(theta[i]);
        s += std::sin(theta[i]);
    }
    return std::atan2(s, c);
}

inline double calculate_order2(const double (&theta)[N])
{
    double c = 0.;
    double s = 0.;
    for (int i = 0; i < N; i++)
    {
        c += std::cos(2. * theta[i]);
        s += std::sin(2. * theta[i]);
    }
    return std::atan2(s, c);
}

inline void coarse_grain(double (&history_density)[cM][cM], double (&history_theta)[cM][cM],
                         const double (&history_theta_sin)[cM][cM], const double (&history_theta_cos)[cM][cM])
{
    for (int i = 0; i < cM; i++)
    {
        for (int j = 0; j < cM; j++)
        {
            if (history_density[j][i] == 0)
            {
                history_theta[j][i] = std::nan("");
            }
            else
            {
                history_theta[j][i] = std::atan2(history_theta_sin[j][i], history_theta_cos[j][i]);
                history_density[j][i] *= density_norm;
            }
        }
    }
}

inline void coarse_grain2(double (&history_density)[cM][cM], double (&history_theta)[cM][cM], double (&history_Omega)[cM][cM],
                          const double (&history_theta_sin)[cM][cM], const double (&history_theta_cos)[cM][cM],
                          const double (&history_Omega_sin)[cM][cM], const double (&history_Omega_cos)[cM][cM])
{
    for (int i = 0; i < cM; i++)
    {
        for (int j = 0; j < cM; j++)
        {
            if (history_density[j][i] == 0)
            {
                history_theta[j][i] = std::nan("");
                history_Omega[j][i] = std::nan("");
            }
            else
            {
                history_theta[j][i] = std::atan2(history_theta_sin[j][i], history_theta_cos[j][i]);
                history_Omega[j][i] = std::atan2(history_Omega_sin[j][i], history_Omega_cos[j][i]);
                history_density[j][i] *= density_norm;
            }
        }
    }
}

inline void save_coarse_grained_data(std::ofstream &ofs, std::string &file_path, const double (&data)[cM][cM])
{
    ofs.open(file_path);
    for (int i = 0; i < cM; i++)
    {
        for (int j = 0; j < cM - 1; j++)
        {
            ofs << data[j][i] << ",";
        }
        ofs << data[cM - 1][i] << "\n";
    }
    ofs.close();
}

inline void save_data(std::ofstream &ofs, std::string &file_path,
                      const double (&x)[N], const double (&y)[N], const double (&theta)[N])
{
    ofs.open(file_path);
    for (int i = 0; i < N; i++)
    {
        ofs << x[i] << ","
            << y[i] << ","
            << theta[i] << "\n";
    }
    ofs.close();
}

inline void save_data2(std::ofstream &ofs, std::string &file_path,
                       const int (&indices)[N], const double (&x)[N], const double (&y)[N], const double (&theta)[N])
{
    ofs.open(file_path);
    for (int i = 0; i < N; i++)
    {
        ofs
            << x[i] << ","
            << y[i] << ","
            << theta[i] << ","
            << indices[i] << "\n";
    }
    ofs.close();
}

inline void save_data_omega(std::ofstream &ofs, std::string &file_path,
                            const int (&indices)[N], const double (&x)[N], const double (&y)[N], const double (&theta)[N], const double (&omega)[N])
{
    ofs.open(file_path);
    for (int i = 0; i < N; i++)
    {
        ofs
            << x[i] << ","
            << y[i] << ","
            << theta[i] << ","
            << omega[i] << ","
            << indices[i] << "\n";
    }
    ofs.close();
}
