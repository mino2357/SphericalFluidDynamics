//
// Flow in a sphere.
// Description: This is a simple code to calculate the velocity field in a sphere.
// 0 <= longitude <= 2pi, - limit_latitudes <= latitude <= limit_latitudes
//
// gnuplot: plot "log.dat" u 1:2:(0.1*$3):(0.1*$4) w vector
// gnuplot> set xr [0:2*pi]
// gnuplot> set yr [-0.25*pi:0.25*pi]

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

constexpr double limit_latitudes = 75.0 * M_PI / 180.0;
constexpr double R = 4.0;
constexpr double T = 1.0;

class FlowInSphere
{
public:
    double radius = 1.0;
    double nu = 0.001;
    double delta_t = 0.00001;
    double delta_longitudes;
    double delta_latitudes;
    size_t N_longitudes;
    size_t N_latitules;
    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> u_longitudes;
    std::vector<std::vector<double>> u_latitudes;
    FlowInSphere(size_t n, size_t m);
    void setVelocity();
    void printVelocity();
    // advection
    void step_advection_eq();
    // pressure
    double pressure(double);
    // max_velocity
    double max_velocity();
};

double FlowInSphere::max_velocity()
{
    double max = 0.0;
    for (size_t i = 0; i < N_longitudes; i++)
    {
        for (size_t j = 0; j < N_latitules; j++)
        {
            double v = sqrt(u_longitudes[i][j]*u_longitudes[i][j] + u_latitudes[i][j]*u_latitudes[i][j]);
            if (v > max)
            {
                max = v;
            }
        }
    }
    return max;
}

double FlowInSphere::pressure(double rho)
{
    return rho * R * T;
}

FlowInSphere::FlowInSphere(size_t n, size_t m)
{
    delta_longitudes = 2.0 * M_PI / (n - 2);
    delta_latitudes = 2.0 * limit_latitudes / (m - 2);
    N_longitudes = n;
    N_latitules = m;
    u_longitudes.resize(n);
    u_latitudes.resize(n);
    rho.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        u_longitudes[i].resize(m);
        u_latitudes[i].resize(m);
        rho[i].resize(m);
    }
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            u_longitudes[i][j] = 0.0;
            u_latitudes[i][j] = 0.0;
            rho[i][j] = 1.0;
        }
    }
}

void FlowInSphere::setVelocity()
{
    for (size_t i = 0; i < N_longitudes; i++)
    {
        for (size_t j = 0; j < N_latitules; j++)
        {
            double alpha = 2.0 * M_PI * i / (N_longitudes - 1);
            double beta = - limit_latitudes + (2.0 * limit_latitudes) * j / (N_latitules - 1);
            u_longitudes[i][j] = 0.0;
            u_latitudes[i][j] = 0.0;
            if (1.0 <= alpha && alpha <= 1.5 && - limit_latitudes / 4.0 <= beta && beta <= limit_latitudes / 4.0)
            {
                u_longitudes[i][j] = 0.0e-1;
                u_latitudes[i][j] = 0.0;
            }
            if (0.0 <= beta)
            {
                u_longitudes[i][j] = 0.0e-1;
                u_latitudes[i][j] = 0.0;
            }
            if (beta <= -0.0)
            {
                u_longitudes[i][j] = 0.0e-1;
                u_latitudes[i][j] = 0.0;
            }
        }
    }
}

void FlowInSphere::printVelocity()
{
    for (size_t i = 0; i < N_longitudes; i++)
    {
        for (size_t j = 0; j < N_latitules; j++)
        {
            double alpha = 2.0 * M_PI * i / (N_longitudes - 2);
            double beta = - limit_latitudes + (2.0 * limit_latitudes) * j / (N_latitules - 2);
            std::cout << alpha * 180.0 / M_PI << " " << beta * 180.0 / M_PI << " " << u_longitudes[i][j] << " " << u_latitudes[i][j]
                        << " " << rho[i][j] << " " << sqrt(u_longitudes[i][j]*u_longitudes[i][j] + u_latitudes[i][j]*u_latitudes[i][j]) << std::endl;
        }
        std::cout << std::endl;
    }
}

void FlowInSphere::step_advection_eq()
{
    // advection equation
    auto u_longitudes_new = u_longitudes;
    auto u_latitudes_new = u_latitudes;
    auto rho_new = rho;
    for (size_t i = 1; i < N_longitudes - 1; i++)
    {
        for (size_t j = 1; j < N_latitules - 1; j++)
        {
            //double alpha = 2.0 * M_PI * i / (N_longitudes - 1);
            double beta = - limit_latitudes + (2.0 * limit_latitudes) * j / (N_latitules - 1);
            auto long_long_long = 0.0;
            auto long_lati_long = 0.0;
            auto lati_long_lati = 0.0;
            auto lati_lati_lati = 0.0;
            if (u_longitudes[i][j] > 0)
            {
                long_long_long = u_longitudes[i][j] * (u_longitudes[i][j] - u_longitudes[i-1][j]) / delta_longitudes;
                long_lati_long = u_longitudes[i][j] * (u_latitudes[i][j] - u_latitudes[i-1][j]) / delta_longitudes;
            }
            else
            {
                long_long_long = u_longitudes[i][j] * (u_longitudes[i+1][j] - u_longitudes[i][j]) / delta_longitudes;
                long_lati_long = u_longitudes[i][j] * (u_latitudes[i+1][j] - u_latitudes[i][j]) / delta_longitudes;
            }
            if(u_latitudes[i][j] > 0)
            {
                lati_long_lati = u_latitudes[i][j] * (u_longitudes[i][j] - u_longitudes[i][j-1]) / delta_latitudes;
                lati_lati_lati = u_latitudes[i][j] * (u_latitudes[i][j] - u_latitudes[i][j-1]) / delta_latitudes;
            }
            else
            {
                lati_long_lati = u_latitudes[i][j] * (u_longitudes[i][j+1] - u_longitudes[i][j]) / delta_latitudes;
                lati_lati_lati = u_latitudes[i][j] * (u_latitudes[i][j+1] - u_latitudes[i][j]) / delta_latitudes;
            }
            u_longitudes_new[i][j] = u_longitudes[i][j] - delta_t * (long_long_long + lati_long_lati)
                                    + delta_t * 2.0 * tan(beta) * u_latitudes[i][j] * u_latitudes[i][j]
                                    - delta_t * (pressure(rho[i+1][j]) - pressure(rho[i-1][j])) / (2.0 * delta_longitudes * rho[i][j] * radius * radius * cos(beta) * cos(beta))
                                    + delta_t * nu / (radius * radius) * ((1.0 / (cos(beta) * cos(beta))) * (u_longitudes[i-1][j] - 2.0 * u_longitudes[i][j] + u_longitudes[i+1][j]) / (delta_longitudes * delta_longitudes)
                                    - tan(beta) * (u_longitudes[i][j+1] - u_longitudes[i][j-1]) / (2.0 * delta_latitudes)
                                    + (u_longitudes[i][j-1] - 2.0 * u_longitudes[i][j] + u_longitudes[i][j+1]) / (delta_latitudes * delta_latitudes));
            u_latitudes_new[i][j] = u_latitudes[i][j] - delta_t * (long_lati_long + lati_lati_lati)
                                    + delta_t * cos(beta) * sin(beta) * u_latitudes[i][j] * u_latitudes[i][j]
                                    - delta_t * (pressure(rho[i][j+1]) - pressure(rho[i][j-1])) / (2.0 * delta_latitudes * rho[i][j] * radius * radius)
                                    + delta_t * nu / (radius * radius) * ((1.0 / (cos(beta) * cos(beta))) * (u_latitudes[i-1][j] - 2.0 * u_latitudes[i][j] + u_latitudes[i+1][j]) / (delta_longitudes * delta_longitudes)
                                    - tan(beta) * (u_latitudes[i][j+1] - u_latitudes[i][j-1]) / (2.0 * delta_latitudes)
                                    + (u_latitudes[i][j-1] - 2.0 * u_latitudes[i][j] + u_latitudes[i][j+1]) / (delta_latitudes * delta_latitudes));
            rho_new[i][j] = rho[i][j] - delta_t * (1.0 / (radius * radius) * (rho[i+1][j] * u_longitudes[i+1][j] - rho[i-1][j] * u_longitudes[i-1][j]) / (2.0 * delta_longitudes)
                                    + (rho[i][j+1] * u_latitudes[i][j+1] - rho[i][j-1] * u_latitudes[i][j-1]) / (2.0 * delta_latitudes));
        }
    }
    u_longitudes = u_longitudes_new;
    u_latitudes = u_latitudes_new;
    rho = rho_new;
    // boundary condition
    for(size_t i = 0; i < N_longitudes; i++)
    {
        u_longitudes[i][0] = 1.0e-1; // u_longitudes[i][N_longitudes-2];
        u_longitudes[i][N_latitules-1] = 1.0e-1; // u_longitudes[i][1];
        u_latitudes[i][0] = 0.0; // u_latitudes[i][N_longitudes-2];
        u_latitudes[i][N_latitules-1] = 0.0; // u_latitudes[i][1];
        rho[i][0] = rho[i][1]; // rho[i][N_longitudes-2];
        rho[i][N_latitules-1] = rho[i][N_latitules-2]; // rho[i][1];
    }
    for (size_t i = 0; i < N_longitudes; i++)
    {
        for (size_t j = 0; j < N_latitules; j++)
        {
            double alpha = 2.0 * M_PI * i / (N_longitudes - 1);
            //double beta = - limit_latitudes + (2.0 * limit_latitudes) * j / (N_latitules - 1);
            u_longitudes[i][0] = 0.2 * sin(2.0 * alpha); // u_longitudes[i][N_longitudes-2];
            u_longitudes[i][N_latitules-1] = 0.2 * cos(2.0 * alpha); // u_longitudes[i][1];
        }
    }
    for(size_t i = 0; i < N_latitules; i++)
    {
        u_longitudes[0][i] = u_longitudes[N_longitudes-2][i];
        u_longitudes[N_longitudes-1][i] = u_longitudes[1][i];
        u_latitudes[0][i] = u_latitudes[N_longitudes-2][i];
        u_latitudes[N_longitudes-1][i] = u_latitudes[1][i];
        rho[0][i] = rho[N_longitudes-2][i];
        rho[N_longitudes-1][i] = rho[1][i];
    }
}

int main()
{
    auto u = FlowInSphere(300, 200);
    u.setVelocity();
    //u.printVelocity();
    //u.step_advection_eq();
    for(auto i = 0; i < 10000000; i++)
    {
        u.step_advection_eq();
        //u.printVelocity();
        if (i % 1000 == 0){
            std::cerr << std::fixed << std::setprecision(15) << i << " " << u.max_velocity() << std::endl;
        }
    }
    u.printVelocity();
    return 0;
}