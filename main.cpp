#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Simple explicit solver for compressible flow on a sphere.
// The grid uses longitude (alpha) and latitude (beta) coordinates.
// Polar caps are removed to avoid singularities at the poles.

namespace {
constexpr double LAT_LIMIT = 75.0 * M_PI / 180.0;  // maximum |latitude|
constexpr double GAS_R = 4.0;                      // gas constant
constexpr double TEMPERATURE = 1.0;                // temperature
}  // namespace

struct Field {
    std::size_t nlon, nlat;
    std::vector<double> data;
    Field() = default;
    Field(std::size_t lon, std::size_t lat, double init = 0.0)
        : nlon(lon), nlat(lat), data(lon * lat, init) {}
    double &operator()(std::size_t i, std::size_t j) {
        return data[i * nlat + j];
    }
    double operator()(std::size_t i, std::size_t j) const {
        return data[i * nlat + j];
    }
};

class SphericalFluidSolver {
  public:
    SphericalFluidSolver(std::size_t nlon, std::size_t nlat, double radius,
                         double viscosity, double dt);

    void initializeSolidRotation(double angular_velocity);
    void step();
    double maxSpeed() const;
    void dump(std::ostream &os) const;

  private:
    std::size_t nlon_;  // including ghost cells
    std::size_t nlat_;
    double radius_;
    double nu_;
    double dt_;
    double dlon_;
    double dlat_;

    Field rho_;
    Field u_lon_;
    Field u_lat_;

    double pressure(double rho) const { return rho * GAS_R * TEMPERATURE; }
    void applyBoundary();
};

SphericalFluidSolver::SphericalFluidSolver(std::size_t nlon, std::size_t nlat,
                                           double radius, double viscosity,
                                           double dt)
    : nlon_(nlon), nlat_(nlat), radius_(radius), nu_(viscosity), dt_(dt) {
    dlon_ = 2.0 * M_PI / (nlon_ - 2);
    dlat_ = 2.0 * LAT_LIMIT / (nlat_ - 2);
    rho_ = Field(nlon_, nlat_, 1.0);
    u_lon_ = Field(nlon_, nlat_, 0.0);
    u_lat_ = Field(nlon_, nlat_, 0.0);
}

void SphericalFluidSolver::initializeSolidRotation(double angular_velocity) {
    for (std::size_t i = 0; i < nlon_; ++i) {
        for (std::size_t j = 0; j < nlat_; ++j) {
            double beta = -LAT_LIMIT + dlat_ * j;
            u_lon_(i, j) = angular_velocity * std::cos(beta);
            u_lat_(i, j) = 0.0;
            rho_(i, j) = 1.0;
        }
    }
    applyBoundary();
}

void SphericalFluidSolver::applyBoundary() {
    // Periodic in longitude
    for (std::size_t j = 0; j < nlat_; ++j) {
        u_lon_(0, j) = u_lon_(nlon_ - 2, j);
        u_lon_(nlon_ - 1, j) = u_lon_(1, j);
        u_lat_(0, j) = u_lat_(nlon_ - 2, j);
        u_lat_(nlon_ - 1, j) = u_lat_(1, j);
        rho_(0, j) = rho_(nlon_ - 2, j);
        rho_(nlon_ - 1, j) = rho_(1, j);
    }
    // Simple Neumann condition near excluded poles
    for (std::size_t i = 0; i < nlon_; ++i) {
        u_lon_(i, 0) = u_lon_(i, 1);
        u_lon_(i, nlat_ - 1) = u_lon_(i, nlat_ - 2);
        u_lat_(i, 0) = 0.0;
        u_lat_(i, nlat_ - 1) = 0.0;
        rho_(i, 0) = rho_(i, 1);
        rho_(i, nlat_ - 1) = rho_(i, nlat_ - 2);
    }
}

void SphericalFluidSolver::step() {
    Field u_lon_new(nlon_, nlat_);
    Field u_lat_new(nlon_, nlat_);
    Field rho_new(nlon_, nlat_);

    for (std::size_t i = 1; i < nlon_ - 1; ++i) {
        for (std::size_t j = 1; j < nlat_ - 1; ++j) {
            double beta = -LAT_LIMIT + dlat_ * j;

            double long_long, long_lati, lati_long, lati_lati;
            if (u_lon_(i, j) > 0) {
                long_long =
                    u_lon_(i, j) * (u_lon_(i, j) - u_lon_(i - 1, j)) / dlon_;
                long_lati =
                    u_lon_(i, j) * (u_lat_(i, j) - u_lat_(i - 1, j)) / dlon_;
            } else {
                long_long =
                    u_lon_(i, j) * (u_lon_(i + 1, j) - u_lon_(i, j)) / dlon_;
                long_lati =
                    u_lon_(i, j) * (u_lat_(i + 1, j) - u_lat_(i, j)) / dlon_;
            }
            if (u_lat_(i, j) > 0) {
                lati_long =
                    u_lat_(i, j) * (u_lon_(i, j) - u_lon_(i, j - 1)) / dlat_;
                lati_lati =
                    u_lat_(i, j) * (u_lat_(i, j) - u_lat_(i, j - 1)) / dlat_;
            } else {
                lati_long =
                    u_lat_(i, j) * (u_lon_(i, j + 1) - u_lon_(i, j)) / dlat_;
                lati_lati =
                    u_lat_(i, j) * (u_lat_(i, j + 1) - u_lat_(i, j)) / dlat_;
            }

            u_lon_new(i, j) =
                u_lon_(i, j) - dt_ * (long_long + lati_long)
                - dt_ * 2.0 * std::tan(beta) * u_lon_(i, j) * u_lat_(i, j)
                - dt_ * (pressure(rho_(i + 1, j)) - pressure(rho_(i - 1, j))) /
                      (2.0 * dlon_ * rho_(i, j) * radius_ * radius_ *
                       std::cos(beta) * std::cos(beta))
                + dt_ * nu_ / (radius_ * radius_) *
                      ((1.0 / (std::cos(beta) * std::cos(beta))) *
                           (u_lon_(i - 1, j) - 2.0 * u_lon_(i, j) +
                            u_lon_(i + 1, j)) /
                           (dlon_ * dlon_)
                       - std::tan(beta) *
                             (u_lon_(i, j + 1) - u_lon_(i, j - 1)) /
                             (2.0 * dlat_)
                       + (u_lon_(i, j - 1) - 2.0 * u_lon_(i, j) +
                          u_lon_(i, j + 1)) /
                             (dlat_ * dlat_));

            u_lat_new(i, j) =
                u_lat_(i, j) - dt_ * (long_lati + lati_lati)
                - dt_ * std::cos(beta) * std::sin(beta) *
                      u_lon_(i, j) * u_lon_(i, j)
                - dt_ * (pressure(rho_(i, j + 1)) - pressure(rho_(i, j - 1))) /
                      (2.0 * dlat_ * rho_(i, j) * radius_ * radius_)
                + dt_ * nu_ / (radius_ * radius_) *
                      ((1.0 / (std::cos(beta) * std::cos(beta))) *
                           (u_lat_(i - 1, j) - 2.0 * u_lat_(i, j) +
                            u_lat_(i + 1, j)) /
                           (dlon_ * dlon_)
                       - std::tan(beta) *
                             (u_lat_(i, j + 1) - u_lat_(i, j - 1)) /
                             (2.0 * dlat_)
                       + (u_lat_(i, j - 1) - 2.0 * u_lat_(i, j) +
                          u_lat_(i, j + 1)) /
                             (dlat_ * dlat_));

            rho_new(i, j) =
                rho_(i, j) -
                dt_ * (1.0 / radius_) *
                    ((rho_(i + 1, j) * u_lon_(i + 1, j) -
                      rho_(i - 1, j) * u_lon_(i - 1, j)) /
                         (2.0 * dlon_)+
                     (rho_(i, j + 1) * u_lat_(i, j + 1) -
                      rho_(i, j - 1) * u_lat_(i, j - 1)) /
                         (2.0 * dlat_));
        }
    }

    u_lon_ = std::move(u_lon_new);
    u_lat_ = std::move(u_lat_new);
    rho_ = std::move(rho_new);
    applyBoundary();
}

double SphericalFluidSolver::maxSpeed() const {
    double max = 0.0;
    for (std::size_t i = 0; i < nlon_; ++i) {
        for (std::size_t j = 0; j < nlat_; ++j) {
            double v =
                std::sqrt(u_lon_(i, j) * u_lon_(i, j) +
                          u_lat_(i, j) * u_lat_(i, j));
            max = std::max(max, v);
        }
    }
    return max;
}

void SphericalFluidSolver::dump(std::ostream &os) const {
    for (std::size_t i = 0; i < nlon_; ++i) {
        for (std::size_t j = 0; j < nlat_; ++j) {
            double alpha = 2.0 * M_PI * i / (nlon_ - 2);
            double beta = -LAT_LIMIT + dlat_ * j;
            os << alpha * 180.0 / M_PI << ' ' << beta * 180.0 / M_PI << ' '
               << u_lon_(i, j) << ' ' << u_lat_(i, j) << ' ' << rho_(i, j)
               << ' '
               << std::sqrt(u_lon_(i, j) * u_lon_(i, j) +
                             u_lat_(i, j) * u_lat_(i, j))
               << '\n';
        }
        os << '\n';
    }
}

int main() {
    SphericalFluidSolver solver(300, 200, 1.0, 0.001, 1.0e-5);
    solver.initializeSolidRotation(0.1);

    for (int step = 0; step < 10000; ++step) {
        solver.step();
        if (step % 1000 == 0) {
            std::cerr << std::fixed << std::setprecision(15) << step << ' '
                      << solver.maxSpeed() << '\n';
        }
    }

    solver.dump(std::cout);
    return 0;
}

