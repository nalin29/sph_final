#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/pick.h"
#include <igl/jet.h>
#include <Eigen/Sparse>
#include <omp.h>

const float G = 10;
const double GRID_WIDTH = 50;
const float GRID_BOTTOM = 0;
bool mouseClicked = false;
Eigen::Vector2d mousePos = Eigen::Vector2d::Zero();

struct SimulationData
{
    int gridn;
    double gridh;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    int nparticles;

    double timestep;
    double viscosity;
    double mouseStrength;
    double gas_const;
    double r;
    double viscLap;
    double k;
    double k_near;
    double rest_density;

    Eigen::MatrixXd u;
    Eigen::VectorXd p;
    Eigen::VectorXd p_near;
    Eigen::MatrixXd markerParticles;
    Eigen::MatrixXd prevMarkerParticles;

    Eigen::MatrixXd force;
    Eigen::VectorXd density;
    Eigen::VectorXd density_near;
};

void initSPH(const unsigned int N, SimulationData &simdata)
{
    int count = 0;
    for (float y = GRID_BOTTOM + GRID_WIDTH/2; y <= 50000; y += simdata.r * 0.5)
    {
        for (float x = 0; x <= GRID_WIDTH; x += simdata.r * 0.5)
        {
            if (count > N -1)
            {
                return;
            }
            simdata.markerParticles.row(count) = Eigen::Vector3d(x, y, 0);
            simdata.prevMarkerParticles.row(count) = simdata.markerParticles.row(count);
            count++;
        }
    }
}

// Initialize a new simulation

void makeGrid(int gridn, SimulationData &result)
{
    result.gridn = gridn;
    result.gridh = 2.0 / double(gridn - 1);
    result.V.resize(4, 3);
    result.V.row(0) = Eigen::Vector3d(-50, 0, 0);
    result.V.row(1) = Eigen::Vector3d(50, 0, 0);
    result.V.row(2) = Eigen::Vector3d(-50, 100, 0);
    result.V.row(3) = Eigen::Vector3d(50, 100, 0);
    gridn = 2;
    result.F.resize(2 * (gridn - 1) * (gridn - 1), 3);
    for (int i = 0; i < gridn - 1; i++)
    {
        for (int j = 0; j < gridn - 1; j++)
        {
            int idx = 2 * (i * (gridn - 1) + j);
            result.F(idx, 0) = i * gridn + j;
            result.F(idx, 1) = i * gridn + (j + 1);
            result.F(idx, 2) = (i + 1) * gridn + (j + 1);
            result.F(idx + 1, 0) = i * gridn + j;
            result.F(idx + 1, 1) = (i + 1) * gridn + (j + 1);
            result.F(idx + 1, 2) = (i + 1) * gridn + j;
        }
    }

    result.u.resize(result.nparticles, 2);
    result.u.setZero();

    result.p.resize(result.nparticles);
    result.p.setZero();

    result.p_near.resize(result.nparticles);
    result.p_near.setZero();

    result.force.resize(result.nparticles, 2);
    result.force.setZero();

    result.density.resize(result.nparticles);
    result.density.setZero();

    result.density_near.resize(result.nparticles);
    result.density_near.setZero();

    result.markerParticles.resize(result.nparticles, 3);
    result.prevMarkerParticles.resize(result.nparticles, 3);

    initSPH(result.nparticles, result);

    result.markerParticles.col(2).setZero();
}

void integrate(SimulationData &simdata)
{
    #pragma omp parallel for 
    for (int i = 0; i < simdata.nparticles; i++)
    {
        simdata.markerParticles.row(i).segment(0, 2) += simdata.force.row(i) * 0.0005;
        simdata.force.row(i) = Eigen::Vector2d(0, -G);
        simdata.u.row(i) = simdata.markerParticles.row(i).segment(0, 2) - simdata.prevMarkerParticles.row(i).segment(0, 2);

        const double max_vel2 = 4.0f;
        const double v2 = simdata.u.row(i).squaredNorm();

        if (v2 > max_vel2)
        {
            simdata.u.row(i) *= 0.5;
        }

        simdata.prevMarkerParticles.row(i) = simdata.markerParticles.row(i);
        simdata.markerParticles.row(i).segment(0, 2) += simdata.u.row(i);

        if (simdata.markerParticles(i, 0) < -GRID_WIDTH)
            simdata.force(i, 0) -= (simdata.markerParticles(i, 0) - -GRID_WIDTH) * 400;
        if (simdata.markerParticles(i, 0) > GRID_WIDTH)
            simdata.force(i, 0) -= (simdata.markerParticles(i, 0) - GRID_WIDTH) * 400;
        if (simdata.markerParticles(i, 1) < GRID_BOTTOM)
            simdata.force(i, 1) -= (simdata.markerParticles(i, 1) - GRID_BOTTOM) * 400;
        if (simdata.markerParticles(i, 1) > 100)
            simdata.force(i, 1) -= (simdata.markerParticles(i, 1) - 100) * 400;

        if (mouseClicked)
        {
            const double mouse_dist = (simdata.markerParticles.row(i).segment(0, 2).transpose() - mousePos).squaredNorm();
            const double max_dist = GRID_WIDTH / 4;

            if (mouse_dist < max_dist * max_dist)
            {
                simdata.force.row(i) -= simdata.mouseStrength * (simdata.markerParticles.row(i).segment(0, 2).transpose() - mousePos) * 2;
            }
        }
    }
}

void computeDensity(SimulationData &simdata)
{
    #pragma omp parallel for 
    for (int i = 0; i < simdata.nparticles; i++)
    {
        simdata.density(i) = 0;
        simdata.density_near(i) = 0;

        double rho = 0;
        double rho_near = 0;

        for (int j = 0; j < simdata.nparticles; j++)
        {
            if (i == j)
                continue;
            Eigen::Vector2d rij = simdata.markerParticles.row(j).segment(0, 2) - simdata.markerParticles.row(i).segment(0, 2);
            double r_len2 = rij.squaredNorm();
            if (r_len2 < simdata.r * simdata.r)
            {
                double r_len = sqrt(r_len2);

                const double q = 1 - (r_len / simdata.r);
                const double q2 = q * q;
                const double q3 = q2 * q;

                rho += q2;
                rho_near += q3;
            }
        }
        simdata.density(i) = rho;
        simdata.density_near(i) = rho_near;
    }
}

void computePressure(SimulationData &simdata)
{
    #pragma omp parallel for 
    for (int i = 0; i < simdata.nparticles; i++)
    {
        simdata.p(i) = simdata.k * (simdata.density(i) - simdata.rest_density);
        simdata.p_near(i) = simdata.k_near * simdata.density_near(i);
    }
}

void pressureForce(SimulationData &simdata)
{
    #pragma omp parallel for 
    for (int i = 0; i < simdata.nparticles; i++)
    {
        for (int j = 0; j < simdata.nparticles; j++)
        {
            if (i == j)
                continue;
            Eigen::Vector2d rij = simdata.markerParticles.row(j).segment(0, 2) - simdata.markerParticles.row(i).segment(0, 2);
            double r_len2 = rij.squaredNorm();
            if (r_len2 < simdata.r * simdata.r)
            {
                double r_len = sqrt(r_len2);
                const double q = 1 - (r_len / simdata.r);
                const double q2 = q * q;

                const double force_p = q * (simdata.p(i) + simdata.p(j)) + q2 * (simdata.p_near(i) + simdata.p_near(j));

                simdata.force.row(i) -= rij.normalized() * force_p;
            }
        }
    }
}

void viscosityForce(SimulationData &simdata)
{
    #pragma omp parallel for 
    for (int i = 0; i < simdata.nparticles; ++i)
    {
        for (int j = 0; j < simdata.nparticles; j++)
        {
            if (i == j)
                continue;
            Eigen::Vector2d rij = simdata.markerParticles.row(j).segment(0, 2) - simdata.markerParticles.row(i).segment(0, 2);
            double r_len2 = rij.squaredNorm();
            if (r_len2 < simdata.r * simdata.r)
            {
                double r_len = sqrt(r_len2);
                simdata.force.row(i) += simdata.viscosity * (simdata.u.row(j) - simdata.u.row(i)) / (simdata.density(j) + simdata.density_near(j)) * simdata.viscLap * (simdata.r - r_len);
            }
        }
    }
}

void simulateOneStep(SimulationData &simdata)
{
    integrate(simdata);
    computeDensity(simdata);
    computePressure(simdata);
    pressureForce(simdata);
    viscosityForce(simdata);
}

int main(int argc, char *argv[])
{
    polyscope::init();

    SimulationData simdata;
    simdata.viscosity = 1000;
    simdata.mouseStrength = 100.0;
    simdata.nparticles = 1500;
    simdata.gas_const = 4;
    simdata.r = 2.5;
    simdata.k = simdata.gas_const;
    simdata.k_near = simdata.k * 10;
    simdata.viscLap = 40.f / (M_PI * pow(simdata.r, 5.f));
    simdata.mouseStrength = 1.0;
    simdata.rest_density = 3;

    makeGrid(20, simdata);

    // Set up rendering

    polyscope::view::style = polyscope::NavigateStyle::Planar;
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;

    auto *pmesh = polyscope::registerSurfaceMesh("Mesh", simdata.V, simdata.F);

    auto *mparticles = polyscope::registerPointCloud("Marker Particles", simdata.markerParticles);
    mparticles->setEnabled(true);
    mparticles->setPointColor({0.2, 1, 0.2});
    mparticles->setPointRadius(0.005);

    // GUI state

    bool isdragging = false;
    float prevx = 0;
    float prevy = 0;

    polyscope::state::userCallback = [&]() -> void
    {
        bool meshdirty = false;

        int oldsize = simdata.gridn;

        if (ImGui::Button("Reset") || ImGui::InputInt("Num Particles", &simdata.nparticles, 100))
        {
            makeGrid(oldsize, simdata);
            meshdirty = true;
        }

        mouseClicked = false;

        if (ImGui::IsMouseDown(0))
        {
            mouseClicked = true;

            ImGuiIO &io = ImGui::GetIO();
            float mousex = io.MousePos.x;
            float mousey = io.MousePos.y;

            float mousexndc = 2.0 * io.MousePos.x / polyscope::view::windowWidth - 1.0;
            float mouseyndc = 1.0 - 2.0 * io.MousePos.y / polyscope::view::windowHeight;
            glm::vec4 mousendc(mousexndc, mouseyndc, 0.0, 1.0);

            glm::vec4 mouseworld = glm::inverse(pmesh->getTransform()) * glm::inverse(polyscope::view::getCameraViewMatrix()) * glm::inverse(polyscope::view::getCameraPerspectiveMatrix()) * mousendc;
            float worldx = mouseworld[0] / mouseworld[3];
            float worldy = mouseworld[1] / mouseworld[3];

            mousePos = Eigen::Vector2d(worldx, worldy);
        }

        ImGui::InputDouble("Viscosity", &simdata.viscosity, 0.1);
        ImGui::InputDouble("Mouse Strength", &simdata.mouseStrength, 0.1);
        ImGui::InputDouble("Rest Density", &simdata.rest_density, 0.1);
        if (ImGui::InputDouble("Radius of Support", &simdata.r, 0.1))
        {
            simdata.viscLap = 40.f / (M_PI * pow(simdata.r, 5.f));
        }

        if (ImGui::InputDouble("Gas Constant", &simdata.gas_const, 0.1))
        {
            simdata.k = simdata.gas_const / 1000.f;
            simdata.k_near = simdata.k * 10;
        }

        simulateOneStep(simdata);

        // Refresh all rendered geometry

        if (meshdirty)
        {
            pmesh = polyscope::registerSurfaceMesh("Mesh", simdata.V, simdata.F);
        }
        mparticles = polyscope::registerPointCloud("Marker Particles", simdata.markerParticles);
    };

    polyscope::show();
}
