#include <fstream>
#include <iostream>
#include <Partio.h>
#include <sys/stat.h>
#include "material_point.h"

template <class T, int dim>
class SimulationDriver
{
public:
    using TV = Eigen::Matrix<T, dim, 1>;

    T dt;
    int point_number;
    material_point<T, dim> mp;

    explicit SimulationDriver(int pnt_num)
    {
        dt = (T)0.001;
        point_number = pnt_num;
        mp = material_point<T, dim>(point_number);
    }

    void advanceOneStep()
    {
        /* clean initialization of grid data */
        std::vector<int> active_nodes;
        std::vector<T> grid_mass = std::vector<T>(mp.grd_node_number, 0);
        std::vector<TV> grid_velocity = std::vector<TV>(mp.grd_node_number, TV::Zero());
        std::vector<TV> grid_force = std::vector<TV>(mp.grd_node_number, TV::Zero());

        /* point to grid */
        mp.trans_pnt_grd(grid_mass, grid_velocity, grid_force, active_nodes, dt);

        /* Boundary conditions */
        mp.boundary_collision(grid_velocity, 3);

        /* grid to point */
        mp.trans_grd_pnt(grid_mass, grid_velocity, grid_force, dt);
    }

    /* dump a bgeo file */
    void dumpBgeo(const std::string &filename)
    {
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH, vH, mH;
        mH = parts->addAttribute("m", Partio::VECTOR, 1);
        posH = parts->addAttribute("position", Partio::VECTOR, 3);
        vH = parts->addAttribute("v", Partio::VECTOR, 3);
        for (int i = 0; i < mp.point_number; i++)
        {
            int idx = parts->addParticle();
            float *m = parts->dataWrite<float>(mH, idx);
            float *p = parts->dataWrite<float>(posH, idx);
            float *v = parts->dataWrite<float>(vH, idx);
            m[0] = (T)mp.point_mass[i];
            for (int j = 0; j < 3; j++)
            {
                p[j] = (T)mp.point_position[i][j];
                v[j] = (T)mp.point_velocity[i][j];
            }
        }
        Partio::write(filename.c_str(), *parts);
        parts->release();
    }

    /* run Euler step like project 1 */
    void run(const int max_frame)
    {
        for (int frame = 1; frame <= max_frame; frame++)
        {
            std::cout << "Frame " << frame << std::endl;
            int N_substeps = (int)(((T)1 / 24) / dt);
            for (int step = 1; step <= N_substeps; step++)
            {
                std::cout << "Step " << step << std::endl;
                advanceOneStep();
            }
            std::string output_folder = "output/";
            mkdir(output_folder.c_str(), 0777);
            std::string filename = output_folder + std::to_string(frame) + ".bgeo";
            dumpBgeo(filename);
            std::cout << std::endl;
        }
    }
};
