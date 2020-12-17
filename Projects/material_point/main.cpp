#include <regex>
#include <random>
#include <chrono>
#include <sstream>
#include <Eigen/Core>

#include "SimulationDriver.h"
#include "mesh_query0.1/mesh_query.h"

/* read the obj file */
void read_obj(std::vector<std::vector<float>> &vertices, std::vector<std::vector<int>> &edges)
{
    std::fstream obj("data/cube.obj");
    std::string obj_line;
    while (std::getline(obj, obj_line))
    {
        /* read vertices with v */
        if (obj_line[0] == 'v')
        {
            std::vector<float> v_buf;
            std::regex space_reg("\\s+"); // split with white space
            std::vector<std::string> split_line(
                std::sregex_token_iterator(obj_line.begin(), obj_line.end(), space_reg, -1),
                std::sregex_token_iterator());
            for (size_t i = 0; i < split_line.size(); i++)
                if (i != 0)
                    v_buf.push_back((std::stof(split_line[i])));
            vertices.push_back(v_buf);
        }

        /* read faces with f */
        if (obj_line[0] == 'f')
        {
            std::vector<int> f_buf;
            std::regex space_reg("\\s+");
            std::vector<std::string> split_line(
                std::sregex_token_iterator(obj_line.begin(), obj_line.end(), space_reg, -1),
                std::sregex_token_iterator());
            for (size_t i = 0; i < split_line.size(); i++)
                if (i != 0)
                    f_buf.push_back((std::stoi(split_line[i])));
            edges.push_back((f_buf));
        }
    }
}

int main()
{
    using T = float;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T, dim, 1>;

    /* read obj file, for mesh query */
    std::vector<std::vector<float>> vertices;
    std::vector<std::vector<int>> edges;
    read_obj(vertices, edges);

    /* parse vertices */
    int num_vertices = (int)vertices.size();
    double vertex_position[num_vertices * 3];
    for (size_t i = 0; i < vertices.size(); i++)
        for (size_t j = 0; j < vertices[i].size(); j++)
            vertex_position[i * 3 + j] = vertices[i][j];

    /* parse face, split rectangle into 2 triangles */
    int num_triangle = edges.size();
    num_triangle = num_triangle * 2;
    int triangle_idx[num_triangle * 3];
    std::vector<std::vector<int>> triangles;
    for (size_t i = 0; i < edges.size(); i++)
    {
        triangle_idx[i * 6] = edges[i][0] - 1;
        triangle_idx[i * 6 + 1] = edges[i][1] - 1;
        triangle_idx[i * 6 + 2] = edges[i][2] - 1;

        triangle_idx[i * 6 + 3] = edges[i][0] - 1;
        triangle_idx[i * 6 + 4] = edges[i][2] - 1;
        triangle_idx[i * 6 + 5] = edges[i][3] - 1;
    }

    MeshObject *mesh_object = construct_mesh_object(num_vertices, &vertex_position[0], num_triangle,
                                                    &triangle_idx[0]);

    /* parameter definition, use 50 grid in [0, 1] region */
    std::vector<TV> point_position;

    /* generate random number in [0.0001, 0.01), use 0.0001 to avoid random give 0 (point overlap with grid node) */
    std::random_device r;
    std::default_random_engine rand(r());
    std::uniform_real_distribution<float> d(0.0001, 0.01);

    for (int grid_x = 0; grid_x < 50; grid_x++) // loop through grid
    {
        for (int grid_y = 0; grid_y < 50; grid_y++)
        {
            for (int grid_z = 0; grid_z < 50; grid_z++)
            {
                T x = 0.02 * grid_x; // transfer grid index to coordinates
                T y = 0.02 * grid_y;
                T z = 0.02 * grid_z;
                for (int i = 0; i < 2; i++) // loop inside gird, sample 8 point
                {
                    for (int j = 0; j < 2; j++)
                    {
                        for (int k = 0; k < 2; k++)
                        {
                            TV point_coordinate;

                            /* coordinate change in 8 sub-grid of size 0.01, with random value */
                            point_coordinate[0] = (T)x + i * 0.01 + d(rand);
                            point_coordinate[1] = (T)y + j * 0.01 + d(rand);
                            point_coordinate[2] = (T)z + k * 0.01 + d(rand);

                            double validation[3];
                            validation[0] = point_coordinate[0];
                            validation[1] = point_coordinate[1];
                            validation[2] = point_coordinate[2];

                            /* check if the point is inside the object */
                            if (point_inside_mesh(&validation[0], mesh_object))
                                point_position.push_back(point_coordinate);
                        }
                    }
                }
            }
        }
    }

    /* start the simulation driver */
    SimulationDriver<T, dim> driver((int)point_position.size());
    std::cout << "===========The simulation driver is started===========" << std::endl;
    std::cout << "======There are " << point_position.size() << " points in this simulation======" << std::endl;
    driver.mp.point_position = point_position;
    driver.run(180);
    return 0;
}
