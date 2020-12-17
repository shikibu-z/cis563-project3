#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

template <class T, int dim>
class material_point
{
public:
    material_point() = default;
    ;
    using TV = Eigen::Matrix<T, dim, 1>;
    using TM = Eigen::Matrix<T, dim, dim>;

    /* parameter definitions */
    TV gravity;
    T youngs_modulus;
    T nu, mu, lambda, rho;

    /* grid definitions */
    T dx;
    T resolution;
    T grd_node_number;

    /* point definitions */
    T point_number;
    T point_volume;
    std::vector<TV> point_position;
    std::vector<T> point_mass;
    std::vector<TV> point_velocity;
    std::vector<TM> point_deform_grad;

    explicit material_point(int pnt_num)
    {
        /* constant initialization */
        gravity[0] = (T)0;
        gravity[1] = (T)-9.8;
        gravity[2] = (T)0;
        youngs_modulus = (T)50000;
        nu = (T)0.2;
        mu = (T)youngs_modulus / 2 * (1 + nu);
        lambda = (T)youngs_modulus * nu / ((1 + nu) * (1 - 2 * nu));
        rho = 1300;

        /* grid initialization */
        dx = (T)0.02;
        resolution = (T)1 / dx + 1;
        grd_node_number = pow(resolution, 3);

        /* point initialization */
        point_number = pnt_num;
        point_volume = pow(dx, 3) / 8;
        point_mass = std::vector<T>(point_number, rho * point_volume);
        point_velocity = std::vector<TV>(point_number, TV::Zero());
        point_deform_grad = std::vector<TM>(point_number, TM::Identity());
    }

    /* The polar svd implementation */
    void polar_svd(TM &u, TM &v, TM &deform_grad)
    {
        /* svd decomposition (ignore sigma, never used) */
        Eigen::JacobiSVD<TM> svd_decomp(deform_grad, Eigen::ComputeThinU | Eigen::ComputeThinV);
        u = svd_decomp.matrixU();
        v = svd_decomp.matrixV();

        /* polar changes */
        if (u.determinant() < 0)
            u.col(2) = -u.col(2);

        if (v.determinant() < 0)
            v.col(2) = -v.col(2);
    }

    /* Computing the kirchoff stress by differentiating Psi */
    void kirchoff_fixed_corotate(TM &stress, TM &deform_grad)
    {
        TM u, v;
        polar_svd(u, v, deform_grad);
        TM rotate = u * v.transpose();
        T jacobi = deform_grad.determinant();
        stress =
            T(2) * mu * (deform_grad - rotate) + lambda * (jacobi - 1) * jacobi * deform_grad.inverse().transpose();
    }

    /* Compute the interpolation weight */
    void inter_weight(TM &w, TM &dw, TV &base, TV &x)
    {
        base = (x.array() - 0.5).floor(); // lower left node
        TV fx = x - base;

        /* quadratic spline weights */
        w.col(0) = 0.5 * (1.5 - fx.array()).pow(2);
        w.col(1) = 0.75 - (fx.array() - 1).pow(2);
        w.col(2) = 0.5 * (fx.array() - 0.5).pow(2);

        dw.col(0) = fx.array() - 1.5;
        dw.col(1) = -2.0 * (fx.array() - 1);
        dw.col(2) = fx.array() - 0.5;

        /* transpose for correct column behavior */
        w.transposeInPlace();
        dw.transposeInPlace();
    }

    /* Transfer point mass and velocity to grid node */
    void trans_pnt_grd(std::vector<T> &grid_mass, std::vector<TV> &grid_velocity, std::vector<TV> &grid_force,
                       std::vector<int> valid_grd_node, T &dt)
    {
        for (int p_idx = 0; p_idx < point_number; p_idx++)
        {
            TV x = point_position[p_idx] / dx;
            TM deform_grad = point_deform_grad[p_idx];
            TM stress, w, dw;
            TV base;

            /* compute stress and weight */
            kirchoff_fixed_corotate(stress, deform_grad);
            inter_weight(w, dw, base, x);

            /* iterate over the 27 related grid node */
            for (int i = 0; i < 3; i++)
            {
                T w_x = w.col(0)[i];
                T dw_x_x = dw.col(0)[i] / dx;
                T node_x = base[0] + i;
                for (int j = 0; j < 3; j++)
                {
                    T w_xy = w_x * w.col(1)[j];
                    T dw_xy_x = dw_x_x * w.col(1)[j];
                    T dw_xy_y = w_x / dx * dw.col(1)[j];
                    T node_y = base[1] + j;
                    for (int k = 0; k < 3; k++)
                    {
                        T w_xyz = w_xy * w.col(2)[k];
                        T dw_xyz_x = dw_xy_x * w.col(2)[k];
                        T dw_xyz_y = dw_xy_y * w.col(2)[k];
                        T dw_xyz_z = w_xy / dx * dw.col(2)[k];
                        T node_z = base[2] + k;

                        TV dw_gradient;
                        dw_gradient[0] = dw_xyz_x;
                        dw_gradient[1] = dw_xyz_y;
                        dw_gradient[2] = dw_xyz_z;

                        int idx = node_x * pow(resolution, 2) + node_y * resolution + node_z;
                        grid_mass[idx] = grid_mass[idx] +
                                         w_xyz * point_mass[p_idx]; // mass
                        grid_velocity[idx] = grid_velocity[idx] +
                                             point_mass[p_idx] * point_velocity[p_idx] * w_xyz; // momentum
                        grid_force[idx] = grid_force[idx] + (-point_volume * stress * deform_grad.transpose()) *
                                                                dw_gradient; // elastic force
                    }
                }
            }
        }

        /* go over grid node for 0 mass */
        for (int i = 0; i < grd_node_number; i++)
        {
            if (grid_mass[i] != 0)
            {
                valid_grd_node.push_back(i);
                grid_velocity[i] = grid_velocity[i] / grid_mass[i];
            }
            else
                grid_velocity[i] = TV::Zero();
        }

        /* add gravity, update grid velocity */
        for (int &iter : valid_grd_node)
        {
            grid_force[iter] = grid_force[iter] + grid_mass[iter] * gravity;
            grid_velocity[iter] = grid_velocity[iter] + dt * grid_force[iter] / grid_mass[iter];
        }
    }

    void boundary_collision(std::vector<TV> &grid_velocity, int ground_layer)
    {
        /* apply boundary condition at the lowest 3 grid on y direction  - hit ground */
        for (int i = 0; i < ground_layer; i++)
            for (int j = 0; j < resolution; j++)
                for (int k = 0; k < resolution; k++)
                    grid_velocity[i * resolution + j * pow(resolution, 2) + k] = TV::Zero();
    }

    void trans_grd_pnt(std::vector<T> &grid_mass, std::vector<TV> &grid_velocity, std::vector<TV> &grid_force, T &dt)
    {
        for (int p_idx = 0; p_idx < point_number; p_idx++)
        {
            /* query values for calculation */
            TV x = point_position[p_idx] / dx;
            TM deform_grad = point_deform_grad[p_idx];
            TM w, dw;
            TV base;
            TM point_velocity_grad = TM::Zero(3, 3);
            TV point_velocity_update = TV::Zero();

            /* compute weight */
            inter_weight(w, dw, base, x);

            /* iterate over the 27 related grid node */
            for (int i = 0; i < 3; i++)
            {
                T w_x = w.col(0)[i];
                T dw_x_x = dw.col(0)[i] / dx;
                T node_x = base[0] + i;
                for (int j = 0; j < 3; j++)
                {
                    T w_xy = w_x * w.col(1)[j];
                    T dw_xy_x = dw_x_x * w.col(1)[j];
                    T dw_xy_y = w_x / dx * dw.col(1)[j];
                    T node_y = base[1] + j;
                    for (int k = 0; k < 3; k++)
                    {
                        T w_xyz = w_xy * w.col(2)[k];
                        T dw_xyz_x = dw_xy_x * w.col(2)[k];
                        T dw_xyz_y = dw_xy_y * w.col(2)[k];
                        T dw_xyz_z = w_xy / dx * dw.col(2)[k];
                        T node_z = base[2] + k;

                        TV dw_gradient;
                        dw_gradient[0] = dw_xyz_x;
                        dw_gradient[1] = dw_xyz_y;
                        dw_gradient[2] = dw_xyz_z;

                        int idx = node_x * pow(resolution, 2) + node_y * resolution + node_z;
                        point_velocity_grad = point_velocity_grad + grid_velocity[idx] *
                                                                        dw_gradient.transpose(); // compute velocity gradient for deform gradient
                        point_velocity_update =
                            point_velocity_update + grid_velocity[idx] * w_xyz; // compute velocity update
                    }
                }
            }

            /* update deformation gradient */
            TM new_deform_grad;
            new_deform_grad = deform_grad * (TM::Identity() + point_velocity_grad * dt);
            point_deform_grad[p_idx] = new_deform_grad;

            /* udpate point position and velocity */
            point_velocity[p_idx] = point_velocity_update;
            point_position[p_idx] = point_position[p_idx] + point_velocity_update * dt;
        }
    }
};
