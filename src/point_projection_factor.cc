//
// Created by lixin04 on 19年12月6日.
//

#include "point_projection_factor.hpp"


Eigen::Matrix2d ProjectionPointInverseDepthFactor::sqrt_info;
double ProjectionPointInverseDepthFactor::sum_t;

ProjectionPointInverseDepthFactor::ProjectionPointInverseDepthFactor(Eigen::Vector3d& ob0, Eigen::Vector3d& obi) :point_obi_(obi), point_ob0_(ob0){}

bool ProjectionPointInverseDepthFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d twc0(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond qwc0(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d twci(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qwci(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

//    Eigen::Vector3d Point(parameters[1][0], parameters[1][1], parameters[1][2]);
    double inverse_depth = parameters[2][0];

    Eigen::Vector3d pts_camera_0 = point_ob0_ / inverse_depth;
    Eigen::Vector3d pts_w = qwc0 * pts_camera_0 + twc0;
    Eigen::Vector3d pts_camera_i = qwci.inverse() * (pts_w - twci);

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual(0) = pts_camera_i(0)/pts_camera_i(2) - point_obi_(0);
    residual(1) = pts_camera_i(1)/pts_camera_i(2) - point_obi_(1);

    residual = sqrt_info * residual;
//    std::cout << residual.transpose() << std::endl;

    if (jacobians)
    {
        Eigen::Matrix3d Rwc0 = qwc0.toRotationMatrix();
        Eigen::Matrix3d Rwci = qwci.toRotationMatrix();
        Eigen::Matrix<double, 2, 3> reduce(2, 3);
        reduce << 1. / pts_camera_i(2), 0, -pts_camera_i(0) / (pts_camera_i(2) * pts_camera_i(2)),
                0, 1. / pts_camera_i(2), -pts_camera_i(1) / (pts_camera_i(2) * pts_camera_i(2));
        reduce = sqrt_info * reduce;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose0(jacobians[0]);

            Eigen::Matrix3d Inde = Eigen::Matrix3d::Identity();
            Eigen::Matrix<double, 3, 6> jaco_pose0;
            jaco_pose0.leftCols<3>() = Rwci.transpose();
            jaco_pose0.rightCols<3>() = -Rwci.transpose()*Rwc0*Utility::skewSymmetric(pts_camera_0);
            jacobian_pose0.leftCols<6>() = reduce * jaco_pose0;
            jacobian_pose0.rightCols<1>().setZero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_posei(jacobians[1]);

            Eigen::Matrix3d Inde = Eigen::Matrix3d::Identity();
            Eigen::Matrix<double, 3, 6> jaco_posei;
            jaco_posei.leftCols<3>() = -Rwci.transpose();
            jaco_posei.rightCols<3>() = Utility::skewSymmetric(Rwci.transpose() * (pts_w - twci));
            jacobian_posei.leftCols<6>() = reduce * jaco_posei;
            jacobian_posei.rightCols<1>().setZero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Vector2d> jacobian_point(jacobians[2]);

            jacobian_point = - reduce * Rwci.transpose()*Rwc0*point_ob0_/(inverse_depth * inverse_depth);

        }

        // Jacobian Check
        if(0 && jacobians[0] && jacobians[1] && jacobians[2])
        {
            const double eps = 1e-6;
            Eigen::Matrix<double, 2, 13> num_jacobian;
            for (int i = 0; i < 13; ++i) {
                Eigen::Vector3d twc0_check(parameters[0][0], parameters[0][1], parameters[0][2]);
                Eigen::Quaterniond qwc0_check(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

                Eigen::Vector3d twci_check(parameters[1][0], parameters[1][1], parameters[1][2]);
                Eigen::Quaterniond qwci_check(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

                double inverse_depth_check = parameters[2][0];

                {
                    int b = i % 3, a = i/3;
                    Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                    if (a == 0)
                        twc0_check = twc0_check + delta;
                    else if (a == 1)
                        qwc0_check = qwc0_check * Utility::deltaQ(delta);
                    else if (a == 2)
                        twci_check = twci_check + delta;
                    else if (a == 3)
                        qwci_check = qwci_check * Utility::deltaQ(delta);
                    else
                        inverse_depth_check += eps;
                }

                Eigen::Vector3d pts_camera_0_check = point_ob0_ / inverse_depth_check;
                Eigen::Vector3d pts_w_check = qwc0_check * pts_camera_0_check + twc0_check;
                Eigen::Vector3d pts_camera_i_check = qwci_check.inverse() * (pts_w_check - twci_check);

                Eigen::Vector2d residual_tmp;
                residual_tmp(0) = pts_camera_i_check(0)/pts_camera_i_check(2) - point_obi_(0);
                residual_tmp(1) = pts_camera_i_check(1)/pts_camera_i_check(2) - point_obi_(1);

                residual_tmp = sqrt_info * residual_tmp;
                num_jacobian.col(i) = (residual_tmp - residual)/eps;
            }

            // check jacobian
            std::cout << "ana = " << std::endl;
            std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jacobians[0]) << std::endl
                      << std::endl;
            std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jacobians[1]) << std::endl
                      << std::endl;
            std::cout << Eigen::Map<Eigen::Vector2d>(jacobians[2]) << std::endl
                      << std::endl;
            std::cout <<"num_jacobian:\n"<< num_jacobian.leftCols(6) <<"\n" << num_jacobian.rightCols(7) <<"\n"<< std::endl;

            std::cout << " ------------------- "<<std::endl;
        }
    }

    return true;
}


double PointOnPlaneFactor::sqrt_info;// = 1.0;
PointOnPlaneFactor::PointOnPlaneFactor(const Eigen::Vector3d& norm, const Eigen::Matrix<double, 3, 2> &tange, const Eigen::Vector3d &_pts_i):
        tangent_base(tange), norm_n(norm), pts_i(_pts_i)
{
};

/*
  parameters[0]:  Twi
  parameters[1]:  Twj
  parameters[2]:  Tbc
  parameters[3]:  line_orth
*/
bool PointOnPlaneFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    double distance = parameters[0][2];
    Eigen::Vector2d tange_n( parameters[0][0], parameters[0][1] );

    Eigen::Vector3d nom_n_now = norm_n + tangent_base * tange_n;

    Eigen::Vector3d Pi(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qi(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    double inv_dep_i = parameters[3][0];

    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;

    residuals[0] = sqrt_info * ( nom_n_now.dot(pts_w) + distance );


//    std::cout << "-----------" << 100 * (nom_n_now.dot(pts_w) + distance) << "---------" << std::endl;

//    std::cout << "sqrt_info:" << sqrt_info << std::endl <<  pts_w << std::endl <<"distance:" << distance << std::endl << "residuals[0]:" << residuals[0] << std::endl;
    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_plan(jacobians[0]);
            Eigen::Vector2d jacob = tangent_base.transpose() * pts_w;
            jacobian_plan(0,0) = jacob(0);
            jacobian_plan(0,1) = jacob(1);
            jacobian_plan(0,2) = 1;

            jacobian_plan =  sqrt_info * jacobian_plan;
        }

        Eigen::Vector3d x  = nom_n_now.transpose()  ;
        Eigen::Vector3d jacobian_word(x(0), x(1), x(2));

        jacobian_word = sqrt_info * jacobian_word;

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose(jacobians[1]);

            jacobian_pose.leftCols<3>() = jacobian_word;
            Eigen::Matrix3d Qi_R = Qi.toRotationMatrix();
            Eigen::Vector3d jacob_Qi = - jacobian_word.transpose() * Qi_R * Utility::skewSymmetric(pts_imu_i);
            jacobian_pose(0, 3) = jacob_Qi(0);
            jacobian_pose(0, 4) = jacob_Qi(1);
            jacobian_pose(0, 5) = jacob_Qi(2);
            jacobian_pose(0, 6) = 0;
        }


        Eigen::Matrix3d Qi_R = Qi.toRotationMatrix();
        Eigen::Vector3d jacobian_imu = jacobian_word.transpose() * Qi_R;


        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_ex(jacobians[2]);

            jacobian_ex.leftCols<3>() = jacobian_imu;
            Eigen::Matrix3d qic_R = qic.toRotationMatrix();
            Eigen::Vector3d jacob_qic = - jacobian_imu.transpose() * qic_R * Utility::skewSymmetric(pts_camera_i);
            jacobian_ex(0, 3) = jacob_qic(0);
            jacobian_ex(0, 4) = jacob_qic(1);
            jacobian_ex(0, 5) = jacob_qic(2);
            jacobian_ex(0, 6) = 0;
        }

        Eigen::Matrix3d qic_R = qic.toRotationMatrix();
        Eigen::Vector3d jacobian_camera = jacobian_imu.transpose() * qic_R;

        if(jacobians[3])
        {
            double jac = jacobian_camera.transpose() * pts_i;
            jacobians[3][0] = - jac / (inv_dep_i * inv_dep_i);
        }



        // Jacobian Check
        if(false)
        {
            const double eps = 1e-6;
            Eigen::Matrix<double, 1, 16> num_jacobian;
            for (int i = 0; i < 16; ++i) {

                Eigen::Vector3d plane_check( parameters[0][0], parameters[0][1], parameters[0][2] );

                Eigen::Vector3d Pi_check(parameters[1][0], parameters[1][1], parameters[1][2]);
                Eigen::Quaterniond Qi_check(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

                Eigen::Vector3d tic_check(parameters[2][0], parameters[2][1], parameters[2][2]);
                Eigen::Quaterniond qic_check(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

                double inv_dep_i_check = parameters[3][0];

//                Eigen::Vector3d point_ckeck(parameters[1][0], parameters[1][1], parameters[1][2]);

                int a = i / 3, b = i % 3;
                Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                if (a == 0)
                    plane_check += delta;
                else if (a == 1)
                    Pi_check += delta;
                else if (a == 2)
                    Qi_check = Qi_check * Utility::deltaQ(delta);
                else if (a == 3)
                    tic_check += delta;
                else if (a == 4)
                    qic_check = qic_check * Utility::deltaQ(delta);
                else if (a == 5)
                    inv_dep_i_check += delta.x();

                double distance_check = plane_check(2);
                Eigen::Vector2d tange_n_check( plane_check(0), plane_check(1) );
                Eigen::Vector3d nom_n_now_check = norm_n + tangent_base * tange_n_check;
                Eigen::Vector3d pts_camera_i_check = pts_i / inv_dep_i_check;
                Eigen::Vector3d pts_imu_i_check = qic_check * pts_camera_i_check + tic_check;
                Eigen::Vector3d pts_w_check = Qi_check * pts_imu_i_check + Pi_check;

                // 误差
                double tmp_residual = sqrt_info * ( nom_n_now_check.dot(pts_w_check) + distance_check );

                num_jacobian[i] = (tmp_residual - residuals[0])/eps;

                std::cout << " point : error " << residuals[0] << ":" << tmp_residual << std::endl;
            }
            std::cout << " analysis:    ";
            for(int i=0;i<3;++i)
                std::cout << jacobians[0][i] << "  ";
            std::cout << std::endl;
            for(int i=0;i<6;++i)
                std::cout << jacobians[1][i] << "  ";
            std::cout << std::endl;
//            for(int i=0;i<6;++i)
//                std::cout << jacobians[2][i] << "  ";
//            std::cout << std::endl;
            std::cout << jacobians[3][0] << std::endl;
            std::cout << " num:         ";
            for(int i=0;i<3;++i)
                std::cout << num_jacobian[i] << "  ";
            std::cout << std::endl;
            for(int i=3;i<9;++i)
                std::cout << num_jacobian[i] << "  ";
            std::cout << std::endl;
            std::cout << num_jacobian[15] << "  "<< std::endl;
            std::cout << std::endl;
            std::cout << " ------------------- "<<std::endl;

        }
    }

    return true;
}


Eigen::Matrix2d PointOnPlaneProjectionFactor::sqrt_info;// = 1.0;
PointOnPlaneProjectionFactor::PointOnPlaneProjectionFactor(
        const Eigen::Vector3d& norm,
        const Eigen::Matrix<double, 3, 2> &tange,
        const Eigen::Vector3d &_ob0,
        const Eigen::Vector3d &_obi):
        tangent_base(tange), norm_n(norm), ob0(_ob0), obi(_obi)
{
};

/*
  parameters[0]:  Twi
  parameters[1]:  Twj
  parameters[2]:  Tbc
  parameters[3]:  line_orth
*/
bool PointOnPlaneProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    double distance = parameters[0][2];
    Eigen::Vector2d tange_n( parameters[0][0], parameters[0][1] );

    Eigen::Vector3d nom_n_now = norm_n + tangent_base * tange_n;

    Eigen::Vector3d P0(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Q0(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d Pi(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond Qi(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    Eigen::Vector3d nom_n_cam0 = Q0.inverse() * nom_n_now;
    double distance_cam0 = distance + ( Q0.inverse() * P0 ).dot(nom_n_cam0);

    double depth_cam = - distance_cam0 / (nom_n_cam0.dot(ob0));

//    double depth_cam = - ( distance + P0.dot(nom_n_now) ) / (nom_n_now.dot( Q0 * ob0 ));
//    std::cout << "depth_cam = " << depth_cam << std::endl;

    Eigen::Vector3d pts_camera_0 = ob0 * depth_cam;
    Eigen::Vector3d pts_w = Q0 * pts_camera_0 + P0;
    Eigen::Vector3d pts_camera_i = Qi.inverse() * (pts_w - Pi);

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual(0) = pts_camera_i(0)/pts_camera_i(2) - obi(0);
    residual(1) = pts_camera_i(1)/pts_camera_i(2) - obi(1);
    residual = sqrt_info * residual;

//    std::cout << residual.transpose() << std::endl;

    if (jacobians)
    {

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_plane(jacobians[0]);
            jacobian_plane.setZero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_posei(jacobians[1]);
            jacobian_posei.setZero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_posej(jacobians[2]);
            jacobian_posej.setZero();

        }

        // Jacobian Check
        if(1)
        {
            const double eps = 1e-6;
            Eigen::Matrix<double, 2, 15> num_jacobian;
            for (int i = 0; i < 15; ++i) {


                double distance_ck = parameters[0][2];
                Eigen::Vector2d tange_n_ck( parameters[0][0], parameters[0][1] );

                Eigen::Vector3d P0_ck(parameters[1][0], parameters[1][1], parameters[1][2]);
                Eigen::Quaterniond Q0_ck(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

                Eigen::Vector3d Pi_ck(parameters[2][0], parameters[2][1], parameters[2][2]);
                Eigen::Quaterniond Qi_ck(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);



                {
                    int b = i % 3, a = i/3;
                    Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;
                    if (a == 0)
                    {
                        if(b==0) tange_n_ck(0) += eps;
                        else if(b==1) tange_n_ck(1) += eps;
                        else if(b==2) distance_ck += eps;
                    }
                    else if (a == 1)
                        P0_ck = P0_ck + delta;
                    else if (a == 2)
                        Q0_ck = Q0_ck * Utility::deltaQ(delta);
                    else if (a == 3)
                        Pi_ck = Pi_ck + delta;
                    else
                        Qi_ck = Qi_ck * Utility::deltaQ(delta);
                }

                Eigen::Vector3d nom_n_now_ck = norm_n + tangent_base * tange_n_ck;

                Eigen::Vector3d nom_n_cam0_ck = Q0_ck.inverse() * nom_n_now_ck;
                double distance_cam0_ck = distance_ck + ( Q0_ck.inverse() * P0_ck ).dot(nom_n_cam0_ck);

                double depth_cam_ck = - distance_cam0_ck / (nom_n_cam0_ck.dot(ob0));

//              double depth_cam = - ( distance + P0.dot(nom_n_now) ) / (nom_n_now.dot( Q0 * ob0 ));
//              std::cout << "depth_cam = " << depth_cam << std::endl;

                Eigen::Vector3d pts_camera_0_ck = ob0 * depth_cam_ck;
                Eigen::Vector3d pts_w_ck = Q0_ck * pts_camera_0_ck + P0_ck;
                Eigen::Vector3d pts_camera_i_ck = Qi_ck.inverse() * (pts_w_ck - Pi_ck);

                Eigen::Vector2d residual_tmp;
                residual_tmp(0) = pts_camera_i_ck(0)/pts_camera_i_ck(2) - obi(0);
                residual_tmp(1) = pts_camera_i_ck(1)/pts_camera_i_ck(2) - obi(1);

                residual_tmp = sqrt_info * residual_tmp;
                num_jacobian.col(i) = (residual_tmp - residual)/eps;
            }

            if (jacobians[0])
            {
                Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_plane(jacobians[0]);
                jacobian_plane = num_jacobian.block(0,0,2,3);
            }

            if (jacobians[1])
            {
                Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_posei(jacobians[1]);
                jacobian_posei.block(0,0,2,6) = num_jacobian.block(0,3,2,6);
            }

            if (jacobians[2])
            {
                Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_posej(jacobians[2]);
                jacobian_posej.block(0,0,2,6) = num_jacobian.block(0,9,2,6);

            }
            // check jacobian
//            std::cout << "ana = " << std::endl;
//            std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jacobians[0]) << std::endl
//                                                                               << std::endl;
//            std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jacobians[1]) << std::endl
//                                                                               << std::endl;
//            std::cout << Eigen::Map<Eigen::Vector2d>(jacobians[2]) << std::endl
//                      << std::endl;
//            std::cout <<"num_jacobian:\n"<< num_jacobian.leftCols(6) <<"\n" << num_jacobian.rightCols(7) <<"\n"<< std::endl;
//
//            std::cout << " ------------------- "<<std::endl;
        }
    }

    return true;
}
