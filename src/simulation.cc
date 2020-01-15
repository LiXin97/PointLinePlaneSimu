//
// Created by lixin04 on 20年1月15日.
//

#include "simulation.hpp"

Simulation::Simulation(double max_nt, double max_nq, double max_pixel_n)
{
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> nt(0., max_nt);
    std::normal_distribution<double> nq(0., max_nq);
    std::normal_distribution<double> pixel_n(0, max_pixel_n);

    generator_ = generator;
    nt_ = nt; nq_ = nq; pixel_n_ = pixel_n;

    constant_pose_num_ = 8;
}

void Simulation::create_keyframe(int frame_num, double radius, int wave_num, double wave_high)
{
    Twcs_.clear();
    for(int n = 0; n < frame_num; ++n )
    {
        // 绕 z轴 旋转
        double theta = n * 2. * M_PI / frame_num;
        double step = n * 1.  / frame_num;
        double wave_theta = n * 2. * M_PI/(frame_num/wave_num);

        Eigen::Matrix4d Twc = Eigen::Matrix4d::Identity();
        Eigen::Vector3d twb = Eigen::Vector3d( radius * cos(theta) , radius * sin(theta), wave_high * std::cos(wave_theta));
        Eigen::Matrix3d Rwb = Eigen::Matrix3d::Identity();

        double theta_ = 0;
        if(theta < M_PI / 2.) theta_ = theta + M_PI/2.;
        else if(theta >= M_PI /2. && theta < M_PI ) theta_ = theta - M_PI*3./2.;
        else if(theta >= M_PI  && theta < M_PI * 3. / 2. ) theta_ = theta+M_PI/2.;
        else if(theta >= M_PI * 3./2. && theta < M_PI * 2 ) theta_ = theta-3.*M_PI/2.;
        Rwb = Eigen::AngleAxisd(theta_, Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(M_PI / 2., Eigen::Vector3d::UnitX());

        Twc.block(0, 0, 3, 3) = Rwb;
        Twc.block(0, 3, 3, 1) = twb;

        Twcs_.push_back(Twc);
        Twcs_true_.push_back(Twc);
    }
}

void Simulation::generate_point_line()
{
    const double distance = 5.;

    std::uniform_real_distribution<double> point_generate( -2., 2. );

    for(int i=0;i<100;++i)
    {
        points_true_.emplace_back( distance, point_generate(generator_), point_generate(generator_) );
    }

    for(int i=0;i<20;++i)
    {
        Eigen::Matrix<double, 3, 2> line;
        line(0,0) = distance;
        line(0,1) = distance;
        line(1,0) = point_generate(generator_);
        line(1,1) = line(1,0);
        line(2,0) = 2.;
        line(2,1) = -2.;
        lines_true_.emplace_back(line);
    }
}

void Simulation::create_point_line_obs(bool add_nose)
{
    for(auto &Point:points_true_)
    {
        std::vector< std::pair<int, Eigen::Vector3d>> obs;
        std::vector< std::pair<int, Eigen::Vector3d>> obs_true;
        for(int i = 0; i< Twcs_.size(); ++i)
        {
            auto Twc = Twcs_[i];
            Eigen::Matrix4d Tcw = Twc.inverse();
            Eigen::Matrix3d Rcw = Tcw.block(0,0,3,3);
            Eigen::Vector3d tcw = Tcw.block(0,3,3,1);

            Eigen::Vector3d ob;

            ob = Rcw * Point + tcw;

            if(ob(2) < 0) continue;
            ob = ob / ob(2);

            Eigen::Vector3d center(0,0,1);
            Eigen::Vector3d ob_cam = ob;
            ob_cam.normalize();
            double fov0 = std::acos(center.dot(ob_cam)); fov0 = fov0 / M_PI * 180.;
            if(fov0 > 60) continue;

            contain_feature_cams_[i] ++;

            obs_true.emplace_back(i, ob);
            if(add_nose && obs_true.size() > 1)
//            if(add_nose )
            {
                Eigen::Vector3d noise ( pixel_n_(generator_),pixel_n_(generator_), 0 );
                ob += noise;
            }
            obs.emplace_back(i, ob);
        }
        if(obs.empty()) continue;
        point_obs_.push_back(obs);
        point_obs_true_.push_back(obs_true);
    }

    for(auto &Line:lines_true_)
    {
        std::vector< std::pair<int, Eigen::Matrix<double, 3, 2>>> obs;
        std::vector< std::pair<int, Eigen::Matrix<double, 3, 2>>> obs_true;
        for(int i = 0; i< Twcs_.size(); ++i)
        {
            auto Twc = Twcs_[i];
            Eigen::Matrix4d Tcw = Twc.inverse();
            Eigen::Matrix3d Rcw = Tcw.block(0,0,3,3);
            Eigen::Vector3d tcw = Tcw.block(0,3,3,1);

            Eigen::Matrix<double, 3, 2> ob;

            ob.block(0,0,3,1) = Rcw * Line.block(0,0,3,1) + tcw;
            ob.block(0,1,3,1) = Rcw * Line.block(0,1,3,1) + tcw;

            if(ob(2,0) < 0 || ob(2,1) < 0) continue;
            ob.col(0) = ob.col(0) / ob.col(0)(2);
            ob.col(1) = ob.col(1) / ob.col(1)(2);

            Eigen::Vector3d center(0,0,1);
            Eigen::Vector3d ob0_cam = ob.col(0);
            Eigen::Vector3d ob1_cam = ob.col(1);
            ob0_cam.normalize();
            ob1_cam.normalize();
            double fov0 = std::acos(center.dot(ob0_cam)); fov0 = fov0 / M_PI * 180.;
            double fov1 = std::acos(center.dot(ob1_cam)); fov1 = fov1 / M_PI * 180.;
//            std::cout << fov0 << " " << fov1 << " " << fov2 << std::endl;
            if(fov0 > 60) continue;
            if(fov1 > 60) continue;

            contain_feature_cams_[i] ++;

            obs_true.emplace_back(i, ob);
//            if(add_nose && obs_true.size() > 1)
            if(add_nose)
            {
                Eigen::Vector3d noise0 ( pixel_n_(generator_),pixel_n_(generator_), 0 );
                Eigen::Vector3d noise1 ( pixel_n_(generator_),pixel_n_(generator_), 0 );

                ob.col(0) += noise0;
                ob.col(1) += noise1;
            }
            obs.emplace_back(i, ob);
        }
        if(obs.empty()) continue;
        line_obs_.push_back(obs);
        line_obs_true_.push_back(obs_true);
    }
}

void Simulation::tri_point()
{
    for(auto &ob:point_obs_)
    {
        // 1. 三角化点
        {
            Eigen::Vector3d point_cam;
            Eigen::MatrixXd A(ob.size() * 2, 4);
            int index = 0;
//            for(auto &ob_in_cam:ob)

            Eigen::Matrix4d Twc0 = Twcs_[ob[0].first];
            for(int i=1;i<ob.size();++i)
            {
                Eigen::Vector3d ob0 = ob[i].second.col(0);

                Eigen::Matrix4d P = Twcs_[ob[i].first].inverse() * Twc0;
                Eigen::Vector3d f = ob0.normalized();
                A.row(index++) = f[0] * P.row(2) - f[2] * P.row(0);
                A.row(index++) = f[1] * P.row(2) - f[2] * P.row(1);
            }

            Eigen::Vector4d svd_V = Eigen::JacobiSVD<Eigen::MatrixXd>(A, Eigen::ComputeThinV).matrixV().rightCols<1>();
            point_cam = svd_V.head(3) / svd_V(3);
            tri_point_inverse_depth_.push_back( 1/point_cam(2) );
        }
    }
}

void Simulation::add_pose_noise()
{
    int constant_block = Twcs_.size() / constant_pose_num_;
    for(int index = 0;index<Twcs_.size(); ++index)
    {
        if(index % constant_block == 0) continue;

        if(contain_feature_cams_[index] < 10) continue;
        auto &Twc = Twcs_[index];

        Eigen::Matrix4d Tcc = Eigen::Matrix4d::Identity();
        Eigen::Vector3d tcc = Eigen::Vector3d(nt_(generator_), nt_(generator_), nt_(generator_));
        Eigen::Matrix3d Rcc = Eigen::Matrix3d::Identity();
        Rcc = Eigen::AngleAxisd(nq_(generator_), Eigen::Vector3d::UnitZ())
              * Eigen::AngleAxisd(nq_(generator_), Eigen::Vector3d::UnitY())
              * Eigen::AngleAxisd(nq_(generator_), Eigen::Vector3d::UnitX());

//        std::cout << tcc.transpose() << std::endl;
        Tcc.block(0, 0, 3, 3) = Rcc;
        Tcc.block(0, 3, 3, 1) = tcc;

        Twc = Twc*Tcc;
    }
}

void Simulation::optimization_point(OptimizationInfo &info, std::vector<Eigen::Matrix4d> &Twcs_opti)
{
    ProjectionPointInverseDepthFactor::sqrt_info = 240 / 1.5 * Eigen::Matrix2d::Identity();
    ceres::Problem problem;

    double para_Feature[tri_point_inverse_depth_.size()][1];
    for(int index = 0; index<tri_point_inverse_depth_.size(); ++index)
    {
        para_Feature[index][0] = tri_point_inverse_depth_[index];
        problem.AddParameterBlock(para_Feature[index], 1);
    }

    int constant_block;
    constant_block = Twcs_.size() / constant_pose_num_;
    double para_Pose[Twcs_.size()][7];
    for(int i=0; i<Twcs_.size(); ++i)
    {
        Eigen::Matrix3d Rot = Twcs_[i].block(0,0,3,3);
        Eigen::Vector3d Tran = Twcs_[i].block(0,3,3,1);
        Eigen::Quaterniond qua(Rot);

        para_Pose[i][0] = Tran(0);
        para_Pose[i][1] = Tran(1);
        para_Pose[i][2] = Tran(2);
        para_Pose[i][3] = qua.x();
        para_Pose[i][4] = qua.y();
        para_Pose[i][5] = qua.z();
        para_Pose[i][6] = qua.w();

        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i], 7, local_parameterization);
        if(i % constant_block == 0 || contain_feature_cams_[i] < 10)
            problem.SetParameterBlockConstant( para_Pose[i] );
    }


    for(int index_featre = 0; index_featre < point_obs_.size(); ++index_featre)
    {
        auto obs_in_cam = point_obs_[index_featre];
        for(int index_cam = 1;index_cam < obs_in_cam.size(); ++index_cam)
        {
            auto obs = obs_in_cam[index_cam];
            ceres::LossFunction* loss_function = nullptr;
            loss_function = new ceres::CauchyLoss(1.);

            ProjectionPointInverseDepthFactor* cost_function = new ProjectionPointInverseDepthFactor( obs_in_cam[0].second, obs.second );

            problem.AddResidualBlock(cost_function, loss_function, para_Pose[obs_in_cam[0].first], para_Pose[obs.first], para_Feature[index_featre]);

        }
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;  // LEVENBERG_MARQUARDT  DOGLEG
//    options.linear_solver_type = ceres::SPARSE_SCHUR; // SPARSE_NORMAL_CHOLESKY  or SPARSE_SCHUR
//    options.max_num_iterations = 10;
    options.minimizer_progress_to_stdout = true;

    TicToc solver_time;
    ceres::Solver::Summary summary;
    ceres::Solve (options, &problem, & summary);
    info.opti_time = solver_time.toc();
    info.opti_step = summary.iterations.size();
    std::cout << summary.FullReport()<<std::endl;

    for(int index = 0; index<tri_point_inverse_depth_.size(); ++index)
    {
        tri_point_inverse_depth_[index] = para_Feature[index][0];
    }

    Twcs_opti = Twcs_;
    for(int i=0; i<Twcs_.size(); ++i)
    {
        Eigen::Quaterniond qua(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]);
        Eigen::Matrix3d rot = qua.toRotationMatrix();
        Eigen::Vector3d tran(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);
        Twcs_opti[i].block(0,3,3,1) = tran;
        Twcs_opti[i].block(0,0,3,3) = rot;
    }
}


void Simulation::optimization_point_on_plane(OptimizationInfo &info, std::vector<Eigen::Matrix4d> &Twcs_opti)
{
    ProjectionPointInverseDepthFactor::sqrt_info = 240 / 1.5 * Eigen::Matrix2d::Identity();
    PointOnPlaneFactor::sqrt_info = 100.0 / 0.5 ;
    ceres::Problem problem;

    double para_Plane[3];
    para_Plane[0] = 0.;
    para_Plane[1] = 0.;
    para_Plane[2] = -5.;

    Eigen::Matrix<double, 3, 2> tangent_base;
    Eigen::Vector3d norm_n( 1., 0., 0. );
    {

        Eigen::Vector3d b, c;
        norm_n.normalize();
        Eigen::Vector3d tmp(0, 1, 0);
        if(norm_n == tmp)
            tmp << 1, 0, 0;
        b = (tmp - norm_n * (norm_n.transpose() * tmp)).normalized();
        c = norm_n.cross(b);
        tangent_base.block<3, 1>(0, 0) = b;
        tangent_base.block<3, 1>(0, 1) = c;
    }

    double para_Feature[tri_point_inverse_depth_.size()][1];
    for(int index = 0; index<tri_point_inverse_depth_.size(); ++index)
    {
        para_Feature[index][0] = tri_point_inverse_depth_[index];
        problem.AddParameterBlock(para_Feature[index], 1);
    }

    int constant_block;
    constant_block = Twcs_.size() / constant_pose_num_;
    double para_Pose[Twcs_.size()][7];
    for(int i=0; i<Twcs_.size(); ++i)
    {
        Eigen::Matrix3d Rot = Twcs_[i].block(0,0,3,3);
        Eigen::Vector3d Tran = Twcs_[i].block(0,3,3,1);
        Eigen::Quaterniond qua(Rot);

        para_Pose[i][0] = Tran(0);
        para_Pose[i][1] = Tran(1);
        para_Pose[i][2] = Tran(2);
        para_Pose[i][3] = qua.x();
        para_Pose[i][4] = qua.y();
        para_Pose[i][5] = qua.z();
        para_Pose[i][6] = qua.w();

        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i], 7, local_parameterization);
        if(i % constant_block == 0 || contain_feature_cams_[i] < 10)
            problem.SetParameterBlockConstant( para_Pose[i] );
    }

    double para_exPose[7];
    {
        Eigen::Matrix3d Rot = Eigen::Matrix3d::Identity();
        Eigen::Vector3d Tran(0.,0.,0.);
        Eigen::Quaterniond qua(Rot);

        para_exPose[0] = Tran(0);
        para_exPose[1] = Tran(1);
        para_exPose[2] = Tran(2);
        para_exPose[3] = qua.x();
        para_exPose[4] = qua.y();
        para_exPose[5] = qua.z();
        para_exPose[6] = qua.w();

        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_exPose, 7, local_parameterization);
        problem.SetParameterBlockConstant( para_exPose );
    }


    for(int index_featre = 0; index_featre < point_obs_.size(); ++index_featre)
    {
        auto obs_in_cam = point_obs_[index_featre];

        {
            ceres::LossFunction* loss = nullptr;
            loss = new ceres::CauchyLoss(1.);

            auto cost_plane = new PointOnPlaneFactor( norm_n, tangent_base, obs_in_cam[0].second );
            problem.AddResidualBlock(cost_plane, loss, para_Plane, para_Pose[obs_in_cam[0].first], para_exPose, para_Feature[index_featre]);
        }

        for(int index_cam = 1;index_cam < obs_in_cam.size(); ++index_cam)
        {
            auto obs = obs_in_cam[index_cam];
            ceres::LossFunction* loss_function = nullptr;
            loss_function = new ceres::CauchyLoss(1.);

            auto cost_function = new ProjectionPointInverseDepthFactor( obs_in_cam[0].second, obs.second );
            problem.AddResidualBlock(cost_function, loss_function, para_Pose[obs_in_cam[0].first], para_Pose[obs.first], para_Feature[index_featre]);

        }
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;  // LEVENBERG_MARQUARDT  DOGLEG
//    options.linear_solver_type = ceres::SPARSE_SCHUR; // SPARSE_NORMAL_CHOLESKY  or SPARSE_SCHUR
//    options.max_num_iterations = 10;
    options.minimizer_progress_to_stdout = true;

    TicToc solver_time;
    ceres::Solver::Summary summary;
    ceres::Solve (options, &problem, & summary);
    info.opti_time = solver_time.toc();
    info.opti_step = summary.iterations.size();
    std::cout << summary.FullReport()<<std::endl;

    for(int index = 0; index<tri_point_inverse_depth_.size(); ++index)
    {
        tri_point_inverse_depth_[index] = para_Feature[index][0];
    }

    Twcs_opti = Twcs_;
    for(int i=0; i<Twcs_.size(); ++i)
    {
        Eigen::Quaterniond qua(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]);
        Eigen::Matrix3d rot = qua.toRotationMatrix();
        Eigen::Vector3d tran(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);
        Twcs_opti[i].block(0,3,3,1) = tran;
        Twcs_opti[i].block(0,0,3,3) = rot;
    }
}

void Simulation::optimization_point_on_plane_projection(OptimizationInfo &info, std::vector<Eigen::Matrix4d> &Twcs_opti)
{
    PointOnPlaneProjectionFactor::sqrt_info = 240 / 1.5 * Eigen::Matrix2d::Identity();
    ceres::Problem problem;

    double para_Plane[3];
    para_Plane[0] = 0.;
    para_Plane[1] = 0.;
    para_Plane[2] = -5.;

    Eigen::Matrix<double, 3, 2> tangent_base;
    Eigen::Vector3d norm_n( 1., 0., 0. );
    {

        Eigen::Vector3d b, c;
        norm_n.normalize();
        Eigen::Vector3d tmp(0, 1, 0);
        if(norm_n == tmp)
            tmp << 1, 0, 0;
        b = (tmp - norm_n * (norm_n.transpose() * tmp)).normalized();
        c = norm_n.cross(b);
        tangent_base.block<3, 1>(0, 0) = b;
        tangent_base.block<3, 1>(0, 1) = c;
    }

    double para_Feature[tri_point_inverse_depth_.size()][1];
    for(int index = 0; index<tri_point_inverse_depth_.size(); ++index)
    {
        para_Feature[index][0] = tri_point_inverse_depth_[index];
        problem.AddParameterBlock(para_Feature[index], 1);
    }

    int constant_block;
    constant_block = Twcs_.size() / constant_pose_num_;
    double para_Pose[Twcs_.size()][7];
    for(int i=0; i<Twcs_.size(); ++i)
    {
        Eigen::Matrix3d Rot = Twcs_[i].block(0,0,3,3);
        Eigen::Vector3d Tran = Twcs_[i].block(0,3,3,1);
        Eigen::Quaterniond qua(Rot);

        para_Pose[i][0] = Tran(0);
        para_Pose[i][1] = Tran(1);
        para_Pose[i][2] = Tran(2);
        para_Pose[i][3] = qua.x();
        para_Pose[i][4] = qua.y();
        para_Pose[i][5] = qua.z();
        para_Pose[i][6] = qua.w();

        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i], 7, local_parameterization);
        if(i % constant_block == 0 || contain_feature_cams_[i] < 10)
            problem.SetParameterBlockConstant( para_Pose[i] );
    }

    double para_exPose[7];
    {
        Eigen::Matrix3d Rot = Eigen::Matrix3d::Identity();
        Eigen::Vector3d Tran(0.,0.,0.);
        Eigen::Quaterniond qua(Rot);

        para_exPose[0] = Tran(0);
        para_exPose[1] = Tran(1);
        para_exPose[2] = Tran(2);
        para_exPose[3] = qua.x();
        para_exPose[4] = qua.y();
        para_exPose[5] = qua.z();
        para_exPose[6] = qua.w();

        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_exPose, 7, local_parameterization);
        problem.SetParameterBlockConstant( para_exPose );
    }


    for(int index_featre = 0; index_featre < point_obs_.size(); ++index_featre)
    {
        auto obs_in_cam = point_obs_[index_featre];

        {
            ceres::LossFunction* loss = nullptr;
            loss = new ceres::CauchyLoss(1.);

            auto cost_plane = new PointOnPlaneFactor( norm_n, tangent_base, obs_in_cam[0].second );
            problem.AddResidualBlock(cost_plane, loss, para_Plane, para_Pose[obs_in_cam[0].first], para_exPose, para_Feature[index_featre]);
        }

        for(int index_cam = 1;index_cam < obs_in_cam.size(); ++index_cam)
        {
            auto obs = obs_in_cam[index_cam];
            ceres::LossFunction* loss_function = nullptr;
            loss_function = new ceres::CauchyLoss(1.);

            auto cost_function = new PointOnPlaneProjectionFactor( norm_n, tangent_base, obs_in_cam[0].second, obs.second );
            problem.AddResidualBlock(cost_function, loss_function,
                    para_Plane,
                    para_Pose[obs_in_cam[0].first],
                    para_Pose[obs.first]);
        }
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;  // LEVENBERG_MARQUARDT  DOGLEG
//    options.linear_solver_type = ceres::SPARSE_SCHUR; // SPARSE_NORMAL_CHOLESKY  or SPARSE_SCHUR
//    options.max_num_iterations = 10;
    options.minimizer_progress_to_stdout = true;

    TicToc solver_time;
    ceres::Solver::Summary summary;
    ceres::Solve (options, &problem, & summary);
    info.opti_time = solver_time.toc();
    info.opti_step = summary.iterations.size();
    std::cout << summary.FullReport()<<std::endl;

    for(int index = 0; index<tri_point_inverse_depth_.size(); ++index)
    {
        tri_point_inverse_depth_[index] = para_Feature[index][0];
    }

    Twcs_opti = Twcs_;
    for(int i=0; i<Twcs_.size(); ++i)
    {
        Eigen::Quaterniond qua(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]);
        Eigen::Matrix3d rot = qua.toRotationMatrix();
        Eigen::Vector3d tran(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);
        Twcs_opti[i].block(0,3,3,1) = tran;
        Twcs_opti[i].block(0,0,3,3) = rot;
    }
}

void Simulation::DrawAllCamera(std::vector<pangolin::OpenGlMatrix>& Ms)
{
    for(auto &M:Ms)
    {
        DrawSigCam(M);
    }
}

void Simulation::DrawSigCam( pangolin::OpenGlMatrix& M )
{
    //相机模型大小：宽度占总宽度比例为0.08
    const float &w = 0.3;
    const float h = w*0.638;
    const float z = w*0.6;

    //百度搜索：glPushMatrix 百度百科
    glPushMatrix();

    //将4*4的矩阵Twc.m右乘一个当前矩阵
    //（由于使用了glPushMatrix函数，因此当前帧矩阵为世界坐标系下的单位矩阵）
    //因为OpenGL中的矩阵为列优先存储，因此实际为Tcw，即相机在世界坐标下的位姿
#ifdef HAVE_GLES
    glMultMatrixf(M.m);
#else
    glMultMatrixd(M.m);
#endif

//    设置绘制图形时线的宽度
    glLineWidth(1.0);
    //设置当前颜色为绿色(相机图标显示为绿色)
    glColor3f(0.0f, 1.0f, 0.0f);
    //用线将下面的顶点两两相连
    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);
    glVertex3f(w, h, z);
    glVertex3f(0, 0, 0);
    glVertex3f(w, -h, z);
    glVertex3f(0, 0, 0);
    glVertex3f(-w, -h, z);
    glVertex3f(0, 0, 0);
    glVertex3f(-w, h, z);

    glVertex3f(w, h, z);
    glVertex3f(w, -h, z);

    glVertex3f(-w, h, z);
    glVertex3f(-w, -h, z);

    glVertex3f(-w, h, z);
    glVertex3f(w, h, z);

    glVertex3f(-w, -h, z);
    glVertex3f(w, -h, z);
    glEnd();


    //双缓存交换缓存以显示图像
//    glutSwapBuffers();

    glPopMatrix();
}

void Simulation::CallCameraTrue(std::vector<pangolin::OpenGlMatrix>& Ms)
{
    {
        Ms.clear();
        for(auto & Twc : Twcs_true_)
        {
            pangolin::OpenGlMatrix M;

            Eigen::Matrix3d Rwc = Twc.block(0,0,3,3);
            Eigen::Vector3d twc = Twc.block(0,3,3,1);

            M.m[0] = Rwc(0,0);
            M.m[1] = Rwc(1,0);
            M.m[2] = Rwc(2,0);
            M.m[3]  = 0.0;

            M.m[4] = Rwc(0,1);
            M.m[5] = Rwc(1,1);
            M.m[6] = Rwc(2,1);
            M.m[7]  = 0.0;

            M.m[8] = Rwc(0,2);
            M.m[9] = Rwc(1,2);
            M.m[10] = Rwc(2,2);
            M.m[11]  = 0.0;

            M.m[12] = twc(0);
            M.m[13] = twc(1);
            M.m[14] = twc(2);
            M.m[15]  = 1.0;

            Ms.push_back(M);
        }
    }
}

void Simulation::DrawTrueLine()
{
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glColor3f(1.0f,0.0f,0.0f);
    for(const auto & Line:lines_true_)
    {
        Eigen::Vector3d line0 = Line.block(0,0,3,1);
        Eigen::Vector3d line1 = Line.block(0,1,3,1);
        glVertex3f(line0(0),line0(1),line0(2));
        glVertex3f(line1(0),line1(1),line1(2));
    }
    glEnd();
}

void Simulation::DrawTriPoint()
{
    for(int i=0;i<point_obs_.size();++i)
    {
        const auto& obs = point_obs_[i];
        const auto& ob = obs.begin();
        const int cam_id = ob->first;
        Eigen::Vector3d point_cam = ob->second;
        point_cam /= tri_point_inverse_depth_[i];

        Eigen::Matrix3d Rwc = Twcs_[cam_id].block(0,0,3,3);
        Eigen::Vector3d twc = Twcs_[cam_id].block(0,3,3,1);

        Eigen::Vector3d point_w = Rwc*point_cam+twc;


        glPointSize(5);
        glBegin(GL_POINTS);
        glColor3f(1.0,1.0,0.0);
        glVertex3d( point_w[0], point_w[1], point_w[2]);
        glEnd();
    }
}

void Simulation::show()
{
    pangolin::CreateWindowAndBind("Simulation",1024,768);

    // 3D Mouse handler requires depth testing to be enabled
    // 启动深度测试，OpenGL只绘制最前面的一层，绘制时检查当前像素前面是否有别的像素，如果别的像素挡住了它，那它就不会绘制
    glEnable(GL_DEPTH_TEST);

    // Issue specific OpenGl we might need
    // 在OpenGL中使用颜色混合
    glEnable(GL_BLEND);
    // 选择混合选项
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // 新建按钮和选择框，第一个参数为按钮的名字，第二个为默认状态，第三个为是否有选择框
    pangolin::CreatePanel("menu").SetBounds(0.0,1.0,0.0,pangolin::Attach::Pix(175));
    pangolin::Var<bool> menuWhite("menu.Show White",false,true);
    pangolin::Var<bool> menuShowCamera("menu.Show Camera",true,true);
    pangolin::Var<bool> menuShowCameraPoint("menu.Show CameraPoint",true,true);
//    pangolin::Var<bool> menuShowCamera("menu.Show Camera",true,true);
//    pangolin::Var<bool> menuShowCamera("menu.Show Camera",true,true);
//    pangolin::Var<bool> menuShowCamera("menu.Show Camera",true,true);
    pangolin::Var<bool> menuShowPoint("menu.Show Point",false,true);
    pangolin::Var<bool> menuShowPointOpti("menu.Show PointOpti",false,true);
    pangolin::Var<bool> menuShowLines("menu.Show Lines",true,true);
    pangolin::Var<bool> menuShowOptiLines("menu.Show Opti Lines",true,true);


    // Define Camera Render Object (for view / scene browsing)
    // 定义相机投影模型：ProjectionMatrix(w, h, fu, fv, u0, v0, zNear, zFar)
    // 定义观测方位向量：观测点位置：(mViewpointX mViewpointY mViewpointZ)
    //                观测目标位置：(0, 0, 0)
    //                观测的方位向量：(0.0,-1.0, 0.0)
    pangolin::OpenGlRenderState s_cam(
            pangolin::ProjectionMatrix(1024,768,500,500,512,389,0.1,1000),
            pangolin::ModelViewLookAt(0,-0.7,-1.8, 0,0,0,0.0,-1.0, 0.0)
    );

    // Add named OpenGL viewport to window and provide 3D Handler
    // 定义显示面板大小，orbslam中有左右两个面板，昨天显示一些按钮，右边显示图形
    // 前两个参数（0.0, 1.0）表明宽度和面板纵向宽度和窗口大小相同
    // 中间两个参数（pangolin::Attach::Pix(175), 1.0）表明右边所有部分用于显示图形
    // 最后一个参数（-1024.0f/768.0f）为显示长宽比
    pangolin::View& d_cam = pangolin::CreateDisplay()
            .SetBounds(0.0, 1.0, pangolin::Attach::Pix(175), 1.0, -1024.0f/768.0f)
            .SetHandler(new pangolin::Handler3D(s_cam));

    pangolin::OpenGlMatrix Twc;
    Twc.SetIdentity();

    std::vector<pangolin::OpenGlMatrix> Twcs;
    Twcs.push_back(Twc);

//    cv::namedWindow("Current Frame");

    bool bFollow = true;
    bool bLocalizationMode = false;

    int Single_rpr_id = 0;
    int Single_line_id = 0;
    bool line_id_change = false;

    int cut_i = 0;

    while( !pangolin::ShouldQuit() )
    {
//        while(!EstOK)
//            cv::waitKey(100);
        // 清除缓冲区中的当前可写的颜色缓冲 和 深度缓冲
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


        d_cam.Activate(s_cam);
        // 步骤3：绘制地图和图像
        // 设置为白色，glClearColor(red, green, blue, alpha），数值范围(0, 1)
        if(menuWhite)
            glClearColor(1.0f,1.0f,1.0f,1.0f);
        else
            glClearColor(0.0f,0.0f,0.0f,0.0f);


        const double lens = 0.5;
        glLineWidth(2.0);
        glBegin(GL_LINES);
        glColor3f(1.0f,0.0f,0.0f);
        glVertex3f(0,0,0);
        glVertex3f(lens,0,0);


        glColor3f(0.0f,1.0f,0.0f);
        glVertex3f(0,0,0);
        glVertex3f(0,lens,0);

        glColor3f(0.0f,0.0f,1.0f);
        glVertex3f(0,0,0);
        glVertex3f(0,0,lens);

        glEnd();


        if(menuShowPoint)
        {
            glPointSize(5);
            glBegin(GL_POINTS);
            glColor3f(1.0,0.0,0.0);
            for(auto p : points_true_)
            {
                glVertex3d( p[0], p[1], p[2]);
            }
            glEnd();
        }

        if(menuShowPointOpti)
            DrawTriPoint();
//        {
//            glPointSize(5);
//            glBegin(GL_POINTS);
//            glColor3f(1.0,0.0,0.0);
//            for(auto p : points_tri_opt_)
//            {
//                glVertex3d( p[0], p[1], p[2]);
//            }
//            glEnd();
//        }

        if(menuShowCamera)
        {
            std::vector<pangolin::OpenGlMatrix> MsTrue;
            CallCameraTrue(MsTrue);
            DrawAllCamera(MsTrue);
        }
//        if(menuShowCameraOpti)
//        {
//            std::vector<pangolin::OpenGlMatrix> Ms_point;
//            CallCameraPoint(Ms_point);
//            DrawAllCamera(Ms_point);
//        }
//
        if(menuShowLines)
            DrawTrueLine();
//
//        if(menuShowOptiLines)
//            DrawOptiiLines();



        pangolin::FinishFrame();
    }
//    std::cerr << "MFK" << std::endl;
    pangolin::DestroyWindow("Simulation");

    pangolin::QuitAll();
}