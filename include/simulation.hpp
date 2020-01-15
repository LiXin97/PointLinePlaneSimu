//
// Created by lixin04 on 20年1月15日.
//
#pragma once

#include <Eigen/Geometry>
#include <pangolin/pangolin.h>
#include <math.h>
#include <vector>
#include <random>
#include "ape.hpp"

#include "point_projection_factor.hpp"
#include "pose_local_parameterization.hpp"

class Simulation
{
public:
    explicit Simulation(double max_nt = 0.1, double max_nq = 1.*M_PI/180., double max_pixel_n = 1./240);

    void generate_point_line();
    void create_keyframe(int frame_num = 150, double radius = 1.5, int wave_num = 10, double wave_high = .75);
    void create_point_line_obs(bool add_nose = true);

    void tri_point();
    void tri_line();

    void optimization_point(OptimizationInfo& info, std::vector<Eigen::Matrix4d>& Twcs_opti);

    void optimization_point_on_plane(OptimizationInfo& info, std::vector<Eigen::Matrix4d>& Twcs_opti);

    void optimization_point_on_plane_projection(OptimizationInfo& info, std::vector<Eigen::Matrix4d>& Twcs_opti);

    void add_pose_noise();

    std::vector<Eigen::Matrix4d> getTwcs_true(){return Twcs_true_;}

    void CallCameraTrue(std::vector<pangolin::OpenGlMatrix>& Ms);
    void DrawAllCamera(std::vector<pangolin::OpenGlMatrix>& Ms);
    void DrawSigCam(pangolin::OpenGlMatrix& M);
    void DrawTrueLine();
    void DrawTriPoint();
    void show();

private:
    std::vector<Eigen::Matrix4d> Twcs_;
    std::vector<Eigen::Matrix4d> Twcs_true_;

    std::vector<Eigen::Vector3d> points_true_;
    std::vector<Eigen::Matrix<double, 3, 2>> lines_true_;

    std::vector<std::vector<std::pair< int,Eigen::Vector3d>>> point_obs_;
    std::vector<std::vector<std::pair< int,Eigen::Vector3d>>> point_obs_true_;

    std::map<int, int> contain_feature_cams_;

    std::vector<std::vector<std::pair<int, Eigen::Matrix<double, 3, 2>>>> line_obs_;
    std::vector<std::vector<std::pair<int, Eigen::Matrix<double, 3, 2>>>> line_obs_true_;

    std::vector<double> tri_point_inverse_depth_;

    std::default_random_engine generator_;
    std::normal_distribution<double> nq_;
    std::normal_distribution<double> nt_;
    std::normal_distribution<double> pixel_n_;

    double constant_pose_num_;
};