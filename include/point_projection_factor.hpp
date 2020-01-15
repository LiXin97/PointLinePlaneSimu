//
// Created by lixin04 on 19年12月6日.
//
#pragma once
#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <ceres/rotation.h>

#include "utility.hpp"
#include "tictoc.hpp"

class ProjectionPointInverseDepthFactor : public ceres::SizedCostFunction<2, 7, 7, 1>
{
public:
    ProjectionPointInverseDepthFactor(Eigen::Vector3d& ob0, Eigen::Vector3d& obi);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);

    static Eigen::Matrix2d sqrt_info;
    static double sum_t;

    Eigen::Vector3d point_obi_;
    Eigen::Vector3d point_ob0_;
};

//class PointOnPlaneFactor : public ceres::SizedCostFunction<2, 7, 7, 1>
//{
//public:
//    PointOnPlaneFactor(Eigen::Vector3d& ob0);
//    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
//    void check(double **parameters);
//
//    static double sqrt_info;
//    static double sum_t;
//
//    Eigen::Vector3d point_obi_;
//    Eigen::Vector3d point_ob0_;
//};

class PointOnPlaneFactor : public ceres::SizedCostFunction<1, 3, 7, 7, 1>
{
public:
    PointOnPlaneFactor(const Eigen::Vector3d& norm, const Eigen::Matrix<double, 3, 2> &tange, const Eigen::Vector3d &_pts_i);

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);
    static double sqrt_info;
    Eigen::Vector3d norm_n;
    Eigen::Vector3d pts_i;
    Eigen::Matrix<double, 3, 2> tangent_base;
};

class PointOnPlaneProjectionFactor : public ceres::SizedCostFunction<2, 3, 7, 7>
{
public:
    PointOnPlaneProjectionFactor(
            const Eigen::Vector3d& norm,
            const Eigen::Matrix<double, 3, 2> &tange,
            const Eigen::Vector3d &_ob0,
            const Eigen::Vector3d& _obi);

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);
    static Eigen::Matrix2d sqrt_info;
    Eigen::Vector3d norm_n;
    Eigen::Vector3d ob0;
    Eigen::Vector3d obi;
    Eigen::Matrix<double, 3, 2> tangent_base;
};