//
// Created by lixin04 on 19年12月6日.
//

#pragma once
#include <Eigen/Geometry>
#include <iostream>
#include <vector>

class OptimizationInfo
{
public:
    double opti_time = 0;
    int opti_step    = 0;
    double error_jaco_time = 0;

    double residual_evaluation_time = 0;
    double jacobian_evaluation_time = 0;
    double linear_solver_time = 0;
    double minimizer_time = 0;
    double total_time = 0;

    void operator += (const OptimizationInfo&T)
    {
        opti_time += T.opti_time;
        opti_step += T.opti_step;
        error_jaco_time += T.error_jaco_time;
        residual_evaluation_time += T.residual_evaluation_time;
        jacobian_evaluation_time += T.jacobian_evaluation_time;
        linear_solver_time += T.linear_solver_time;
        minimizer_time += T.minimizer_time;
        total_time += T.total_time;
    }

    void operator /= (const double T)
    {
        opti_time /= T;
        opti_step /= T;
        error_jaco_time /= T;
        residual_evaluation_time /= T;
        jacobian_evaluation_time /= T;
        linear_solver_time /= T;
        minimizer_time /= T;
        total_time /= T;
    }

    void coutinfo()
    {
        std::cout << "----------------------------------------- "  << std::endl;
        std::cout << "averg all time: " << opti_time << std::endl;
        std::cout << "averg error jaco time: " << error_jaco_time << std::endl;
        std::cout << "averg step: " << opti_step << std::endl;
        std::cout << "averg residual_evaluation_time: " << residual_evaluation_time << std::endl;
        std::cout << "averg jacobian_evaluation_time: " << jacobian_evaluation_time << std::endl;
        std::cout << "averg linear_solver_time: " << linear_solver_time << std::endl;
        std::cout << "averg minimizer_time: " << minimizer_time << std::endl;
        std::cout << "averg total_time: " << total_time << std::endl;
        std::cout << "----------------------------------------- "  << std::endl;
    }
};


class AbsoluteError{
public:

    void operator += (const AbsoluteError&T)
    {
        rmse += T.rmse;
        mean += T.mean;
        max += T.max;
        min += T.min;
    }

    void operator /= (const double T)
    {
        rmse /= T;
        mean /= T;
        max /= T;
        min /= T;
    }

    AbsoluteError(double _rmse, double _mean, double _max, double _min):rmse(_rmse), mean(_mean), max(_max), min(_min){}
    AbsoluteError():rmse(0), mean(0), max(0), min(0){}
    AbsoluteError(const AbsoluteError& T) = default;

//private:
    double rmse = 0;
    double mean = 0;
    double max  = 0;
    double min  = 0;
};


class APE{
public:
    void operator += (const APE&T)
    {
        tranlational_error += T.tranlational_error;
        rotational_error += T.rotational_error;
    }

    void operator /= (const double T)
    {
        tranlational_error /= T;
        rotational_error /= T;
    }

    void coutrpe()
    {
        std::cout << "----------------------------------------- "  << std::endl;
        std::cout << "translational_error.rmse: " << tranlational_error.rmse << std::endl;
        std::cout << "translational_error.mean: " << tranlational_error.mean << std::endl;
        std::cout << "translational_error.max: " << tranlational_error.max << std::endl;
        std::cout << "translational_error.min: " << tranlational_error.min << std::endl;


        std::cout << "rotational_error.rmse: " << rotational_error.rmse <<std::endl;
        std::cout << "rotational_error.mean: " << rotational_error.mean << std::endl;
        std::cout << "rotational_error.max: " << rotational_error.max << std::endl;
        std::cout << "rotational_error.min: " << rotational_error.min << std::endl;
        std::cout << "----------------------------------------- "  << std::endl;
    }

    APE() = default;
    APE(AbsoluteError &_tranlational_error, AbsoluteError &_rotational_error):tranlational_error(_tranlational_error), rotational_error(_rotational_error){}

//private:
    AbsoluteError tranlational_error;
    AbsoluteError rotational_error;

};

inline double computeAngle(const Eigen::Matrix3d& T) {
    //    # an invitation to 3-d vision, p 27
    double rotation_trace = T.trace();
    return acos(std::min(double(1), std::max(double(-1), (rotation_trace - 1)/2)));
}


APE evaluateTrajectory(std::vector<Eigen::Matrix4d>& Twcs_gt, std::vector<Eigen::Matrix4d>& Twcs_esti);
