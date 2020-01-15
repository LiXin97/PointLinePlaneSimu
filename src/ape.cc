//
// Created by lixin04 on 19年12月6日.
//

#include "ape.hpp"


APE evaluateTrajectory(std::vector<Eigen::Matrix4d>& Twcs_gt, std::vector<Eigen::Matrix4d>& Twcs_esti)
{
    std::vector<Eigen::Matrix4d> diffTs(Twcs_esti.size());
    for(int i=0;i<Twcs_gt.size();++i)
    {
        diffTs[i] = Twcs_esti[i] * Twcs_gt[i].inverse();
    }

    double sq_sum = 0.f, sum = 0.f, rot_sq_sum = 0.f, rot_sum = 0.f;
    double min = 1e3f, max = -1e3f, rot_min = 1e3f, rot_max = -1e3f;
    for(auto & diffT : diffTs) {
        const double trans = diffT.block(0,3,3,1).norm();
        const double& rot = computeAngle(diffT.block(0,0,3,3));
        sq_sum += trans*trans;
        rot_sq_sum += rot*rot;
        sum += trans;
        rot_sum += rot;
        if(trans > max)
            max = trans;
        if(trans < min)
            min = trans;
        if(rot > rot_max)
            rot_max = rot;
        if(rot < rot_min)
            rot_min = rot;
    }

//    std::cout << max << " - " << rot_max * 180.0 / M_PI << std::endl;

    double rmse_tran = sqrt(sq_sum/diffTs.size());
    double mean_tran = sum/diffTs.size();
    double max_tran = max;
    double min_tran = min;
    AbsoluteError tran(rmse_tran, mean_tran, max_tran, min_tran);

    double rmse_rot = sqrt(rot_sq_sum/diffTs.size()) * 180.0 / M_PI;
    double mean_rot = rot_sum/diffTs.size() * 180.0 / M_PI;
    double max_rot = rot_max * 180.0 / M_PI;
    double min_rot = rot_min * 180.0 / M_PI;
    AbsoluteError rotation(rmse_rot, mean_rot, max_rot, min_rot);

    APE ape(tran, rotation);
//    rpe.coutrpe();
    return ape;
}