//
// Created by lixin04 on 20年1月15日.
//

#include "simulation.hpp"

int main()
{

    APE test0, test1, test2, test3, test4, test5, test6;
    OptimizationInfo info0, info1, info2, info3, info4, info5, info6;

    int loop_num = 10;
    for(int loop = 0; loop < loop_num; ++loop)
    {
        Simulation simu;

        simu.create_keyframe();
        simu.generate_point_line();
        simu.create_point_line_obs();
        simu.add_pose_noise();
        simu.tri_point();

        OptimizationInfo tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;

        std::vector<Eigen::Matrix4d> Twcs_true = simu.getTwcs_true();

//        return 0;

        std::vector<Eigen::Matrix4d> Twcs_point;
        simu.optimization_point(tmp0, Twcs_point);
        test0 += evaluateTrajectory(Twcs_true, Twcs_point);
        info0 += tmp0;

        std::vector<Eigen::Matrix4d> Twcs_point_on_plane;
        simu.optimization_point_on_plane(tmp1, Twcs_point_on_plane);
        test1 += evaluateTrajectory(Twcs_true, Twcs_point_on_plane);
        info1+= tmp1;

        std::vector<Eigen::Matrix4d> Twcs_point_on_plane_projection;
        simu.optimization_point_on_plane_projection(tmp2, Twcs_point_on_plane_projection);
        test2 += evaluateTrajectory(Twcs_true, Twcs_point_on_plane_projection);
        info2 += tmp2;

//        return 0;

//        simu.show();

    }

    test0 /= loop_num;
    test1 /= loop_num;
    test2 /= loop_num;
    test3 /= loop_num;
    test4 /= loop_num;
    test5 /= loop_num;
    test6 /= loop_num;

    info0 /= loop_num;
    info1 /= loop_num;
    info2 /= loop_num;
    info3 /= loop_num;
    info4 /= loop_num;
    info5 /= loop_num;
    info6 /= loop_num;
    std::cout << "point: " << std::endl;
    test0.coutrpe();
    info0.coutinfo();
    std::cout << "point on plan: " << std::endl;
    test1.coutrpe();
    info1.coutinfo();
    std::cout << "point on plan projection: " << std::endl;
    test2.coutrpe();
    info2.coutinfo();
//    std::cout << "line lixin4: " << std::endl;
//    test3.coutrpe();
//    info3.coutinfo();
//    std::cout << "line in cam: " << std::endl;
//    test4.coutrpe();
//    info4.coutinfo();
//    std::cout << "line in cam qua: " << std::endl;
//    test5.coutrpe();
//    info5.coutinfo();
//    std::cout << "line zhaoliang qua: " << std::endl;
//    test6.coutrpe();
//    info6.coutinfo();

    return 0;
}