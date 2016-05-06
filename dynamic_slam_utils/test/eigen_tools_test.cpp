// Bring in my package's API, which is what I'm testing
#include <dynamic_slam_utils/eigen_tools.h>
// Bring in gtest
#include <gtest/gtest.h>

using namespace eigt;

double EPSILON = 0.001;

TEST(EigenTools, createTransFromPosesSimple)
{
  try {
    pose2d_t first_pose;
    first_pose << 0, 0, 0;
    pose2d_t second_pose;
    second_pose << 10, 10, 0;
    transform2d_t res_trans;
    res_trans.matrix() << 1, 0, 10, 0, 1, 10, 0, 0, 1;
    Eigen::Matrix3d calc = transBtwPoses(first_pose, second_pose).matrix();
    EXPECT_TRUE((res_trans.matrix() - calc).cwiseAbs().sum() < EPSILON);
    // ADD_FAILURE() <<res_trans.matrix()<< "\n\n"<< calc;
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, createTransFromPosesAdvance)
{
  try {
    pose2d_t first_pose;
    first_pose << 10, 10, 0;
    pose2d_t second_pose;
    second_pose << 20, 20, 0;
    transform2d_t res_trans;
    res_trans.matrix() << 1, 0, 10, 0, 1, 10, 0, 0, 1;
    Eigen::Matrix3d calc = transBtwPoses(first_pose, second_pose).matrix();
    EXPECT_TRUE((res_trans.matrix() - calc).cwiseAbs().sum() < EPSILON);
    // ADD_FAILURE() <<res_trans.matrix()<< "\n\n"<< calc;
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, createTransFromPosesAngle)
{
  try {
    pose2d_t first_pose;
    first_pose << 0, 0, 0;
    pose2d_t second_pose;
    second_pose << 0, 0, M_PI;
    transform2d_t res_trans;
    res_trans.matrix() << -1, 0, 0, 0, -1, 0, 0, 0, 1;
    Eigen::Matrix3d calc = transBtwPoses(first_pose, second_pose).matrix();
    EXPECT_TRUE((res_trans.matrix() - calc).cwiseAbs().sum() < EPSILON);
    // ADD_FAILURE() <<res_trans.matrix()<< "\n\n"<< calc;
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, CycleTest)
{
  try {
    pose2d_t first_pose;
    first_pose << 10, 20, 0;
    pose2d_t second_pose;
    second_pose << -10, 30, -M_PI;
    transform2d_t calc_trans = transBtwPoses(first_pose, second_pose);
    pose2d_t calc_sec_pose = transformPose(first_pose,calc_trans);
    std::cout<<getAngle(calc_trans)<<"\n";
    std::cout<<getPoseFromTransform(calc_trans)<<"\n";
    std::cout<<"target pose: " <<second_pose.transpose()<<" calc target pose: " << calc_sec_pose.transpose()<< "\n";
    EXPECT_TRUE((calc_sec_pose - second_pose).cwiseAbs().sum() < EPSILON);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, CycleTest2)
{
  try {
    pose2d_t first_pose;
    first_pose << -10, 5, 0;
    pose2d_t second_pose;
    second_pose << 10, 6, 2;
    transform2d_t calc_trans = transBtwPoses(first_pose, second_pose);
    pose2d_t calc_sec_pose = transformPose(first_pose,calc_trans);
    std::cout<<getAngle(calc_trans)<<"\n";
    std::cout<<getPoseFromTransform(calc_trans)<<"\n";
    std::cout<<"target pose: " <<second_pose.transpose()<<" calc target pose: " << calc_sec_pose.transpose()<< "\n";
    EXPECT_TRUE((calc_sec_pose - second_pose).cwiseAbs().sum() < EPSILON);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, CycleTest3)
{
  try {
    pose2d_t first_pose;
    first_pose << -10, 5, 1;
    pose2d_t second_pose;
    second_pose << 10, 6, M_PI;
    transform2d_t calc_trans = transBtwPoses(first_pose, second_pose);
    pose2d_t calc_sec_pose = transformPose(first_pose,calc_trans);
    std::cout<<getAngle(calc_trans)<<"\n";
    std::cout<<getPoseFromTransform(calc_trans)<<"\n";
    std::cout<<"target pose: " <<second_pose.transpose()<<" calc target pose: " << calc_sec_pose.transpose()<< "\n";
    EXPECT_TRUE((calc_sec_pose - second_pose).cwiseAbs().sum() < EPSILON);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, complexInterfaceTest)
{
  try {
    pose2d_t source_pose;
    pose2d_t target_pose;
    transform2d_t trans_btw_poses;
    pose2d_t trans_btw_poses_short;
    double angle_diff;
    double displacement;

    source_pose << -10, 5, 1.2;
    target_pose << 10, 6, 2;
    angle_diff = 0.8;
    displacement = std::sqrt(std::pow(source_pose(0) - target_pose(0), 2) +
                             std::pow(source_pose(1) - target_pose(1), 2));
    trans_btw_poses_short << 6.31512, 19.0031, 0.8;
    transform2d_t calc_trans = transBtwPoses(source_pose, target_pose);
    double calc_angle = getAngle(calc_trans);
    double calc_angle_diff = getAngleDiffrence(source_pose,target_pose);
    double calc_displacement = getDisplacement(calc_trans);
    pose2d_t calc_pose = getPoseFromTransform(calc_trans);
    pose2d_t calc_trans_pose = transformPose(source_pose, calc_trans);
    std::cout<< calc_trans.matrix() << "\n";
    std::cout<<"angle: " <<calc_angle << "\n";
    std::cout<<"angle diffrence: " <<calc_angle_diff << "\n";
    std::cout<<"displacement: " <<calc_displacement << "\n";
    std::cout<<"pose from transform: " <<calc_pose << "\n";
     std::cout<<"transformed source pose: " <<calc_trans_pose.transpose() <<
              "target pose reference: " << target_pose.transpose()<<"\n";

    EXPECT_TRUE((trans_btw_poses_short - calc_pose)
                    .cwiseAbs()
                    .sum() < EPSILON);

    EXPECT_TRUE(std::abs(angle_diff - calc_angle) < EPSILON);
    EXPECT_TRUE(std::abs(displacement - calc_displacement) < EPSILON);
    EXPECT_TRUE((target_pose - calc_trans_pose).cwiseAbs().sum() <
        EPSILON);

  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, normalizeAngle){
  EXPECT_LT(normalizeAngle(M_PI) - M_PI,EPSILON);
  EXPECT_LT(normalizeAngle(3*M_PI) - M_PI,EPSILON);
  EXPECT_LT(normalizeAngle(2*M_PI) ,EPSILON);
  EXPECT_LT(normalizeAngle(3*(M_PI/2)) - M_PI + (M_PI/2) ,EPSILON);

  EXPECT_LT(normalizeAngle(-M_PI) + M_PI,EPSILON);
  EXPECT_LT(normalizeAngle(-2*M_PI),EPSILON);
  EXPECT_LT(normalizeAngle(-3*M_PI) + M_PI,EPSILON);
}


TEST(EigenTools, getTransFromPoseTest)
{
  try {
    pose2d_t first_pose;
    first_pose << -10, 5, 0;
    pose2d_t second_pose;
    second_pose << 10, 6, 2;
    transform2d_t calc_trans = transBtwPoses(first_pose, second_pose);
    pose2d_t calc_pose = getPoseFromTransform(calc_trans);
    EXPECT_LT((calc_trans.matrix() - getTransFromPose(calc_pose).matrix()).cwiseAbs().sum(),EPSILON);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
