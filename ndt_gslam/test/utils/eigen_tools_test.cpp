// Bring in my package's API, which is what I'm testing
#include <ndt_gslam/utils/eigen_tools.h>
// Bring in gtest
#include <gtest/gtest.h>

using namespace eigt;

double EPSILON = 0.001;

TEST(EigenTools, createTransFromPosesSimple)
{
  try {
    pose2d_t<double> first_pose;
    first_pose << 0, 0, 0;
    pose2d_t<double> second_pose;
    second_pose << 10, 10, 0;
    transform2d_t<double> res_trans;
    res_trans.matrix() << 1, 0, 10, 0, 1, 10, 0, 0, 1;
    auto calc = transBtwPoses(first_pose, second_pose).matrix();
    EXPECT_TRUE((res_trans.matrix() - calc).cwiseAbs().sum() < EPSILON);
    // ADD_FAILURE() <<res_trans.matrix()<< "\n\n"<< calc;
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, createTransFromPosesAdvance)
{
  try {
    pose2d_t<double> first_pose;
    first_pose << 10, 10, 0;
    pose2d_t<double> second_pose;
    second_pose << 20, 20, 0;
    transform2d_t<double> res_trans;
    res_trans.matrix() << 1, 0, 10, 0, 1, 10, 0, 0, 1;
    auto calc = transBtwPoses(first_pose, second_pose).matrix();
    EXPECT_TRUE((res_trans.matrix() - calc).cwiseAbs().sum() < EPSILON);
    // ADD_FAILURE() <<res_trans.matrix()<< "\n\n"<< calc;
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, createTransFromPosesAngle)
{
  try {
    pose2d_t<double> first_pose;
    first_pose << 0, 0, 0;
    pose2d_t<double> second_pose;
    second_pose << 0, 0, M_PI;
    transform2d_t<double> res_trans;
    res_trans.matrix() << -1, 0, 0, 0, -1, 0, 0, 0, 1;
    auto calc = transBtwPoses(first_pose, second_pose).matrix();
    EXPECT_LT((res_trans.matrix() - calc).cwiseAbs().sum(), EPSILON);
    // ADD_FAILURE() <<res_trans.matrix()<< "\n\n"<< calc;
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}
TEST(EigenTools, transBtwPoses)
{
  pose2d_t<double> first_pose;
  first_pose << 0, 0, 0;
  pose2d_t<double> second_pose;
  second_pose << 29.1929, 130.67, 0.833175;
  transform2d_t<double> res_trans;
  res_trans.matrix() = eigt::vecToMat2d<double>(second_pose);
  auto calc = transBtwPoses(first_pose, second_pose);
  std::cout << getPoseFromTransform(calc).transpose() << std::endl;
  std::cout << getPoseFromTransform(res_trans).transpose() << std::endl;

  EXPECT_LT((res_trans.matrix() - calc.matrix()).cwiseAbs().sum(), EPSILON);
  EXPECT_LT((eigt::vecToMat2d<double>(second_pose).matrix() -
             getTransFromPose(second_pose).matrix())
                .cwiseAbs()
                .sum(),
            EPSILON);
  EXPECT_LT((calc.matrix() - eigt::getTransFromPose(second_pose).matrix())
                .cwiseAbs()
                .sum(),
            EPSILON);
  EXPECT_LT(
      (eigt::transformPose(first_pose, calc) - second_pose).cwiseAbs().sum(),
      EPSILON);
}

TEST(EigenTools, transBtwPoses2)
{
  pose2d_t<double> first_pose;
  first_pose << 29.1929, 130.67, 0.833175;
  pose2d_t<double> second_pose;
  second_pose << 0, 0, 0;
  transform2d_t<double> res_trans;
  res_trans.matrix() = eigt::vecToMat2d<double>(first_pose).inverse();
  auto calc = transBtwPoses(first_pose, second_pose);
  std::cout << getPoseFromTransform(calc).transpose() << std::endl;
  std::cout << getPoseFromTransform(res_trans).transpose() << std::endl;
  EXPECT_LT(
      (eigt::transformPose(first_pose, calc) - second_pose).cwiseAbs().sum(),
      EPSILON);
}

TEST(EigenTools, CycleTest)
{
  try {
    pose2d_t<double> first_pose;
    first_pose << 10, 20, 0;
    pose2d_t<double> second_pose;
    second_pose << -10, 30, -M_PI;
    auto calc_trans = transBtwPoses(first_pose, second_pose);
    auto calc_sec_pose = transformPose(first_pose, calc_trans);
    std::cout << getAngle(calc_trans) << "\n";
    std::cout << getPoseFromTransform(calc_trans) << "\n";
    std::cout << "target pose: " << second_pose.transpose()
              << " calc target pose: " << calc_sec_pose.transpose() << "\n";
    EXPECT_TRUE((calc_sec_pose - second_pose).cwiseAbs().sum() < EPSILON);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, CycleTest2)
{
  try {
    pose2d_t<double> first_pose;
    first_pose << -10, 5, 0;
    pose2d_t<double> second_pose;
    second_pose << 10, 6, 2;
    auto calc_trans = transBtwPoses(first_pose, second_pose);
    auto calc_sec_pose = transformPose(first_pose, calc_trans);
    std::cout << getAngle(calc_trans) << "\n";
    std::cout << getPoseFromTransform(calc_trans) << "\n";
    std::cout << "target pose: " << second_pose.transpose()
              << " calc target pose: " << calc_sec_pose.transpose() << "\n";
    EXPECT_TRUE((calc_sec_pose - second_pose).cwiseAbs().sum() < EPSILON);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, CycleTest3)
{
  try {
    pose2d_t<double> first_pose;
    first_pose << -10, 5, 1;
    pose2d_t<double> second_pose;
    second_pose << 10, 6, M_PI;
    auto calc_trans = transBtwPoses(first_pose, second_pose);
    auto calc_sec_pose = transformPose(first_pose, calc_trans);
    std::cout << getAngle(calc_trans) << "\n";
    std::cout << getPoseFromTransform(calc_trans) << "\n";
    std::cout << "target pose: " << second_pose.transpose()
              << " calc target pose: " << calc_sec_pose.transpose() << "\n";
    EXPECT_TRUE((calc_sec_pose - second_pose).cwiseAbs().sum() < EPSILON);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, complexInterfaceTest)
{
  try {
    pose2d_t<double> source_pose;
    pose2d_t<double> target_pose;
    transform2d_t<double> trans_btw_poses;
    pose2d_t<double> trans_btw_poses_short;
    double angle_diff;
    double displacement;

    source_pose << -10, 5, 1.2;
    target_pose << 10, 6, 2;
    angle_diff = 0.8;
    displacement = std::sqrt(std::pow(source_pose(0) - target_pose(0), 2) +
                             std::pow(source_pose(1) - target_pose(1), 2));
    trans_btw_poses_short << 20, 1, 0.8;
    auto calc_trans = transBtwPoses(source_pose, target_pose);
    auto calc_angle = getAngle(calc_trans);
    auto calc_angle_diff = getAngleDiffrence(source_pose, target_pose);
    auto calc_displacement = getDisplacement(calc_trans);
    auto calc_pose = getPoseFromTransform(calc_trans);
    auto calc_trans_pose = transformPose(source_pose, calc_trans);
    std::cout << calc_trans.matrix() << "\n";
    std::cout << "angle: " << calc_angle << "\n";
    std::cout << "angle diffrence: " << calc_angle_diff << "\n";
    std::cout << "displacement: " << calc_displacement << "\n";
    std::cout << "pose from transform: " << calc_pose << "\n";
    std::cout << "transformed source pose: " << calc_trans_pose.transpose()
              << "target pose reference: " << target_pose.transpose() << "\n";

    EXPECT_LT((trans_btw_poses_short - calc_pose).cwiseAbs().sum(), EPSILON);

    EXPECT_TRUE(std::abs(angle_diff - calc_angle) < EPSILON);
    EXPECT_TRUE(std::abs(displacement - calc_displacement) < EPSILON);
    EXPECT_TRUE((target_pose - calc_trans_pose).cwiseAbs().sum() < EPSILON);

  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, normalizeAngle)
{
  EXPECT_LT(normalizeAngle(M_PI) - M_PI, EPSILON);
  EXPECT_LT(normalizeAngle(3 * M_PI) - M_PI, EPSILON);
  EXPECT_LT(normalizeAngle(2 * M_PI), EPSILON);
  EXPECT_LT(normalizeAngle(3 * (M_PI / 2)) - M_PI + (M_PI / 2), EPSILON);

  EXPECT_LT(normalizeAngle(-M_PI) + M_PI, EPSILON);
  EXPECT_LT(normalizeAngle(-2 * M_PI), EPSILON);
  EXPECT_LT(normalizeAngle(-3 * M_PI) + M_PI, EPSILON);
}

TEST(EigenTools, getTransFromPoseTest)
{
  try {
    pose2d_t<double> first_pose;
    first_pose << -10, 5, 0;
    pose2d_t<double> second_pose;
    second_pose << 10, 6, 2;
    auto calc_trans = transBtwPoses(first_pose, second_pose);
    auto calc_pose = getPoseFromTransform(calc_trans);
    EXPECT_LT((calc_trans.matrix() - getTransFromPose(calc_pose).matrix())
                  .cwiseAbs()
                  .sum(),
              EPSILON);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, transBtwFrames)
{
  pose2d_t<double> first_pose;
  first_pose << 2, 2, M_PI / 2;
  pose2d_t<double> second_pose;
  second_pose << 0, 0, 0;
  Eigen::Vector2d point(0, 0), point1(0, 1), point2(1, 0), point3(1, 1);
  Eigen::Vector2d point_r(2, 2), point1_r(1, 2), point2_r(2, 3), point3_r(1, 3);
  transform2d_t<double> res;
  res = transBtwFrames(first_pose, second_pose);
  std::cout << res * point << std::endl;
  std::cout << res * point1 << std::endl;
  std::cout << res * point2 << std::endl;
  std::cout << res * point3 << std::endl;
  EXPECT_LT((res * point - point_r).cwiseAbs().sum(), EPSILON);
  EXPECT_LT((res * point1 - point1_r).cwiseAbs().sum(), EPSILON);
  EXPECT_LT((res * point2 - point2_r).cwiseAbs().sum(), EPSILON);
  EXPECT_LT((res * point3 - point3_r).cwiseAbs().sum(), EPSILON);
}

TEST(EigenTools, transBtwFrames2)
{
  pose2d_t<double> first_pose;
  first_pose << 2, 2, M_PI / 4;
  pose2d_t<double> second_pose;
  second_pose << 0, 0, M_PI / 2;
  Eigen::Vector2d point(0, 0), point1(1, 0), point2(0, 1);
  Eigen::Vector2d point_r(2, -2), point1_r(2.7, -2.7), point2_r(2.7, -1.3);
  transform2d_t<double> res;
  res = transBtwFrames(first_pose, second_pose);
  std::cout << res * point << std::endl;
  std::cout << res * point1 << std::endl;
  std::cout << res * point2 << std::endl;
  EXPECT_LT((res * point - point_r).cwiseAbs().sum(), EPSILON * 100);
  EXPECT_LT((res * point1 - point1_r).cwiseAbs().sum(), EPSILON * 100);
  EXPECT_LT((res * point2 - point2_r).cwiseAbs().sum(), EPSILON * 100);
}

TEST(EigenTools, transBtwFrames3)
{
  pose2d_t<double> first_pose;
  first_pose << 2, 2, M_PI / 4;
  pose2d_t<double> second_pose;
  second_pose << -1, 3, M_PI / 2;
  Eigen::Vector2d point(0, 0), point1(1, 0), point2(0, 1);
  Eigen::Vector2d point_r(-1, -3), point1_r(-0.3, -3.7), point2_r(-0.29, -2.29);
  transform2d_t<double> res;
  res = transBtwFrames(first_pose, second_pose);
  std::cout << res * point << std::endl;
  std::cout << res * point1 << std::endl;
  std::cout << res * point2 << std::endl;
  EXPECT_LT((res * point - point_r).cwiseAbs().sum(), EPSILON * 100);
  EXPECT_LT((res * point1 - point1_r).cwiseAbs().sum(), EPSILON * 100);
  EXPECT_LT((res * point2 - point2_r).cwiseAbs().sum(), EPSILON * 100);
}

TEST(EigenTools, transBtwFrames4)
{
  pose2d_t<double> first_pose;
  first_pose << 2, 2, 0;
  pose2d_t<double> second_pose;
  second_pose << -1, 3, 0;
  Eigen::Vector2d point(0, 0), point1(1, 0), point2(0, 1);
  Eigen::Vector2d point_r(3, -1), point1_r(4, -1), point2_r(3, 0);
  transform2d_t<double> res;
  res = transBtwFrames(first_pose, second_pose);
  std::cout << res * point << std::endl;
  std::cout << res * point1 << std::endl;
  std::cout << res * point2 << std::endl;
  EXPECT_LT((res * point - point_r).cwiseAbs().sum(), EPSILON * 100);
  EXPECT_LT((res * point1 - point1_r).cwiseAbs().sum(), EPSILON * 100);
  EXPECT_LT((res * point2 - point2_r).cwiseAbs().sum(), EPSILON * 100);
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
