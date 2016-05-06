#ifndef DYNAMIC_SLAM_UTILS_EIGEN_TOOLS
#define DYNAMIC_SLAM_UTILS_EIGEN_TOOLS

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592
#endif

namespace eigt
{

typedef Eigen::Transform<double, 2, 1> transform2d_t;


typedef Eigen::Matrix<double, 2, 1> point2d_t;


typedef  Eigen::Matrix<double, 3, 1> pose2d_t;
// ************************* DECLARATION *************************************

/// Calculates transformation between inputed poses.

transform2d_t transBtwPoses(const pose2d_t &from, const pose2d_t &to);

// Transforms pose in global frame to its new posistion based on transform

pose2d_t transformPose(const pose2d_t &pose,
                          const transform2d_t &trans);

// maps from transformation matrix 3x3 to vector
// [delta_x,delta_y, detla_angle]

pose2d_t getPoseFromTransform(const transform2d_t &trans);

// returns angle which is rotational part of transformation

double getAngle(const transform2d_t &trans);

// returns square distance from translation part of transform

double getDisplacement(const transform2d_t &trans);

// returns absolute angular diffrence between poses. Selects shorter angle.

double getAngleDiffrence(const pose2d_t &from, const pose2d_t &to);

// Maps 4x4 transformation matrix to 3x3 transformation matrix

transform2d_t convertToTransform(const Eigen::Matrix<double, 4, 4> &trans);

// MAPS from 3x3 transformation matrix to 4x4 Transformation matrix

Eigen::Matrix<double, 4, 4> convertFromTransform(const transform2d_t &trans);

// Maps transformation encoded in pose vector [delta_x,delta_y,delta_angle] to
// transformation matrix

eigt::transform2d_t getTransFromPose(const pose2d_t &trans);


double normalizeAngle(double angle);

}
#endif
