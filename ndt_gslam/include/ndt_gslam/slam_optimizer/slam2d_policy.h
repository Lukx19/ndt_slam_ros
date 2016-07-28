#ifndef NDT_GSLAM_SLAM2D_POLICY
#define NDT_GSLAM_SLAM2D_POLICY

#include <math.h>
#include <ndt_gslam/utils/eigen_tools.h>
#include <Eigen/Dense>
struct Slam2d_Policy {
private:
  // typedef Eigen::Matrix<double, 3, 3> TransformMatrixRaw;
public:
  // change sizes to multiples of 2 for use of vector instructions
  typedef size_t Id;
  typedef Eigen::Matrix<double, 3, 3> JacobianMatrix;
  typedef Eigen::Matrix<double, 3, 3> InformMatrix;
  typedef Eigen::Matrix<double, 3, 3> CovarMatrix;
  typedef Eigen::Matrix<double, 3, 1> ErrorVector;
  typedef std::pair<JacobianMatrix, JacobianMatrix> JacobianPair;
  typedef Eigen::Matrix<double, 3, 1> Pose;
  typedef Eigen::Transform<double, 2, Eigen::TransformTraits::Affine>
      TransformMatrix;

  static const size_t block_size = 3;  // change to 6 for 3D world

  static TransformMatrix vecToTransMat(const Pose &p)
  {
    return eigt::getTransFromPose(p);
  }

  static Pose transMatToVec(const TransformMatrix &trans)
  {
    return eigt::getPoseFromTransform(trans);
  }

  static JacobianPair calcJacobianBlocks(const Pose &xi, const Pose &xj,
                                         const Pose &trans)
  {
    JacobianMatrix aij = JacobianMatrix::Ones();
    JacobianMatrix bij = JacobianMatrix::Ones();
    double angle = trans(2) + xi(2);
    auto dist = xj - xi;
    double co = cos(angle);
    double si = sin(angle);
    aij(0, 0) = -co;
    aij(0, 1) = -si;
    aij(0, 2) = -si * dist(0) + co * dist(1);
    aij(1, 0) = si;
    aij(1, 1) = -co;
    aij(1, 2) = -co * dist(0) - si * dist(1);
    aij(2, 0) = 0;
    aij(2, 1) = 0;
    aij(2, 2) = -1;

    bij(0, 0) = co;
    bij(0, 1) = si;
    bij(0, 2) = 0;
    bij(1, 0) = -si;
    bij(1, 1) = co;
    bij(1, 2) = 0;
    bij(2, 0) = 0;
    bij(2, 1) = 0;
    bij(2, 2) = 1;

    return std::make_pair(std::move(aij), std::move(bij));
  }

  static ErrorVector calcError(const Pose &x, const Pose &y, const Pose &trans)
  {
    TransformMatrix xi = vecToTransMat(x);
    TransformMatrix xj = vecToTransMat(y);
    TransformMatrix zij = vecToTransMat(trans);
    return transMatToVec(zij.inverse() * (xi.inverse() * xj));
  }

  static Pose addPoses(const Pose &a, const Pose &b)
  {
    Pose ret;
    ret = a + b;
    ret(2) = eigt::normalizeAngle(a(2) + b(2));
    return ret;
  }
};

#endif
