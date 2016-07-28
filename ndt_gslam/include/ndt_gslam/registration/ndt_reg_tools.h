#ifndef NDT_GSLAM_NDT_REG_TOOLS
#define NDT_GSLAM_NDT_REG_TOOLS
#include <Eigen/Dense>

namespace pcl
{
namespace ndt_reg
{
/** \brief Class to store vector value and first and second derivatives
* (grad vector and hessian matrix), so they can be returned easily from
* functions
*/
template <unsigned N = 3, typename T = double>
struct ScoreAndDerivatives {
  ScoreAndDerivatives() : hessian_(), gradient_(), value_()
  {
    hessian_.setZero();
    gradient_.setZero();
    value_ = 0;
  }

  Eigen::Matrix<T, N, N> hessian_;
  Eigen::Matrix<T, N, 1> gradient_;
  T value_;

  static ScoreAndDerivatives<N, T> Zero()
  {
    ScoreAndDerivatives<N, T> r;
    r.hessian_ = Eigen::Matrix<T, N, N>::Zero();
    r.gradient_ = Eigen::Matrix<T, N, 1>::Zero();
    r.value_ = 0;
    return r;
  }

  ScoreAndDerivatives<N, T> &operator+=(ScoreAndDerivatives<N, T> const &r)
  {
    hessian_ += r.hessian_;
    gradient_ += r.gradient_;
    value_ += r.value_;
    return *this;
  }
};
template <unsigned N, typename T>
ScoreAndDerivatives<N, T> operator+(const ScoreAndDerivatives<N, T> &lhs,
                                    const ScoreAndDerivatives<N, T> &rhs)
{
  ScoreAndDerivatives<N, T> ret;
  ret += lhs;
  ret += rhs;
  return ret;
}

//////////////////////////////////// PARAMETERS FOR FITTING OF COVARIANCE

struct FittingParams {
  double gauss_d1_;
  double gauss_d2_;
  double gauss_d2__half_;

  FittingParams(double outlier_ratio, double resolution)
    : gauss_d1_(0), gauss_d2_(0), gauss_d2__half_(0)
  {
    calcParams(outlier_ratio, resolution);
  }

private:
  void calcParams(double outlier_ratio, double resolution)
  {
    double gauss_c1, gauss_c2, gauss_d3;
    // Initializes the guassian fitting parameters (eq. 6.8) [Magnusson 2009]
    gauss_c1 = 10.0 * (1 - outlier_ratio);
    gauss_c2 = outlier_ratio / pow(resolution, 2);
    gauss_d3 = -log(gauss_c2);
    gauss_d1_ = -log(gauss_c1 + gauss_c2) - gauss_d3;
    gauss_d2_ = -2 * log((-log(gauss_c1 * exp(-0.5) + gauss_c2) - gauss_d3) /
                         gauss_d1_);
    gauss_d2__half_ = gauss_d2_ / 2;
  }
};

struct JacobianHessianDerivatives {
  Eigen::Matrix<double, 3, 3> Jest;
  Eigen::Matrix<double, 9, 3> Hest;
  Eigen::Matrix<double, 3, 9> Zest;
  Eigen::Matrix<double, 9, 9> ZHest;

  JacobianHessianDerivatives()
  {
    setZero();
  }
  void setZero()
  {
    Jest.setZero();
    Hest.setZero();
    Zest.setZero();
    ZHest.setZero();
  }
  static JacobianHessianDerivatives Zero()
  {
    JacobianHessianDerivatives res;
    res.setZero();
    return res;
  }
};

bool updateIntervalMT(double &a_l, double &f_l, double &g_l, double &a_u,
                      double &f_u, double &g_u, double a_t, double f_t,
                      double g_t)
{
  // Case U1 in Update Algorithm and Case a in Modified Update Algorithm
  // [More,
  // Thuente 1994]
  if (f_t > f_l) {
    a_u = a_t;
    f_u = f_t;
    g_u = g_t;
    return (false);
  }
  // Case U2 in Update Algorithm and Case b in Modified Update Algorithm
  // [More,
  // Thuente 1994]
  else if (g_t * (a_l - a_t) > 0) {
    a_l = a_t;
    f_l = f_t;
    g_l = g_t;
    return (false);
  }
  // Case U3 in Update Algorithm and Case c in Modified Update Algorithm
  // [More,
  // Thuente 1994]
  else if (g_t * (a_l - a_t) < 0) {
    a_u = a_l;
    f_u = f_l;
    g_u = g_l;

    a_l = a_t;
    f_l = f_t;
    g_l = g_t;
    return (false);
  }
  // Interval Converged
  else
    return (true);
}

double trialValueSelectionMT(double a_l, double f_l, double g_l, double a_u,
                             double f_u, double g_u, double a_t, double f_t,
                             double g_t)
{
  // Case 1 in Trial Value Selection [More, Thuente 1994]
  if (f_t > f_l) {
    // Calculate the minimizer of the cubic that interpolates f_l, f_t, g_l
    // and
    // g_t
    // Equation 2.4.52 [Sun, Yuan 2006]
    double z = 3 * (f_t - f_l) / (a_t - a_l) - g_t - g_l;
    double w = std::sqrt(z * z - g_t * g_l);
    // Equation 2.4.56 [Sun, Yuan 2006]
    double a_c = a_l + (a_t - a_l) * (w - g_l - z) / (g_t - g_l + 2 * w);

    // Calculate the minimizer of the quadratic that interpolates f_l, f_t and
    // g_l
    // Equation 2.4.2 [Sun, Yuan 2006]
    double a_q =
        a_l - 0.5 * (a_l - a_t) * g_l / (g_l - (f_l - f_t) / (a_l - a_t));

    if (std::fabs(a_c - a_l) < std::fabs(a_q - a_l))
      return (a_c);
    else
      return (0.5 * (a_q + a_c));
  }
  // Case 2 in Trial Value Selection [More, Thuente 1994]
  else if (g_t * g_l < 0) {
    // Calculate the minimizer of the cubic that interpolates f_l, f_t, g_l
    // and
    // g_t
    // Equation 2.4.52 [Sun, Yuan 2006]
    double z = 3 * (f_t - f_l) / (a_t - a_l) - g_t - g_l;
    double w = std::sqrt(z * z - g_t * g_l);
    // Equation 2.4.56 [Sun, Yuan 2006]
    double a_c = a_l + (a_t - a_l) * (w - g_l - z) / (g_t - g_l + 2 * w);

    // Calculate the minimizer of the quadratic that interpolates f_l, g_l and
    // g_t
    // Equation 2.4.5 [Sun, Yuan 2006]
    double a_s = a_l - (a_l - a_t) / (g_l - g_t) * g_l;

    if (std::fabs(a_c - a_t) >= std::fabs(a_s - a_t))
      return (a_c);
    else
      return (a_s);
  }
  // Case 3 in Trial Value Selection [More, Thuente 1994]
  else if (std::fabs(g_t) <= std::fabs(g_l)) {
    // Calculate the minimizer of the cubic that interpolates f_l, f_t, g_l
    // and
    // g_t
    // Equation 2.4.52 [Sun, Yuan 2006]
    double z = 3 * (f_t - f_l) / (a_t - a_l) - g_t - g_l;
    double w = std::sqrt(z * z - g_t * g_l);
    double a_c = a_l + (a_t - a_l) * (w - g_l - z) / (g_t - g_l + 2 * w);

    // Calculate the minimizer of the quadratic that interpolates g_l and g_t
    // Equation 2.4.5 [Sun, Yuan 2006]
    double a_s = a_l - (a_l - a_t) / (g_l - g_t) * g_l;

    double a_t_next;

    if (std::fabs(a_c - a_t) < std::fabs(a_s - a_t))
      a_t_next = a_c;
    else
      a_t_next = a_s;

    if (a_t > a_l)
      return (std::min(a_t + 0.66 * (a_u - a_t), a_t_next));
    else
      return (std::max(a_t + 0.66 * (a_u - a_t), a_t_next));
  }
  // Case 4 in Trial Value Selection [More, Thuente 1994]
  else {
    // Calculate the minimizer of the cubic that interpolates f_u, f_t, g_u
    // and
    // g_t
    // Equation 2.4.52 [Sun, Yuan 2006]
    double z = 3 * (f_t - f_u) / (a_t - a_u) - g_t - g_u;
    double w = std::sqrt(z * z - g_t * g_u);
    // Equation 2.4.56 [Sun, Yuan 2006]
    return (a_u + (a_t - a_u) * (w - g_u - z) / (g_t - g_u + 2 * w));
  }
}
/** \brief Auxilary function used to determin endpoints of More-Thuente
 * interval.
  * \note \f$ \psi(\alpha) \f$ in Equation 1.6 (Moore, Thuente 1994)
  * \param[in] a the step length, \f$ \alpha \f$ in More-Thuente (1994)
  * \param[in] f_a function value at step length a, \f$ \phi(\alpha) \f$ in
 * More-Thuente (1994)
  * \param[in] f_0 initial function value, \f$ \phi(0) \f$ in Moore-Thuente
 * (1994)
  * \param[in] g_0 initial function gradiant, \f$ \phi'(0) \f$ in More-Thuente
 * (1994)
  * \param[in] mu the step length, constant \f$ \mu \f$ in Equation 1.1 [More,
 * Thuente 1994]
  * \return sufficent decrease value
  */
double auxilaryFunction_PsiMT(double a, double f_a, double f_0, double g_0,
                              double mu = 1.e-4)
{
  return (f_a - f_0 - mu * g_0 * a);
}

/** \brief Auxilary function derivative used to determin endpoints of
 * More-Thuente interval.
  * \note \f$ \psi'(\alpha) \f$, derivative of Equation 1.6 (Moore, Thuente
 * 1994)
  * \param[in] g_a function gradient at step length a, \f$ \phi'(\alpha) \f$
 * in More-Thuente (1994)
  * \param[in] g_0 initial function gradiant, \f$ \phi'(0) \f$ in More-Thuente
 * (1994)
  * \param[in] mu the step length, constant \f$ \mu \f$ in Equation 1.1 [More,
 * Thuente 1994]
  * \return sufficent decrease derivative
  */
double auxilaryFunction_dPsiMT(double g_a, double g_0, double mu = 1.e-4)
{
  return (g_a - mu * g_0);
}

/**
 * @brief      Converts vector representation of the transformation into matrix
 *             representation
 *
 * @param[in]  trans  The transformation
 *
 * @tparam     T      Scalar of the matrix
 *
 * @return     transformation matrix
 */
template <typename T = float>
Eigen::Matrix<T, 4, 4> vecToMat(const Eigen::Vector3d &trans)
{
  Eigen::Matrix<T, 4, 4> trans_mat = Eigen::Matrix<T, 4, 4>::Identity();

  trans_mat.block(0, 0, 3, 3).matrix() =
      Eigen::Matrix<T, 3, 3>(Eigen::AngleAxis<T>(
          static_cast<T>(trans(2)), Eigen::Matrix<T, 3, 1>::UnitZ()));

  trans_mat.block(0, 3, 3, 1).matrix() = Eigen::Matrix<T, 3, 1>(
      static_cast<T>(trans(0)), static_cast<T>(trans(1)), 0.0);

  return trans_mat;
}

/**
 * @brief      Convert matrix representation of translation into vector
 *             representation of the pose
 *
 * @param[in]  trans  The transformation
 *
 * @tparam     T      Scalar of the matrix
 *
 * @return     The pose
 */
template <typename T = float>
Eigen::Vector3d matToVec(const Eigen::Matrix<T, 4, 4> &trans)
{
  Eigen::Vector3d vec;
  Eigen::Transform<T, 3, Eigen::Affine, Eigen::ColMajor> trans_mat(trans);
  Eigen::Matrix<T, 3, 1> translation = trans_mat.translation();
  vec << translation(0), translation(1),
      std::atan2(trans_mat.rotation()(1, 0), trans_mat.rotation()(0, 0));
  return vec;
}

}  // end of ndt_reg
}  // end of pcl

#endif
