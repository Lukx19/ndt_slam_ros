#ifndef NDT_GSLAM_NDT2D
#define NDT_GSLAM_NDT2D

#include <ndt_gslam/ndt/ndt_cell.h>
#include <ndt_gslam/ndt/ndt_grid2d.h>
#include <ndt_gslam/registration/ndt_reg_tools.h>
#include <pcl/filters/voxel_grid_covariance.h>
#include <pcl/registration/registration.h>
#include <ros/ros.h>
#include <Eigen/Dense>

namespace pcl
{
template <typename PointSource, typename PointTarget,
          typename CellType = slamuk::NDTCell>
class NormalDistributionsTransform2DEx
    : public Registration<PointSource, PointTarget>
{
protected:
  typedef typename Registration<PointSource, PointTarget>::PointCloudSource
      PclSource;
  typedef typename PclSource::Ptr PclSourcePtr;
  typedef typename PclSource::ConstPtr PclSourceConstPtr;

  typedef typename Registration<PointSource, PointTarget>::PointCloudTarget
      PclTarget;
  typedef typename PclTarget::Ptr PclTargetPtr;
  typedef typename PclTarget::ConstPtr PclTargetConstPtr;

  typedef PointIndices::Ptr PointIndicesPtr;
  typedef PointIndices::ConstPtr PointIndicesConstPtr;

  /** \brief Typename of searchable voxel grid containing mean and covariance.*/
  typedef slamuk::NDTGrid2D<CellType, PointTarget> GridTarget;
  typedef slamuk::NDTGrid2D<CellType, PointSource> GridSource;
  typedef typename GridTarget::ConstPtr GridTargetConstPtr;
  typedef typename GridSource::ConstPtr GridSourceConstPtr;

public:
  typedef boost::shared_ptr<
      NormalDistributionsTransform2DEx<PointSource, PointTarget>>
      Ptr;
  typedef boost::shared_ptr<
      const NormalDistributionsTransform2DEx<PointSource, PointTarget>>
      ConstPtr;
  typedef Eigen::Vector3d VectorTrans;

  NormalDistributionsTransform2DEx();

  /** \brief Empty destructor */
  virtual ~NormalDistributionsTransform2DEx()
  {
  }

  virtual void setInputTarget(const PclTargetConstPtr &cloud)
  {
    Registration<PointSource, PointTarget>::setInputTarget(cloud);
    GridSource *tmp = new GridSource();
    tmp->setCellSize(cell_sizes_.back());
    tmp->initializeSimple(*cloud);
    target_grid_ = GridTargetConstPtr(tmp);
  }

  virtual void setInputTarget(const GridTargetConstPtr &grid)
  {
    target_grid_ = grid;
    PclTargetConstPtr pcl = grid->getMeans();
    Registration<PointSource, PointTarget>::setInputTarget(pcl);
  }

  inline void setCellSize(float cell_size)
  {
    if (initCellSizes(cell_size)) {
      initParams();
    }
  }

  inline float getCellSize() const
  {
    return (cell_sizes_.back());
  }

  inline void setNumLayers(size_t num)
  {
    layer_count_ = num;
    initCellSizes(getCellSize());
    initParams();
  }

  inline size_t getNumLayers()
  {
    return layer_count_;
  }
  /** \brief Get the point cloud outlier ratio.
    * \return outlier ratio
    */
  inline double getOulierRatio() const
  {
    return (outlier_ratio_);
  }

  /** \brief Set/change the point cloud outlier ratio.
    * \param[in] outlier_ratio outlier ratio
    */
  inline void setOulierRatio(double outlier_ratio)
  {
    outlier_ratio_ = outlier_ratio;
    initParams();
  }
  /** \brief Get the newton line search maximum step length.
  * \return maximum step length
  */
  inline double getStepSize() const
  {
    return (step_size_);
  }

  /** \brief Set/change the newton line search maximum step length.
    * \param[in] step_size maximum step length
    */
  inline void setStepSize(double step_size)
  {
    step_size_ = step_size;
  }
  /** \brief Get the registration alignment probability.
    * \return transformation probability
    */
  inline double getTransformationProbability() const
  {
    return (trans_probability_);
  }

  /** \brief Get the number of iterations required to calculate alignment.
    * \return final number of iterations
    */
  inline int getFinalNumIteration() const
  {
    return (nr_iterations_);
  }

  inline Eigen::Matrix3d getCovariance() const
  {
    return covariance_;
  }

  inline Eigen::Matrix3d getInformMatrix() const
  {
    return inform_matrix_;
  }

protected:
  using Registration<PointSource, PointTarget>::reg_name_;
  using Registration<PointSource, PointTarget>::getClassName;
  using Registration<PointSource, PointTarget>::input_;
  using Registration<PointSource, PointTarget>::indices_;
  using Registration<PointSource, PointTarget>::target_;
  using Registration<PointSource, PointTarget>::nr_iterations_;
  using Registration<PointSource, PointTarget>::max_iterations_;
  using Registration<PointSource, PointTarget>::previous_transformation_;
  using Registration<PointSource, PointTarget>::final_transformation_;
  using Registration<PointSource, PointTarget>::transformation_;
  using Registration<PointSource, PointTarget>::transformation_epsilon_;
  using Registration<PointSource, PointTarget>::converged_;
  using Registration<PointSource, PointTarget>::corr_dist_threshold_;
  using Registration<PointSource, PointTarget>::inlier_threshold_;
  using Registration<PointSource, PointTarget>::update_visualizer_;

  /** \brief The side length of voxels. */
  std::vector<float> cell_sizes_;

  /** \brief The maximum step length. */
  double step_size_;

  /** \brief The ratio of outliers of points w.r.t. a normal distribution,
   * Equation 6.7 [Magnusson 2009]. */
  double outlier_ratio_;

  std::vector<ndt_reg::FittingParams> params_;

  /** \brief The probability score of the transform applied to the input cloud,
   * Equation 6.9 and 6.10 [Magnusson 2009]. */
  double trans_probability_;

  Eigen::Matrix3d covariance_;
  Eigen::Matrix3d inform_matrix_;

  size_t layer_count_;
  GridTargetConstPtr target_grid_;

  /** \brief Estimate the transformation and returns the transformed source
   * (input) as output.
    * \param[out] output the resultant input transfomed point cloud dataset
    */
  virtual void computeTransformation(PclSource &output)
  {
    computeTransformation(output, Eigen::Matrix4f::Identity());
  }

  /** \brief Estimate the transformation and returns the transformed source
   * (input) as output.
    * \param[out] output the resultant input transfomed point cloud dataset
    * \param[in] guess the initial gross estimation of the transformation
    */
  virtual void computeTransformation(PclSource &output,
                                     const Eigen::Matrix4f &guess);

  virtual bool computeSingleGrid(PclSource &output,
                                 const GridTarget &target_grid,
                                 const Eigen::Matrix4f &guess,
                                 const ndt_reg::FittingParams &param);

  virtual inline void initParams()
  {
    params_.clear();
    for (size_t i = 0; i < cell_sizes_.size(); ++i) {
      params_.push_back(
          ndt_reg::FittingParams(outlier_ratio_, 1 / cell_sizes_[i]));
    }
  }

  virtual inline bool initCellSizes(float base_size)
  {
    if (cell_sizes_.size() > 0 && cell_sizes_.back() == base_size &&
        cell_sizes_.size() == layer_count_)
      return false;

    cell_sizes_.clear();
    for (int i = layer_count_ - 1; i >= 0; --i) {
      cell_sizes_.push_back(base_size * std::pow(2, i));
    }
    return true;
  }

  /** \brief Estimate score, hessian and gradient for selected source point
   * cloud based on given transformation.
    * \param[in] ndt_grid NDT target grid used as registration target.
    * \param[in] source Source point cloud used to register agains NDT target
   * grid
    * \param[in] trans Transformation used for projecting points of source to
   * grid cells of NDT grid
    * \return Estimated parameters score, gradient, hessian.
    */
  virtual ndt_reg::ScoreAndDerivatives<3, double>
  calcScore(const GridTarget &ndt_grid, const PclSource &source,
            const Eigen::Vector3d &trans, const ndt_reg::FittingParams &param);
  /** \brief Estimate score, hessian and gradient for selected point.
    * \param[in] cell NDT target grid cell for selected point.
    * \param[in] trans_pt Transformed source cloud point.
    * \param[in] pt Source cloud point befor transformation.
    * \param[in] cos_theta Cosine of transformation angle.
    * \param[in] sin_theta Sine of transformation angle.
    * \return Estimated parameters score, gradient, hessian.
    */
  virtual ndt_reg::ScoreAndDerivatives<3, double>
  calcPointScore(const PointSource &pt, const CellType *cell,
                 const PointSource &trans_pt, double cos_theta,
                 double sin_theta, const ndt_reg::FittingParams &param) const;

  // linear search methods
  double computeStepLengthMT(
      const Eigen::Matrix<double, 3, 1> &x,
      Eigen::Matrix<double, 3, 1> &step_dir, double step_init, double step_max,
      double step_min, const GridTarget &ndt_grid,
      const ndt_reg::ScoreAndDerivatives<3, double> &score,
      const PclSource &input, const ndt_reg::FittingParams &param);
};

template <typename PointSource, typename PointTarget, typename CellType>
NormalDistributionsTransform2DEx<PointSource, PointTarget,
                                 CellType>::NormalDistributionsTransform2DEx()
  : step_size_(0.05)
  , outlier_ratio_(0.30)
  , trans_probability_()
  , layer_count_(4)
{
  nr_iterations_ = 0;
  max_iterations_ = 35;
  transformation_epsilon_ = 0.1;
  converged_ = false;
  initCellSizes(0.25);
  initParams();
}
/////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
void NormalDistributionsTransform2DEx<PointSource, PointTarget, CellType>::
    computeTransformation(PclSource &output, const Eigen::Matrix4f &guess)
{
  PclSource out;
  ROS_DEBUG_STREAM("[NDT2D]: guess:" << ndt_reg::matToVec(guess).transpose());
  converged_ = false;
  if (!target_grid_) {
    ROS_ERROR_STREAM("[NDT2D]: target cloud or grid not set");
    return;
  }
  if (!input_) {
    ROS_ERROR_STREAM("[NDT2D]: source cloud not set");
    return;
  }
  transformation_ = guess;
  for (size_t i = 0; i < layer_count_; ++i) {
    if (!computeSingleGrid(out, target_grid_->createCoarserGrid(cell_sizes_[i]),
                           transformation_, params_[i])) {
      converged_ = false;
      return;
    }
  }
  ROS_DEBUG_STREAM("[NDT2D]: final trans:"
                   << ndt_reg::matToVec(transformation_).transpose());
  transformPointCloud(*input_, output, transformation_);
  final_transformation_ = transformation_;
  converged_ = true;
}

////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
bool NormalDistributionsTransform2DEx<PointSource, PointTarget, CellType>::
    computeSingleGrid(PclSource &output, const GridTarget &target_grid,
                      const Eigen::Matrix4f &guess,
                      const ndt_reg::FittingParams &param)
{
  nr_iterations_ = 0;
  converged_ = false;
  // Initialise final transformation to the guessed one
  final_transformation_ = guess;
  previous_transformation_ = final_transformation_;
  // variables needed for calculation
  Eigen::Vector3d xytheta_p = ndt_reg::matToVec(guess);
  Eigen::Vector3d delta_xytheta_p;
  Eigen::Matrix4f p;
  double delta_p_norm;
  ndt_reg::ScoreAndDerivatives<> score;
  while (!converged_) {
    score = calcScore(target_grid, *input_, xytheta_p, param);
    // Solve for decent direction using newton method, line 23 in Algorithm 2
    // [Magnusson 2009]
    Eigen::JacobiSVD<Eigen::Matrix3d> sv(
        score.hessian_, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // Negative for maximization as opposed to minimization
    delta_xytheta_p = sv.solve(-score.gradient_);

    // Calculate step length with guarnteed sufficient decrease [More, Thuente
    // 1994]
    delta_p_norm = delta_xytheta_p.norm();
    if (delta_p_norm == 0 || delta_p_norm != delta_p_norm) {
      trans_probability_ =
          score.value_ / static_cast<double>(input_->points.size());
      converged_ = delta_p_norm == delta_p_norm;
      covariance_.setIdentity();
      inform_matrix_.setIdentity();
      ROS_ERROR_STREAM(
          "[NDT2D]: Not enough overlap. Probability: " << trans_probability_);
      return false;
    }
    delta_xytheta_p.normalize();
    delta_p_norm = computeStepLengthMT(xytheta_p, delta_xytheta_p, delta_p_norm,
                                       step_size_, transformation_epsilon_ / 2,
                                       target_grid, score, output, param);
    // delta_p_norm = 0.8
    delta_xytheta_p *= delta_p_norm;
    xytheta_p += delta_xytheta_p;
    p = ndt_reg::vecToMat<float>(xytheta_p);

    // convergence testing
    if (nr_iterations_ > max_iterations_ ||
        (nr_iterations_ &&
         (std::abs(delta_p_norm) < transformation_epsilon_))) {
      converged_ = true;
    }

    ++nr_iterations_;
    previous_transformation_ = transformation_;
    transformation_ = p;
    trans_probability_ =
        score.value_ / static_cast<double>(input_->points.size());
    ROS_DEBUG_STREAM("[NDT2D]: Step: "
                     << delta_p_norm
                     << " Delta: " << delta_xytheta_p.transpose()
                     << " Score: " << score.value_
                     << "probability of match: " << trans_probability_);
    // Update Visualizer (untested)
    if (update_visualizer_ != 0) {
      transformPointCloud(*input_, output, transformation_);
      update_visualizer_(output, std::vector<int>(), *target_,
                         std::vector<int>());
    }
  }
  final_transformation_ = transformation_;
  transformPointCloud(*input_, output, transformation_);
  covariance_ = score.hessian_;
  inform_matrix_ = score.hessian_.inverse();
  return true;
}

/////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
ndt_reg::ScoreAndDerivatives<3, double>
NormalDistributionsTransform2DEx<PointSource, PointTarget, CellType>::calcScore(
    const GridTarget &ndt_grid, const PclSource &source,
    const Eigen::Vector3d &trans, const ndt_reg::FittingParams &param)
{
  PclSource trans_source;
  ndt_reg::ScoreAndDerivatives<3, double> res;
  double cos_theta, sin_theta;
  // prepare angles for jaccobian and hessian derivatives usage
  // simplify small angles
  if (std::abs(trans(2)) < 10e-5) {
    cos_theta = 1.0;
    sin_theta = 0.0;
  } else {
    cos_theta = std::cos(trans(2));
    sin_theta = std::sin(trans(2));
  }
  res = ndt_reg::ScoreAndDerivatives<3, double>::Zero();
  transformPointCloud(source, trans_source, ndt_reg::vecToMat(trans));

  // calculate score for each point
  typedef typename GridTarget::CellPtrVector TargetVec;
  for (size_t i = 0; i < trans_source.size(); ++i) {
    Eigen::Vector2d pt(trans_source[i].x, trans_source[i].y);

    TargetVec neighbors = ndt_grid.getKNearestNeighbors(pt, 4);

    for (size_t nb = 0; nb < neighbors.size(); nb++) {
      res += calcPointScore(source[i], neighbors[nb], trans_source[i],
                            cos_theta, sin_theta, param);
    }
    // ROS_DEBUG_STREAM("Number of neighbours for point: "<<count);
  }
  return res;
}

////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
ndt_reg::ScoreAndDerivatives<3, double> NormalDistributionsTransform2DEx<
    PointSource, PointTarget,
    CellType>::calcPointScore(const PointSource &pt, const CellType *cell,
                              const PointSource &trans_pt, double cos_theta,
                              double sin_theta,
                              const ndt_reg::FittingParams &param) const
{
  ndt_reg::ScoreAndDerivatives<3, double> res;
  Eigen::Matrix<double, 2, 3> jacobian;
  Eigen::Vector2d hess_derivative, hess_derivatives2;
  double x = pt.x;
  double y = pt.y;
  double co = cos_theta;
  double si = sin_theta;

  jacobian << 1, 0, -x * si - y * co, 0, 1, x * co - y * si;
  hess_derivative << -x * co + y * si, -x * si - y * co;
  hess_derivatives2 << 0, 0;

  Eigen::Matrix2d icov = cell->getICov().block(0, 0, 2, 2);
  Eigen::Vector2d mean = cell->getMean().head(2);
  Eigen::Vector2d trans_point(trans_pt.x, trans_pt.y);
  Eigen::Vector2d pt_mean_diff = trans_point - mean;
  Eigen::Vector2d diff_icov = pt_mean_diff.transpose() * icov;

  double e_x_cov_x =
      std::exp(-param.gauss_d2__half_ * diff_icov.dot(pt_mean_diff));
  double d2_e_x_cov_x = param.gauss_d2_ * e_x_cov_x;
  double d1d2_e_x_cov_x = param.gauss_d1_ * d2_e_x_cov_x;

  if (d2_e_x_cov_x > 1 || d2_e_x_cov_x < 0 || d2_e_x_cov_x != d2_e_x_cov_x) {
    ROS_DEBUG_STREAM("[NDT2D]:Wrong probability density function (PDF)");
    return res.Zero();
  }
  // calculate score (eq. 6.9) [Magnusson 2009]
  res.value_ = -param.gauss_d1_ * e_x_cov_x;

  Eigen::Vector2d hess_part;
  for (long i = 0; i < 3; ++i) {
    // calculate gradient (eq. 6.12) [Magnusson 2009]
    res.gradient_(i) = diff_icov.dot(jacobian.col(i)) * d1d2_e_x_cov_x;
    for (long j = 0; j < 3; ++j) {
      // calculate hessian (eq. 6.13) [Magnusson 2009]
      if (i == 2 && j == 2) {
        hess_part = hess_derivatives2;
      } else {
        hess_part = hess_derivative;
      }
      double jacc_i = diff_icov.dot(jacobian.col(i));
      double jacc_j = diff_icov.dot(jacobian.col(j));
      res.hessian_(i, j) =
          d1d2_e_x_cov_x *
          (-param.gauss_d2_ * jacc_i * jacc_j + diff_icov.dot(hess_part) +
           jacobian.col(j).dot(icov * jacobian.col(i)));
    }
  }
  return res;
}

///////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
double NormalDistributionsTransform2DEx<PointSource, PointTarget, CellType>::
    computeStepLengthMT(const Eigen::Matrix<double, 3, 1> &x,
                        Eigen::Matrix<double, 3, 1> &step_dir, double step_init,
                        double step_max, double step_min,
                        const GridTarget &ndt_grid,
                        const ndt_reg::ScoreAndDerivatives<3, double> &score,
                        const PclSource &input,
                        const ndt_reg::FittingParams &param)
{
  Eigen::Matrix4f transformation;
  PclSource trans_cloud;
  // transformation.setIdentity();
  ndt_reg::ScoreAndDerivatives<3, double> score_vals = score;
  // Set the value of phi(0), Equation 1.3 [More, Thuente 1994]
  double phi_0 = -score_vals.value_;
  // Set the value of phi'(0), Equation 1.3 [More, Thuente 1994]
  double d_phi_0 = -(score_vals.gradient_.dot(step_dir));

  Eigen::Matrix<double, 3, 1> x_t;

  if (d_phi_0 >= 0) {
    // Not a decent direction
    if (d_phi_0 == 0)
      return 0;
    else {
      // Reverse step direction and calculate optimal step.
      d_phi_0 *= -1;
      step_dir *= -1;
    }
  }

  // The Search Algorithm for T(mu) [More, Thuente 1994]

  int max_step_iterations = 10;
  int step_iterations = 0;

  // Sufficient decreace constant, Equation 1.1 [More, Thuete 1994]
  double mu = 1.e-4;
  // Curvature condition constant, Equation 1.2 [More, Thuete 1994]
  double nu = 0.9;

  // Initial endpoints of Interval I,
  double a_l = 0, a_u = 0;

  // Auxiliary function psi is used until I is determined ot be a closed
  // interval, Equation 2.1 [More, Thuente 1994]
  double f_l = ndt_reg::auxilaryFunction_PsiMT(a_l, phi_0, phi_0, d_phi_0, mu);
  double g_l = ndt_reg::auxilaryFunction_dPsiMT(d_phi_0, d_phi_0, mu);

  double f_u = ndt_reg::auxilaryFunction_PsiMT(a_u, phi_0, phi_0, d_phi_0, mu);
  double g_u = ndt_reg::auxilaryFunction_dPsiMT(d_phi_0, d_phi_0, mu);

  // Check used to allow More-Thuente step length calculation to be skipped by
  // making step_min == step_max
  bool interval_converged = (step_max - step_min) > 0, open_interval = true;

  double a_t = step_init;
  a_t = std::min(a_t, step_max);
  a_t = std::max(a_t, step_min);

  x_t = x + step_dir * a_t;

  // transformation = ndt_reg::vecToMat(x_t);

  // Updates score, gradient and hessian.  Hessian calculation is unessisary but
  // testing showed that most step calculations use the
  // initial step suggestion and recalculation the reusable portions of the
  // hessian would intail more computation time.
  score_vals = calcScore(ndt_grid, input, x_t, param);

  // Calculate phi(alpha_t)
  double phi_t = -score_vals.value_;
  // Calculate phi'(alpha_t)
  double d_phi_t = -(score_vals.gradient_.dot(step_dir));

  // Calculate psi(alpha_t)
  double psi_t =
      ndt_reg::auxilaryFunction_PsiMT(a_t, phi_t, phi_0, d_phi_0, mu);
  // Calculate psi'(alpha_t)
  double d_psi_t = ndt_reg::auxilaryFunction_dPsiMT(d_phi_t, d_phi_0, mu);

  // Iterate until max number of iterations, interval convergance or a value
  // satisfies the sufficient decrease, Equation 1.1, and curvature condition,
  // Equation 1.2 [More, Thuente 1994]
  while (!interval_converged && step_iterations < max_step_iterations &&
         !(psi_t <= 0 /*Sufficient Decrease*/ &&
           d_phi_t <= -nu * d_phi_0 /*Curvature Condition*/)) {
    // Use auxilary function if interval I is not closed
    if (open_interval) {
      a_t = ndt_reg::trialValueSelectionMT(a_l, f_l, g_l, a_u, f_u, g_u, a_t,
                                           psi_t, d_psi_t);
    } else {
      a_t = ndt_reg::trialValueSelectionMT(a_l, f_l, g_l, a_u, f_u, g_u, a_t,
                                           phi_t, d_phi_t);
    }

    a_t = std::min(a_t, step_max);
    a_t = std::max(a_t, step_min);

    x_t = x + step_dir * a_t;
    // transformation = ndt_reg::vecToMat(x_t);

    score_vals = calcScore(ndt_grid, input, x_t, param);

    // Calculate phi(alpha_t+)
    phi_t = -score_vals.value_;
    // Calculate phi'(alpha_t+)
    d_phi_t = -(score_vals.gradient_.dot(step_dir));

    // Calculate psi(alpha_t+)
    psi_t = ndt_reg::auxilaryFunction_PsiMT(a_t, phi_t, phi_0, d_phi_0, mu);
    // Calculate psi'(alpha_t+)
    d_psi_t = ndt_reg::auxilaryFunction_dPsiMT(d_phi_t, d_phi_0, mu);

    // Check if I is now a closed interval
    if (open_interval && (psi_t <= 0 && d_psi_t >= 0)) {
      open_interval = false;

      // Converts f_l and g_l from psi to phi
      f_l = f_l + phi_0 - mu * d_phi_0 * a_l;
      g_l = g_l + mu * d_phi_0;

      // Converts f_u and g_u from psi to phi
      f_u = f_u + phi_0 - mu * d_phi_0 * a_u;
      g_u = g_u + mu * d_phi_0;
    }

    if (open_interval) {
      // Update interval end points using Updating Algorithm [More, Thuente
      // 1994]
      interval_converged = ndt_reg::updateIntervalMT(a_l, f_l, g_l, a_u, f_u,
                                                     g_u, a_t, psi_t, d_psi_t);
    } else {
      // Update interval end points using Modified Updating Algorithm [More,
      // Thuente 1994]
      interval_converged = ndt_reg::updateIntervalMT(a_l, f_l, g_l, a_u, f_u,
                                                     g_u, a_t, phi_t, d_phi_t);
    }

    step_iterations++;
  }

  // If inner loop was run then hessian needs to be calculated.
  // Hessian is unnessisary for step length determination but gradients are
  // required
  // so derivative and transform data is stored for the next iteration.
  // if (step_iterations)
  //   computeHessian (hessian, trans_cloud, x_t);

  return (a_t);
}

}  // end of namespace

#endif
