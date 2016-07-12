#ifndef GRAPH_SLAM_UK_D2D_NDT2D
#define GRAPH_SLAM_UK_D2D_NDT2D

#include <graph_slam_uk/ndt/cell_policy2d.h>
#include <graph_slam_uk/ndt/ndt_cell.h>
#include <graph_slam_uk/ndt/ndt_grid2d.h>
#include <graph_slam_uk/utils/covariance_inverse.h>
#include <pcl/common/time.h>
#include <pcl/filters/voxel_grid_covariance.h>
#include <pcl/registration/registration.h>
#include <ros/ros.h>
#include <Eigen/Dense>

namespace pcl
{
namespace d2d_ndt2d
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
}  // end of namespace d2d_ndt2d

template <typename PointSource, typename PointTarget,
          typename CellType = slamuk::NDTCell<slamuk::CellPolicy2d>>
class D2DNormalDistributionsTransform2D
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

  /** \brief Typename of searchable voxel grid containing mean and covariance.
   */
  typedef slamuk::NDTGrid2D<CellType, PointTarget> GridTarget;
  typedef slamuk::NDTGrid2D<CellType, PointSource> GridSource;
  typedef typename GridTarget::ConstPtr GridTargetConstPtr;
  typedef typename GridSource::ConstPtr GridSourceConstPtr;

  // /** \brief Typename of const pointer to searchable voxel grid. */
  // typedef const GridTarget *TargetGridConstPtr;
  // * \brief Typename of const pointer to searchable voxel grid leaf.
  // typedef typename GridTarget::LeafConstPtr TargetGridLeafConstPtr;

public:
  typedef boost::shared_ptr<
      D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>>
      Ptr;
  typedef boost::shared_ptr<const D2DNormalDistributionsTransform2D<
      PointSource, PointTarget, CellType>>
      ConstPtr;
  typedef Eigen::Vector3d VectorTrans;
  /** \brief Constructor.
    * Sets \ref outlier_ratio_ to 0.35, \ref step_size_ to 0.05 and \ref
   * resolution_ to 1.0
    */
  D2DNormalDistributionsTransform2D();

  /** \brief Empty destructor */
  virtual ~D2DNormalDistributionsTransform2D()
  {
  }

  inline void setNumLayers(size_t num)
  {
    layer_count_ = num;
    initCellSizes();
    initParams();
  }

  inline size_t getNumLayers()
  {
    return layer_count_;
  }
  /** \brief Set/change the voxel grid cell size for the finnest grid. Other
   * grids
   * will have higher length of cell (coarse)
    * \param[in] cel_size cell size in meters
    */
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

  void enableMultithreading(unsigned int thread_count)
  {
    threads_ = thread_count;
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

  virtual void setInputTarget(const PclTargetConstPtr &cloud)
  {
    Registration<PointSource, PointTarget>::setInputTarget(cloud);
    GridSource *tmp = new GridSource();
    tmp->setCellSize(cell_sizes_.back());
    tmp->initializeSimple(*cloud);
    target_grid_ = GridTargetConstPtr(tmp);
  }

  virtual void setInputSource(const PclSourceConstPtr &cloud)
  {
    Registration<PointSource, PointTarget>::setInputSource(cloud);
    GridSource *tmp = new GridSource();
    tmp->setCellSize(cell_sizes_.back());
    tmp->initializeSimple(*cloud);
    source_grid_ = GridSourceConstPtr(tmp);
  }

  // using Registration<PointSource, PointTarget>::setInputSource;
  virtual void setInputSource(const GridSourceConstPtr &grid)
  {
    source_grid_ = grid;
    PclSourceConstPtr pcl = grid->getMeans();
    Registration<PointSource, PointTarget>::setInputSource(pcl);
    // setCellSize(source_grid_->getCellSize());
  }

  // using Registration<PointSource, PointTarget>::setInputTarget;
  virtual void setInputTarget(const GridTargetConstPtr &grid)
  {
    target_grid_ = grid;
    PclTargetConstPtr pcl = grid->getMeans();
    Registration<PointSource, PointTarget>::setInputTarget(pcl);
    // setCellSize(target_grid_->getCellSize());
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
  using Registration<PointSource, PointTarget>::target_cloud_updated_;
  using Registration<PointSource, PointTarget>::source_cloud_updated_;

  /** \brief The side length of voxels. */
  std::vector<float> cell_sizes_;

  /** \brief The maximum step length. */
  double step_size_;

  /** \brief The ratio of outliers of points w.r.t. a normal distribution,
   * Equation 6.7 [Magnusson 2009]. */
  double outlier_ratio_;

  /** \brief The normalization constants used fit the point distribution to a
   * normal d  // std::vector<GridTarget> target_cells_;istribution, Equation
   * 6.8 [Magnusson 2009]. */
  std::vector<d2d_ndt2d::FittingParams> params_;

  /** \brief The probability score of the transform applied to the input cloud,
   * Equation 6.9 and 6.10 [Magnusson 2009]. */
  double trans_probability_;

  Eigen::Matrix3d covariance_;
  Eigen::Matrix3d inform_matrix_;

  size_t layer_count_;
  bool target_grid_updated_;
  bool source_grid_updated_;
  GridSourceConstPtr source_grid_;
  GridTargetConstPtr target_grid_;

  unsigned int threads_;

  /** \brief Initialize fitting parameters for normal distrubution in cells for
   * every resolution.
  */
  virtual inline void initParams()
  {
    params_.clear();
    for (size_t i = 0; i < cell_sizes_.size(); ++i) {
      params_.push_back(
          d2d_ndt2d::FittingParams(outlier_ratio_, 1 / cell_sizes_[i]));
    }
  }
  /** \brief Initialize cell_sizes. First grid will have cells of base_size
  * length.base_size
  * Other layers have cell sizes in multiples of 2. e.g.
  * 0.25, 0.5, 1, 2
   * Numer of layers is chosen based on layer_count parameter. Cell sizes
  * are sorted from coarsest grid to finest grid size
  * \return Returns true if made any changes to cell_sizes_
  */
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

  /** \brief Estimate the transformation and returns the transformed source
   * (input) as output.
    * \param[out] output the resultant input transfomed point cloud dataset
    * \param[in] guess the initial gross estimation of the transformation
    */
  virtual void computeTransformation(PclSource &output,
                                     const Eigen::Matrix4f &guess);

  /** \brief Estimate the transformation and returns the transformed source
   * (input) as output.
    * \param[out] output the resultant input transfomed point cloud dataset
    */
  virtual void computeTransformation(PclSource &output)
  {
    computeTransformation(output, Eigen::Matrix4f::Identity());
  }

  virtual bool computeSingleGrid(const GridSource &source_grid,
                                 const Eigen::Matrix4f &guess,
                                 const GridTarget &target_grid,
                                 const d2d_ndt2d::FittingParams &param,
                                 Eigen::Matrix4f &trans);

  virtual d2d_ndt2d::ScoreAndDerivatives<3, double>
  calcScore(const d2d_ndt2d::FittingParams &param, const GridSource &sourceNDT,
            const Eigen::Vector3d &trans, const GridTarget &targetNDT,
            bool calc_hessian);

  virtual void computeDerivatives(const Eigen::Vector3d &x,
                                  const Eigen::Matrix3d &cov,
                                  d2d_ndt2d::JacobianHessianDerivatives &data,
                                  bool calc_hessian);

  virtual d2d_ndt2d::ScoreAndDerivatives<3, double>
  calcSourceCellScore(const Eigen::Vector3d &mean_source,
                      const Eigen::Matrix3d &cov_source, const CellType *cell_t,
                      const d2d_ndt2d::JacobianHessianDerivatives &deriv,
                      const d2d_ndt2d::FittingParams &param, bool calc_hessian);

  // linear search methods//////////////////////////////////////////
  virtual double computeStepLengthMT(
      const Eigen::Matrix<double, 3, 1> &x,
      Eigen::Matrix<double, 3, 1> &step_dir, double step_init, double step_max,
      double step_min, const GridSource &source_grid,
      const d2d_ndt2d::ScoreAndDerivatives<3, double> &score,
      const GridTarget &target_grid, const d2d_ndt2d::FittingParams &param);

  bool updateIntervalMT(double &a_l, double &f_l, double &g_l, double &a_u,
                        double &f_u, double &g_u, double a_t, double f_t,
                        double g_t) const;

  double trialValueSelectionMT(double a_l, double f_l, double g_l, double a_u,
                               double f_u, double g_u, double a_t, double f_t,
                               double g_t) const;
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
  inline double auxilaryFunction_PsiMT(double a, double f_a, double f_0,
                                       double g_0, double mu = 1.e-4) const
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
  inline double auxilaryFunction_dPsiMT(double g_a, double g_0,
                                        double mu = 1.e-4) const
  {
    return (g_a - mu * g_0);
  }
  template <typename T = float>
  Eigen::Matrix<T, 4, 4> vecToMat(const Eigen::Vector3d &trans) const;
  template <typename T = float>
  Eigen::Vector3d matToVec(const Eigen::Matrix<T, 4, 4> &trans) const;
};
////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
D2DNormalDistributionsTransform2D<PointSource, PointTarget,
                                  CellType>::D2DNormalDistributionsTransform2D()
  : step_size_(0.01)
  , outlier_ratio_(0.99)
  , trans_probability_()
  , layer_count_(4)
  , target_grid_updated_(false)
  , source_grid_updated_(false)
  , threads_(2)
{
  nr_iterations_ = 0;
  max_iterations_ = 35;
  transformation_epsilon_ = 0.1;
  converged_ = false;
  target_cloud_updated_ = false;
  source_cloud_updated_ = false;
  initCellSizes(0.25);
  initParams();
}

////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
void D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::
    computeTransformation(PclSource &output, const Eigen::Matrix4f &guess)
{
  ROS_DEBUG_STREAM("[D2D_NDT2D]: guess:" << matToVec(guess).transpose());
  Eigen::Matrix4f trans = guess;
  converged_ = false;
  if (!target_grid_) {
    ROS_ERROR_STREAM("[D2D_NDT2D]: target cloud or grid not set");
    return;
  }
  if (!source_grid_) {
    ROS_ERROR_STREAM("[D2D_NDT2D]: source cloud or grid not set");
    return;
  }
  for (size_t i = 0; i < layer_count_; ++i) {
    if (!computeSingleGrid(source_grid_->createCoarserGrid(cell_sizes_[i]),
                           trans,
                           target_grid_->createCoarserGrid(cell_sizes_[i]),
                           params_[i], trans)) {
      converged_ = false;
      return;
    }
  }
  ROS_DEBUG_STREAM("[D2D_NDT2D]: final trans:" << matToVec(trans).transpose());
  transformPointCloud(*input_, output, trans);
  final_transformation_ = trans;
  converged_ = true;
}
////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
bool D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::
    computeSingleGrid(const GridSource &source_grid,
                      const Eigen::Matrix4f &guess,
                      const GridTarget &target_grid,
                      const d2d_ndt2d::FittingParams &param,
                      Eigen::Matrix4f &trans)
{
  nr_iterations_ = 0;
  converged_ = false;
  // Initialise final transformation to the guessed one
  final_transformation_ = guess;
  previous_transformation_ = final_transformation_;
  // variables needed for calculation
  Eigen::Vector3d xytheta_p = matToVec(guess);
  Eigen::Vector3d delta_xytheta_p;
  Eigen::Matrix4f p;
  double delta_p_norm;
  d2d_ndt2d::ScoreAndDerivatives<3, double> score;
  while (!converged_) {
    score = calcScore(param, source_grid, xytheta_p, target_grid, true);

    // Solve for decent direction using newton method
    Eigen::JacobiSVD<Eigen::Matrix3d> sv(
        score.hessian_, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // Negative for maximization as opposed to minimization
    delta_xytheta_p = sv.solve(-score.gradient_);

    delta_p_norm = delta_xytheta_p.norm();
    if (delta_p_norm == 0 || delta_p_norm != delta_p_norm) {
      trans_probability_ =
          score.value_ / static_cast<double>(input_->points.size());
      converged_ = delta_p_norm == delta_p_norm;
      covariance_.setIdentity();
      inform_matrix_.setIdentity();
      ROS_ERROR_STREAM("[D2D_NDT2D]:Not enough overlap. Probability: "
                       << trans_probability_);
      return false;
    }
    delta_xytheta_p.normalize();
    // Calculate step length with guarnteed sufficient decrease [More, Thuente
    // 1994]
    delta_p_norm = computeStepLengthMT(xytheta_p, delta_xytheta_p, delta_p_norm,
                                       step_size_, transformation_epsilon_ / 2,
                                       source_grid, score, target_grid, param);
    // delta_p_norm = 0.8
    delta_xytheta_p *= delta_p_norm;
    xytheta_p += delta_xytheta_p;
    p = vecToMat(xytheta_p);

    ++nr_iterations_;
    previous_transformation_ = transformation_;
    transformation_ = p;
    trans_probability_ =
        score.value_ / static_cast<double>(input_->points.size());
    // ROS_DEBUG_STREAM("[D2D_NDT2D]: Step: "
    //                  << delta_p_norm
    //                  << " Delta: " << delta_xytheta_p.transpose()
    //                  << " Score: " << score.value_
    //                  << " probability of match: " << trans_probability_
    //                  << " current transformation: \n" << xytheta_p);
    // convergence testing
    if (nr_iterations_ >= max_iterations_ ||
        (nr_iterations_ &&
         (std::abs(delta_p_norm) < transformation_epsilon_))) {
      converged_ = true;
    }

    // // Update Visualizer (untested)
    // if (update_visualizer_ != 0) {
    //   transformPointCloud(*input_, output, transformation_);
    //   update_visualizer_(output, std::vector<int>(), *target_,
    //                      std::vector<int>());
    // }
  }
  trans = p;
  // final_transformation_ = transformation_;
  covariance_ = score.hessian_;
  inform_matrix_ = score.hessian_.inverse();
  return true;
}
////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
d2d_ndt2d::ScoreAndDerivatives<3, double>
D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::calcScore(
    const d2d_ndt2d::FittingParams &param, const GridSource &sourceNDT,
    const Eigen::Vector3d &trans, const GridTarget &targetNDT,
    bool calc_hessian)
{
  typedef typename GridSource::CellPtrVector SourceVec;
  typedef typename GridTarget::CellPtrVector TargetVec;
  typedef d2d_ndt2d::ScoreAndDerivatives<3, double> ReturnVals;

  ReturnVals res;

  Eigen::Transform<double, 3, Eigen::Affine, Eigen::ColMajor> trans_mat;
  trans_mat.matrix() = vecToMat<double>(trans);

  SourceVec source_cells = sourceNDT.getGaussianCells();
  // std::cout << "ndt cells: " << source_cells.size() << std::endl;
  std::vector<ReturnVals> omp_ret;
  omp_ret.resize(threads_);
  for (size_t i = 0; i < threads_; ++i) {
    omp_ret.push_back(ReturnVals());
  }

#pragma omp parallel num_threads(threads_)
  {
#pragma omp for
    for (size_t cell_id = 0; cell_id < source_cells.size(); ++cell_id) {
      int thread_id = omp_get_thread_num();
      Eigen::Vector3d mean_source;
      Eigen::Matrix3d cov_source;
      ReturnVals local_ret;
      // TRANSFORMATION OF SOURCE GRID
      Eigen::Vector3d cell_mean;
      Eigen::Matrix3d cell_cov;
      cell_mean = source_cells[cell_id]->getMean();
      cell_cov = source_cells[cell_id]->getCov();
      mean_source = trans_mat * cell_mean;
      cov_source =
          trans_mat.rotation() * cell_cov * trans_mat.rotation().transpose();
      // compute derivatives of score function
      d2d_ndt2d::JacobianHessianDerivatives partial_derivatives;
      computeDerivatives(mean_source, cov_source, partial_derivatives,
                         calc_hessian);
      TargetVec neighbors =
          targetNDT.getKNearestNeighbors(mean_source.head(2), 2);
      // std::cout << "neighbors: " << neighbors.size() << std::endl;
      for (size_t i = 0; i < neighbors.size(); ++i) {
        local_ret +=
            calcSourceCellScore(mean_source, cov_source, neighbors[i],
                                partial_derivatives, param, calc_hessian);
        // std::cout << "\ncovar:\n" << neighbors[i]->getCov() << std::endl
        //           << neighbors[i]->getICov()
        //           << "\nmean: " << neighbors[i]->getMean().transpose()
        //           << "\n points: " << neighbors[i]->points()
        //           << "\n----------\n";
      }
      omp_ret[thread_id] += local_ret;
    }
  }
  for (size_t i = 0; i < omp_ret.size(); ++i) {
    res += omp_ret[i];
  }
  return res;
}

////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
void D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::
    computeDerivatives(const Eigen::Vector3d &x, const Eigen::Matrix3d &cov,
                       d2d_ndt2d::JacobianHessianDerivatives &data,
                       bool calc_hessian)
{
  data.setZero();
  data.Jest.block<2, 2>(0, 0).setIdentity();
  data.Jest(0, 2) = -x(1);
  data.Jest(1, 2) = x(0);

  //_Zest
  data.Zest.block<3, 3>(0, 6) << -2 * cov(0, 1), -cov(1, 1) + cov(0, 0),
      -cov(1, 2), -cov(1, 1) + cov(0, 0), 2 * cov(0, 1), cov(0, 2), -cov(1, 2),
      cov(0, 2), 0;

  if (calc_hessian) {
    data.Hest.block<3, 1>(6, 2) << -x(0), -x(1), 0;
    data.ZHest.block<3, 3>(6, 6) << 2 * cov(1, 1) - 2 * cov(0, 0),
        -4 * cov(0, 1), -cov(0, 2), -4 * cov(0, 1),
        2 * cov(0, 0) - 2 * cov(1, 1), -cov(1, 2), -cov(0, 2), -cov(1, 2), 0;
  }
}
////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
d2d_ndt2d::ScoreAndDerivatives<3, double>
D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::
    calcSourceCellScore(const Eigen::Vector3d &mean_source,
                        const Eigen::Matrix3d &cov_source,
                        const CellType *cell_t,
                        const d2d_ndt2d::JacobianHessianDerivatives &deriv,
                        const d2d_ndt2d::FittingParams &param,
                        bool calc_hessian)
{
  d2d_ndt2d::ScoreAndDerivatives<3, double> res;

  // declaration
  Eigen::Vector3d diff_mean;
  Eigen::Matrix3d cov_sum, icov;
  cov_sum.setIdentity();
  double det = 0;
  bool exists = false;
  double dist;
  // vars for gradient
  Eigen::Matrix<double, 3, 1> xtBJ, xtBZBx, Q;
  // vars for hessian
  Eigen::Matrix<double, 3, 3> xtBZBJ, xtBH, xtBZBZBx, xtBZhBx;
  Eigen::Matrix<double, 1, 3> TMP1, xtB;

  xtBJ.setZero();
  xtBZBx.setZero();
  Q.setZero();
  xtBZBJ.setZero();
  xtBH.setZero();
  xtBZBZBx.setZero();
  xtBZhBx.setZero();
  TMP1.setZero();
  xtB.setZero();

  diff_mean = (mean_source - cell_t->getMean());
  cov_sum = (cell_t->getCov() + cov_source);

  cov_sum.computeInverseAndDetWithCheck(icov, det, exists);
  if (!exists) {
    return res.Zero();
  }
  dist = (diff_mean).dot(icov * (diff_mean));
  if (dist * 0 != 0) {
    return res.Zero();
  }
  res.value_ = -param.gauss_d1_ * std::exp(-param.gauss_d2__half_ * dist);

  xtB = diff_mean.transpose() * icov;
  xtBJ = xtB * deriv.Jest;

  TMP1 = xtB * deriv.Zest.block<3, 3>(0, 6) * icov;
  xtBZBx(2) = TMP1 * diff_mean;
  if (calc_hessian) {
    xtBZBJ.col(2) = (TMP1 * deriv.Jest).transpose();
    for (unsigned int j = 0; j < 3; j++) {
      xtBH(2, j) = xtB * deriv.Hest.block<3, 1>(6, j);
      xtBZBZBx(2, j) =
          TMP1 * deriv.Zest.block<3, 3>(0, 3 * j) * icov * diff_mean;
      xtBZhBx(2, j) =
          xtB * deriv.ZHest.block<3, 3>(6, 3 * j) * icov * diff_mean;
    }
  }
  Q = 2 * xtBJ - xtBZBx;
  double factor = -(param.gauss_d2__half_) * res.value_;
  res.gradient_ += Q * factor;

  if (calc_hessian) {
    res.hessian_ += factor * (2 * deriv.Jest.transpose() * icov * deriv.Jest +
                              2 * xtBH - xtBZhBx - 2 * xtBZBJ.transpose() -
                              2 * xtBZBJ + xtBZBZBx + xtBZBZBx.transpose() -
                              param.gauss_d2__half_ * Q * Q.transpose());
  }

  return res;
}

///////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
double D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::
    computeStepLengthMT(const Eigen::Matrix<double, 3, 1> &x,
                        Eigen::Matrix<double, 3, 1> &step_dir, double step_init,
                        double step_max, double step_min,
                        const GridSource &source_grid,
                        const d2d_ndt2d::ScoreAndDerivatives<3, double> &score,
                        const GridTarget &target_grid,
                        const d2d_ndt2d::FittingParams &param)
{
  Eigen::Matrix4f transformation;
  PclSource trans_cloud;
  // transformation.setIdentity();
  d2d_ndt2d::ScoreAndDerivatives<3, double> score_vals = score;
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
  double f_l = auxilaryFunction_PsiMT(a_l, phi_0, phi_0, d_phi_0, mu);
  double g_l = auxilaryFunction_dPsiMT(d_phi_0, d_phi_0, mu);

  double f_u = auxilaryFunction_PsiMT(a_u, phi_0, phi_0, d_phi_0, mu);
  double g_u = auxilaryFunction_dPsiMT(d_phi_0, d_phi_0, mu);

  // Check used to allow More-Thuente step length calculation to be skipped by
  // making step_min == step_max
  bool interval_converged = (step_max - step_min) > 0, open_interval = true;

  double a_t = step_init;
  a_t = std::min(a_t, step_max);
  a_t = std::max(a_t, step_min);

  x_t = x + step_dir * a_t;

  // transformation = vecToMat(x_t);

  // Updates score, gradient and hessian.  Hessian calculation is unessisary
  // but
  // testing showed that most step calculations use the
  // initial step suggestion and recalculation the reusable portions of the
  // hessian would intail more computation time.
  score_vals = calcScore(param, source_grid, x_t, target_grid, false);

  // Calculate phi(alpha_t)
  double phi_t = -score_vals.value_;
  // Calculate phi'(alpha_t)
  double d_phi_t = -(score_vals.gradient_.dot(step_dir));

  // Calculate psi(alpha_t)
  double psi_t = auxilaryFunction_PsiMT(a_t, phi_t, phi_0, d_phi_0, mu);
  // Calculate psi'(alpha_t)
  double d_psi_t = auxilaryFunction_dPsiMT(d_phi_t, d_phi_0, mu);

  // Iterate until max number of iterations, interval convergance or a value
  // satisfies the sufficient decrease, Equation 1.1, and curvature condition,
  // Equation 1.2 [More, Thuente 1994]
  while (!interval_converged && step_iterations < max_step_iterations &&
         !(psi_t <= 0 /*Sufficient Decrease*/ &&
           d_phi_t <= -nu * d_phi_0 /*Curvature Condition*/)) {
    // Use auxilary function if interval I is not closed
    if (open_interval) {
      a_t = trialValueSelectionMT(a_l, f_l, g_l, a_u, f_u, g_u, a_t, psi_t,
                                  d_psi_t);
    } else {
      a_t = trialValueSelectionMT(a_l, f_l, g_l, a_u, f_u, g_u, a_t, phi_t,
                                  d_phi_t);
    }

    a_t = std::min(a_t, step_max);
    a_t = std::max(a_t, step_min);

    x_t = x + step_dir * a_t;
    // transformation = vecToMat(x_t);

    score_vals = calcScore(param, source_grid, x_t, target_grid, false);

    // Calculate phi(alpha_t+)
    phi_t = -score_vals.value_;
    // Calculate phi'(alpha_t+)
    d_phi_t = -(score_vals.gradient_.dot(step_dir));

    // Calculate psi(alpha_t+)
    psi_t = auxilaryFunction_PsiMT(a_t, phi_t, phi_0, d_phi_0, mu);
    // Calculate psi'(alpha_t+)
    d_psi_t = auxilaryFunction_dPsiMT(d_phi_t, d_phi_0, mu);

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
      interval_converged =
          updateIntervalMT(a_l, f_l, g_l, a_u, f_u, g_u, a_t, psi_t, d_psi_t);
    } else {
      // Update interval end points using Modified Updating Algorithm [More,
      // Thuente 1994]
      interval_converged =
          updateIntervalMT(a_l, f_l, g_l, a_u, f_u, g_u, a_t, phi_t, d_phi_t);
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
bool D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::
    updateIntervalMT(double &a_l, double &f_l, double &g_l, double &a_u,
                     double &f_u, double &g_u, double a_t, double f_t,
                     double g_t) const
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
double D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::
    trialValueSelectionMT(double a_l, double f_l, double g_l, double a_u,
                          double f_u, double g_u, double a_t, double f_t,
                          double g_t) const
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

template <typename PointSource, typename PointTarget, typename CellType>
template <typename T>
Eigen::Matrix<T, 4, 4>
D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::vecToMat(
    const Eigen::Vector3d &trans) const
{
  Eigen::Matrix<T, 4, 4> trans_mat = Eigen::Matrix<T, 4, 4>::Identity();

  trans_mat.block(0, 0, 3, 3).matrix() =
      Eigen::Matrix<T, 3, 3>(Eigen::AngleAxis<T>(
          static_cast<T>(trans(2)), Eigen::Matrix<T, 3, 1>::UnitZ()));

  trans_mat.block(0, 3, 3, 1).matrix() = Eigen::Matrix<T, 3, 1>(
      static_cast<T>(trans(0)), static_cast<T>(trans(1)), 0.0);

  return trans_mat;
}

template <typename PointSource, typename PointTarget, typename CellType>
template <typename T>
Eigen::Vector3d
D2DNormalDistributionsTransform2D<PointSource, PointTarget, CellType>::matToVec(
    const Eigen::Matrix<T, 4, 4> &trans) const
{
  Eigen::Vector3d vec;
  Eigen::Transform<T, 3, Eigen::Affine, Eigen::ColMajor> trans_mat(trans);
  Eigen::Matrix<T, 3, 1> translation = trans_mat.translation();
  vec << translation(0), translation(1),
      std::atan2(trans_mat.rotation()(1, 0), trans_mat.rotation()(0, 0));
  return vec;
}

}  // end of pcl namespace

#endif