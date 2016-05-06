#ifndef NDT_SCANMATCHING_NDT2D
#define NDT_SCANMATCHING_NDT2D

#include <ros/ros.h>
#include <pcl/registration/registration.h>
#include <pcl/filters/voxel_grid_covariance.h>
#include <Eigen/Dense>

namespace pcl
{

  namespace ndt2d
  {
      /** \brief Class to store vector value and first and second derivatives
        * (grad vector and hessian matrix), so they can be returned easily from
        * functions
        */
      template<unsigned N=3, typename T=double>
      struct ScoreAndDerivatives
      {
        ScoreAndDerivatives () : hessian (), gradient (), value () {}

        Eigen::Matrix<T, N, N> hessian;
        Eigen::Matrix<T, N, 1>    gradient;
        T value;

        static ScoreAndDerivatives<N,T>
        Zero ()
        {
          ScoreAndDerivatives<N,T> r;
          r.hessian = Eigen::Matrix<T, N, N>::Zero ();
          r.gradient    = Eigen::Matrix<T, N, 1>::Zero ();
          r.value   = 0;
          return r;
        }

        ScoreAndDerivatives<N,T>&
        operator+= (ScoreAndDerivatives<N,T> const& r)
        {
          hessian += r.hessian;
          gradient    += r.gradient;
          value   += r.value;
          return *this;
        }
      };
      template<unsigned N, typename T>
      ScoreAndDerivatives<N,T> operator+(const ScoreAndDerivatives<N,T> & lhs,const ScoreAndDerivatives<N,T> & rhs)
      {
          ScoreAndDerivatives<N,T> ret;
          ret+=lhs;
          ret+=rhs;
          return ret;
      }
  } // end of namespace ndt2d

template <typename PointSource, typename PointTarget>
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

  /** \brief Typename of searchable voxel grid containing mean and covariance.
   */
  typedef VoxelGridCovariance<PointTarget> TargetGrid;
  /** \brief Typename of pointer to searchable voxel grid. */
  typedef TargetGrid* TargetGridPtr;
  /** \brief Typename of const pointer to searchable voxel grid. */
  typedef const TargetGrid *TargetGridConstPtr;
  /** \brief Typename of const pointer to searchable voxel grid leaf. */
  typedef typename TargetGrid::LeafConstPtr TargetGridLeafConstPtr;

public:
  typedef boost::shared_ptr<
      NormalDistributionsTransform2DEx<PointSource, PointTarget> > Ptr;
  typedef boost::shared_ptr<
      const NormalDistributionsTransform2DEx<PointSource, PointTarget> > ConstPtr;
  typedef Eigen::Vector3d VectorTrans;
  /** \brief Constructor.
    * Sets \ref outlier_ratio_ to 0.35, \ref step_size_ to 0.05 and \ref
   * resolution_ to 1.0
    */
  NormalDistributionsTransform2DEx();

  /** \brief Empty destructor */
  virtual ~NormalDistributionsTransform2DEx()
  {
  }

  /** \brief Provide a pointer to the input target (e.g., the point cloud that
   * we want to align the input source to).
    * \param[in] cloud the input point cloud target
    */
  inline void setInputTarget(const PclSourceConstPtr &cloud)
  {
    Registration<PointSource, PointTarget>::setInputTarget(cloud);
    initGrid();
  }

  /** \brief Set/change the voxel grid resolution.
    * \param[in] resolution side length of voxels
    */
  inline void setResolution(float resolution)
  {
    // Prevents unnessary voxel initiations
    if (resolution_ != resolution) {
      resolution_ = resolution;
      if (input_)
        initGrid();
      initParams();
    }
  }

  /** \brief Get voxel grid resolution.
    * \return side length of voxels
    */
  inline float getResolution() const
  {
    return (resolution_);
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

  /** \brief The voxel grid generated from target cloud containing point means and covariances. */
  TargetGrid target_cells_;
  /** \brief The side length of voxels. */
  float resolution_;

  /** \brief The maximum step length. */
  double step_size_;

  /** \brief The ratio of outliers of points w.r.t. a normal distribution,
   * Equation 6.7 [Magnusson 2009]. */
  double outlier_ratio_;

  /** \brief The normalization constants used fit the point distribution to a
   * normal distribution, Equation 6.8 [Magnusson 2009]. */
  double gauss_d1_, gauss_d2_, gauss_d2__half_;

  /** \brief The probability score of the transform applied to the input cloud,
   * Equation 6.9 and 6.10 [Magnusson 2009]. */
  double trans_probability_;

  Eigen::Matrix3d covariance_;
  Eigen::Matrix3d inform_matrix_;

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
                                     const Eigen::Matrix4f &guess)
  {
    computeSingleGrid(output,target_cells_,guess);
  }

  virtual bool computeSingleGrid(PclSource & output,TargetGrid & target_grid,const Eigen::Matrix4f &guess);
  /** \brief Initiate covariance voxel structure. */
  virtual void inline initGrid()
  {
    target_cells_.setLeafSize(resolution_, resolution_, resolution_);
    target_cells_.setInputCloud(target_);
    // Initiate voxel structure.
    target_cells_.filter(true);
  }

  virtual inline void initParams(){
    double gauss_c1, gauss_c2, gauss_d3;
    // Initializes the guassian fitting parameters (eq. 6.8) [Magnusson 2009]
    gauss_c1 = 10.0 * (1 - outlier_ratio_);
    gauss_c2 = outlier_ratio_ / pow (resolution_, 2);
    gauss_d3 = -log (gauss_c2);
    gauss_d1_ = -log ( gauss_c1 + gauss_c2 ) - gauss_d3;
    gauss_d2_ = -2 * log ((-log ( gauss_c1 * exp ( -0.5 ) + gauss_c2 ) - gauss_d3) / gauss_d1_);
    gauss_d2__half_ = gauss_d2_ /2;
  }
  /** \brief Estimate score, hessian and gradient for selected source point cloud based on given transformation.
    * \param[in] ndt_grid NDT target grid used as registration target.
    * \param[in] source Source point cloud used to register agains NDT target grid
    * \param[in] trans Transformation used for projecting points of source to grid cells of NDT grid
    * \return Estimated parameters score, gradient, hessian.
    */
  virtual ndt2d::ScoreAndDerivatives<3, double>
  calcScore(TargetGrid &ndt_grid,
            const PclSource &source, const Eigen::Vector3d &trans);
  /** \brief Estimate score, hessian and gradient for selected point.
    * \param[in] cell NDT target grid cell for selected point.
    * \param[in] trans_pt Transformed source cloud point.
    * \param[in] pt Source cloud point befor transformation.
    * \param[in] cos_theta Cosine of transformation angle.
    * \param[in] sin_theta Sine of transformation angle.
    * \return Estimated parameters score, gradient, hessian.
    */
  virtual ndt2d::ScoreAndDerivatives<3, double> calcPointScore(
      const PointSource &pt, const typename TargetGrid::LeafConstPtr &cell,
      const PointSource &trans_pt, double cos_theta, double sin_theta) const;

  // linear search methods
    double computeStepLengthMT(const Eigen::Matrix<double, 3, 1> &x,
                               Eigen::Matrix<double, 3, 1> &step_dir,
                               double step_init, double step_max, double step_min,
                               TargetGrid &ndt_grid,
                               const ndt2d::ScoreAndDerivatives<3, double> &score,
                               const PclSource &input);
    bool updateIntervalMT(double &a_l, double &f_l, double &g_l,
                                   double &a_u, double &f_u, double &g_u,
                                   double a_t, double f_t, double g_t) const;
    double trialValueSelectionMT(double a_l, double f_l, double g_l, double a_u,
                                 double f_u, double g_u, double a_t, double f_t,
                                 double g_t) const;
    /** \brief Auxilary function used to determin endpoints of More-Thuente interval.
      * \note \f$ \psi(\alpha) \f$ in Equation 1.6 (Moore, Thuente 1994)
      * \param[in] a the step length, \f$ \alpha \f$ in More-Thuente (1994)
      * \param[in] f_a function value at step length a, \f$ \phi(\alpha) \f$ in More-Thuente (1994)
      * \param[in] f_0 initial function value, \f$ \phi(0) \f$ in Moore-Thuente (1994)
      * \param[in] g_0 initial function gradiant, \f$ \phi'(0) \f$ in More-Thuente (1994)
      * \param[in] mu the step length, constant \f$ \mu \f$ in Equation 1.1 [More, Thuente 1994]
      * \return sufficent decrease value
      */
    inline double
    auxilaryFunction_PsiMT (double a, double f_a, double f_0, double g_0, double mu = 1.e-4)const
    {
      return (f_a - f_0 - mu * g_0 * a);
    }

    /** \brief Auxilary function derivative used to determin endpoints of More-Thuente interval.
      * \note \f$ \psi'(\alpha) \f$, derivative of Equation 1.6 (Moore, Thuente 1994)
      * \param[in] g_a function gradient at step length a, \f$ \phi'(\alpha) \f$ in More-Thuente (1994)
      * \param[in] g_0 initial function gradiant, \f$ \phi'(0) \f$ in More-Thuente (1994)
      * \param[in] mu the step length, constant \f$ \mu \f$ in Equation 1.1 [More, Thuente 1994]
      * \return sufficent decrease derivative
      */
    inline double
    auxilaryFunction_dPsiMT (double g_a, double g_0, double mu = 1.e-4) const
    {
      return (g_a - mu * g_0);
    }

  Eigen::Matrix4f vecToMat(const Eigen::Vector3d & trans) const;
  Eigen::Vector3d matToVec(const Eigen::Matrix4f & trans) const;
};

template <typename PointSource, typename PointTarget>
NormalDistributionsTransform2DEx<PointSource,PointTarget>::NormalDistributionsTransform2DEx()
  : target_cells_()
  , resolution_ (1.0f)
  , step_size_ (0.1)
  , outlier_ratio_ (0.55)
  , gauss_d1_ ()
  , gauss_d2_ ()
  , gauss_d2__half_()
  , trans_probability_ ()
{
  nr_iterations_ = 0;
  max_iterations_ = 35;
  transformation_epsilon_ = 0.1;
  converged_ = false;
  initParams();
}

template <typename PointSource, typename PointTarget>
bool NormalDistributionsTransform2DEx<PointSource,PointTarget>::computeSingleGrid(
    PclSource &output, TargetGrid &target_grid,
    const Eigen::Matrix4f &guess)
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
    ndt2d::ScoreAndDerivatives<> score;
    while (!converged_) {

      score = calcScore(target_grid,*input_,xytheta_p);
      // Solve for decent direction using newton method, line 23 in Algorithm 2 [Magnusson 2009]
      Eigen::JacobiSVD<Eigen::Matrix3d > sv (score.hessian, Eigen::ComputeFullU | Eigen::ComputeFullV);
      // Negative for maximization as opposed to minimization
      delta_xytheta_p = sv.solve (-score.gradient);

      //Calculate step length with guarnteed sufficient decrease [More, Thuente 1994]
      delta_p_norm = delta_xytheta_p.norm ();
      if (delta_p_norm == 0 || delta_p_norm != delta_p_norm)
      {
        trans_probability_ = score.value/ static_cast<double> (input_->points.size ());
        converged_ = delta_p_norm == delta_p_norm;
        covariance_.setIdentity();
        inform_matrix_.setIdentity();
        ROS_ERROR_STREAM("Not enough overlap. Probability: "<<trans_probability_);
        return false;
      }
      delta_xytheta_p.normalize();
      delta_p_norm = computeStepLengthMT (xytheta_p, delta_xytheta_p, delta_p_norm, step_size_, transformation_epsilon_ / 2, target_grid,score,output);
      //delta_p_norm = 0.8
      delta_xytheta_p *= delta_p_norm;
      xytheta_p+= delta_xytheta_p;
      p = vecToMat(xytheta_p);

      // convergence testing
      if (nr_iterations_ > max_iterations_ ||
          (nr_iterations_ &&
           (std::abs(delta_p_norm) < transformation_epsilon_))) {
        converged_ = true;
      }

      ++nr_iterations_;
      previous_transformation_ = transformation_;
      transformation_ = p;
      trans_probability_ = score.value / static_cast<double> (input_->points.size ());
      ROS_DEBUG_STREAM("[NDT2D]: Step: " << delta_p_norm << " Delta: "
                                         << delta_xytheta_p.transpose()
                                         << " Score: " << score.value
          << "probability of match: "<< trans_probability_);
      // Update Visualizer (untested)
      if (update_visualizer_ != 0) {
        transformPointCloud(*input_, output, transformation_);
        update_visualizer_(output, std::vector<int>(), *target_,
                           std::vector<int>());
      }
    }
    final_transformation_ = transformation_;
    transformPointCloud(*input_, output, transformation_);
    covariance_ = score.hessian;
    inform_matrix_ = score.hessian.inverse();
    return true;
}

/////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget>
ndt2d::ScoreAndDerivatives<3, double>
NormalDistributionsTransform2DEx<PointSource, PointTarget>::calcScore(
    TargetGrid &ndt_grid, const PclSource &source,
    const Eigen::Vector3d &trans)
{
  PclSource trans_source;
  ndt2d::ScoreAndDerivatives<3, double> res;
  double cos_theta, sin_theta;
  // prepare angles for jaccobian and hessian derivatives usage
  // simplify small angles
  if (std::abs(trans(2)) < 10e-5)
  {
    cos_theta = 1.0;
    sin_theta = 0.0;
  }
  else
  {
    cos_theta = std::cos (trans(2));
    sin_theta = std::sin (trans(2));
  }
  res = ndt2d::ScoreAndDerivatives<3, double>::Zero ();
  transformPointCloud (source, trans_source, vecToMat(trans));
  // calculate score for each point
  typedef std::vector<typename TargetGrid::LeafConstPtr> NeighborsVector;
  for (size_t i = 0; i < trans_source.size (); ++i){
    NeighborsVector neighborhood;
    std::vector<float> distances;
    ndt_grid.radiusSearch (trans_source[i], resolution_, neighborhood, distances);
    size_t count  =0;
    for (typename NeighborsVector::iterator neighborhood_it = neighborhood.begin();
         neighborhood_it != neighborhood.end(); neighborhood_it++)
    {
      count++;
      res += calcPointScore(source[i],*neighborhood_it,trans_source[i], cos_theta, sin_theta);
    }
    //ROS_DEBUG_STREAM("Number of neighbours for point: "<<count);
  }
  return res;
}

////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget>
ndt2d::ScoreAndDerivatives<3, double>
NormalDistributionsTransform2DEx<PointSource, PointTarget>::calcPointScore(
    const PointSource &pt,const typename TargetGrid::LeafConstPtr &cell,  const PointSource & trans_pt,double cos_theta,
    double sin_theta) const
{
  ndt2d::ScoreAndDerivatives<3, double> res;
  Eigen::Matrix<double, 2, 3> jacobian;
  Eigen::Vector2d hess_derivative,hess_derivatives2;
  double x = pt.x;
  double y = pt.y;
  double co = cos_theta;
  double si = sin_theta;
  jacobian << 1, 0, -x * si - y * co,
              0, 1,  x * co - y * si;
  hess_derivative << -x * co + y * si,
                     -x * si - y * co;
  hess_derivatives2 <<0,0;
  Eigen::Matrix2d icov = cell->getInverseCov().block(0,0,2,2);
  Eigen::Vector2d mean = cell->getMean().head(2);
  Eigen::Vector2d trans_point(trans_pt.x,trans_pt.y);
  Eigen::Vector2d pt_mean_diff = trans_point - mean;
  Eigen::Vector2d diff_icov = pt_mean_diff.transpose() * icov;
  double e_x_cov_x = std::exp(-gauss_d2__half_ * diff_icov.dot(pt_mean_diff));
  double d2_e_x_cov_x = gauss_d2_ * e_x_cov_x;
  double d1d2_e_x_cov_x = gauss_d1_ * d2_e_x_cov_x;

  if(d2_e_x_cov_x > 1 || d2_e_x_cov_x <0 || d2_e_x_cov_x != d2_e_x_cov_x)
  {
    ROS_DEBUG_STREAM("Wrong probability density function (PDF)");
    return res.Zero();
  }
  //calculate score (eq. 6.9) [Magnusson 2009]
  res.value = -gauss_d1_ * e_x_cov_x;

  Eigen::Vector2d hess_part;
  for(long i = 0; i< 3;++i){
       // calculate gradient (eq. 6.12) [Magnusson 2009]
      res.gradient(i) = diff_icov.dot(jacobian.col(i)) * d1d2_e_x_cov_x;
    for(long j = 0; j< 3 ;++j){
      //calculate hessian (eq. 6.13) [Magnusson 2009]
      if(i == 2 && j == 2){
        hess_part = hess_derivatives2;
      }else{
        hess_part = hess_derivative;
      }
      double jacc_i = diff_icov.dot(jacobian.col(i));
      double jacc_j = diff_icov.dot(jacobian.col(j));
      res.hessian(i, j) =d1d2_e_x_cov_x *
          (-gauss_d2_ * jacc_i * jacc_j + diff_icov.dot(hess_part) +
           jacobian.col(j).dot(icov * jacobian.col(i)));
    }
  }
  return res;
}

///////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget>
double
NormalDistributionsTransform2DEx<PointSource, PointTarget>::computeStepLengthMT(
    const Eigen::Matrix<double, 3, 1> &x,
    Eigen::Matrix<double, 3, 1> &step_dir,
    double step_init,
    double step_max,
    double step_min,
    TargetGrid &ndt_grid,
    const ndt2d::ScoreAndDerivatives<3, double> &score,
    const PclSource &input)
{
  Eigen::Matrix4f transformation;
  PclSource trans_cloud;
  //transformation.setIdentity();
  ndt2d::ScoreAndDerivatives<3, double> score_vals = score;
  // Set the value of phi(0), Equation 1.3 [More, Thuente 1994]
  double phi_0 = -score_vals.value;
  // Set the value of phi'(0), Equation 1.3 [More, Thuente 1994]
  double d_phi_0 = -(score_vals.gradient.dot (step_dir));

  Eigen::Matrix<double, 3, 1>  x_t;

  if (d_phi_0 >= 0)
  {
    // Not a decent direction
    if (d_phi_0 == 0)
      return 0;
    else
    {
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

  // Auxiliary function psi is used until I is determined ot be a closed interval, Equation 2.1 [More, Thuente 1994]
  double f_l = auxilaryFunction_PsiMT (a_l, phi_0, phi_0, d_phi_0, mu);
  double g_l = auxilaryFunction_dPsiMT (d_phi_0, d_phi_0, mu);

  double f_u = auxilaryFunction_PsiMT (a_u, phi_0, phi_0, d_phi_0, mu);
  double g_u = auxilaryFunction_dPsiMT (d_phi_0, d_phi_0, mu);

  // Check used to allow More-Thuente step length calculation to be skipped by making step_min == step_max
  bool interval_converged = (step_max - step_min) > 0, open_interval = true;

  double a_t = step_init;
  a_t = std::min (a_t, step_max);
  a_t = std::max (a_t, step_min);

  x_t = x + step_dir * a_t;

  //transformation = vecToMat(x_t);

  // Updates score, gradient and hessian.  Hessian calculation is unessisary but testing showed that most step calculations use the
  // initial step suggestion and recalculation the reusable portions of the hessian would intail more computation time.
  score_vals = calcScore(ndt_grid,input,x_t);

  // Calculate phi(alpha_t)
  double phi_t = -score_vals.value;
  // Calculate phi'(alpha_t)
  double d_phi_t = -(score_vals.gradient.dot (step_dir));

  // Calculate psi(alpha_t)
  double psi_t = auxilaryFunction_PsiMT (a_t, phi_t, phi_0, d_phi_0, mu);
  // Calculate psi'(alpha_t)
  double d_psi_t = auxilaryFunction_dPsiMT (d_phi_t, d_phi_0, mu);

  // Iterate until max number of iterations, interval convergance or a value satisfies the sufficient decrease, Equation 1.1, and curvature condition, Equation 1.2 [More, Thuente 1994]
  while (!interval_converged && step_iterations < max_step_iterations && !(psi_t <= 0 /*Sufficient Decrease*/ && d_phi_t <= -nu * d_phi_0 /*Curvature Condition*/))
  {
    // Use auxilary function if interval I is not closed
    if (open_interval)
    {
      a_t = trialValueSelectionMT (a_l, f_l, g_l,
                                   a_u, f_u, g_u,
                                   a_t, psi_t, d_psi_t);
    }
    else
    {
      a_t = trialValueSelectionMT (a_l, f_l, g_l,
                                   a_u, f_u, g_u,
                                   a_t, phi_t, d_phi_t);
    }

    a_t = std::min (a_t, step_max);
    a_t = std::max (a_t, step_min);

    x_t = x + step_dir * a_t;
    //transformation = vecToMat(x_t);

    score_vals = calcScore(ndt_grid,input,x_t);

    // Calculate phi(alpha_t+)
    phi_t = -score_vals.value;
    // Calculate phi'(alpha_t+)
    d_phi_t = -(score_vals.gradient.dot (step_dir));

    // Calculate psi(alpha_t+)
    psi_t = auxilaryFunction_PsiMT (a_t, phi_t, phi_0, d_phi_0, mu);
    // Calculate psi'(alpha_t+)
    d_psi_t = auxilaryFunction_dPsiMT (d_phi_t, d_phi_0, mu);

    // Check if I is now a closed interval
    if (open_interval && (psi_t <= 0 && d_psi_t >= 0))
    {
      open_interval = false;

      // Converts f_l and g_l from psi to phi
      f_l = f_l + phi_0 - mu * d_phi_0 * a_l;
      g_l = g_l + mu * d_phi_0;

      // Converts f_u and g_u from psi to phi
      f_u = f_u + phi_0 - mu * d_phi_0 * a_u;
      g_u = g_u + mu * d_phi_0;
    }

    if (open_interval)
    {
      // Update interval end points using Updating Algorithm [More, Thuente 1994]
      interval_converged = updateIntervalMT (a_l, f_l, g_l,
                                             a_u, f_u, g_u,
                                             a_t, psi_t, d_psi_t);
    }
    else
    {
      // Update interval end points using Modified Updating Algorithm [More, Thuente 1994]
      interval_converged = updateIntervalMT (a_l, f_l, g_l,
                                             a_u, f_u, g_u,
                                             a_t, phi_t, d_phi_t);
    }

    step_iterations++;
  }

  // If inner loop was run then hessian needs to be calculated.
  // Hessian is unnessisary for step length determination but gradients are required
  // so derivative and transform data is stored for the next iteration.
  // if (step_iterations)
  //   computeHessian (hessian, trans_cloud, x_t);

  return (a_t);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget>
bool
NormalDistributionsTransform2DEx<PointSource, PointTarget>::updateIntervalMT(
    double &a_l, double &f_l, double &g_l, double &a_u, double &f_u,
    double &g_u, double a_t, double f_t, double g_t) const
{
  // Case U1 in Update Algorithm and Case a in Modified Update Algorithm [More, Thuente 1994]
  if (f_t > f_l)
  {
    a_u = a_t;
    f_u = f_t;
    g_u = g_t;
    return (false);
  }
  // Case U2 in Update Algorithm and Case b in Modified Update Algorithm [More, Thuente 1994]
  else
  if (g_t * (a_l - a_t) > 0)
  {
    a_l = a_t;
    f_l = f_t;
    g_l = g_t;
    return (false);
  }
  // Case U3 in Update Algorithm and Case c in Modified Update Algorithm [More, Thuente 1994]
  else
  if (g_t * (a_l - a_t) < 0)
  {
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
template <typename PointSource, typename PointTarget>
double NormalDistributionsTransform2DEx<
    PointSource, PointTarget>::trialValueSelectionMT(double a_l, double f_l,
                                                     double g_l, double a_u,
                                                     double f_u, double g_u,
                                                     double a_t, double f_t,
                                                     double g_t) const
{
  // Case 1 in Trial Value Selection [More, Thuente 1994]
  if (f_t > f_l)
  {
    // Calculate the minimizer of the cubic that interpolates f_l, f_t, g_l and g_t
    // Equation 2.4.52 [Sun, Yuan 2006]
    double z = 3 * (f_t - f_l) / (a_t - a_l) - g_t - g_l;
    double w = std::sqrt (z * z - g_t * g_l);
    // Equation 2.4.56 [Sun, Yuan 2006]
    double a_c = a_l + (a_t - a_l) * (w - g_l - z) / (g_t - g_l + 2 * w);

    // Calculate the minimizer of the quadratic that interpolates f_l, f_t and g_l
    // Equation 2.4.2 [Sun, Yuan 2006]
    double a_q = a_l - 0.5 * (a_l - a_t) * g_l / (g_l - (f_l - f_t) / (a_l - a_t));

    if (std::fabs (a_c - a_l) < std::fabs (a_q - a_l))
      return (a_c);
    else
      return (0.5 * (a_q + a_c));
  }
  // Case 2 in Trial Value Selection [More, Thuente 1994]
  else
  if (g_t * g_l < 0)
  {
    // Calculate the minimizer of the cubic that interpolates f_l, f_t, g_l and g_t
    // Equation 2.4.52 [Sun, Yuan 2006]
    double z = 3 * (f_t - f_l) / (a_t - a_l) - g_t - g_l;
    double w = std::sqrt (z * z - g_t * g_l);
    // Equation 2.4.56 [Sun, Yuan 2006]
    double a_c = a_l + (a_t - a_l) * (w - g_l - z) / (g_t - g_l + 2 * w);

    // Calculate the minimizer of the quadratic that interpolates f_l, g_l and g_t
    // Equation 2.4.5 [Sun, Yuan 2006]
    double a_s = a_l - (a_l - a_t) / (g_l - g_t) * g_l;

    if (std::fabs (a_c - a_t) >= std::fabs (a_s - a_t))
      return (a_c);
    else
      return (a_s);
  }
  // Case 3 in Trial Value Selection [More, Thuente 1994]
  else
  if (std::fabs (g_t) <= std::fabs (g_l))
  {
    // Calculate the minimizer of the cubic that interpolates f_l, f_t, g_l and g_t
    // Equation 2.4.52 [Sun, Yuan 2006]
    double z = 3 * (f_t - f_l) / (a_t - a_l) - g_t - g_l;
    double w = std::sqrt (z * z - g_t * g_l);
    double a_c = a_l + (a_t - a_l) * (w - g_l - z) / (g_t - g_l + 2 * w);

    // Calculate the minimizer of the quadratic that interpolates g_l and g_t
    // Equation 2.4.5 [Sun, Yuan 2006]
    double a_s = a_l - (a_l - a_t) / (g_l - g_t) * g_l;

    double a_t_next;

    if (std::fabs (a_c - a_t) < std::fabs (a_s - a_t))
      a_t_next = a_c;
    else
      a_t_next = a_s;

    if (a_t > a_l)
      return (std::min (a_t + 0.66 * (a_u - a_t), a_t_next));
    else
      return (std::max (a_t + 0.66 * (a_u - a_t), a_t_next));
  }
  // Case 4 in Trial Value Selection [More, Thuente 1994]
  else
  {
    // Calculate the minimizer of the cubic that interpolates f_u, f_t, g_u and g_t
    // Equation 2.4.52 [Sun, Yuan 2006]
    double z = 3 * (f_t - f_u) / (a_t - a_u) - g_t - g_u;
    double w = std::sqrt (z * z - g_t * g_u);
    // Equation 2.4.56 [Sun, Yuan 2006]
    return (a_u + (a_t - a_u) * (w - g_u - z) / (g_t - g_u + 2 * w));
  }
}


template <typename PointSource, typename PointTarget>
Eigen::Matrix4f NormalDistributionsTransform2DEx<
    PointSource, PointTarget>::vecToMat(const Eigen::Vector3d &trans) const
{
  Eigen::Matrix4f trans_mat = Eigen::Matrix4f::Identity();

  trans_mat.block<3, 3>(0, 0).matrix() = Eigen::Matrix3f(Eigen::AngleAxisf(
      static_cast<float>(trans(2)), Eigen::Vector3f::UnitZ()));

  trans_mat.block<3, 1>(0, 3).matrix() = Eigen::Vector3f(
      static_cast<float>(trans(0)), static_cast<float>(trans(1)), 0.0f);

  return trans_mat;
}

template <typename PointSource, typename PointTarget>
Eigen::Vector3d NormalDistributionsTransform2DEx<
    PointSource, PointTarget>::matToVec(const Eigen::Matrix4f &trans) const
{
  Eigen::Vector3d vec;
  Eigen::Affine3f trans_mat(trans);
  Eigen::Vector3f translation = trans_mat.translation();
  vec << translation(0), translation(1),
      std::atan2(trans_mat.rotation()(1, 0), trans_mat.rotation()(0, 0));
  return vec;
}

} // end of namespace

#endif
