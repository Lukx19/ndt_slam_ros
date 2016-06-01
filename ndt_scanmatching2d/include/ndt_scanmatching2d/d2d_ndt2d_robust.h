#ifndef NDT_SCANMATCHING2D_D2D_NDT2D_ROBUST
#define NDT_SCANMATCHING2D_D2D_NDT2D_ROBUST

#include <ndt_scanmatching2d/correlative_estimation_tools.h>
#include <ndt_scanmatching2d/d2d_ndt2d.h>
#include <ndt_scanmatching2d/correlative_estimation2d.h>
#include <pcl/registration/registration.h>

namespace pcl
{
template <typename PointSource, typename PointTarget>
class D2DNormalDistributionsTransform2DRobust
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

public:
  typedef boost::shared_ptr<
      D2DNormalDistributionsTransform2DRobust<PointSource, PointTarget>> Ptr;
  typedef boost::shared_ptr<const D2DNormalDistributionsTransform2DRobust<
      PointSource, PointTarget>> ConstPtr;
  typedef Eigen::Vector3d VectorTrans;
  /** \brief Constructor.
    * Sets \ref outlier_ratio_ to 0.35, \ref step_size_ to 0.05 and \ref
   * resolution_ to 1.0
    */
  D2DNormalDistributionsTransform2DRobust();

  /** \brief Empty destructor */
  virtual ~D2DNormalDistributionsTransform2DRobust()
  {
  }

  virtual void setInputSource(const PclSourceConstPtr &cloud)
  {
    d2d_.setInputSource(cloud);
    corr_est_.setInputSource(cloud);
    Registration<PointSource, PointTarget>::setInputSource(cloud);
  }

  virtual void setInputTarget(const PclTargetConstPtr &cloud)
  {
    Registration<PointSource, PointTarget>::setInputTarget(cloud);
    d2d_.setInputTarget(cloud);
    corr_est_.setInputTarget(cloud);
  }

  inline void setNumLayers(size_t num)
  {
    d2d_.setNumLayers(num);
  }

  inline size_t getNumLayers()
  {
    return d2d_.getNumLayers();
  }
  /** \brief Set/change the voxel grid resolution for largest grid(finest).
   * Other grids
   * will have smaller resolution
    * \param[in] resolution side length of voxels
    */
  inline void setResolution(float resolution)
  {
    d2d_.setResolution(resolution);
    resolution_ = resolution;
  }

  /** \brief Get voxel grid resolution.
    * \return side length of the moast coarse
    */
  inline float getResolution() const
  {
    return resolution_;
  }

  /** \brief Get the newton line search maximum step length.
    * \return maximum step length
    */
  inline double getStepSize() const
  {
    return d2d_.getStepSize();
  }

  /** \brief Set/change the newton line search maximum step length.
    * \param[in] step_size maximum step length
    */
  inline void setStepSize(double step_size)
  {
    d2d_.setStepSize(step_size);
  }

  /** \brief Get the point cloud outlier ratio.
    * \return outlier ratio
    */
  inline double getOulierRatio() const
  {
    return d2d_.getOulierRatio();
  }

  /** \brief Set/change the point cloud outlier ratio.
    * \param[in] outlier_ratio outlier ratio
    */
  inline void setOulierRatio(double outlier_ratio)
  {
    d2d_.setOulierRatio(outlier_ratio);
  }

  /** \brief Get the registration alignment probability.
    * \return transformation probability
    */
  inline float getTransformationProbability() const
  {
    return trans_probability_;
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
    return d2d_.getCovariance();
  }

  inline Eigen::Matrix3d getInformMatrix() const
  {
    return d2d_.getInformMatrix();
  }

protected:
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

  D2DNormalDistributionsTransform2D<PointSource, PointTarget> d2d_;
  CorrelativeEstimation<PointSource, PointTarget> corr_est_;

  float resolution_;
  float trans_probability_;

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

  virtual bool proofTransform(const Eigen::Matrix4f &trans);
};

//////////////////////////////IMPLEMENTATION
template <typename PointSource, typename PointTarget>
D2DNormalDistributionsTransform2DRobust<
    PointSource, PointTarget>::D2DNormalDistributionsTransform2DRobust()
  : resolution_(2), trans_probability_(0)
{
  d2d_.setMaximumIterations(max_iterations_);
}
//////////////////
template <typename PointSource, typename PointTarget>
void D2DNormalDistributionsTransform2DRobust<PointSource, PointTarget>::
    computeTransformation(PclSource &output, const Eigen::Matrix4f &guess)
{
  corr_est_.align(output, guess);
  if (!corr_est_.hasConverged()) {
    converged_ = false;
    return;
  }
  d2d_.align(output, corr_est_.getFinalTransformation());
  if (!d2d_.hasConverged() || !proofTransform(d2d_.getFinalTransformation())) {
    converged_ = false;
    return;
  }

  // output data
  converged_ = true;
  final_transformation_ = d2d_.getFinalTransformation();
  transformPointCloud(*input_, output, final_transformation_);
}

template <typename PointSource, typename PointTarget>
bool D2DNormalDistributionsTransform2DRobust<
    PointSource, PointTarget>::proofTransform(const Eigen::Matrix4f &trans)
{
  ml_corr::LookUpTable<PointTarget> proof_grid;
  proof_grid.initGrid(*target_,1/resolution_,0.5);
  PclSource output;
  transformPointCloud(*input_, output, trans);
  double score = proof_grid.getScore(output);
  ROS_DEBUG_STREAM("proofer score: "<< score);
  if(score < 0.2)
    return false;
  else
    return true;

}
}  // end of namespace pcl

#endif
