#ifndef NDT_GSLAM_D2D_NDT2D_ROBUST
#define NDT_GSLAM_D2D_NDT2D_ROBUST

#include <ndt_gslam/registration/correlative_estimation2d.h>
#include <ndt_gslam/registration/correlative_estimation_tools.h>
#include <ndt_gslam/registration/d2d_ndt2d.h>
#include <pcl/registration/registration.h>

#include <ndt_gslam/ndt/ndt_cell.h>
#include <ndt_gslam/ndt/ndt_grid2d.h>
#include <exception>

namespace pcl
{
template <typename PointSource, typename PointTarget,
          typename CellType = slamuk::NDTCell>
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

  typedef slamuk::NDTGrid2D<CellType, PointTarget> GridTarget;
  typedef slamuk::NDTGrid2D<CellType, PointSource> GridSource;
  typedef typename GridTarget::ConstPtr GridTargetConstPtr;
  typedef typename GridSource::ConstPtr GridSourceConstPtr;

public:
  typedef boost::shared_ptr<
      D2DNormalDistributionsTransform2DRobust<PointSource, PointTarget>>
      Ptr;
  typedef boost::shared_ptr<
      const D2DNormalDistributionsTransform2DRobust<PointSource, PointTarget>>
      ConstPtr;
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

  virtual void setInputSource(const GridSourceConstPtr &grid)
  {
    d2d_.setInputSource(grid);
    PclSourceConstPtr pcl = grid->getMeans();
    corr_est_.setInputSource(pcl);
    Registration<PointSource, PointTarget>::setInputSource(pcl);
  }

  virtual void setInputTarget(const GridTargetConstPtr &grid)
  {
    PclTargetConstPtr pcl = grid->getMeans();
    Registration<PointSource, PointTarget>::setInputTarget(pcl);
    d2d_.setInputTarget(grid);
    corr_est_.setInputTarget(pcl);
    cell_size_ = grid->getCellSize();
    corr_est_.setCoarseStep(cell_size_);
  }

  void setRejectionLimit(float limit)
  {
    if (limit < 1 && limit > 0)
      rejection_limit_ = limit;
    else
      throw std::invalid_argument("Rejection limit in robust scan-matching out "
                                  "of bounds [0,1]");
  }
  void setNumLayers(size_t num)
  {
    d2d_.setNumLayers(num);
  }

  size_t getNumLayers()
  {
    return d2d_.getNumLayers();
  }
  /** \brief Set/change the voxel grid cell size for largest grid(finest).
   * Other grids
   * will have smaller longer cell size (coarser)
    * \param[in] base_size side length of voxels
    */
  void setCellSize(float base_size)
  {
    d2d_.setCellSize(base_size);
    cell_size_ = base_size;
  }

  /** \brief Get voxel grid resolution.
    * \return side length of the moast coarse
    */
  float getCellSize() const
  {
    return cell_size_;
  }

  /** \brief Get the newton line search maximum step length.
    * \return maximum step length
    */
  double getStepSize() const
  {
    return d2d_.getStepSize();
  }

  /** \brief Set/change the newton line search maximum step length.
    * \param[in] step_size maximum step length
    */
  void setStepSize(double step_size)
  {
    d2d_.setStepSize(step_size);
  }

  double getOutlierRatio() const
  {
    return d2d_.getOutlierRatio();
  }

  void setOulierRatio(double step_size)
  {
    d2d_.setOulierRatio(step_size);
  }

  float getAlignmentQuality() const
  {
    return alignment_quality_;
  }
  /** \brief Get the number of iterations required to calculate alignment.
    * \return final number of iterations
    */
  int getFinalNumIteration() const
  {
    return d2d_.getFinalNumIteration();
  }

  Eigen::Matrix3d getCovariance() const
  {
    return d2d_.getCovariance();
  }

  Eigen::Matrix3d getInformMatrix() const
  {
    return d2d_.getInformMatrix();
  }

  void setTranslationRange(float range)
  {
    corr_est_.setTranslationRange(range);
  }
  void setRotationRange(float range)
  {
    corr_est_.setRotationRange(range);
  }
  void enableMultithreading(unsigned int thread_count)
  {
    d2d_.enableMultithreading(thread_count);
    corr_est_.enableMultithreading(thread_count);
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

  float cell_size_;
  float alignment_quality_;
  float rejection_limit_;

  /** \brief Estimate the transformation and returns the transformed source
   * (input) as output.
    * \param[out] output the resultant input transformed point cloud dataset
    */
  virtual void computeTransformation(PclSource &output)
  {
    computeTransformation(output, Eigen::Matrix4f::Identity());
  }

  /** \brief Estimate the transformation and returns the transformed source
   * (input) as output.
    * \param[out] output the resultant input transformed point cloud dataset
    * \param[in] guess the initial gross estimation of the transformation
    */
  virtual void computeTransformation(PclSource &output,
                                     const Eigen::Matrix4f &guess);

  virtual double proofTransform(const Eigen::Matrix4f &trans);
};

//////////////////////////////IMPLEMENTATION
template <typename PointSource, typename PointTarget, typename CellType>
D2DNormalDistributionsTransform2DRobust<PointSource, PointTarget, CellType>::
    D2DNormalDistributionsTransform2DRobust()
  : cell_size_(0.25), alignment_quality_(0), rejection_limit_(0.4)
{
  d2d_.setCellSize(cell_size_);
  d2d_.setMaximumIterations(10);
}
//////////////////
template <typename PointSource, typename PointTarget, typename CellType>
void D2DNormalDistributionsTransform2DRobust<
    PointSource, PointTarget,
    CellType>::computeTransformation(PclSource &output,
                                     const Eigen::Matrix4f &guess)
{
  // standard D2D match try for good guess estimations
  d2d_.align(output, guess);
  double score = proofTransform(d2d_.getFinalTransformation());
  Eigen::Matrix4f first_trans;
  // test if first d2d has return bad result
  if (!(d2d_.hasConverged() && score > 0.7)) {
    // bad result -> robust alignment needed
    first_trans = d2d_.getFinalTransformation();
    corr_est_.align(output, guess);
    if (!corr_est_.hasConverged()) {
      converged_ = false;
      alignment_quality_ = 0;
      final_transformation_.setIdentity();
      return;
    }
    // second d2d -> precise alignment
    d2d_.align(output, corr_est_.getFinalTransformation());
    if (!d2d_.hasConverged()) {  //||
      //! proofTransform(d2d_.getFinalTransformation()))
      //{
      converged_ = false;
      final_transformation_ = corr_est_.getFinalTransformation();
      alignment_quality_ = 0;
      return;
    }
  }
  // score result
  double score2 = proofTransform(d2d_.getFinalTransformation());
  alignment_quality_ = score;
  // robust alignment still not good enough
  if (score2 < rejection_limit_) {
    // Maybe at least first d2d got some reasonable result
    if (score > 0.6) {
      converged_ = true;
      final_transformation_ = first_trans;
    } else {
      // everything is bad probably not the same place
      converged_ = false;
      final_transformation_ = d2d_.getFinalTransformation();
    }
  } else {
    // we got good result by robust algorithm
    converged_ = true;
    alignment_quality_ = score2;
    final_transformation_ = d2d_.getFinalTransformation();
  }
  // output cloud transform
  transformPointCloud(*input_, output, final_transformation_);
}

/////////////////////
template <typename PointSource, typename PointTarget, typename CellType>
double D2DNormalDistributionsTransform2DRobust<
    PointSource, PointTarget,
    CellType>::proofTransform(const Eigen::Matrix4f &trans)
{
  ml_corr::LookUpTable<PointTarget> proof_grid;
  proof_grid.initGrid(*target_, cell_size_, 0.5);
  PclSource output;
  transformPointCloud(*input_, output, trans);
  double score = proof_grid.getScore(output);
  ROS_DEBUG_STREAM("[D2DRobust]:proofer score: " << score);
  return score;
}
}  // end of namespace pcl

#endif
