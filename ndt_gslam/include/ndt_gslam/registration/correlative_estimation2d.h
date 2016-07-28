#ifndef NDT_GSLAM_CORRELATIVE_ESTIMATION
#define NDT_GSLAM_CORRELATIVE_ESTIMATION

#include <ndt_gslam/registration/correlative_estimation_tools.h>
#include <ndt_gslam/utils/eigen_tools.h>
#include <pcl/common/time.h>
#include <pcl/filters/voxel_grid_covariance.h>
#include <pcl/registration/registration.h>
#include <ros/ros.h>
#include <Eigen/Dense>

namespace pcl
{
template <typename PointSource, typename PointTarget>
class CorrelativeEstimation : public Registration<PointSource, PointTarget>
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
  typedef VoxelGridCovariance<PointSource> SourceGrid;

  /** \brief Typename of const pointer to searchable voxel grid. */
  typedef const TargetGrid *TargetGridConstPtr;
  /** \brief Typename of const pointer to searchable voxel grid leaf. */
  typedef typename TargetGrid::LeafConstPtr TargetGridLeafConstPtr;

public:
  typedef boost::shared_ptr<CorrelativeEstimation<PointSource, PointTarget>> Ptr;
  typedef boost::shared_ptr<const CorrelativeEstimation<PointSource, PointTarget>>
      ConstPtr;
  typedef Eigen::Vector3d VectorTrans;

  CorrelativeEstimation();

  /**
   * @brief      Empty destructor
   */
  virtual ~CorrelativeEstimation()
  {
  }

  /**
   * @brief      Sets the cell size of the coarse grid.
   *
   * @param[in]  step  The step
   */
  inline void setCoarseStep(float step)
  {
    coarse_step_ = step;
  }

  /**
   * @brief      Gets Sets the cell size of the coarse grid.
   *
   * @return     The coarse step.
   */
  inline float getCoarseStep() const
  {
    return coarse_step_;
  }

  /**
   * @brief      Sets the translation range.
   *
   * @param[in]  range  The range
   */
  void setTranslationRange(float range)
  {
    translation_range_ = range;
  }
  /**
   * @brief      Sets the rotation range.
   *
   * @param[in]  range  The range
   */
  void setRotationRange(float range)
  {
    rotation_range_ = range;
  }
  /**
   * @brief      Enables the multithreading.
   *
   * @param[in]  thread_count  The thread count
   */
  void enableMultithreading(unsigned int thread_count)
  {
    threads_ = thread_count;
  }

protected:
  using Registration<PointSource, PointTarget>::input_;
  using Registration<PointSource, PointTarget>::target_;
  using Registration<PointSource, PointTarget>::final_transformation_;
  using Registration<PointSource, PointTarget>::converged_;

  float coarse_step_;
  float coarse_rot_step_;
  float translation_range_;  // search from [-range,range]
  float rotation_range_;     // maximum is [-Pi,Pi]
  unsigned int threads_;

  ml_corr::LookUpTable<PointTarget> coarse_lookup_;
  ml_corr::LookUpTable<PointTarget> fine_lookup_;

  /**
   * @brief      Estimate the transformation and returns the transformed source
   *             (input) as output.
   *
   * @param[out] output  the resultant input transfomed point cloud dataset
   */
  virtual void computeTransformation(PclSource &output)
  {
    computeTransformation(output, Eigen::Matrix4f::Identity());
  }

  /**
    * @brief      Estimate the transformation and returns the transformed source
    *             (input) as output.
    *
    * @param[out] output  the resultant input transfomed point cloud dataset
    * @param[in]  guess   the initial gross estimation of the transformation
    */
  virtual void computeTransformation(PclSource &output,
                                     const Eigen::Matrix4f &guess);

  template <typename T = float, typename In>
  Eigen::Matrix<T, 4, 4> vecToMat(const Eigen::Matrix<In, 3, 1> &trans) const;
};

template <typename PointSource, typename PointTarget>
CorrelativeEstimation<PointSource, PointTarget>::CorrelativeEstimation()
  : coarse_step_(0.5f)
  , coarse_rot_step_(0.2f)
  , translation_range_(4.0f)
  , rotation_range_(1.6f)
  , threads_(2)
{
  converged_ = false;
}

template <typename PointSource, typename PointTarget>
void CorrelativeEstimation<PointSource, PointTarget>::computeTransformation(
    PclSource &output, const Eigen::Matrix4f &guess)
{
  pcl::PointCloud<PointSource> guess_pcl;  // pcl with points transformed by
                                           // initial guess
  transformPointCloud(*input_, guess_pcl, guess);
  // initialize lookup tables
  coarse_lookup_.initGrid(*target_, coarse_step_, 0.5f);
  ROS_DEBUG_STREAM("[CorrelativeEstimation]: grids initialized");
  //////////////////
  ml_corr::SearchVoxel voxel2;
  voxel2.transform_ << 0, 0, 0;
  voxel2.score_ = coarse_lookup_.getScore(guess_pcl);
  ROS_DEBUG_STREAM("[CorrelativeEstimation]: Coarse Voxel guess = trans: "
                   << voxel2.transform_.transpose()
                   << " score: " << voxel2.score_);

  voxel2.transform_ << 0, 0, 0;
  voxel2.score_ = coarse_lookup_.getScore(*target_);
  ROS_DEBUG_STREAM("[CorrelativeEstimation]: Coarse Voxel target= trans: "
                   << voxel2.transform_.transpose()
                   << " score: " << voxel2.score_);
  ROS_DEBUG_STREAM("[CorrelativeEstimation]: maximum on target: "
                   << coarse_lookup_.getMaxScore());
  // prepare array for all results
  std::vector<ml_corr::SearchVoxel> search_voxels;
  size_t elements = static_cast<size_t>(
      std::ceil(std::pow((2 * translation_range_) / coarse_step_, 2) *
                (rotation_range_ / coarse_rot_step_)));

  search_voxels.reserve(elements);

  // preparre storage for thrread based result voxels
  std::vector<std::vector<ml_corr::SearchVoxel>> search_voxels_thread(threads_);
  for (size_t i = 0; i < threads_; ++i) {
    search_voxels_thread.push_back(std::vector<ml_corr::SearchVoxel>());
  }

  // prepare all rotations
  std::vector<float> rotations;
  for (float theta = -rotation_range_; theta < rotation_range_;
       theta += coarse_rot_step_) {
    rotations.push_back(theta);
  }
  {
    // pcl::ScopeTime t_init("[CorrelativeEstimation]: total calculation
    // time:");

    pcl::PointCloud<PointSource> rot_pcl;  // cloud with points after applied
                                           // rotation
    std::vector<ml_corr::IndexPoint> index_pcl;  // rotated points to indexes
    std::vector<ml_corr::IndexPoint> transl_index_pcl;  // translated rotated
                                                        // indexes

#pragma omp parallel for private(rot_pcl, index_pcl,                           \
                                 transl_index_pcl) num_threads(threads_)
    for (size_t theta_id = 0; theta_id < rotations.size(); ++theta_id) {
      // prepare storage for calculated voxels
      std::vector<ml_corr::SearchVoxel> local_voxels;
      local_voxels.reserve(static_cast<size_t>(
          std::ceil(std::pow(2 * translation_range_ / coarse_step_, 2))));
      // rotate point cloud
      ml_corr::rotatePointCloud(guess_pcl, rot_pcl, rotations[theta_id]);
      index_pcl = coarse_lookup_.toIndexes(rot_pcl);
      // iterate over all possible translation x values
      for (float dx = -translation_range_; dx < translation_range_;
           dx += coarse_step_) {
        // iterate over all possible y translation values
        for (float dy = -translation_range_; dy < translation_range_;
             dy += coarse_step_) {
          coarse_lookup_.transformIndexes(index_pcl, transl_index_pcl, dx, dy);
          ml_corr::SearchVoxel voxel;
          voxel.transform_ << dx, dy, rotations[theta_id];
          voxel.score_ = coarse_lookup_.getScore(transl_index_pcl);
          local_voxels.push_back(std::move(voxel));
        }
      }
      int thread_id = omp_get_thread_num();
      std::move(local_voxels.begin(), local_voxels.end(),
                std::back_inserter(search_voxels_thread[thread_id]));
    }

    // end of parallel world
    // merge parallel results
    for (size_t i = 0; i < threads_; ++i) {
      std::copy(search_voxels_thread[i].begin(), search_voxels_thread[i].end(),
                std::back_inserter(search_voxels));
    }
  }  // end of measuring scope
  // find voxel with best score = best chance for good match
  std::sort(search_voxels.begin(), search_voxels.end());
  ml_corr::SearchVoxel best_voxel;

  Eigen::Vector3f best_trans = search_voxels.back().transform_;
  ROS_DEBUG_STREAM(
      "[CorrelativeEstimation]: final_delta_trans:" << best_trans.transpose());

  // output data
  converged_ = true;
  final_transformation_ = vecToMat(best_trans);  // * guess;
  transformPointCloud(*input_, output, final_transformation_);
}

/////////////////////////////////
template <typename PointSource, typename PointTarget>
template <typename T, typename In>
Eigen::Matrix<T, 4, 4> CorrelativeEstimation<PointSource, PointTarget>::vecToMat(
    const Eigen::Matrix<In, 3, 1> &trans) const
{
  Eigen::Matrix<T, 4, 4> trans_mat = Eigen::Matrix<T, 4, 4>::Identity();

  trans_mat.block(0, 0, 3, 3).matrix() =
      Eigen::Matrix<T, 3, 3>(Eigen::AngleAxis<T>(
          static_cast<T>(trans(2)), Eigen::Matrix<T, 3, 1>::UnitZ()));

  trans_mat.block(0, 3, 3, 1).matrix() = Eigen::Matrix<T, 3, 1>(
      static_cast<T>(trans(0)), static_cast<T>(trans(1)), 0.0);

  return trans_mat;
}

}  // end of pcl namespace

#endif
