#ifndef GRAPH_SLAM_UK_GRAPH_SLAM2D_ALG
#define GRAPH_SLAM_UK_GRAPH_SLAM2D_ALG

#include <Eigen/Dense>
#include <dynamic_slam_utils/eigen_tools.h>
#include <graph_slam_uk/graph_slam_interfaces.h>
#include <graph_slam_uk/ndt_scanmatcher.h>
#include <graph_slam_uk/ndt_grid2d_holder.h>
#include <dynamic_slam_utils/covariance_wrapper.h>
#include <opencv/cv.h>

#include <ndt_scanmatching2d/d2d_ndt.h>

namespace slamuk
{
class NdtSlamAlgortihm
{
  typedef eigt::transform2d_t<double> Transform;
  typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;
  typedef NDTGrid2D<NDTCell<CellPolicy2d>, pcl::PointXYZ> FrameType;
  typedef NDTGrid2DHolder FrameTypeHolder;

public:
  NdtSlamAlgortihm(IGraphOptimalizer2d<FrameType> &opt_engine);
  inline Eigen::Matrix3d update(const Eigen::Matrix3d &motion,
                                const Eigen::Matrix3d &covariance,
                                const PointCloud &pcl);
  cv::Mat &getOccupMap()
  {
    map_.getOccupancyMap();
  }
  bool initialized()
  {
    return initialized_;
  }
  void setRunWindowRadius(float radius)
  {
    win_radius_ = radius;
  }

protected:
  bool initialized_ = false;
  float win_radius_;
  size_t last_node_id_;
  Transform position_;
  Transform postion_win_;     // position of running window in respect to world
                              // frame
  Transform transform_temp_;  // transform of current temp frame
  CovarianceWrapper covar_temp_;
  NDTMapper<FrameType> map_;
  FrameType running_window_;
  FrameType frame_temp_;

  IGraphOptimalizer2d<FrameTypeHolder> *opt_engine_;
  // matcher used only for incremental matching against running window
  pcl::D2DNormalDistributionsTransform2D inc_matcher_;
};

Eigen::Matrix3d NdtSlamAlgortihm::update(const Eigen::Matrix3d &motion,
                                         const Eigen::Matrix3d &covariance,
                                         const PointCloud &pcl)
{
  if (!initialized_) {
    running_window_.setOrigin(Pose(0, 0, 0));
    frame_temp_.setOrigin(Pose(0, 0, 0));
    position_.setIdentity();

    initialized_ = true;
    // enters first scan as it is to running window
    running_window_.enalrge(win_radius_, win_radius_, win_radius_, win_radius_);
    running_window_.initialize(pcl);
    frame_temp_.initialize(pcl);
    return position_.matrix();
  }
  Transform position_odom = position_ * motion;
  FrameType local;
  PointCloud pcl_out;
  local.initialize(pcl);
  local.setOrigin(eigt::getPoseFromTransform(position_odom));
  inc_matcher_.setTarget(&running_window_);
  inc_matcher_.setSource(&local);
  inc_matcher_.align(pcl_out);
  if (inc_matcher_.hasConverged()) {
    // prepare transformation from successful scan registration
    Transform trans_mat =
        eigt::convertToTransform(inc_matcher_.getFinalTransformation());
    Pose trans_vec = eigt::getPoseFromTransform(trans_mat);
    // update running window only based on incremental scan-matching
    running_window_.move(tarns_vec.head(2));
    running_window_.mergeInTraced(pcl, eigt::getPoseFromTransform(position_),
                                  false);
    // increasing transformation and covariance of local temp frame
    transform_temp_ = transform_temp_ * trans_mat_;
    covar_temp_.addToCovar(trans_vec, covariance);
    if (movedEnough(transfom_temp, covar_temp_)) {
      // robot moved out of reasonable bounds for single temp frame
      FrameType *frame = map_.addFrame(std::move(frame_temp_));
      frame_temp_ = FrameType();
      frame_temp_.setOrigin(position_);
      frame_temp_.initialize(pcl);
      FrameTypeHolder fram_wrap(frame);
      last_node_id_ = opt_engine_->addPose(frame->getOrigin(), frame_wrap);
      opt_engine_->addLastConstrain(eigt::getPoseFromTransform(transform_temp_),
                                    covar_temp_.covar_.inverse());
      covar_temp_ = CovarianceWrapper();
      transform_temp_.setIdentity();
      // find loop closures and optimize graph
      if (opt_engine_->tryLoopClose()) {
        std::vector<size_t> changed_nodes;
        // optimize graph
        // changed_nodes = opt_engine_->optimizeIterationaly();
        // update frames in map
        for (const size_t &id : changed_nodes) {
          map_.updateFrame(&opt_engine_->getPoseData(id).getDataPointer(),
                           opt_engine_->getPoseLocation(id));
        }
        map_.recalc()
      }
      // update position based on information from pose graph
      position_ =
          eigt::getTransFromPose(opt_engine_->getPoseLocation(last_node_id_));
    } else {
      frame_temp_.mergeInTraced(pcl, eigt::getPoseFromTransform(position_),
                                true);
      position_ = position_ * trans_mat;
    }
  } else {
    // scan-matching failed, use odometry
    // do not merge in scan - avoiding problems whit bad local map
    position_ = position_odom;
  }
  return position_.matrix();
}

// double getDoubleTime()
// {
//   struct timeval time;
//   gettimeofday(&time, NULL);
//   return time.tv_sec + time.tv_usec * 1e-6;
// }
}  // end of namespace slamuk

#endif
