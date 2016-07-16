#ifndef GRAPH_SLAM_UK_GRAPH_SLAM2D_ALG
#define GRAPH_SLAM_UK_GRAPH_SLAM2D_ALG

#include <graph_slam_uk/ndt_grid2d_holder.h>
#include <graph_slam_uk/slam_optimizer/graph_slam_interfaces.h>
#include <graph_slam_uk/slam_optimizer/ndt_scanmatcher.h>
#include <graph_slam_uk/slam_optimizer/slam2d.h>
#include <graph_slam_uk/utils/covariance_wrapper.h>
#include <graph_slam_uk/utils/eigen_tools.h>
#include <opencv/cv.h>
#include <ros/ros.h>
#include <Eigen/Dense>

#include <graph_slam_uk/ndt/cell_policy2d.h>
#include <graph_slam_uk/ndt/ndt_cell.h>
#include <graph_slam_uk/ndt/ndt_mapper.h>
#include <graph_slam_uk/ndt/output_msgs.h>
#include <graph_slam_uk/registration/d2d_ndt2d.h>
#include <graph_slam_uk/slam_algorithm_interface.h>
#include <nav_msgs/OccupancyGrid.h>

namespace slamuk
{
class NdtSlamAlgortihm : public ISlamAlgorithm
{
  typedef eigt::pose2d_t<double> Pose;
  typedef ISlamAlgorithm::PointType PointType;
  typedef NDTCell<CellPolicy2d> CellType;
  typedef NDTGrid2D<NDTCell<CellPolicy2d>, PointType> FrameType;
  typedef FrameType::Ptr FrameTypePtr;
  typedef NDTGrid2DHolder<NDTCell<CellPolicy2d>, PointType> FrameTypeHolder;

public:
  typedef ISlamAlgorithm::Transform Transform;
  typedef ISlamAlgorithm::Covar Covar;
  typedef ISlamAlgorithm::PointCloud PointCloud;

public:
  NdtSlamAlgortihm();
  virtual inline Transform update(const Transform &motion,
                                  const Covar &covariance,
                                  const PointCloud &pcl,
                                  const ros::Time &update_time) override;
  virtual nav_msgs::OccupancyGrid
  calcOccupancyGrid(std::string &world_frame_id) const override
  {
    auto map = map_.calcOccupancyGridMsg();
    map.header.frame_id = world_frame_id;
    return map;
  }

  virtual visualization_msgs::MarkerArray
  getGraphSerialized(std::string world_frame_id) const override
  {
    return opt_engine_->getGraphSerialized(world_frame_id);
  }
  virtual graph_slam_uk::NDTMapMsg
  getWindowMap(std::string &world_frame_id) const
  {
    return toNdtMapMsg(running_window_->serialize(), world_frame_id);
  }

  virtual void setRunWindowRadius(float radius) override
  {
    win_radius_ = radius;
  }

protected:
  bool initialized_;
  bool first_node_added_;
  float win_radius_;

  size_t last_node_id_;
  CovarianceWrapper last_covar_;
  Transform last_trans_;

  Transform position_;
  Transform position_cumul_;
  Transform transform_temp_;  // transform of current temp frame
  CovarianceWrapper covar_temp_;
  NDTMapper<CellType, PointType> map_;
  FrameTypePtr running_window_;
  FrameTypePtr frame_temp_;

  std::unique_ptr<IScanmatcher2d<FrameTypeHolder>> loop_scanmatcher_;
  std::unique_ptr<IGraphOptimalizer2d<FrameTypeHolder>> opt_engine_;
  // matcher used only for incremental matching against running window
  pcl::D2DNormalDistributionsTransform2D<PointType, PointType, CellType>
      inc_matcher_;

  double max_uncertanity_in_window_;
  double max_euclidean_distance_traveled_;

  ros::Time window_update_time_;
  unsigned int window_seq_;

  inline bool movedEnough(const Transform &trans,
                          const CovarianceWrapper &covar) const;
  inline graph_slam_uk::NDTMapMsg
  toNdtMapMsg(const NDTGridMsg &msg, const std::string &fixed_frame) const;
};
NdtSlamAlgortihm::NdtSlamAlgortihm()
  : initialized_(false)
  , first_node_added_(false)
  , win_radius_(30.0f)
  , last_node_id_(0)
  , position_(Transform::Identity())
  , position_cumul_(Transform::Identity())
  , transform_temp_(Transform::Identity())
  , covar_temp_()
  , map_()
  , running_window_(new FrameType())
  , frame_temp_(new FrameType())
  , loop_scanmatcher_(new NdtScanmatcher<FrameTypeHolder>())
  , opt_engine_(new Slam2D<FrameTypeHolder>(*loop_scanmatcher_))
  , inc_matcher_()
  , max_uncertanity_in_window_(10.0)
  , max_euclidean_distance_traveled_(5.0)
{
}

NdtSlamAlgortihm::Transform
NdtSlamAlgortihm::update(const Transform &motion, const Covar &covariance,
                         const PointCloud &pcl, const ros::Time &update_time)
{
  // PROCESS OF ADDING TO POSE GRAPH
  // 1. program is integrating measurements, transformations and covariances
  // 2. robot moved enough
  // 3 local map is added to mapper.
  // 4 a)  integrated covariance is copied to last known covariance var
  //   b) integrated transform is  copied to last known transform var
  // 5. first node is added to pose graph with current aggregated values (temp
  // variables)
  // 6. repeat 1.
  // 7. execute 2. 3.
  // 8. add new node to graph
  // 9. add edge between last two nodes in graph with values saved in step 3
  // 10. execute 4. and repeat process
  if (!initialized_) {
    running_window_->setOrigin(FrameType::Pose(0, 0, 0));
    frame_temp_->setOrigin(FrameType::Pose(0, 0, 0));
    position_.setIdentity();

    initialized_ = true;
    // enters first scan as it is to running window
    running_window_->enlarge(win_radius_, win_radius_, win_radius_,
                             win_radius_);
    running_window_->initialize(pcl);
    frame_temp_->initialize(pcl);
    window_update_time_ = update_time;
    window_seq_ = 0;
    return position_;
  }
  Transform position_odom = position_ * motion;
  FrameTypePtr local(new FrameType());
  PointCloud pcl_out;
  local->setOrigin(running_window_->getOrigin());
  local->initializeSimple(pcl);

  inc_matcher_.setInputTarget(running_window_);
  inc_matcher_.setInputSource(local);
  inc_matcher_.align(
      pcl_out,
      eigt::convertFromTransform(position_cumul_ * motion).cast<float>());
  ROS_INFO_STREAM("[SLAM ALGORITHM]: incremental scanamtching converged:"
                  << inc_matcher_.hasConverged());
  if (inc_matcher_.hasConverged()) {
    // prepare transformation from successful scan registration
    Transform trans_mat = eigt::convertToTransform<double>(
        inc_matcher_.getFinalTransformation().cast<double>());
    // updating cumulative position. This position meassures transformation is
    // relative to running window origin
    position_cumul_ = trans_mat;
    local->transform(position_cumul_);
    // merge in new data to running window
    running_window_->mergeInTraced(*local, true, false);
    // move running window only horizontaly or verticaly if needed
    position_cumul_ = running_window_->move(position_cumul_);

    // info for output msg
    window_update_time_ = update_time;
    ++window_seq_;

    // increasing transformation and covariance of local temp frame
    Transform diff_trans = eigt::transBtwPoses(
        eigt::getPoseFromTransform(position_), local->getOrigin());
    std::cout << "diff trans:"
              << eigt::getPoseFromTransform(diff_trans).transpose()
              << std::endl;
    transform_temp_ = transform_temp_ * diff_trans;
    covar_temp_.addToCovar(covariance, diff_trans);
    if (movedEnough(transform_temp_, covar_temp_)) {
      ROS_INFO_STREAM("[SLAM ALGORITHM]: creating new frame");
      // robot moved out of reasonable bounds for single temp frame
      map_.addFrame(frame_temp_, update_time);
      FrameTypeHolder frame_wrap(frame_temp_);
      last_node_id_ =
          opt_engine_->addPose(frame_temp_->getOrigin(), frame_wrap);
      // add edge measurement only after first node was inserted
      if (first_node_added_) {
        opt_engine_->addLastConstrain(eigt::getPoseFromTransform(last_trans_),
                                      last_covar_.covar_.inverse());
      }
      last_covar_ = covar_temp_;
      last_trans_ = transform_temp_;
      first_node_added_ = true;
      // initialize new grid for next round
      frame_temp_ = FrameTypePtr(new FrameType());
      frame_temp_->setOrigin(local->getOrigin());
      frame_temp_->mergeInTraced(*local, false, true);
      covar_temp_ = CovarianceWrapper();
      transform_temp_.setIdentity();
      // find loop closures and optimize graph
      // if (opt_engine_->tryLoopClose()) {
      //   std::vector<size_t> changed_nodes;
      //   // optimize graph
      //   if (opt_engine_->optimalizeIterationaly()) {
      //     // optimization updates origins in grid through FrameHolder
      //     // interface
      //     map_.recalc(update_time);
      //   }
      // }
      // update position based on information from pose graph plus last addition
      // which is currently in next tem frae...so it is not considered by pose
      // graph yet
      std::cout << "last nodes pose in graph: "
                << opt_engine_->getPoseLocation(last_node_id_).transpose()
                << std::endl;
      position_ =
          eigt::getTransFromPose(opt_engine_->getPoseLocation(last_node_id_)) *
          diff_trans;
    } else {
      // robot is still in bounds of reasonable error from incremental
      // scan matching. We don't want to add to many nodes into pose graph
      frame_temp_->mergeInTraced(*local, true, true);
      position_ = eigt::getTransFromPose(local->getOrigin());
    }
  } else {
    ROS_INFO_STREAM("[SLAM ALGORITHM]: unsucessful scanmatching-> using "
                    "odometry");
    // scan-matching failed, use odometry
    // do not merge in scan - avoiding problems whit bad local map
    position_ = position_odom;
  }
  std::cout << "transform_temp: "
            << eigt::getPoseFromTransform(transform_temp_).transpose()
            << std::endl;
  std::cout << "position: " << eigt::getPoseFromTransform(position_).transpose()
            << std::endl;
  return position_;
}

bool NdtSlamAlgortihm::movedEnough(const Transform &trans,
                                   const CovarianceWrapper &covar) const
{
  double uncertainty = (covar.covar_.determinant());
  double euclidean_dist = trans.translation().norm();
  ROS_DEBUG_STREAM("[NDT_SLAM_ALGORITHM]:uncertainity: "
                   << uncertainty << "  moved_distance: " << euclidean_dist);
  if (uncertainty > max_uncertanity_in_window_ ||
      euclidean_dist > max_euclidean_distance_traveled_)
    return true;
  else
    return false;
}

graph_slam_uk::NDTMapMsg NdtSlamAlgortihm::toNdtMapMsg(
    const NDTGridMsg &msg, const std::string &fixed_frame) const
{
  graph_slam_uk::NDTMapMsg ros_msg;
  ros_msg.header.frame_id = fixed_frame;
  ros_msg.header.stamp = window_update_time_;
  ros_msg.header.seq = window_seq_;
  ros_msg.x_cen = msg.origin_(0);
  ros_msg.y_cen = msg.origin_(1);
  ros_msg.z_cen = msg.origin_(2);
  ros_msg.x_size = msg.size_(0);
  ros_msg.y_size = msg.size_(1);
  ros_msg.z_size = msg.size_(2);
  ros_msg.x_cell_size = msg.cell_sizes_(0);
  ros_msg.y_cell_size = msg.cell_sizes_(1);
  ros_msg.z_cell_size = msg.cell_sizes_(2);
  for (auto &&cell : msg.cells_) {
    graph_slam_uk::NDTCellMsg cl;
    cl.mean_x = cell.mean_(0);
    cl.mean_y = cell.mean_(1);
    cl.mean_z = cell.mean_(2);
    cl.occupancy = cell.occupancy_;
    cl.N = cell.points_;
    for (size_t i = 0; i < cell.cov_.rows(); ++i) {
      for (size_t j = 0; j < cell.cov_.cols(); ++j) {
        cl.cov_matrix.push_back(cell.cov_(i, j));
      }
    }
    ros_msg.cells.push_back(std::move(cl));
  }
  return ros_msg;
}
}  // end of namespace slamuk

#endif
