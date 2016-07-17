#ifndef GRAPH_SLAM_UK_GRAPH_SLAM2D_ALG
#define GRAPH_SLAM_UK_GRAPH_SLAM2D_ALG

#include <graph_slam_uk/ndt_grid2d_holder.h>
#include <graph_slam_uk/slam_optimizer/graph_slam_interfaces.h>
#include <graph_slam_uk/slam_optimizer/ndt_scanmatcher.h>
#include <graph_slam_uk/slam_optimizer/slam2d.h>
#include <graph_slam_uk/utils/covariance_wrapper.h>
#include <graph_slam_uk/utils/eigen_tools.h>
#include <graph_slam_uk/utils/point_cloud_tools.h>
#include <opencv/cv.h>
#include <ros/ros.h>
#include <Eigen/Dense>

#include <graph_slam_uk/ndt/cell_policy2d.h>
#include <graph_slam_uk/ndt/ndt_cell.h>
#include <graph_slam_uk/ndt/ndt_mapper.h>
#include <graph_slam_uk/ndt/output_msgs.h>
#include <graph_slam_uk/registration/d2d_ndt2d.h>
#include <graph_slam_uk/registration/ndt2d.h>
#include <graph_slam_uk/slam_algorithm_interface.h>
#include <nav_msgs/OccupancyGrid.h>

namespace slamuk
{
class NdtSlamAlgortihm : public ISlamAlgorithm
{
  typedef ISlamAlgorithm::PointType PointType;
  typedef NDTCell<CellPolicy2d> CellType;
  typedef NDTGrid2D<NDTCell<CellPolicy2d>, PointType> FrameType;
  typedef FrameType::Ptr FrameTypePtr;
  typedef NDTGrid2DHolder<NDTCell<CellPolicy2d>, PointType> FrameTypeHolder;
  typedef NdtScanmatcher<FrameTypeHolder> LoopScanmatcher;

public:
  typedef ISlamAlgorithm::Transform Transform;
  typedef ISlamAlgorithm::Covar Covar;
  typedef ISlamAlgorithm::PointCloud PointCloud;
  typedef ISlamAlgorithm::Pose Pose;

public:
  NdtSlamAlgortihm();
  virtual inline Pose update(const Transform &motion, const Covar &covariance,
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
  // PARAMETERS
  float win_radius_;
  double min_distance_;
  double min_rotation_;
  double max_uncertanity_in_window_;
  double max_euclidean_distance_traveled_;
  float cell_size_;

  size_t last_node_id_;
  CovarianceWrapper last_covar_;
  Transform last_trans_;
  Transform unused_trans_;  // transformation of unused odometry
  Pose addition_;

  Pose position_;
  Transform position_cumul_;
  Transform transform_temp_;  // transform of current temp frame
  Pose last_matcher_pose_;
  CovarianceWrapper covar_temp_;
  NDTMapper<CellType, PointType> map_;
  FrameTypePtr running_window_;
  FrameTypePtr frame_temp_;

  std::unique_ptr<IGraphOptimalizer2d<FrameTypeHolder>> opt_engine_;
  // matcher used only for incremental matching against running window
  pcl::D2DNormalDistributionsTransform2D<PointType, PointType, CellType>
      inc_matcher_;
  // pcl::NormalDistributionsTransform2DEx<PointType, PointType, CellType>
  //     inc_matcher_;

  unsigned int window_seq_;
  ros::Time window_update_time_;
  std::vector<FrameTypePtr> frames_;

  inline void updateGraph(const Transform &increase, const Covar &covariance,
                          const PointCloud &new_scan,
                          const ros::Time &update_time);
  inline bool movedEnough(const Transform &trans) const;
  inline bool movedEnough(const Transform &trans,
                          const CovarianceWrapper &covar) const;
  inline graph_slam_uk::NDTMapMsg
  toNdtMapMsg(const NDTGridMsg &msg, const std::string &fixed_frame) const;
};
NdtSlamAlgortihm::NdtSlamAlgortihm()
  : initialized_(false)
  , first_node_added_(false)
  , win_radius_(15.0f)
  , min_distance_(0.1)
  , min_rotation_(0.1)
  , max_uncertanity_in_window_(1000000.0)
  , max_euclidean_distance_traveled_(3.0)
  , cell_size_(0.25)
  , last_node_id_(0)
  , last_trans_(Transform::Identity())
  , unused_trans_(Transform::Identity())
  , position_(Pose::Zero())
  , position_cumul_(Transform::Identity())
  , transform_temp_(Transform::Identity())
  , last_matcher_pose_(Pose::Zero())
  , covar_temp_()
  , map_()
  , running_window_(new FrameType())
  , frame_temp_(new FrameType())
  , opt_engine_(new Slam2D<FrameTypeHolder>(
        std::unique_ptr<LoopScanmatcher>(new LoopScanmatcher())))
  , inc_matcher_()
  , window_seq_(0)
{
  inc_matcher_.setStepSize(0.01);
  inc_matcher_.setOulierRatio(0.99);
}

NdtSlamAlgortihm::Pose NdtSlamAlgortihm::update(const Transform &motion,
                                                const Covar &covariance,
                                                const PointCloud &pcl,
                                                const ros::Time &update_time)
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
    running_window_->setCellSize(cell_size_);
    frame_temp_->setOrigin(FrameType::Pose(0, 0, 0));
    frame_temp_->setCellSize(cell_size_);
    position_.setZero();
    last_matcher_pose_.setZero();
    initialized_ = true;
    // enters first scan as it is to running window
    running_window_->enlarge(-win_radius_, -win_radius_, win_radius_,
                             win_radius_);
    running_window_->initialize(pcl);
    frame_temp_->initialize(pcl);
    window_update_time_ = update_time;
    window_seq_ = 0;
    addition_.setZero();
    return position_;
  }
  // unused_trans_ = unused_trans_ * motion;
  // if (!movedEnough(unused_trans_)) {
  //   return position_ * unused_trans_;
  // }
  // // ROS_INFO_STREAM("[SLAM_ALGORITHM]: Robot moved enough: "
  //               << eigt::getPoseFromTransform(unused_trans_).transpose());
  FrameTypePtr local = FrameTypePtr(new FrameType());
  PointCloud::Ptr pcl_out(new PointCloud());
  local->setCellSize(cell_size_);
  local->setOrigin(running_window_->getOrigin());
  local->initializeSimple(pcl);

  inc_matcher_.setInputTarget(running_window_);
  inc_matcher_.setInputSource(local);
  // inc_matcher_.setInputSource(local->getMeans());
  inc_matcher_.align(*pcl_out,
                     eigt::convertFromTransform(position_cumul_).cast<float>());
  //
  ROS_INFO_STREAM("[SLAM ALGORITHM]: incremental scanamtching converged:"
                  << inc_matcher_.hasConverged());
  if (inc_matcher_.hasConverged()) {
    // prepare transformation from successful scan registration
    Transform registration = eigt::convertToTransform<double>(
        inc_matcher_.getFinalTransformation().cast<double>());

    // GRAPH SLAM STUFF////////////////////////////////////////////
    Pose matcher_pose =
        eigt::transformPose(running_window_->getOrigin(), registration);

    // calculate pose increase with respect to last known pose in worls of
    // scanmatching
    Transform increase = eigt::transBtwPoses(last_matcher_pose_, matcher_pose);
    last_matcher_pose_ = matcher_pose;
    // std::cout << "diff trans:"
    //           << eigt::getPoseFromTransform(increase).transpose() <<
    //           std::endl;
    // std::cout << "matcher_pose: " << matcher_pose.transpose() << std::endl;
    // update and optimize graph
    updateGraph(increase, Covar::Identity() * 0.0001, pcl, update_time);

    // RUNNING WINDOW STUFF ///////////////////////////////////
    // updating cumulative position. This position meassures
    // transformation is
    // relative to running window origin
    position_cumul_ = registration;
    local->transform(registration);
    // merge in new data to running window
    running_window_->mergeInTraced(*local, true, false);
    // move running window only horizontaly or verticaly if needed
    position_cumul_ = running_window_->move(registration);

    // info for output msg
    window_update_time_ = update_time;
    ++window_seq_;

  } else {
    ROS_INFO_STREAM("[SLAM ALGORITHM]: unsucessful scanmatching-> using "
                    "odometry");
    // scan-matching failed, use odometry
    // do not merge in scan - avoiding problems whit bad local map
    position_ = eigt::transformPose(position_, motion);
  }
  // unused_trans_.setIdentity();
  return position_;
}

void NdtSlamAlgortihm::updateGraph(const Transform &increase,
                                   const Covar &covariance,
                                   const PointCloud &new_scan,
                                   const ros::Time &update_time)
{
  transform_temp_ = eigt::getTransFromPose(eigt::transformPose(
      eigt::getPoseFromTransform(transform_temp_), increase));
  covar_temp_.addToCovar(covariance, increase);

  if (movedEnough(transform_temp_, covar_temp_)) {
    ROS_INFO_STREAM("[SLAM ALGORITHM]: creating new frame");
    // robot moved out of reasonable bounds for single temp
    // frame
    frames_.emplace_back(new FrameType());
    frames_.back().swap(frame_temp_);
    map_.addFrame(frames_.back(), update_time);
    FrameTypeHolder frame_wrap(frames_.back());

    last_node_id_ =
        opt_engine_->addPose(frames_.back()->getOrigin(), frame_wrap);
    // add edge measurement only after first node was inserted
    if (first_node_added_) {
      opt_engine_->addLastConstrain(eigt::getPoseFromTransform(last_trans_),
                                    last_covar_.covar_.inverse());
    }
    last_covar_ = covar_temp_;
    last_trans_ = transform_temp_;
    covar_temp_ = CovarianceWrapper();
    transform_temp_.setIdentity();
    first_node_added_ = true;
    pcl::visualizePcl<PointType>(running_window_->getMeans());
    // find loop closures and optimize graph
    Pose pre_opt_pose = opt_engine_->getPoseLocation(last_node_id_);
    if (opt_engine_->tryLoopClose()) {
      // optimize graph
      if (opt_engine_->optimalizeIterationaly()) {
        // optimization updates origins in grid through FrameHolder
        // interface
        // map_.recalc(update_time);
      }
    }
    // initialize new grid for next round
    Pose after_opt_pose = opt_engine_->getPoseLocation(last_node_id_);
    after_opt_pose = Pose(10, 0, 2);
    addition_ = Pose(10, 0, 0);
    // std::cout << "robot: " << robot.transpose() << std::endl;
    // std::cout << "farem: " << new_scan->getOrigin().transpose() <<
    // std::endl;
    // std::cout << "farem2: " << running_window_->getOrigin().transpose()
    //           << std::endl;
    // position_ = Pose(10, 0, 2);
    position_ = eigt::transformPose(pre_opt_pose, last_trans_);
    // position_ = eigt::transformPose(
    //     position_, eigt::transBtwPoses(after_opt_pose, pre_opt_pose));
    PointCloud pcl_out;
    pcl::transformPointCloud(
        new_scan, pcl_out,
        eigt::convertFromTransform(eigt::getTransFromPose(position_)));
    pcl::visualizePcl<PointType>(frames_.back()->getMeans());

    frame_temp_->setOrigin(position_);
    frame_temp_->setCellSize(cell_size_);
    // frame_temp_->mergeInTraced(*new_scan, true, true);
    frame_temp_->mergeInTraced(pcl_out, Eigen::Vector3d(0, 0, 0), true);
    // std::cout << "last nodes pose in graph: "
    //           << opt_engine_->getPoseLocation(last_node_id_).transpose()
    //           << std::endl;
    pcl::visualizePcl<PointType>(frame_temp_->getMeans());

  } else {
    // robot is still in bounds of reasonable error from
    // incremental
    // scan matching. We don't want to add to many nodes into
    // pose graph
    // frame_temp_->mergeInTraced(*new_scan, true, true);

    position_ = eigt::transformPose(position_, increase);
    PointCloud pcl_out;
    pcl::transformPointCloud(
        new_scan, pcl_out,
        eigt::convertFromTransform(eigt::getTransFromPose(position_)));
    frame_temp_->mergeInTraced(pcl_out, Eigen::Vector3d(0, 0, 0), true);
    // if (std::abs(addition_.sum()) > 0)
    //   pcl::visualizePcl<PointType>(frame_temp_->getMeans());
  }

  std::cout << "position: " << position_.transpose() << std::endl << std::endl;
}

bool NdtSlamAlgortihm::movedEnough(const Transform &trans,
                                   const CovarianceWrapper &covar) const
{
  double uncertainty = (covar.covar_.determinant());
  std::cout << uncertainty << std::endl;

  double euclidean_dist = std::abs(eigt::getDisplacement(trans));
  std::cout << euclidean_dist << std::endl;
  ROS_DEBUG_STREAM("[NDT_SLAM_ALGORITHM]:uncertainity: "
                   << uncertainty << "  moved_distance: " << euclidean_dist);
  if (/*uncertainty > max_uncertanity_in_window_ ||*/
      euclidean_dist > max_euclidean_distance_traveled_)
    return true;
  else
    return false;
}

bool NdtSlamAlgortihm::movedEnough(const Transform &trans) const
{
  double rotation = std::abs(eigt::getAngle(trans));
  double translation = std::abs(eigt::getDisplacement(trans));
  if (rotation > min_rotation_ || translation > min_distance_)
    return true;
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
    for (int i = 0; i < cell.cov_.rows(); ++i) {
      for (int j = 0; j < cell.cov_.cols(); ++j) {
        cl.cov_matrix.push_back(cell.cov_(i, j));
      }
    }
    ros_msg.cells.push_back(std::move(cl));
  }
  return ros_msg;
}
}  // end of namespace slamuk

#endif
