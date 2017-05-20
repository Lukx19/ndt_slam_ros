#ifndef NDT_GSLAM_GRAPH_SLAM2D_ALG
#define NDT_GSLAM_GRAPH_SLAM2D_ALG

#include <ndt_gslam/ndt_grid2d_holder.h>
#include <ndt_gslam/slam_optimizer/graph_slam_interfaces.h>
#include <ndt_gslam/slam_optimizer/ndt_scanmatcher.h>
#include <ndt_gslam/slam_optimizer/slam2d.h>
#include <ndt_gslam/utils/covariance_wrapper.h>
#include <ndt_gslam/utils/eigen_tools.h>
#include <ndt_gslam/utils/point_cloud_tools.h>
#include <opencv/cv.h>
#include <ros/ros.h>
#include <Eigen/Dense>

#include <nav_msgs/OccupancyGrid.h>
#include <ndt_gslam/ndt/ndt_cell.h>
#include <ndt_gslam/ndt/ndt_mapper.h>
#include <ndt_gslam/ndt/output_msgs.h>
#include <ndt_gslam/registration/d2d_ndt2d.h>
#include <ndt_gslam/registration/ndt2d.h>
#include <ndt_gslam/slam_algorithm_interface.h>
#include <pcl_ros/point_cloud.h>

#include <atomic>
#include <deque>
#include <mutex>
#include <thread>

namespace slamuk
{
class NdtSlamAlgortihm : public ISlamAlgorithm
{
  typedef ISlamAlgorithm::PointType PointType;
  typedef NDTCell CellType;
  typedef NDTGrid2D<NDTCell, PointType> FrameType;
  typedef FrameType::Ptr FrameTypePtr;
  typedef NDTGrid2DHolder<NDTCell, PointType> FrameTypeHolder;
  typedef NdtScanmatcher<FrameTypeHolder> LoopScanmatcher;

public:
  typedef ISlamAlgorithm::Transform Transform;
  typedef ISlamAlgorithm::Covar Covar;
  typedef ISlamAlgorithm::PointCloud PointCloud;
  typedef ISlamAlgorithm::Pose Pose;

private:
  struct BufferedFrame {
    FrameTypePtr frame_;
    ros::Time stamp_;
    Transform tf_last_to_this_;
    CovarianceWrapper cov_last_to_this_;
  };

public:
  NdtSlamAlgortihm(ros::NodeHandle &nh_private);
  ~NdtSlamAlgortihm();

  virtual inline Pose update(const Transform &motion, const Covar &covariance,
                             const PointCloud &pcl,
                             const ros::Time &update_time) override;
  virtual nav_msgs::OccupancyGrid
  getOccupancyGrid(const std::string &world_frame_id) override;

  virtual visualization_msgs::MarkerArray
  getGraphSerialized(const std::string &world_frame_id) override
  {
    return opt_engine_->getGraphSerialized(world_frame_id);
  }

  virtual ndt_gslam::NDTMapMsg
  getNDTMap(const std::string &world_frame_id) override;

protected:
  bool initialized_;
  // PARAMETERS
  double min_distance_;
  double min_rotation_;
  double max_uncertanity_in_window_;
  double max_euclidean_distance_traveled_;
  float cell_size_;

  size_t last_node_id_;
  Transform unused_trans_;  // transformation of unused odometry

  Transform position_;
  Transform position_cumul_;
  Transform transform_temp_;  // transform of current temp frame
  Pose last_matcher_pose_;
  CovarianceWrapper covar_temp_;
  NDTMapper<CellType, PointType> map_;
  FrameTypePtr frame_temp_;

  std::unique_ptr<IGraphOptimalizer2d<FrameTypeHolder>> opt_engine_;

  unsigned int window_seq_;
  ros::Time window_update_time_;
  std::vector<FrameTypePtr> frames_;

  std::deque<BufferedFrame> built_frames_;
  std::mutex built_frames_mutex_;
  std::mutex position_mutex_;
  std::mutex map_mutex_;
  std::atomic<bool> exit_;
  std::thread backend_thread_;

  inline void updateGraph(const Transform &increase, const Covar &covariance,
                          const PointCloud &new_scan,
                          const ros::Time &update_time);
  inline bool movedEnough(const Transform &trans) const;
  inline bool movedEnough(const Transform &trans,
                          const CovarianceWrapper &covar) const;
  inline ndt_gslam::NDTMapMsg toNdtMapMsg(const NDTGridMsg &msg,
                                          const std::string &fixed_frame) const;

  void runBackend();
};
NdtSlamAlgortihm::NdtSlamAlgortihm(ros::NodeHandle &nh_private)
  : initialized_(false)
  , min_distance_(0.2)
  , min_rotation_(0.1)
  , max_uncertanity_in_window_(1000000.0)
  , max_euclidean_distance_traveled_(1.5)
  , cell_size_(0.25)
  , last_node_id_(0)
  , unused_trans_(Transform::Identity())
  , position_(Transform::Identity())
  , position_cumul_(Transform::Identity())
  , transform_temp_(Transform::Identity())
  , last_matcher_pose_(Pose::Zero())
  , covar_temp_()
  , map_(cell_size_)
  , frame_temp_()
  , opt_engine_(new Slam2D<FrameTypeHolder>(
        std::unique_ptr<LoopScanmatcher>(new LoopScanmatcher())))
  , window_seq_(0)
  , exit_(false)
  , backend_thread_(&NdtSlamAlgortihm::runBackend, this)
{
  max_euclidean_distance_traveled_ =
      static_cast<float>(nh_private.param<double>("node_gen_distance", 2.0));

  float loop_max_dist =
      static_cast<float>(nh_private.param<double>("loop_max_distance", 30.0));
  opt_engine_->setLoopGenerationMaxDist(loop_max_dist);

  float loop_min_dist =
      static_cast<float>(nh_private.param<double>("loop_min_distance", 14.0));
  opt_engine_->setLoopGenerationMinDist(loop_min_dist);

  float loop_score_threshold =
      static_cast<float>(nh_private.param<double>("loop_score_threshold", 0.6));
  opt_engine_->setLoopRegistrationScore(loop_score_threshold);

  cell_size_ = nh_private.param<double>("cell_size", 0.25);

  frame_temp_.reset(new FrameType(cell_size_, FrameType::Pose(0, 0, 0)));
  map_ = NDTMapper<CellType, PointType>(cell_size_);
}

NdtSlamAlgortihm::~NdtSlamAlgortihm()
{
  exit_.store(true);
  backend_thread_.join();
}

NdtSlamAlgortihm::Pose NdtSlamAlgortihm::update(const Transform &motion,
                                                const Covar &covariance,
                                                const PointCloud &pcl,
                                                const ros::Time &update_time)
{
  if (pcl.empty()) {
    ROS_WARN_STREAM("[SLAM ALGORITHM]: Input point cloud is empty.");
    position_ = position_ * motion;
    return eigt::getPoseFromTransform(position_ * motion);
  }

  if (!initialized_) {
    initialized_ = true;
    // enters first scan as it is to running window
    frame_temp_->initializeSimple(pcl);
    window_update_time_ = update_time;
    window_seq_ = 0;
    ROS_INFO_STREAM("[SLAM ALGORITHM]: NDT Slam algorithm initialized!");
    return eigt::getPoseFromTransform(position_);
  }
  transform_temp_ = transform_temp_ * motion;
  covar_temp_.addToCovar(covariance, motion);

  {
    std::lock_guard<std::mutex> lock(position_mutex_);
    position_ = position_ * motion;
  }

  if (movedEnough(transform_temp_, covar_temp_)) {
    ROS_INFO_STREAM("[SLAM ALGORITHM]: creating new frame");
    {
      std::lock_guard<std::mutex> lock(built_frames_mutex_);
      BufferedFrame frame;
      frame.frame_.reset(new FrameType(cell_size_));
      frame.stamp_ = update_time;
      frame.frame_.swap(frame_temp_);
      frame.cov_last_to_this_ = covar_temp_;
      frame.tf_last_to_this_ = transform_temp_;
      built_frames_.emplace_back(std::move(frame));
    }
    covar_temp_ = CovarianceWrapper();
    transform_temp_.setIdentity();
    frame_temp_->initialize(pcl);
    {
      std::lock_guard<std::mutex> lock(position_mutex_);
      frame_temp_->setOrigin(eigt::getPoseFromTransform(position_));
    }

  } else {
    PointCloud pcl_out;
    pcl::transformPointCloud(pcl, pcl_out,
                             eigt::convertFromTransform(transform_temp_));
    frame_temp_->mergeInTraced(pcl_out, frame_temp_->getOrigin(), true);
  }
  return eigt::getPoseFromTransform(position_);
}

nav_msgs::OccupancyGrid
NdtSlamAlgortihm::getOccupancyGrid(const std::string &world_frame_id)
{
  std::lock_guard<std::mutex> lock(map_mutex_);
  auto map = map_.calcOccupancyGridMsg();
  map.header.frame_id = world_frame_id;
  map.header.stamp = ros::Time::now();
  return map;
}

ndt_gslam::NDTMapMsg
NdtSlamAlgortihm::getNDTMap(const std::string &world_frame_id)
{
  return toNdtMapMsg(map_.getNDTMap()->serialize(), world_frame_id);
}

bool NdtSlamAlgortihm::movedEnough(const Transform &trans,
                                   const CovarianceWrapper &covar) const
{
  double uncertainty = (covar.covar_.determinant());
  // std::cout << uncertainty << std::endl;

  double euclidean_dist = std::abs(eigt::getDisplacement(trans));
  // std::cout << euclidean_dist << std::endl;
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

ndt_gslam::NDTMapMsg NdtSlamAlgortihm::toNdtMapMsg(
    const NDTGridMsg &msg, const std::string &fixed_frame) const
{
  ndt_gslam::NDTMapMsg ros_msg;
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
    ndt_gslam::NDTCellMsg cl;
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

void NdtSlamAlgortihm::runBackend()
{
  bool first_node_added = false;
  size_t last_node_id = 0;
  while (!exit_.load()) {
    if (built_frames_.size() > 0) {
      ROS_INFO_STREAM("NDTSLAM_ALGORITHM: frame used. Remaining:"
                      << built_frames_.size() - 1);
      BufferedFrame fr;
      {
        std::lock_guard<std::mutex> lock(built_frames_mutex_);
        fr = built_frames_.front();
        built_frames_.pop_front();
      }
      frames_.emplace_back(fr.frame_);
      {
        std::lock_guard<std::mutex> lock(map_mutex_);
        map_.addFrame(frames_.back(), fr.stamp_);
      }

      FrameTypeHolder frame_wrap(frames_.back());

      last_node_id =
          opt_engine_->addPose(frames_.back()->getOrigin(), frame_wrap);
      // add edge measurement only after first node was inserted
      if (first_node_added) {
        opt_engine_->addLastConstrain(
            eigt::getPoseFromTransform(fr.tf_last_to_this_),
            fr.cov_last_to_this_.covar_.inverse());
      }
      first_node_added = true;
      // find loop closures and optimize graph
      if (opt_engine_->tryLoopClose()) {
        // optimize graph
        //        if (opt_engine_->optimalizeIterationaly()) {
        //          // optimization updates origins in grid through FrameHolder
        //          // interface
        //          {
        //            std::lock_guard<std::mutex> lock(map_mutex_);
        //            map_.recalc(fr.stamp_);
        //          }
        //          {
        //            std::lock_guard<std::mutex> lock(position_mutex_);
        //            position_ = eigt::getTransFromPose(
        //                opt_engine_->getPoseLocation(last_node_id));
        //          }
        //        }
      }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }
}
}  // end of namespace slamuk

#endif
