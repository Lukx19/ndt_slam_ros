#include <graph_slam_uk/pcl_ndt_scanmatcher.h>

using namespace slamuk;


MatchResult PclNdtScanmatcher::match(const pcl_constptr_t &source,
                                     const pcl_constptr_t &target,
                                     const Eigen::Matrix3d &initial_guess)
{
    pcl_matcher_.setInputSource (source);
    pcl_matcher_.setInputTarget(target);
    // Set initial alignment estimate found using robot odometry.
    Eigen::Matrix<double, 4, 4>init_guess = eigt::convertFromTransform(eigt::transform2d_t<double>(initial_guess));
     // Calculating required rigid transform to align the input cloud to the target cloud.
     pcl_ptr_t output_cloud (new pcl_t());
     pcl_matcher_.align (*output_cloud,init_guess.cast<float>());

     ROS_INFO_STREAM("PCL_NDT2D:Normal Distributions Transform has converged:" << pcl_matcher_.hasConverged ()
          << " Probability: " << pcl_matcher_.getTransformationProbability());
     MatchResult res;
     res.success_ = pcl_matcher_.hasConverged ();
     res.score_ =  pcl_matcher_.getFitnessScore();
     res.inform_ =pcl_matcher_.getInformMatrix();
     res.transform_ = eigt::convertToTransform<double>(
                      pcl_matcher_.getFinalTransformation().cast<double>());
     return res;
}
void PclNdtScanmatcher::setGridStep(double step)
{
    pcl_matcher_.setResolution(static_cast<float>(1/step));
}
void PclNdtScanmatcher::setMaxRange(double range)
{
}
void PclNdtScanmatcher::setTransformationEpsilon(double epsilon)
{
    pcl_matcher_.setTransformationEpsilon(epsilon);
}
