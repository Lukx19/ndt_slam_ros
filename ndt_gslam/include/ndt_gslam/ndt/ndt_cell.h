#ifndef NDT_GSLAM_NDT_CELL2D
#define NDT_GSLAM_NDT_CELL2D

#include <ndt_gslam/ndt/output_msgs.h>
#include <Eigen/Dense>
#include <iostream>
#include <ostream>

namespace slamuk
{
class NDTCell
{
  constexpr static float MAX_OCCUPANCY = 1;
  constexpr static float MIN_OCCUPANCY = -1;
  constexpr static float LOG_LIKE_OCCUPANCY = 0.405465108;
  constexpr static float SENSOR_NOISE = 0.01;
  constexpr static size_t MAX_POINTS = 1000000000;

  constexpr static float NO_GAUSSIAN_OCCUP_UPDATE = -0.2;

public:
  typedef typename Eigen::Vector3d Vector;
  typedef typename Eigen::Matrix3d Matrix;
  typedef Eigen::Transform<double, 3, Eigen::TransformTraits::Affine> Transform;

public:
  NDTCell();

  /**
   * @brief      {Apply RCU NDT edge merging with occupancy update }
   *
   * @param[in]  other  The other cell
   *
   * @return     {Merged cells}
   */
  NDTCell &operator+=(const NDTCell &other);

  /**
   * @brief      Determines if it has Gaussian distribution inside.
   *
   * @return     True if has Gaussian, False otherwise.
   */
  bool hasGaussian() const
  {
    return gaussian_;
  }
  /**
   * @brief      { Returns number of points used in the cell. }
   *
   * @return     { Number of points. }
   */
  size_t points() const
  {
    return points_;
  }
  /**
   * @brief      Gets cell centroid.
   *
   * @return     The centroid.
   */
  const Vector &getCentroid() const
  {
    return centroid_;
  }

  /**
   * @brief      Sets the centroid.
   *
   * @param[in]  centroid  The centroid
   */
  void setCentroid(const Vector &centroid)
  {
    centroid_ = centroid;
  }

  /**
   * @brief      Gets the mean.
   *
   * @return     The mean.
   */
  const Vector &getMean() const
  {
    return mean_;
  }

  /**
   * @brief      Sets the mean.
   *
   * @param[in]  mean  The mean
   */
  void setMean(const Vector &mean)
  {
    mean_ = mean;
  }
  /**
   * @brief      Gets the covariance.
   *
   * @return     The covariance.
   */
  const Matrix &getCov() const
  {
    return cov_;
  }
  /**
   * @brief      Gets inverse covariance.
   *
   * @return     The inverse covariance.
   */
  const Matrix &getICov() const
  {
    return icov_;
  }
  /**
   * @brief      Gets the occupancy.
   *
   * @return     The occupancy scaled to range[1,100] and -1 for unoccupied.
   */
  int8_t getOccupancy() const
  {
    if (occup_ >= MAX_OCCUPANCY)
      return 100;
    if (occup_ < MIN_OCCUPANCY)
      return 0;
    return static_cast<int8_t>(std::floor(occup_ * 100));
  }
  /**
   * @brief      Gets the raw occupancy.
   *
   * @return     The occupancy in original range [-1,1].
   */
  float getOccupancyRaw() const
  {
    return occup_;
  }
  /**
   * @brief      Sets the occupancy.
   *
   * @param[in]  occup  Values in range [1,100] and -1 for unoccupied.
   */
  void setOccupancy(int8_t occup)
  {
    occup_ = occup;
  }

  /**
   * @brief      Adds a point.
   *
   * @param[in]  pt    The point
   */
  void addPoint(const Vector &pt)
  {
    points_vec_.push_back(pt);
  }
  /**
   * @brief      Calculates the Gaussian.
   */
  void computeGaussian();

  /**
   * @brief      {Update occupancy based on ray-tracing starting from start and
   *             ending in end. If ray passes through normal distribution in
   *             this cell. This cell will receive lower occupancy}
   *
   * @param[in]  start       The start of ray-tracing
   * @param[in]  end         The end of ray-tracing
   * @param[in]  new_points  The number of points in NDT cell at the end of the
   *                         ray. This behaves as weight in computation.
   */
  void updateOccupancy(const Vector &start, const Vector &end,
                       size_t new_points);

  /**
   * @brief      { Transforms normal distribution inside this cell }
   *
   * @param[in]  trans  The transformation
   */
  void transform(const Transform &trans);

  /**
   * @brief      { Returns and cell representation in message format }
   *
   * @return     { Message representing normal distribution of this cell }
   */
  NDTCellMsg serialize() const;

  inline friend std::ostream &operator<<(std::ostream &os, const NDTCell &cell);

  /**
   * @brief      Returns a string representation of the object.
   *
   * @return     String representation of the object.
   */
  std::string toString() const
  {
    std::stringstream os;
    os << "occup: " << occup_ << "gaussian: " << gaussian_
       << " mean: " << mean_.transpose() << " points: " << points_ << std::endl;
    return os.str();
  }

private:
  Vector mean_;
  Matrix cov_;
  Matrix icov_;
  float occup_;
  size_t points_;
  bool gaussian_;
  Vector centroid_;
  std::vector<Vector> points_vec_;

  void updateOccupancy(float occup);
  void rescaleCovar();
  double calcMaxLikelihoodOnLine(const Vector &start, const Vector &end,
                                 Vector &pt) const;
};

std::ostream &operator<<(std::ostream &os, const NDTCell &cell)
{
  os << std::floor(cell.occup_);
  return os;
}
}  // end of slamuk namsace
#endif
