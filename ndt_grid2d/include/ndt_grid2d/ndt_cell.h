#ifndef NDT_GRID2D_NDT_CELL2D
#define NDT_GRID2D_NDT_CELL2D

#include <Eigen/Dense>

namespace slamuk
{
template <class Policy>
class NDTCell
{
  typedef typename Policy::Vector Vector;
  typedef typename Policy::Matrix Matrix;

public:
  NDTCell()
  {
  }

  NDTCell &operator+=(const NDTCell &other);
  bool hasGaussian() const
  {
    return gaussian_;
  }
  size_t points() const
  {
    return points_;
  }
  const Eigen::Vector2d &getMean() const
  {
    return mean_;
  }
  const Eigen::Matrix2d &getCov() const
  {
    return cov_;
  }
  const Eigen::Matrix2d &getICov() const
  {
    return icov_;
  }
  int8_t getOccupancy() const
  {
    if (occup_ > Policy::max_occupancy)
      return 100;
    if (occup_ < 0)
      return 0;
    return static_cast<int8_t>(occup_);
  }

  void addPoint(const Vector &pt)
  {
    points_vec_.push_back(pt);
  }
  void computeGaussian();
  void updateOccupancy(const Vector &start, const Vector &end,
                       size_t new_points);

  void transform(const typename Policy::Transform &trans);

protected:
  Vector mean_;
  Matrix cov_;
  Matrix icov_;
  Vector evals_;
  Matrix evecs_;
  float occup_;
  size_t points_;
  bool gaussian_;
  std::vector<Vector> points_vec_;

  void updateOccupancy(float occup);
  void rescaleCovar();
  double calcMaxLikelihoodOnLine(const Vector &start, const Vector &end,
                                 Vector &pt) const;
};

template <class Policy>
NDTCell<Policy> &NDTCell<Policy>::operator+=(const NDTCell<Policy> &other)
{
  // invalid cell with too little points
  if (other.points() <= 2)
    return *this;
  if (!other.hasGaussian())
    return *this;
  if (this->hasGaussian()) {
    Vector m_sum1 = mean_ * static_cast<double>(points_);
    Vector m_sum2 = other.mean_ * static_cast<double>(other.points_);

    Matrix c_sum1 = cov_ * static_cast<double>(points_ - 1);
    Matrix c_sum2 = other.cov_ * static_cast<double>(other.points_ - 1);

    size_t points_sum = other.points_ + this->points_;

    mean_ = (m_sum2 + m_sum1) / points_sum;

    double w1 = static_cast<double>(this->points_) /
                static_cast<double>(other.points_ * points_sum);
    double w2 =
        static_cast<double>(other.points_) / static_cast<double>(this->points_);
    Matrix c_sum3 =
        c_sum1 + c_sum2 +
        w1 * (w2 * m_sum1 - m_sum2) * (w2 * m_sum1 - m_sum2).transpose();
    points_ = points_sum;
    cov_ = (1.0 / static_cast<double>(points_)) * c_sum3;
    float occup_addition =
        static_cast<float>(other.points_) * Policy::log_like_occ_;
    updateOccupancy(occup_addition);
    if (points_ > Policy::max_points)
      points_ = Policy::max_points;
    if (occup_ < 0)
      gaussian_ = false;
    else
      rescaleCovar();
  } else {
    this->operator=(other);
    gaussian_ = true;
  }

  return *this;
}

template <class Policy>
void NDTCell<Policy>::computeGaussian()
{
  // update occupancy probability
  float occup_addition =
      static_cast<float>(points_vec_.size()) * Policy::log_like_occ_;
  updateOccupancy(occup_addition);
  if (occup_ <= 0) {
    gaussian_ = false;
    return;
  }
  // compute gaussian of new points
  Vector mean_sum2;
  for (int i = 0; i < Policy::dim_; ++i) {
    mean_sum2(i) = 0;
  }
  for (auto &&pt : points_vec_) {
    mean_sum2 += pt;
  }
  Vector mean2 = mean_sum2 / points_vec_.size();

  Eigen::MatrixXd mp;
  mp.resize(points_vec_.size(), Policy::dim_);
  for (int i = 0; i < points_vec_.size(); ++i) {
    mp.col(i) = points_vec_[i] - mean2;
  }
  Matrix cov_sum2 = mp.transpose() * mp;
  Matrix cov2 = cov_sum2 / (points_vec_.size() - 1);

  if (!gaussian_) {
    // cell do not have any gaussian information calculated
    mean_ = mean2;
    cov_ = cov2;
    points_ = points_vec_.size();
    rescaleCovar();
  } else {
    // previously calculated gaussian needs to be updated from points added
    Vector mean_sum1 = mean_ * points_;
    Matrix cov_sum1 = cov_ * (points_ - 1);
    double points1 = static_cast<double>(points_);
    double points2 = static_cast<double>(points_vec_.size());
    double w1 = points1 / (points2 * (points1 + points2));
    double w2 = points2 / points1;
    Vector mean_sum_comb = mean_sum1 + mean_sum2;
    Matrix cov_sum_comb = cov_sum1 + cov_sum2 +
                          w1 * ((w2 * mean_sum1 - mean_sum2) *
                                (w2 * mean_sum1 - mean_sum2).transpose());

    points_ += points_vec_.size();
    // restrict number of points in cell
    if (Policy::max_points > 0 && Policy::max_points_ < points_) {
      mean_sum_comb *= (static_cast<double>(Policy::max_points_) /
                        static_cast<double>(points_));
      cov_sum_comb *= (static_cast<double>(Policy::max_points_ - 1) /
                       static_cast<double>(points_ - 1));
      points_ = Policy::max_points_;
    }
    mean_ = mean_sum_comb / points_;
    cov_ = cov_sum_comb / (points_ - 1);
    this->rescaleCovar();
  }
  points_vec_.clear();
}

template <class Policy>
void NDTCell<Policy>::updateOccupancy(const Vector &start, const Vector &end,
                                      size_t new_points)
{
  if (!gaussian_) {
    updateOccupancy(-0.85 * new_points);
    gaussian_ = false;
  }
  Vector pt_out;
  double lik = calcMaxLikelihoodOnLine(start, end, pt_out);
  double l2_target = (pt_out - end).norm();

  double dist = (start - pt_out).norm();
  double sigma_dist =
      0.5 * (dist / 30.0);  /// test for distance based sensor noise
  double snoise = sigma_dist + Policy::sensor_noise_;
  /// This is the probability of max lik point being endpoint
  double thr = exp(-0.5 * (l2_target * l2_target) / (snoise * snoise));
  lik = lik * (1.0 - thr);
  if (lik < 0.3)
    return;
  lik = 0.1 * lik + 0.5;  /// Evidence value for empty - alpha * p(x);
  double logoddlik = log((1.0 - lik) / (lik));
  updateOccupancy(new_points * logoddlik);
  if (occup_ <= 0) {
    gaussian_ = false;
  }
}

template <class Policy>
void NDTCell<Policy>::transform(const typename Policy::Transform &trans)
{
}
// PROTECTED///////////

template <class Policy>
void NDTCell<Policy>::updateOccupancy(float occup)
{
  // -1 unoccupied, [0,100] - occupied
  occup_ += occup;
  if (occup_ > Policy::max_occupancy_)
    occup_ = Policy::max_occupancy_;
  if (occup_ < -Policy::max_occupancy_)
    occup_ = -Policy::max_occupancy_;
}

template <class Policy>
void NDTCell<Policy>::rescaleCovar()
{
  Eigen::SelfAdjointEigenSolver<Matrix> solver(cov_);
  evecs_ = solver.eigenvectors().real();
  evals_ = solver.eigenvalues().real();
  double eval_factor = 100;
  if (evals_.sum() <= 0) {
    gaussian_ = false;
  } else {
    gaussian_ = true;

    bool recalc = false;
    // guard against near singular matrices::
    int id_max;
    double max_eval = evals_.maxCoeff(&id_max);
    // test if every eigenval is big enough
    for (int i = 0; i < Policy::dim_; ++i) {
      if (max_eval > evals_(i) * eval_factor) {
        evals_(i) = evals_(id_max) / eval_factor;
        recalc = true;
      }
    }
    Matrix lamda;
    lamda = evals_.asDiagonal();
    if (recalc) {
      cov_ = evecs_ * lamda * (evecs_.transpose());
    }
    // compute inverse covariance
    icov_ = evecs_ * (lamda.inverse()) * (evecs_.transpose());
  }
}

template <class Policy>
double NDTCell<Policy>::calcMaxLikelihoodOnLine(const Vector &start,
                                                const Vector &end,
                                                Vector &pt) const
{
  Vector l = (end - start) / (end - start).norm();
  Vector a = icov_ * l;
  Vector b = (end - mean_);
  double sigma = a.cwiseProduct(l).sum();
  if (sigma == 0)
    return 1.0;
  // maximalization of parameter t
  double t = a.cwiseProduct(b).sum() / sigma;

  pt = l * t + end;  // spot on line with maximal likelihood with respect to
                     // gaussian in this cell

  double likelihood = (pt - mean_).dot(icov_ * (pt - mean_));
  if (std::isnan(likelihood))
    return -1;
  return std::exp(-likelihood / 2);
}

}  // end of slamuk namsace
#endif
