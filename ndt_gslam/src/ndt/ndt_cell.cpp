#include <ndt_gslam/ndt/ndt_cell.h>
using namespace slamuk;

NDTCell::NDTCell()
  : mean_(Vector::Zero())
  , cov_(Matrix::Identity())
  , icov_(Matrix::Identity())
  , occup_(0)
  , points_(0)
  , gaussian_(false)
{
}

NDTCell &NDTCell::operator+=(const NDTCell &other)
{
  if (other.occup_ < 0) {
    return *this;
  }
  // invalid cell with too little points
  if (!other.hasGaussian()) {
    //    //    std::cout << "mergig cell with little points \n";
    //    std::copy(other.points_vec_.begin(), other.points_vec_.end(),
    //              std::back_inserter(points_vec_));
    //    computeGaussian();
    return *this;
  }

  if (this->hasGaussian()) {
    auto res = combineGaussians(points_, mean_, cov_, other.points_,
                                other.mean_, other.cov_);
    points_ = std::get<0>(res);
    mean_ = std::get<1>(res);
    cov_ = std::get<2>(res);
    float occup_addition =
        static_cast<float>(other.points_) * LOG_LIKE_OCCUPANCY;
    updateOccupancy(occup_addition);
    if (points_ > MAX_POINTS)
      points_ = MAX_POINTS;
    if (occup_ < 0)
      gaussian_ = false;
    else
      rescaleCovar();
  } else {
    // no previous data. Use calculated Gaussian from other cell
    this->operator=(other);
    gaussian_ = true;
  }

  return *this;
}

void NDTCell::computeGaussian()
{
  // update occupancy probability
  float occup_addition =
      static_cast<float>(points_vec_.size()) * LOG_LIKE_OCCUPANCY;
  updateOccupancy(occup_addition);
  if (occup_ < 0) {
    gaussian_ = false;
    return;
  }
  // compute gaussian of new points
  Vector mean_add = Vector::Zero();
  Matrix cov_add = Matrix::Identity();
  for (auto &&pt : points_vec_) {
    mean_add += pt;
    cov_add += (pt * pt.transpose());
  }
  Vector mean2 = mean_add / static_cast<double>(points_vec_.size());
  if (points_vec_.size() < MIN_POINTS) {
    if (!gaussian_)
      mean_ = mean2;
    return;
  }

  // // simgle pass covariance calculation
  Matrix cov_temp = (cov_add - 2 * (mean_add * mean2.transpose())) /
                        static_cast<double>(points_vec_.size()) +
                    mean2 * mean2.transpose();
  double norm = static_cast<double>(points_vec_.size() - 1) /
                static_cast<double>(points_vec_.size());
  Matrix cov2 = cov_temp * norm;

  if (!gaussian_) {
    // cell do not have any gaussian information calculated
    this->mean_ = mean2;
    this->cov_ = cov2;
    this->points_ = points_vec_.size();
    rescaleCovar();
  } else {
    // previously calculated gaussian needs to be updated from points added
    auto res =
        combineGaussians(points_, mean_, cov_, points_vec_.size(), mean2, cov2);
    points_ = std::get<0>(res);
    mean_ = std::get<1>(res);
    cov_ = std::get<2>(res);
    if (points_ > MAX_POINTS)
      points_ = MAX_POINTS;
    rescaleCovar();
  }
  points_vec_.clear();
}

void NDTCell::updateOccupancy(const Vector &start, const Vector &end,
                              size_t new_points)
{
  if (!gaussian_) {
    updateOccupancy(NO_GAUSSIAN_OCCUP_UPDATE * new_points);
    if (occup_ < 0)
      gaussian_ = false;
    // std::cout << "\nNOTHING\n";
    return;
  }
  // std::cout << "UpdateOccupancy: start: " << start.transpose()
  //           << " end: " << end.transpose() << "adding pts: " << new_points
  // << "cell mean: " << mean_.transpose() << std::endl;
  if ((mean_ - end).cwiseAbs().sum() < 0.01) {
    // we are at the end keep occupancy
    updateOccupancy(MAX_OCCUPANCY);
    // std::cout << "END\n";
    return;
  }

  Vector pt_out;
  double lik = calcMaxLikelihoodOnLine(start, end, pt_out);
  // std::cout << "lik: " << lik << " point: " << pt_out.transpose() <<
  // std::endl;
  double l2_target = (pt_out - end).norm();

  // double dist = (start - pt_out).norm();
  // double sigma_dist =
  //     0.5 * (dist / 30.0);  /// test for distance based sensor noise
  // double snoise = sigma_dist + SENSOR_NOISE;
  double snoise = SENSOR_NOISE;
  /// This is the probability of max lik point being endpoint
  double thr = exp(-0.5 * (l2_target * l2_target) / (snoise * snoise));
  lik = lik * (1.0 - thr);
  // std::cout << "threashold: " << thr << std::endl;
  // if (lik < 0.45) {
  //   // cell is observed as empty
  //   // gaussian_ = false;
  //   // updateOccupancy(0);
  //   return;
  // }
  lik = 0.1 * lik + 0.5;  /// scaling to have lik in range [0,0.5]
  double logodd_lik = log((1.0 - lik) / lik);
  // std::cout << "\nnew occup: " << new_points * logodd_lik << std::endl
  //           << std::endl;
  updateOccupancy(new_points * logodd_lik);
  if (occup_ <= 0) {
    gaussian_ = false;
  }
}

void NDTCell::transform(const Transform &trans)
{
  mean_ = trans * mean_;
  if (gaussian_) {
    cov_ = trans.rotation() * cov_ * trans.rotation().transpose();
    rescaleCovar();
  }
  for (auto &&pt : points_vec_) {
    pt = trans * pt;
  }
}

NDTCellMsg NDTCell::serialize() const
{
  NDTCellMsg msg;
  msg.mean_ = mean_;
  msg.cov_.setIdentity();
  msg.cov_ = cov_;
  msg.occupancy_ = occup_;
  msg.points_ = points_;
  return msg;
}

NDTCell::IntensityCloud NDTCell::sample(size_t samples, float std) const
{
  IntensityCloud cloud;

  if (hasGaussian()) {
    typedef eigt::NormalRandomVariable<float, 3> Generator;
    Generator generator(mean_.cast<float>(), cov_.cast<float>());
    cloud.reserve(samples);
    for (size_t i = 0; i < samples; ++i) {
      Generator::Vector sample = generator(std);
      Eigen::Vector4f point;
      point << sample(0), sample(1), sample(2), occup_;
      cloud.push_back(std::move(point));
    }
  }
  return cloud;
}

// PRIVATE///////////

void NDTCell::updateOccupancy(float occup)
{
  // -1 unoccupied, [0,1] - occupied
  occup_ += occup;
  if (occup_ > MAX_OCCUPANCY)
    occup_ = MAX_OCCUPANCY;
  if (occup_ < MIN_OCCUPANCY)
    occup_ = MIN_OCCUPANCY;
}

void NDTCell::rescaleCovar()
{
  Eigen::SelfAdjointEigenSolver<Matrix> solver(cov_);
  Matrix evacs = solver.eigenvectors();
  Matrix evals = solver.eigenvalues().asDiagonal();
  if (evals(0, 0) < 0 || evals(1, 1) < 0 || evals(2, 2) <= 0)
    gaussian_ = false;

  double min_eval = 0.01 * evals(2, 2);
  if (evals(0, 0) < min_eval)
    evals(0, 0) = min_eval;
  if (evals(1, 1) < min_eval)
    evals(1, 1) = min_eval;
  cov_ = evacs * evals * evacs.transpose();
  icov_ = cov_.inverse();
  if (icov_.maxCoeff() == std::numeric_limits<float>::infinity() ||
      icov_.minCoeff() == -std::numeric_limits<float>::infinity())
    gaussian_ = false;
  gaussian_ = true;
}

double NDTCell::calcMaxLikelihoodOnLine(const Vector &start, const Vector &end,
                                        Vector &pt) const
{
  Vector l = (end - start) / (end - start).norm();
  Vector a = icov_ * l;
  Vector b = (end - mean_);
  double sigma = a.cwiseProduct(l).sum();
  if (sigma == 0)
    return 1.0;
  // maximalization of parameter t
  double t = -a.cwiseProduct(b).sum() / sigma;

  pt = l * t + end;  // spot on line with maximal likelihood with respect to
                     // gaussian in this cell

  if (!gaussian_)
    return -1;
  double likelihood = (pt - mean_).dot(icov_ * (pt - mean_));
  if (std::isnan(likelihood))
    return -1;
  return std::exp(-likelihood / 2);
}

std::tuple<size_t, NDTCell::Vector, NDTCell::Matrix> NDTCell::combineGaussians(
    size_t points1, const NDTCell::Vector &mean1, const NDTCell::Matrix &cov1,
    size_t points2, const NDTCell::Vector &mean2,
    const NDTCell::Matrix &cov2) const
{
  Vector m_sum1 = mean1 * static_cast<double>(points1);
  Vector m_sum2 = mean2 * static_cast<double>(points2);

  Matrix c_sum1 = cov1 * static_cast<double>(points1 - 1);
  Matrix c_sum2 = cov2 * static_cast<double>(points2 - 1);

  size_t points_sum = points2 + points1;

  double w1 =
      static_cast<double>(points1) / static_cast<double>(points2 * points_sum);
  double w2 = static_cast<double>(points2) / static_cast<double>(points1);

  Matrix c_sum3 =
      c_sum1 + c_sum2 +
      w1 * (w2 * m_sum1 - m_sum2) * (w2 * m_sum1 - m_sum2).transpose();

  return std::make_tuple(points1 + points2,
                         (m_sum2 + m_sum1) / static_cast<double>(points_sum),
                         (1.0 / static_cast<double>(points_sum - 1)) * c_sum3);

  //  if (points_ > MAX_POINTS) {
  //    mean_sum_comb *=
  //        (static_cast<double>(MAX_POINTS) / static_cast<double>(points_));
  //    cov_sum_comb *= (static_cast<double>(MAX_POINTS - 1) /
  //                     static_cast<double>(points_ - 1));
  //    points_ = MAX_POINTS;
  //  }
  //  mean_ = mean_sum_comb / points_;
  //  cov_ = cov_sum_comb / (points_ - 1);
}
