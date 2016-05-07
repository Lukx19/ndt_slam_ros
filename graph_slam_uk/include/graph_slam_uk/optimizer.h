#ifndef GRAPH_SLAM_UK_OPTIMIZER
#define GRAPH_SLAM_UK_OPTIMIZER

#include <graph_slam_uk/pose_graph.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <exception>
#include <iostream>
#include <chrono>
#include <ros/ros.h>

namespace slamuk
{
template <typename P, typename T>
class Optimizer
{
  typedef Eigen::Triplet<double> Triple;
  typedef std::vector<Triple> BlockBuffer;
  typedef typename P::JacobianMatrix HBlock;  // hessian small block type
public:
  Optimizer()
  {
  }
  bool optimizeGraph(Graph<P, T> &g, double epsilon, size_t iterations);
  double calcTotalError(const Graph<P, T> &g) const;
  double calcMaxUpdate(Eigen::VectorXf &deltax);

protected:
  void addToHessian(BlockBuffer &buff, HBlock &hessian, size_t i, size_t j);
};

template <typename P, typename T>
bool Optimizer<P, T>::optimizeGraph(Graph<P, T> &g, double epsilon,
                                    size_t iterations)
{
  std::chrono::time_point<std::chrono::system_clock> s_upd, e_upd, s_hess,
      e_hess, s_calc, e_calc, s_inside, e_inside;
  std::chrono::duration<double> elapsed_seconds;
  double total_time_inside = 0;

  typedef Eigen::SparseMatrix<double> H;  // big sparse hessian matrix type

  size_t block_dim = static_cast<int>(P::block_size);
  int hessian_size = static_cast<int>(g.nodeCount() * P::block_size);

  H hessian(hessian_size, hessian_size);
  Eigen::VectorXd b;
  b.resize(hessian_size);
  std::vector<Triple> feeder;
  feeder.reserve(g.edgeCount() * block_dim * block_dim * 2 * 2 +
                 1000);  // optimalize this allocation
  typename P::ErrorVector err;
  typename P::JacobianPair jacobians;

  HBlock hessian_ii, hessian_ij, hessian_jj, hessian_ji, A_ij_inf, B_ij_inf;
  bool is_analyzed = false;
  Eigen::SimplicialLDLT<H> solver;
  std::vector<size_t> priors = {0};
  for (size_t iter = 0; iter < iterations; ++iter) {
    // initzialing all structures to be empty
    hessian.setZero();
    b.setZero();
    feeder.clear();

    s_inside = std::chrono::system_clock::now();
    for (auto it = g.beginEdge(); it != g.endEdge(); ++it) {
      if (it->getState() != Edge<P, T>::State::ACTIVE)
        continue;

      size_t i = it->getFrom()->getId();
      size_t j = it->getTo()->getId();

      bool i_prior =
          std::find(priors.begin(), priors.end(), i) == priors.end() ? false :
                                                                       true;
      bool j_prior =
          std::find(priors.begin(), priors.end(), j) == priors.end() ? false :
                                                                       true;

      if (i_prior && j_prior)
        continue;
      // calculate error function

      err = std::move(it->getError());

      // calculate jaccobians

      jacobians = move(it->getJacobianBlocks());
      // pre-calculate just (Jij)transposed*information matrix

      A_ij_inf = jacobians.first.transpose() * (it->getInformationMatrix());
      B_ij_inf = jacobians.second.transpose() * (it->getInformationMatrix());

      if (i_prior) {
        hessian_jj = B_ij_inf * jacobians.second;
        // add hessian_jj to matrix with all hessians
        addToHessian(feeder, hessian_jj, j, j);
        b.segment(j * block_dim, block_dim) +=
            (B_ij_inf * err);  //.segment(0,block_dim);
      } else if (j_prior) {
        hessian_ii = A_ij_inf * jacobians.first;
        // add hessian_ii to matrix with all hessians
        addToHessian(feeder, hessian_ii, i, i);
        b.segment(i * block_dim, block_dim) +=
            (A_ij_inf * err);  //.segment(0,block_dim);
      } else {
        hessian_ii = A_ij_inf * jacobians.first;
        hessian_ij = A_ij_inf * jacobians.second;
        hessian_jj = B_ij_inf * jacobians.second;
        hessian_ji = hessian_ij.transpose();

        // add hessian_ii to matrix with all hessians
        addToHessian(feeder, hessian_ii, i, i);
        // add hessian_jj to matrix with all hessians
        addToHessian(feeder, hessian_jj, j, j);
        // add hessian_ij to matrix with all hessians
        addToHessian(feeder, hessian_ij, i, j);
        // add hessian_ji to matrix with all hessians
        addToHessian(feeder, hessian_ji, j, i);

        // calculate b_i and b_j which are two blocks in whole b vector.
        // and add this parts to b vector
        b.segment(i * block_dim, block_dim) +=
            (A_ij_inf * err);  //.segment(0,block_dim);
        b.segment(j * block_dim, block_dim) +=
            (B_ij_inf * err);  //.segment(0,block_dim);
      }
    }
    e_inside = std::chrono::system_clock::now();
    // fix first node in hessian by adding identity matrix to first block of
    // hessian
    HBlock identity;
    identity.setIdentity();
    addToHessian(feeder, identity, 0, 0);

    // solve the system
    s_hess = std::chrono::system_clock::now();
    hessian.setFromTriplets(feeder.begin(), feeder.end());
    e_hess = std::chrono::system_clock::now();
    std::cout << "HESIAN created" << std::endl;
    // std::cout<<hessian.cols()<<std::endl;
    //  std::cout<<hessian.nonZeros()<<std::endl;
    std::cout << "feeder size" << feeder.size() << std::endl;
    std::cout << "Sparcity in percents: "
              << (double)hessian.nonZeros() /
                     (double)(hessian.rows() * hessian.cols()) * 100
              << std::endl;
    // Eigen::SparseQR<H,Eigen::COLAMDOrdering<int>> solver;

    //  std::cout<<b<<std::endl;

    // Eigen::ConjugateGradient<H> solver;
    // Eigen::BiCGSTAB<H> solver;
    s_calc = std::chrono::system_clock::now();

    // if(!is_analyzed){
    //    solver.analyzePattern(hessian);
    //    is_analyzed = true;
    //}

    solver.compute(hessian);
    if (solver.info() != Eigen::Success) {
      ROS_ERROR("Error in solving linear system-hessian ");
      if (iter > 0)
        return true;
      return false;
    }
    // std::cout<<"HESIAN computed"<<std::endl;
    Eigen::VectorXd deltax = solver.solve(-b);
    if (solver.info() != Eigen::Success) {
      ROS_ERROR("Error in solving linear system-b vector");
      if (iter > 0)
        return true;
      return false;
    }
    // updating all node's positions based on result from solver
    // std::cout<<deltax.segment(0,3)<<std::endl;
    e_calc = std::chrono::system_clock::now();
    s_upd = std::chrono::system_clock::now();
    size_t nd = 0;
    for (auto it = g.beginNode(); it != g.endNode(); ++it) {
      typename P::Pose pose = P::Pose::Zero();
      // pose.segment(0,block_dim)=deltax.segment(nd,block_dim);
      it->addToPose(std::move(deltax.segment(nd, block_dim)));
      nd += block_dim;
    }
    e_upd = std::chrono::system_clock::now();
    // end calculation if change is too small;
    if (deltax.norm() < epsilon)
      break;
    std::cout << "Total error:" << calcTotalError(g) << std::endl;
    elapsed_seconds = (e_hess - s_hess);
    std::cout << "Total Hessian build time: " << elapsed_seconds.count()
              << std::endl;
    elapsed_seconds = (e_inside - s_inside);
    std::cout << "Total Inside  time: " << elapsed_seconds.count() << std::endl;
    elapsed_seconds = (e_calc - s_calc);
    std::cout << "Total solve time: " << elapsed_seconds.count() << std::endl;
    elapsed_seconds = (e_upd - s_upd);
    std::cout << "Total update time: " << elapsed_seconds.count() << std::endl;
    std::cout << std::endl;
  }
  std::cout << "Calc done" << std::endl;
  return true;
}

template <typename P, typename T>
double Optimizer<P, T>::calcTotalError(const Graph<P, T> &g) const
{
  double total_err = 0;
  for (auto it = g.cbeginEdge(); it != g.cendEdge(); ++it) {
    typename P::ErrorVector err = it->getError();
    total_err += err.transpose() * it->getInformationMatrix() * err;
    // std::cout<<"err"<<it->getError()<<std::endl;
  }

  return total_err;
}

template <typename P, typename T>
double Optimizer<P, T>::calcMaxUpdate(Eigen::VectorXf &deltax)
{
  size_t block_dim = static_cast<int>(P::block_size);
  std::vector<double> displace;
  displace.reserve(deltax.size() / block_dim);
  for (size_t i = 0; i < deltax.size() / block_dim; i += block_dim) {
    displace.push_back(deltax.segment(i, block_dim).lpNorm<1>());
  }
  return *std::max_element(displace.begin(), displace.end());
}

template <typename P, typename T>
void Optimizer<P, T>::addToHessian(BlockBuffer &buff, HBlock &hessian, size_t i,
                                   size_t j)
{
  // add hessian to  BlockBuffer buff with specific possition i , j
  size_t block_dim = static_cast<int>(P::block_size);
  for (size_t row = 0; row < block_dim; ++row) {
    for (size_t col = 0; col < block_dim; ++col) {
      buff.emplace_back(
          Triple(i * block_dim + row, j * block_dim + col, hessian(row, col)));
    }
  }
}
} // end of namespace slamuk
#endif
