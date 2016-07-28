#ifndef EDGE_SE3SWITCHABLE_PRIOR_H_
#define EDGE_SE3SWITCHABLE_PRIOR_H_

#include <g2o/core/base_unary_edge.h>
#include <ndt_gslam/slam_optimizer/switchable_constraints/vertex_switchLinear.h>

class EdgeSwitchPrior : public g2o::BaseUnaryEdge<1, double, VertexSwitchLinear>
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  EdgeSwitchPrior(){};

  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  void linearizeOplus();
  void computeError();
};
#endif