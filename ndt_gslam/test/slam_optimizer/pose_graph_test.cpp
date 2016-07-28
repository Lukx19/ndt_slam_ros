#include <gtest/gtest.h>
#include <ndt_gslam/slam_optimizer/pose_graph.h>
#include <ndt_gslam/slam_optimizer/slam2d_policy.h>
using namespace slamuk;

typedef Slam2d_Policy P;
typedef std::vector<int> obj_t;
typedef Edge<Slam2d_Policy, obj_t> edge_t;
typedef std::unique_ptr<edge_t> edge_ptr_t;
typedef Node<Slam2d_Policy, obj_t> node_t;
typedef std::unique_ptr<node_t> node_ptr_t;
typedef Graph<Slam2d_Policy, obj_t> graph_t;

void prepareGraph(graph_t &graph)
{
  obj_t obj;
  P::Pose p;
  P::InformMatrix inform;
  for (size_t i = 0; i < 20; ++i) {
    graph.addNode(node_t(p, obj));
  }
  for (size_t i = 0; i < 19; ++i) {
    graph.addEdge(edge_t(&graph.getNode(i), &graph.getNode(i + 1), p, inform));
  }
}

void prepareGraphIncrementaly(graph_t &graph)
{
  obj_t obj;
  P::Pose p;
  P::InformMatrix inform;
  graph.addNode(node_t(p, obj));
  for (size_t i = 0; i < 19; ++i) {
    graph.addNode(node_t(p, obj));
    graph.addEdge(edge_t(&graph.getNode(i), &graph.getNode(i + 1), p, inform));
  }
}

TEST(PoseGraph, buildIncGraph)
{
  try {
    graph_t graph;
    prepareGraphIncrementaly(graph);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(PoseGraph, addingNodesAndEdges)
{
  try {
    graph_t graph;
    P::Pose p;
    P::InformMatrix inform;
    obj_t obj;
    for (size_t i = 0; i < 20; ++i) {
      size_t id = graph.addNode(node_t(p, obj));
      // size_t id = graph.addNode(node_t(p,obj));
      EXPECT_TRUE(id == i);
    }
    for (size_t i = 0; i < 19; ++i) {
      size_t id = graph.addEdge(
          edge_t(&graph.getNode(i), &graph.getNode(i + 1), p, inform));
      EXPECT_TRUE(id == i);
    }

    EXPECT_TRUE(graph.nodeCount() == 20);
    EXPECT_TRUE(graph.edgeCount() == 19);
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(PoseGraph, removingEdgesAndAddBack)
{
  try {
    P::Pose p;
    P::InformMatrix inform;
    graph_t graph;
    prepareGraph(graph);
    EXPECT_TRUE(graph.removeEdge(10));
    EXPECT_TRUE(graph.removeEdge(11));
    EXPECT_FALSE(graph.removeEdge(200));
    EXPECT_FALSE(graph.removeEdge(10));
    EXPECT_EQ(11, graph.addEdge(
                      edge_t(&graph.getNode(1), &graph.getNode(2), p, inform)));
    EXPECT_EQ(10, graph.addEdge(
                      edge_t(&graph.getNode(1), &graph.getNode(2), p, inform)));
    EXPECT_EQ(19, graph.addEdge(
                      edge_t(&graph.getNode(1), &graph.getNode(2), p, inform)));
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(PoseGraph, checkCorrectEdgeRoutes)
{
  try {
    P::Pose p;
    P::InformMatrix inform;
    graph_t graph;
    prepareGraphIncrementaly(graph);
    for (size_t i = 0; i < 19; ++i) {
      EXPECT_EQ(i, graph.getEdge(i).getFrom()->getId());
      EXPECT_EQ(i + 1, graph.getEdge(i).getTo()->getId());
    }
    size_t i = 0;
    for (auto it = graph.cbeginEdge(); it != graph.cendEdge(); ++it) {
      EXPECT_EQ(i, it->getFrom()->getId());
      EXPECT_EQ(i + 1, it->getTo()->getId());
      ++i;
    }
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(PoseGraph, iterators)
{
  graph_t graph;
  prepareGraphIncrementaly(graph);
  size_t i = 0;
  std::vector<size_t> real_idx = {0,  1,  2,  3,  5,  6,  9,
                                  10, 11, 12, 13, 14, 15, 16};
  EXPECT_TRUE(graph.removeEdge(4));
  EXPECT_TRUE(graph.removeEdge(7));
  EXPECT_TRUE(graph.removeEdge(8));
  EXPECT_TRUE(graph.removeEdge(17));
  EXPECT_TRUE(graph.removeEdge(18));
  for (auto it = graph.cbeginEdge(); it != graph.cendEdge(); ++it) {
    EXPECT_EQ(real_idx[i], it->getFrom()->getId());
    ++i;
  }
  i = 0;
  for (auto it = graph.beginEdge(); it != graph.endEdge(); ++it) {
    EXPECT_EQ(real_idx[i], it->getFrom()->getId());
    ++i;
  }
  for (auto it = graph.beginEdge(); it != graph.endEdge(); ++it) {
    it->setId(9999);
    EXPECT_EQ(9999, it->getId());
  }
}

TEST(PoseGraph, inputOutputEdges)
{
  graph_t graph;
  obj_t obj;
  P::Pose p;
  P::InformMatrix inform;
  prepareGraphIncrementaly(graph);
  graph.addEdge(edge_t(&graph.getNode(1), &graph.getNode(3), p, inform));
  graph.addEdge(edge_t(&graph.getNode(3), &graph.getNode(1), p, inform));
  // 1->2->3->4
  // 1->3
  // 3->1
  EXPECT_EQ(2, graph.getNode(3).getEdgesIn().size());
  EXPECT_EQ(2, graph.getNode(3).getEdgesOut().size());

  EXPECT_EQ(1, graph.getNode(1).getEdgesOut()[0]->getId());
  EXPECT_EQ(2, graph.getNode(1).getEdgesOut()[0]->getTo()->getId());
  EXPECT_EQ(1, graph.getNode(1).getEdgesOut()[0]->getFrom()->getId());

  EXPECT_EQ(0, graph.getNode(1).getEdgesIn()[0]->getId());
  EXPECT_EQ(0, graph.getNode(1).getEdgesIn()[0]->getFrom()->getId());
  EXPECT_EQ(1, graph.getNode(1).getEdgesIn()[0]->getTo()->getId());

  EXPECT_EQ(3, graph.getNode(3).getEdgesOut()[0]->getFrom()->getId());
  EXPECT_EQ(4, graph.getNode(3).getEdgesOut()[0]->getTo()->getId());
  EXPECT_EQ(1, graph.getNode(3).getEdgesOut()[1]->getTo()->getId());

  EXPECT_EQ(2, graph.getNode(3).getEdgesIn()[0]->getFrom()->getId());
  EXPECT_EQ(1, graph.getNode(3).getEdgesIn()[1]->getFrom()->getId());
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
