#include <fstream>
#include <iostream>
#include <graph_slam_uk/slam_optimizer/dataset_parser.h>
#include <graph_slam_uk/slam_optimizer/gauss_newton_optimalizer2d.h>
#include <graph_slam_uk/slam_optimizer/pcl_ndt_scanmatcher.h>
#include <string>
#include <vector>

using namespace slamuk;
// commind line arguments in this order:
// input file path
// number of iterations
// output file path

int main(int argc, char **argv)
{
  std::cout << "starting parser!\n";
  std::fstream in, out;
  std::vector<std::string> arg(argv, argv + argc);
  if (arg.size() < 3) {
    std::cout << "not enought parameters: dataset_in, iteration, dataset_out\n";
    return 0;
  }
  in.open(arg[1], std::fstream::in);
  out.open(arg[3], std::fstream::out);
  std::string line;
  std::vector<std::string> parts;

  PclNdtScanmatcher scanmatcher;
  GaussNewtonOptimalize2d<pcl::PointCloud<pcl::PointXYZ>::Ptr> optimalizer(
      scanmatcher);
  DatasetIO<pcl::PointCloud<pcl::PointXYZ>::Ptr> d(optimalizer);
  d.loadGraphToro(in);
  optimalizer.setMaxIterations(static_cast<size_t>(std::stoi(arg[2])));

  std::cout << "Graph initialized" << std::endl;
  std::cout << d.loadedNodes() << std::endl;
  std::cout << d.loadedEdges() << std::endl;
  std::chrono::time_point<std::chrono::system_clock> s_calc, e_calc;
  std::chrono::duration<double> elapsed_seconds;

  s_calc = std::chrono::system_clock::now();
  optimalizer.optimalize();
  e_calc = std::chrono::system_clock::now();

  elapsed_seconds = (e_calc - s_calc);
  std::cout << "Total algorithm run: " << elapsed_seconds.count() << std::endl;

  d.saveGraphG2o(out);
  out.close();
  std::ofstream out2;
  out2.open("pgraph.dot");
  optimalizer.getGraphSerialized(out2);
  out2.close();

  return 0;
}
