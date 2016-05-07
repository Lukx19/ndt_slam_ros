#ifndef GRAPH_SLAM_UK_DATASET_PARSER
#define GRAPH_SLAM_UK_DATASET_PARSER

#include <vector>
#include <graph_slam_uk/graph_slam_interfaces.h>

namespace slamuk
{
template <typename T>
class DatasetIO
{
public:
  DatasetIO(IGraphOptimalizer2d<T>& optimalizer) : opt_engine_(&optimalizer)
  {
  }
  void loadGraphToro(std::istream& inp);
  void loadGraphG2o(std::istream& inp);
  void saveGraphTORO(std::ostream& out);
  void saveGraphG2o(std::ostream& out);

  size_t loadedNodes() const
  {
    return node_ids_.size();
  }
  size_t loadedEdges() const
  {
    return edge_ids_.size();
  }

private:
  IGraphOptimalizer2d<T>* opt_engine_;
  T obj_;
  std::vector<size_t> edge_ids_;
  std::vector<size_t> node_ids_;
  void tokenize(const std::string& str, std::vector<std::string>& tokens,
                const std::string& delimiters = " ");
};

template <typename T>
void DatasetIO<T>::loadGraphToro(std::istream& inp)
{
  std::string line;
  std::vector<std::string> parts;
  edge_ids_.clear();
  node_ids_.clear();
  while (std::getline(inp, line)) {
    tokenize(line, parts, " ");
    if (parts[0] == "VERTEX2") {
      Eigen::Vector3d p;
      p << std::stof(parts[2]), std::stof(parts[3]), std::stof(parts[4]);
      node_ids_.push_back(opt_engine_->addPose(p, obj_));
    } else if (parts[0] == "EDGE2") {
      Eigen::Vector3d p;
      Eigen::Matrix3d inf;
      p << std::stof(parts[3]), std::stof(parts[4]), std::stof(parts[5]);
      inf << std::stof(parts[6]), std::stof(parts[7]), std::stof(parts[10]),
          std::stof(parts[7]), std::stof(parts[8]), std::stof(parts[11]),
          std::stof(parts[10]), std::stof(parts[11]), std::stof(parts[9]);
      edge_ids_.push_back(opt_engine_->addConstrain(
          std::stoi(parts[1]), std::stoi(parts[2]), p, inf.inverse()));
    }

    parts.clear();
    line.clear();
  }
}

template <typename T>
void DatasetIO<T>::loadGraphG2o(std::istream& inp)
{
  std::string line;
  std::vector<std::string> parts;
  edge_ids_.clear();
  node_ids_.clear();
  while (std::getline(inp, line)) {
    tokenize(line, parts, " ");
    if (parts[0] == "VERTEX2_SE2") {
      Eigen::Vector3d p;
      p << std::stof(parts[2]), std::stof(parts[3]), std::stof(parts[4]);
      node_ids_.push_back(opt_engine_->addPose(p, obj_));
    } else if (parts[0] == "EDGE2_SE2") {
      Eigen::Vector3d p;
      Eigen::Matrix3d inf;
      p << std::stof(parts[3]), std::stof(parts[4]), std::stof(parts[5]);
      inf << std::stof(parts[6]), std::stof(parts[7]), std::stof(parts[8]),
          std::stof(parts[7]), std::stof(parts[9]), std::stof(parts[10]),
          std::stof(parts[8]), std::stof(parts[10]), std::stof(parts[11]);
      edge_ids_.push_back(opt_engine_->addConstrain(
          std::stoi(parts[1]), std::stoi(parts[2]), p, inf.inverse()));
    }

    parts.clear();
    line.clear();
  }
}

template <typename T>
void DatasetIO<T>::saveGraphG2o(std::ostream& out)
{
  for (auto id : node_ids_) {
    Eigen::Vector3d p = opt_engine_->getPoseLocation(id);
    p.transpose();
    out << "VERTEX_SE2 " << id << " ";
    out << p(0) << " " << p(1) << " " << p(2) << " ";
    out << std::endl;
  }

  for (size_t id : edge_ids_) {
    Eigen::Vector3d p = opt_engine_->getConstrainTransform(id);
    Eigen::Matrix3d inf = opt_engine_->getConstrainInformMat(id);
    p.transpose();
    std::pair<size_t, size_t> nodes = opt_engine_->getConstrainPoses(id);
    out << "EDGE_SE2 " << nodes.first << " " << nodes.second << " ";
    out << p(0) << " " << p(1) << " " << p(2) << " ";
    out << inf(0, 0) << " " << inf(0, 1) << " " << inf(0, 2) << " " << inf(1, 1)
        << " " << inf(1, 2) << " " << inf(2, 2);
    out << std::endl;
  }
}

template <typename T>
void DatasetIO<T>::tokenize(const std::string& str,
                            std::vector<std::string>& tokens,
                            const std::string& delimiters)
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}
}
#endif
