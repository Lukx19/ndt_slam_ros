#ifndef NDT_GSLAM_UTILS_STRING_TOOLS
#define NDT_GSLAM_UTILS_STRING_TOOLS

#include <string>
#include <vector>

namespace slamuk
{
std::vector<std::string> split(const std::string& data,
                               const std::string& token);

// IMPLEMENTATION
std::vector<std::string> split(const std::string& data,
                               const std::string& token)
{
  std::vector<std::string> output;
  size_t pos = std::string::npos;
  std::string temp = data;
  do {
    pos = temp.find(token);
    output.push_back(temp.substr(0, pos));
    if (std::string::npos != pos)
      temp = temp.substr(pos + token.size());
  } while (std::string::npos != pos);
  return output;
}
}

#endif
