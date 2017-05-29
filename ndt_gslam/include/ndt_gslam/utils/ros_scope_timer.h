#ifndef ROS_SCOPE_TIMER_H
#define ROS_SCOPE_TIMER_H
#include <ros/ros.h>
#include <chrono>
#include <iostream>
namespace ros
{
class ScopeTimer
{
public:
  ScopeTimer(const std::string &desc)
    : start_(std::chrono::high_resolution_clock::now()), desc_(desc)
  {
  }
  ~ScopeTimer()
  {
    auto end = std::chrono::high_resolution_clock::now();
    long ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start_)
            .count();
    long sec =
        std::chrono::duration_cast<std::chrono::seconds>(end - start_).count();
    long min = sec / 60;
    long h = min / 60;
    min = min - h * 60;
    sec = sec - min * 60 - h * 3600;
    ROS_DEBUG_STREAM(desc_ << ": took " << ms << "ms which is:" << h << "h "
                           << min << "m " << sec << "s ");
  }

private:
  std::chrono::high_resolution_clock::time_point start_;
  std::string desc_;
};
}
#endif  // ROS_SCOPE_TIMER_H
