#ifndef NDT_GRID2D_CELL_POLICY2D
#define NDT_GRID2D_CELL_POLICY2D

#include <Eigen/Dense>

class CellPolicy2d{
public:
    typedef Eigen::Vector2d Vector;
    typedef Eigen::Matrix2d Matrix;
    typedef Eigen::Matrix3d Transform;

    size_t dim_ = 2;
    size_t max_points_ = 1e9;
    float max_occupancy_ = 255.0f;
    float sensor_noise_ = 0.01f;
    double log_like_occ_ = std::log(0.6f / (1.0f - 0.6f));

};


#endif