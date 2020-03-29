// Copyright 2018-present, Simone Fontana
// Distributed under the GNU GPL 3.0 License (https://www.gnu.org/licenses/gpl-3.0.html)

#ifndef PSO_REGISTRATION_UTILITIES_HPP_
#define PSO_REGISTRATION_UTILITIES_HPP_

#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "pcl/kdtree/kdtree.h"
#include "pcl/point_cloud.h"

#include "pso_registration/probabilistic_weights.hpp"

namespace pso_registration {
inline double pi() { return std::atan(1) * 4; }

inline double l1_distance(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2,
                          pcl::KdTree<pcl::PointXYZ>::Ptr cloud2_tree) {
  double sum = 0;
  for (std::size_t i = 0; i < cloud1->size(); i++) {
    std::vector<int> neighbours;
    std::vector<float> distances;
    neighbours.reserve(1);
    distances.reserve(1);
    cloud2_tree->nearestKSearch(*cloud1, i, 1, neighbours, distances);
    auto closest = cloud2->at(neighbours[0]);
    double diff_x = cloud1->at(i).x - closest.x, diff_y = cloud1->at(i).y - closest.y,
           diff_z = cloud1->at(i).z - closest.z;
    sum += abs(diff_x) + abs(diff_y) + abs(diff_z);
  }
  return sum;
}

inline double l2_distance(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2,
                          pcl::KdTree<pcl::PointXYZ>::Ptr cloud2_tree) {
  double sum = 0;
  int num_neig = 10;
  prob_point_cloud_registration::ProbabilisticWeights weight_updater(std::numeric_limits<double>::infinity(), 3,
                                                                     num_neig);
  Eigen::SparseMatrix<double, Eigen::RowMajor> data_association(cloud1->size(), cloud2->size());
  std::vector<Eigen::Triplet<double>> tripletList;
  std::vector<double> squared_errors;
  squared_errors.reserve(data_association.cols() * num_neig);
  for (std::size_t i = 0; i < cloud1->size(); i++) {
    std::vector<int> neighbours;
    std::vector<float> distances;
    cloud2_tree->nearestKSearch(*cloud1, i, num_neig, neighbours, distances);
    int k = 0;
    for (int j : neighbours) {
      tripletList.push_back(Eigen::Triplet<double>(i, j, distances[k]));
      squared_errors.push_back(distances[k] * distances[k]);
      k++;
    }
  }
  data_association.setFromTriplets(tripletList.begin(), tripletList.end());
  data_association.makeCompressed();
  Eigen::SparseMatrix<double, Eigen::RowMajor> weights =
      weight_updater.updateWeights(data_association, squared_errors);
  for (int i = 0, k = 0; i < weights.outerSize(); ++i) {
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(weights, i); it; ++it) {
      sum += it.value() * squared_errors[k];
      k++;
    }
  }
  return sum;
}

inline double robust_l2_distance(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2,
                                 pcl::KdTree<pcl::PointXYZ>::Ptr cloud2_tree) {
  double sum = 0;
  double median_distance;
  std::vector<double> all_distances;
  for (std::size_t i = 0; i < cloud1->size(); i++) {
    std::vector<int> neighbours;
    std::vector<float> distances;
    neighbours.reserve(1);
    distances.reserve(1);
    cloud2_tree->nearestKSearch(*cloud1, i, 1, neighbours, distances);
    all_distances.push_back(distances[0]);
  }

  std::sort(all_distances.begin(), all_distances.end());
  if (all_distances.size() % 2 != 0) {
    median_distance = all_distances[(all_distances.size() + 1) / 2];
  } else {
    median_distance = (all_distances[all_distances.size() / 2] + all_distances[(all_distances.size() / 2) + 1]) / 2.0;
  }
  int num_filtered = 0;
  for (auto it = all_distances.begin(); it != all_distances.end(); it++) {
    if (*it <= median_distance * 3 && *it >= median_distance / 3) {
      sum += (*it);
      num_filtered++;
    }
  }
  if (num_filtered < 10) {
    return std::numeric_limits<double>::max();
  }
  return sum;
}

inline double robust_l2_distance(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2,
                                 pcl::KdTree<pcl::PointXYZ>::Ptr cloud2_tree, double factor) {
  double sum = 0;
  double median_distance;
  std::vector<double> all_distances;
  for (std::size_t i = 0; i < cloud1->size(); i++) {
    std::vector<int> neighbours;
    std::vector<float> distances;
    neighbours.reserve(1);
    distances.reserve(1);
    cloud2_tree->nearestKSearch(*cloud1, i, 1, neighbours, distances);
    all_distances.push_back(distances[0]);
  }

  std::sort(all_distances.begin(), all_distances.end());
  if (all_distances.size() % 2 != 0) {
    median_distance = all_distances[(all_distances.size() + 1) / 2];
  } else {
    median_distance = (all_distances[all_distances.size() / 2] + all_distances[(all_distances.size() / 2) + 1]) / 2.0;
  }
  int num_filtered = 0;
  for (auto it = all_distances.begin(); it != all_distances.end(); it++) {
    if (*it <= median_distance * factor && *it >= median_distance / factor) {
      sum += (*it);
      num_filtered++;
    }
  }
  if (num_filtered < 10) {
    return std::numeric_limits<double>::max();
  }
  return sum;
}

inline double robust_normalized_l2_distance(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1,
                                            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2,
                                            pcl::KdTree<pcl::PointXYZ>::Ptr cloud2_tree) {
  double sum = 0;
  double median_distance;
  std::vector<double> all_distances;
  for (std::size_t i = 0; i < cloud1->size(); i++) {
    std::vector<int> neighbours;
    std::vector<float> distances;
    neighbours.reserve(1);
    distances.reserve(1);
    cloud2_tree->nearestKSearch(*cloud1, i, 1, neighbours, distances);
    all_distances.push_back(distances[0]);
  }

  std::sort(all_distances.begin(), all_distances.end());
  if (all_distances.size() % 2 != 0) {
    median_distance = all_distances[(all_distances.size() + 1) / 2];
  } else {
    median_distance = (all_distances[all_distances.size() / 2] + all_distances[(all_distances.size() / 2) + 1]) / 2.0;
  }
  int num_filtered = 0;
  for (auto it = all_distances.begin(); it != all_distances.end(); it++) {
    if (*it <= median_distance * 3 && *it >= median_distance / 3) {
      sum += (*it);
      num_filtered++;
    }
  }
  if (num_filtered < (cloud2->size() * 0.01)) {
    return std::numeric_limits<double>::max();
  }
  return sum / num_filtered;
}

inline Eigen::Quaterniond euler2Quaternion(const double roll, const double pitch, const double yaw) {
  Eigen::AngleAxisd rollAngle(roll, Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd pitchAngle(pitch, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd yawAngle(yaw, Eigen::Vector3d::UnitZ());

  Eigen::Quaterniond q = yawAngle * pitchAngle * rollAngle;
  return q;
}

}  // namespace pso_registration

#endif  // PSO_REGISTRATION_UTILITIES_HPP_
