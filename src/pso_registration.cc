// Copyright 2018-present, Simone Fontana
// Distributed under the GNU GPL 3.0 License (https://www.gnu.org/licenses/gpl-3.0.html)

#include <Eigen/Core>
#include <memory>
#include <string>

#include "pcl/common/transforms.h"
#include "pcl/filters/random_sample.h"
#include "pcl/filters/voxel_grid.h"
#include "pcl/io/pcd_io.h"
#include "pcl/kdtree/kdtree_flann.h"
#include "pcl/point_types.h"
#include "pcl/visualization/pcl_visualizer.h"
#include "spdlog/spdlog.h"
#include "tclap/CmdLine.h"

#include "pso_registration/metric.hpp"
#include "pso_registration/particle.hpp"
#include "pso_registration/swarm.hpp"
#include "pso_registration/utilities.hpp"

typedef pcl::PointXYZ PointType;

using pso_registration::Particle;
using pso_registration::Swarm;

std::string matrix_to_string(const Eigen::MatrixXd &matrix) {
  std::stringstream ret;
  for (int i = 0; i < matrix.size(); i++) {
    ret << matrix(i) << " ";
  }
  return ret.str();
}

int main(int argc, char **argv) {
  std::string source_file_name;
  std::string target_file_name;
  std::string ground_truth_name;
  std::string metric;
  bool display = false;
  float source_filter_size, target_filter_size, random_perc;
  Eigen::Matrix4d initial_transformation = Eigen::Matrix4d::Identity();

  int num_part = 0, num_gen = 0;

  try {
    TCLAP::CmdLine cmd("PSO Point Clouds Registration", ' ', "1.0");
    TCLAP::UnlabeledValueArg<std::string> source_file_name_arg("source_file_name", "The path of the source point cloud",
                                                               true, "source_cloud.pcd", "string", cmd);
    TCLAP::UnlabeledValueArg<std::string> target_file_name_arg(
        "target_file_namedd", "The path of the target point cloud", true, "target_cloud.pcd", "string", cmd);

    TCLAP::UnlabeledMultiArg<double> transformation_arg("trans",
                                                        "The elements of the "
                                                        "4x4 transformation "
                                                        "matrix (without "
                                                        "the  0 0 0 1 "
                                                        "line) to apply to "
                                                        "the source cloud "
                                                        "[m1..12]",
                                                        false,
                                                        "|m1 m2 m3 m4|\n|m5 "
                                                        "m6 m7 m8|\n|m9 m10 "
                                                        "m11 m12|\n|0 0 0 1|",
                                                        cmd);
    TCLAP::ValueArg<int> num_part_arg("p", "num_part", "The number of particles of the swarm", false, 50, "int", cmd);
    TCLAP::ValueArg<int> num_gen_arg("e", "num_it", "The number of iterations (generations) of the algorithm", false,
                                     1000, "int", cmd);
    TCLAP::ValueArg<std::string> metric_arg("m", "metric",
                                            "The metric to use. One of: l1, l2, robust_l2, normalized_robust_l2", false,
                                            "l1", "string", cmd);
    TCLAP::ValueArg<float> source_filter_arg(
        "s", "source_filter_size", "The leaf size of the voxel filter of the source cloud", false, 0, "float", cmd);
    TCLAP::ValueArg<float> target_filter_arg(
        "t", "target_filter_size", "The leaf size of the voxel filter of the target cloud", false, 0, "float", cmd);
    TCLAP::ValueArg<float> perc_filter_arg(
        "r", "random_perc",
        "The percentage of points (of the source and target cloud) to use for the alignment (randomly sampled).", false,
        100, "float", cmd);
    TCLAP::SwitchArg verbose_arg("v", "verbose", "Verbosity", cmd, false);
    TCLAP::SwitchArg display_arg("d", "display", "Display registration progress", cmd, false);
    cmd.parse(argc, argv);

    source_file_name = source_file_name_arg.getValue();
    target_file_name = target_file_name_arg.getValue();
    num_part = num_part_arg.getValue();
    num_gen = num_gen_arg.getValue();
    metric = metric_arg.getValue();
    source_filter_size = source_filter_arg.getValue();
    target_filter_size = target_filter_arg.getValue();
    display = display_arg.getValue();
    random_perc = perc_filter_arg.getValue();
    if (verbose_arg.getValue()) {
      spdlog::set_level(spdlog::level::debug);
    }
    if (transformation_arg.isSet()) {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
          initial_transformation(i, j) = transformation_arg.getValue()[i * 4 + j];
        }
      }
    }
    spdlog::debug("Using {} as initial transformation", matrix_to_string(initial_transformation));
  } catch (TCLAP::ArgException &e) {
    spdlog::critical("error: {0} for arg {1}", e.error(), e.argId());
    exit(EXIT_FAILURE);
  }

  spdlog::debug("Loading source point cloud from {}", source_file_name);

  pcl::PointCloud<PointType>::Ptr source_cloud(new pcl::PointCloud<PointType>);
  if (pcl::io::loadPCDFile(source_file_name, *source_cloud) == -1) {
    spdlog::critical("Could not load source cloud {}, closing...", source_file_name);

    exit(EXIT_FAILURE);
  }

  spdlog::debug("Loading target point cloud from {}", target_file_name);
  pcl::PointCloud<PointType>::Ptr target_cloud(new pcl::PointCloud<PointType>);
  if (pcl::io::loadPCDFile<PointType>(target_file_name, *target_cloud) == -1) {
    spdlog::critical("Could not load source cloud {}, closing...", target_file_name);
    exit(EXIT_FAILURE);
  }

  std::vector<int> nan_points;
  pcl::removeNaNFromPointCloud(*source_cloud, *source_cloud, nan_points);
  pcl::removeNaNFromPointCloud(*target_cloud, *target_cloud, nan_points);

  pcl::PointCloud<PointType>::Ptr moved_filtered_source_cloud(new pcl::PointCloud<PointType>);
  pcl::copyPointCloud(*source_cloud, *moved_filtered_source_cloud);
  pcl::transformPointCloud(*moved_filtered_source_cloud, *moved_filtered_source_cloud, initial_transformation);

  double initial_error = point_cloud_registration_benchmark::calculate_error(source_cloud, moved_filtered_source_cloud);

  pcl::VoxelGrid<pcl::PointXYZ> voxel_filter;
  if (source_filter_size > 0) {
    spdlog::debug("Downsampling source cloud with a leaf size of {}", source_filter_size);
    voxel_filter.setInputCloud(moved_filtered_source_cloud);
    voxel_filter.setLeafSize(source_filter_size, source_filter_size, source_filter_size);
    voxel_filter.filter(*moved_filtered_source_cloud);
  }
  if (target_filter_size > 0) {
    spdlog::debug("Downsampling target cloud with a leaf size of {}", target_filter_size);
    voxel_filter.setInputCloud(target_cloud);
    voxel_filter.setLeafSize(target_filter_size, target_filter_size, target_filter_size);
    voxel_filter.filter(*target_cloud);
  }

  if (random_perc != 100) {
    spdlog::debug("Randomly using {}% of points", random_perc);
    pcl::RandomSample<PointType> sample(true);
    sample.setInputCloud(moved_filtered_source_cloud);
    sample.setSample(moved_filtered_source_cloud->size() * (random_perc / 100));
    sample.filter(*moved_filtered_source_cloud);

    sample.setInputCloud(target_cloud);
    sample.setSample(target_cloud->size() * random_perc / 100);
    sample.filter(*target_cloud);
  }

  double (*score_function)(pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr,
                           pcl::KdTree<PointType>::Ptr);
  if (metric == "l2") {
    spdlog::debug("Using l2 distance");
    score_function = &pso_registration::l2_distance;
  } else if (metric == "robust_l2") {
    spdlog::debug("Using robust l2 distance");
    score_function = &pso_registration::robust_l2_distance;
  } else if (metric == "robust_normalized_l2") {
    spdlog::debug("Using robust normalized l2 distance");
    score_function = &pso_registration::robust_normalized_l2_distance;
  } else {
    spdlog::debug("Using l1 distance");
    score_function = &pso_registration::l1_distance;
  }

  pcl::visualization::PCLVisualizer::Ptr viewer;
  if (display) {
    viewer.reset(new pcl::visualization::PCLVisualizer("PSO Registration"));
    viewer->setBackgroundColor(255, 255, 255);

    pcl::visualization::PointCloudColorHandlerCustom<PointType> target_color(target_cloud, 0, 255, 0);
    viewer->addPointCloud<PointType>(target_cloud, target_color, "target");

    pcl::visualization::PointCloudColorHandlerCustom<PointType> moved_source_color(moved_filtered_source_cloud, 0, 0,
                                                                                   255);
    viewer->addPointCloud<PointType>(moved_filtered_source_cloud, moved_source_color, "source");

    pcl::visualization::PointCloudColorHandlerCustom<PointType> gt_source_color(source_cloud, 255, 0, 0);
    viewer->addPointCloud<PointType>(source_cloud, gt_source_color, "gt_source");
  }
  Swarm swarm;

  pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr target_tree (new pcl::KdTreeFLANN<pcl::PointXYZ>);
  target_tree->setEpsilon(3.6);
  target_tree->setInputCloud(target_cloud);
  for (int i = 0; i < num_part; i++) {
    swarm.add_particle(Particle(moved_filtered_source_cloud, target_cloud, target_tree, i, score_function));
  }
  Particle best;
  swarm.init();
  spdlog::debug("{0:<10} {1:>10}", "Generation", "Best (id: tx ty tz angle azimuth elevation --> score)");
  for (int i = 0; i < num_gen; i++) {
    std::ostringstream os;
    swarm.evolve();
    best = swarm.getBest();
    os << best;
    spdlog::debug("{0:<10} {1:>10}", i, os.str());
    if (display) {
      viewer->updatePointCloudPose("source", best.getTransformation().cast<float>());
      viewer->spinOnce(1);
    }
  }

  auto estimated_transformation = swarm.getBest().getTransformation();
  pcl::transformPointCloud(*source_cloud, *moved_filtered_source_cloud, estimated_transformation);
  double error = point_cloud_registration_benchmark::calculate_error(source_cloud, moved_filtered_source_cloud);

  std::cout.precision(15);
  std::cout << initial_error << ", " << error << ", ";
  Eigen::Matrix4d matrix = estimated_transformation.matrix();
  for (int i = 0; i < matrix.size(); i++) {
    std::cout << matrix(i) << " ";
  }
  return 0;
}
