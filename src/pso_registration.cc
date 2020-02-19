#include <stdlib.h>

#include <Eigen/Core>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "boost/make_shared.hpp"
#include "pcl/common/transforms.h"
#include "pcl/conversions.h"
#include "pcl/filters/filter.h"
#include "pcl/filters/voxel_grid.h"
#include "pcl/io/pcd_io.h"
#include "pcl/point_types.h"
#include "pcl/visualization/pcl_visualizer.h"
#include "pso_registration/particle.hpp"
#include "pso_registration/swarm.hpp"
#include "pso_registration/utilities.hpp"
#include "pso_registration/metric.hpp"
#include "spdlog/spdlog.h"
#include "tclap/CmdLine.h"

typedef pcl::PointXYZ PointType;

using pso_registration::calculateMSE;
using pso_registration::Particle;
using pso_registration::Swarm;

int main(int argc, char **argv) {
  std::string source_file_name;
  std::string target_file_name;
  std::string ground_truth_name;
  std::string metric;
  bool verbose = false;
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
    TCLAP::ValueArg<std::string> metric_arg("m", "metric", "The metric to use. One of: l1, l2, robust_l2, normalized_robust_l2", false,
                                     "l1", "string", cmd);
    TCLAP::SwitchArg verbose_arg("v", "verbose", "Verbosity", cmd, false);
    cmd.parse(argc, argv);

    source_file_name = source_file_name_arg.getValue();
    target_file_name = target_file_name_arg.getValue();
    num_part = num_part_arg.getValue();
    num_gen = num_gen_arg.getValue();
    metric = metric_arg.getValue();
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
  } catch (TCLAP::ArgException &e) {
    spdlog::critical("error: {0} for arg {1}", e.error(), e.argId());
    exit(EXIT_FAILURE);
  }

  spdlog::debug("Loading source point cloud from {}", source_file_name);

  pcl::PointCloud<PointType>::Ptr source_cloud = boost::make_shared<pcl::PointCloud<PointType>>();
  if (pcl::io::loadPCDFile(source_file_name, *source_cloud) == -1) {
    if (verbose) {
      spdlog::critical("Could not load source cloud {}, closing...", source_file_name);
    }
    exit(EXIT_FAILURE);
  }

  spdlog::debug("Loading target point cloud from {}", target_file_name);
  pcl::PointCloud<PointType>::Ptr target_cloud = boost::make_shared<pcl::PointCloud<PointType>>();
  if (pcl::io::loadPCDFile<PointType>(target_file_name, *target_cloud) == -1) {
    spdlog::critical("Could not load source cloud {}, closing...", target_file_name);
    exit(EXIT_FAILURE);
  }

  std::vector<int> nan_points;
  pcl::removeNaNFromPointCloud(*source_cloud, *source_cloud, nan_points);
  pcl::removeNaNFromPointCloud(*target_cloud, *target_cloud, nan_points);
  
  pcl::PointCloud<PointType>::Ptr  moved_source;
  moved_source.reset(new pcl::PointCloud<PointType>);
  pcl::transformPointCloud(*source_cloud, *moved_source, initial_transformation);
  double initial_error = point_cloud_registration_benchmark::calculate_error(source_cloud, moved_source);

  double (*score_function)(pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr);
  if(metric == "l2"){
    spdlog::debug("Using l2 distance");
    score_function = &pso_registration::l2_distance;
  }else if(metric =="robust_l2"){
    spdlog::debug("Using robust l2 distance");
    score_function = &pso_registration::robust_l2_distance;
  }else if(metric == "robust_normalized_l2"){
    spdlog::debug("Using robust normalized l2 distance");
    score_function = &pso_registration::robust_normalized_l2_distance;
  }else{
    spdlog::debug("Using l1 distance");
    score_function = &pso_registration::l1_distance;
  }

  Swarm swarm;
  for (int i = 0; i < num_part; i++) {
    swarm.add_particle(Particle(moved_source, target_cloud, i, score_function));
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
  }

  auto estimated_transformation = swarm.getBest().getTransformation();
  pcl::transformPointCloud(*moved_source, *moved_source, estimated_transformation);
  double error = point_cloud_registration_benchmark::calculate_error(source_cloud, moved_source);
  
  std::cout.precision(15);
  std::cout << initial_error << ", " << error << ", ";
  Eigen::Matrix4d matrix = estimated_transformation.matrix();
  for (int i = 0; i < matrix.size(); i++) {
    std::cout << matrix(i) << " ";
  }
  return 0;
}
