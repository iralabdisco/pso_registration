# pso_registration
This project aims at developing a technique for global point clouds registration using Particle Swarm Optimization. The software is meant to be used in conjunction with any other fine registration technique to align two generic point clouds without specifying any initial guess on the transformation. The two point clouds do not need to be already roughly aligned. The suggested usage is:
 - Obtain an initial guess using pso_registration
 - If the alignment obtained with pso_registration is not accurate enough, refine with ICP, G-ICP, ecc...

### Compiling ###
Dependencies:
 - [Boost](https://www.boost.org/)
 - [TCLap](http://tclap.sourceforge.net/)
 - [PCL](http://pointclouds.org/)
 - [spdlog](https://github.com/gabime/spdlog)
 - A c++14 compiler
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

### Usage ###
   `./pso_registration  [-v] [-g <string>] [-e <int>] [-p <int>] [-t
                       <float>] [-s <float>] [--] [--version] [-h] <source_cloud>
                       <target_cloud>`
                       
Where: 

   `-v,  --verbose` 
      Verbosity. Is set, display the current status of the alignment.

   `-g <string>,  --ground_truth <string>`
   The path of the ground truth for the source cloud, if available. It is used only to calculate the error of the obtained  alignment w.r.t. the ground truth.

   `-e <int>,  --num_it <int>`
     The number of iterations (generations) of particle swarm optimization algorithm

   `-p <int>,  --num_part <int>`
     The number of particles of the swarm

   `-t <float>,  --target_filter_size <float>`
     The leaf size of the voxel filter of the target cloud. The target point cloud can be downsampled using a [VoxelGrid Filter](http://pointclouds.org/documentation/tutorials/voxel_grid.php). Use a filter size of 0 if downsampling is not required.

   `-s <float>,  --source_filter_size <float>`
     The leaf size of the voxel filter of the source cloud. The source point cloud can be downsampled using a [VoxelGrid Filter](http://pointclouds.org/documentation/tutorials/voxel_grid.php). Use a filter size of 0 if downsampling is not required.

   `--,  --ignore_rest`
     Ignores the rest of the labeled arguments following this flag.

   `--version`
     Displays version information and exits.

   `-h,  --help`
     Displays usage information and exits.

   `<source_cloud>
     (required)`  The path of the source point cloud

   `<target_cloud>
     (required)`  The path of the target point cloud

The default values for the number of particles (50) and the number of generations (1000) should work for most registrations. The number of generations can be lowered to improve the speed of the algorithm.
