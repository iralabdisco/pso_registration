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
USAGE: 

   ./pso_registration  [-d] [-v] [-r <float>] [-t <float>] [-s <float>] [-m
                       <string>] [-e <int>] [-p <int>] [--] [--version]
                       [-h] <string> <string> <|m1 m2 m3 m4|

                       |m5 m6 m7 m8|

                       |m9 m10 m11 m12|

                       |0 0 0 1|> ...


Where: 

   -d,  --display
     Display registration progress

   -v,  --verbose
     Verbosity

   -r <float>,  --random_perc <float>
     The percentage of points (of the source and target cloud) to use for
     the alignment (randomly sampled).

   -t <float>,  --target_filter_size <float>
     The leaf size of the voxel filter of the target cloud

   -s <float>,  --source_filter_size <float>
     The leaf size of the voxel filter of the source cloud

   -m <string>,  --metric <string>
     The metric to use. One of: l1, l2, robust_l2, normalized_robust_l2

   -e <int>,  --num_it <int>
     The number of iterations (generations) of the algorithm

   -p <int>,  --num_part <int>
     The number of particles of the swarm

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <string>
     (required)  The path of the source point cloud

   <string>
     (required)  The path of the target point cloud

   <|m1 m2 m3 m4|

      |m5 m6 m7 m8|

      |m9 m10 m11 m12|

      |0 0 0 1|> 
     The elements of the 4x4 transformation matrix (without the  0 0 0 1
     line) to apply to the source cloud before the registration [m1..12]
