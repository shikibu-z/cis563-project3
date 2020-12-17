# CIS 563 Project 3 - Material Point Method
**This project requries specific a C++ environment from [CISPBA](https://github.com/cffjiang/cispba)**
```
sudo apt install libsuitesparse-dev libxcursor-dev libxinerama-dev libopenblas-dev
sudo apt install make cmake g++ libeigen3-dev gfortran libmetis-dev libopenvdb-dev
sudo apt install libboost-all-dev libilmbase-dev libopenexr-dev libtbb2 libtbb-dev
sudo apt install libz-dev clang-format-6.0 xorg-dev libglu1-mesa-dev
```
**Make sure to check the above commands for required software packages.**
```
cis-563-project3
├── cis563-project3.patch
├── CMake
│   ├── DownloadProject.cmake
│   ├── DownloadProject.CMakeLists.cmake.in
│   ├── DownloadProject.LICENSE
│   ├── FindCapnProto.cmake
│   ├── FindEigen3.cmake
│   ├── FindHalf.cmake
│   ├── FindMetis.cmake
│   ├── FindOpenVDB.cmake
│   ├── FindTBB.cmake
│   ├── FindZeroMQ.cmake
│   ├── FindZeroMQPP.cmake
│   └── UseLATEX.cmake
├── CMakeLists.txt
├── Deps
│   ├── CMakeLists.txt
│   └── partio.patch
├── Makefile
├── Projects
│   ├── CMakeLists.txt
│   └── material_point
│       ├── CMakeLists.txt
│       ├── data
│       │   └── cube.obj
│       ├── main.cpp
│       ├── material_point.h
│       ├── mesh_query0.1
│       │   ├── bounding_box.h
│       │   ├── bounding_box_tree.cpp
│       │   ├── bounding_box_tree.h
│       │   ├── Makefile
│       │   ├── mesh_query.cpp
│       │   ├── mesh_query.h
│       │   ├── predicates.cpp
│       │   ├── predicates.h
│       │   ├── README
│       │   ├── util.h
│       │   └── vec.h
│       └── SimulationDriver.h
├── README.md
├── report.pdf
└── Scripts
    ├── Makefile.in
    └── valgrind.supp

```
+ Please make sure you are at C++ 11 or above.
+ This project requries mainly `Eigen` and `Partio` library.
+ Detailed code layout and explainations could be found at `report.pdf`.
+ Type `make` at the top level directory will compile. 
+ Type `./material_point` at `/Projects/material_point/` will run the simulation, `.bgeo` outputs could be found at `Projects/material_point/output/`.
+ Based on different operating systems and setups, some source files might needed to be changed. Please refer to the `cis563-project3.patch` file for more information. Patching this file will enable correct compilation from `macOS 10.15.7` to `Ubuntu 20.04`. **If you are already setup correctly following the above instruction on a Linux machine, you don't need to worry about this file. This directory is already patched for you**.
