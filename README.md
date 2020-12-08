# CIS 563 Project 3 - Material Point Method
**This project requries specific CMake environment from [CISPBA](https://github.com/cffjiang/cispba).**
```
cis-563-project3
├── CMake
│   ├── DownloadProject.CMakeLists.cmake.in
│   ├── DownloadProject.LICENSE
│   ├── DownloadProject.cmake
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
│       ├── SimulationDriver.h
│       ├── data
│       │   └── cube.obj
│       ├── main.cpp
│       ├── material_point.h
│       └── mesh_query0.1
│           ├── Makefile
│           ├── README
│           ├── bounding_box.h
│           ├── bounding_box_tree.cpp
│           ├── bounding_box_tree.h
│           ├── mesh_query.cpp
│           ├── mesh_query.h
│           ├── predicates.cpp
│           ├── predicates.h
│           ├── util.h
│           └── vec.h
├── README.md
├── Scripts
│   ├── Makefile.in
│   └── valgrind.supp
├── cube.mov
└── report.pdf
```
+ Please make sure you are at C++ 11 or above.
+ This project requrie `Eigen` and `Partio` library.
+ Detailed code layout and explainations could be found at `/report.pdf`.
+ Type `make` at the top level directory will compile. 
+ Type `./material_point` at `/Project/material_point/` will run the simulation.