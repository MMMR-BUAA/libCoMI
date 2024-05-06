**CoMI** is the C++ library of center of inertia and mass.   
**CoMI** need Eigen and pinocchio.  
**CoMI** is compiled from MMMR laboratory of Beihang university.  

## Table of contents 

  - [CoMI files architecture](#CoMI-files-architecture)
  - [build and install](#build-and-install)
  - [test](#test)
  - [Usage](#Usage)
  - [Documentation](#documentation)
  - [Acknowledgments](#acknowledgments)

## files architecture
—CoMI: the lib of center of mass and inertia.  
—doc: the theory and introduction of CoMI.  
—test: the test progrom of CoMI.  

## build and install
Make sure that Eigen and pinocchio are installed.
In the CoMI/build folder, build CoMI:
```bibtex
cmake ..
```
and install CoMI:
```bibtex
sudo make install
```

## test
The test model is a two-link mechanism.
In the test/build folder, build test:
```bibtex
cmake ..
```
and install CoMI:
```bibtex
sudo make install
```
running test programmer:
```bibtex
./pinocchio_test
```

## Usage
**CoMI** follows the template concept.
You can add the Makefile and include the head in your project to use **CoMI**.
In the make file, add:
```bibtex
find_package(CoMI REQUIRED)
include_directories(${CoMI_INCLUDE_DIRS})
```
In the source file, add:
```bibtex
#include "CoMI/CoI/centerofinertia.hpp"
```
And you can use the template class or function in your project.

## Acknowledgments
The development of **CoMI** is actively supported by the National Natural Sceince Foundation of China (Grant No.).




