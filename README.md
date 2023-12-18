# README #

#### BLAS and LAPACK ####
```bash
$ sudo apt-get install libblas-dev liblapack-dev libopenblas-base libopenblas-dev
```

#### OPEN MP #####
```bash
$ sudo apt-get install libomp-dev
```


#### MPI ####

[Download](https://www.open-mpi.org/)

```bash
$ mkdir $OPENMPIBUILDDIR
$ cd $MPIDIR
$ ./configure --prefix=$OPENMPIBUILDDIR
$ make (all)
$ make install
```

execute programm with
mpirun -np NumProcessors programm



#### suitesparse: ####

needs blas and lapack and openblas

##### option 1 : install
```bash
sudo apt-get install libsuitesparse-dev
```

##### option 2: Include in source/externalLib (recommended):
unzip to SuiteSparse
cd SuiteSparse
(if IPOPT is used: delete folder SuiteSparseDirectory/metis-5)
make

note:
in /usr/include/suitesparse, the cmake-file in cmake/findSuiteSparse helps to find the library


#### eigen ####

##### option 1 : install
sudo apt-get install libeigen3-dev

##### option 2: Include in source/externalLib (recommended)


#### boost#: ####

##### option 1: Install from archive manager: 
UNIX: sudo apt-get install libboost-all-dev


##### option 2: build from source

download from boost.org/users/download
extract to $BOOSTDIR, e.g. $BOOSTDIR = /home/user/prog/boost_1_75_0

UNIX:
Open Terminal and navigate to $BOOSTDIR
Run ./bootstrap.sh

To install, the options are listed by using
./b2 --help

Consruct build and install directory:
mkdir boost-build
mkdir boost-install

Run
./b2 --build-dir="$BOOSTDIR/boost-build" --prefix="$BOOSTDIR/boost-install" --build-type=minimal --layout=system toolset=gcc install

(layout=versioned produces another directory structure and libraries are name with version numbers. It seems that cmake cannot find the versioned libraries?)
(build-type=complete installs all libraries, but minimal should be sufficient. Note that minimal is platform-dependent.)


Run ./b2 install --prefix=PREFIX where PREFIX is a directory where you want Boost.Build to be installed.



WINDOWS: 
open shell-script bootstrap

On the command line, go to the root of the unpacked tree.
Run ./bootstrap.sh
Run ./b2 install --prefix=PREFIX where PREFIX is a directory where you want Boost.Build to be installed.


#### VTK: ####
[Download] (http://www.vtk.org/Wiki/VTK/Configure_and_Build)

Requirement: OpenGL
sudo apt-get install freeglut3-dev

```bash
unzip source code into directory VTK
cd VTK
mkdir build
cd build
ccmake ../
configure (c)
set CMAKE_Build_Type to Release
configure (c)
generate (g)
make
```

##### to run VTK on server:
install xvfb: sudo apt-get install xvfb
run programm with nohup xvfb-run -a --server-args "-screen 0 1920x1080x24" ./programm &>/dev/null& 2>&1 

##### info:
                Legacy (.vtk)                        New (XML)
structured      STRUCTURED_POINTS                    .vti (vtkImageData)
                STRUCTURED_GRID                      .vts
unstructured    RECTLINEAR_GRID                      .vtr
                POLYDATA                             .vtp
                UNSTRUCTURED_GRID                    .vtu


#### IPOPT  ####

requires fortran compiler, e.g. 
sudo apt-get install gfortran

0) unzip the IPOPT-package

1) Install the ThirdParty-packages :
$ cd $IPOPTDIR/ThirdParty/Blas
$ ./get.Blas
$ cd ../Lapack
$ ./get.Lapack
$ cd ../ASL
$ ./get.ASL 
$ cd ../Mumps
$ ./get.Mumps
$ cd ../Metis
$ ./get.Metis

To compile the HSL code as part of IPOPT: 
unpack the archive
then move and rename the resulting directory to $IPOPTDIR/ThirdParty/HSL/coinhsl

2) 
$ cd $IPOPTDIR
$ mkdir build
$ cd build
$ ../configure
$ make 
$ make test
$ make install


#### qt ####
https://www.qt.io/download-qt-installer

or
sudo apt-get install qt5-default

#### gcc #####

check version with

```bash
$ gcc --version
```
to intall version 10:
sudo apt-get install gcc-10
sudo apt-get install g++-10

sudo update-alternatives --remove-all gcc

sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 --slave /usr/bin/g++ g++ /usr/bin/g++-10

sudo update-alternatives --config gcc

#### Cmake ####

check version with

```bash
$ cmake --version
```

for C++17 at least version 3.8 is required
[Download](https://cmake.org/download/)

the installation requires openssl

```bash
sudo apt-get install libssl-dev
```

to guarantee that ccmake is istalled:
sudo apt-get install libncurses5-dev libncursesw5-dev

```bash
$ cd $CMAKEDIR
$ bash ./configure   or  ./bootstrap
(to install cmake-gui: ./configure --qt-gui )
$ make
$ make install
```


#### CMAKE-Setup: ####

if external library is installed in source/externalLibs, it should automatically work

else use

* EIGEN_INCLUDE_DIR:PATH=$EIGENDIR
* SUITESPARSE_LIBRARY_DIR:PATH=$SuiteSparseDir/lib
* SUITESPARSE_CONFIG_LIB:FILEPATH=$SuiteSparseDir/lib/libsuitesparseconfig.so
* CHOLMOD_INCLUDE_DIR:PATH=$SuiteSparseDir/include
* IPOPT_INCLUDE_DIR:PATH=$IPOPTDIR/build/include/coin
* IPOPT_LIBRARY:FILEPATH=$IPOPTDIR/build/lib/libipopt.so
* VTK_DIR:PATH=$VTKDIR/build
* MPI_INCLUDE_DIR=$OPENMPIBUILDDIR/include

if boost is installed from source:
* PESOPT_BOOST_DIR:PATH=$BOOSTDIR/boost-install

#### Build documentation with DOXYGEN #### 

1) install:
sudo apt-get install doxygen 
2) generate:
cd source
doxygen doxygen.conf


#### use debugger gdb #### 

```bash
$ gdb ./programm 
$ run arguments
$ bt
```


### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact


### Usage ###
The main executable is projects/dkt/paramSurfaces/ShellDeformIsometrySemiNonLinAdaptive

Parameter files are in /parameters/dkt (Plate, Half-Cylinder, Saddle)
