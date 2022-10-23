# CMAP-DeepBKZ

CMAP-DeepBKZ solver \[[Tat+21b](#ref.cmapdeepbkz)\] that is a configuration for parallelizing lattice-based reduction based on CMAP-LAP \[[Tat+21a](#ref.cmaplap)\]. CMAP-LAP is a C++ & MPI parallel framework for the lattice problems \[[Tat+21a](#ref.cmaplap)\].

These frameworks and solvers are distributed under the LGPL version 3 or later, in accordance with the license of the UG framework. Commercial licenses are available through <licenses@zib.de>.
In addition, some programs independent of the UG framework are distributed under the MIT License.

<br>

# Table of contents
- [CMAP-DeepBKZ](#cmap-deepbkz)
- [Table of contents](#table-of-contents)
- [Setup](#setup)
- [Compiling](#compiling)
- [Usage](#usage)
- [Bibliography](#bibliography)

<br>

# Setup

We provide Dockerfile for building virtual environments and making binaries. If you do not use docker, you set up third-party libraries using `setup.sh` or install them manually.

## Setup on Docker

If you use docker to set up, you just have to run the following building command. The docker image will be created from ubuntu:20.04. Then you run the docker image.

```
docker build -t cmapdeepbkz .
```


## Setup Manually

If you do not use docker, you have to install third-party libraries as follows.

- CMake (version >= 3.18)
- NTL
- Eigen
- Boost (version == 1.75)

**Note** We have checked our framework and sovlers with CMake(v3.22.2), NTL(v11.5.1), Eigen(v3.4.0) and Boost(v1.75).
The commands to install these libraries are as follows.
```bash
# install cmake (version 3.22.3)
wget https://github.com/Kitware/CMake/releases/download/v3.22.3/cmake-3.22.3.tar.gz
tar -xf cmake-3.22.3.tar.gz
cd cmake-3.22.3
./bootstrap [--prefix=PREFIX] [--parallel=PARALLEL]
make & make install

# install NTL (version 11.5.1)
wget https://libntl.org/ntl-11.5.1.tar.gz
tar -xf ntl-11.5.1.tar.gz
cd ntl-11.5.1/src
./configure [PREFIX=PREFIX]
make & make install

# install Eigen (version 3.4.0)
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xf eigen-3.4.0.tar.gz
cd eigen-3.4.0
mkdir build
cd build
cmake .. [-DCMAKE_INSTALL_PREFIX=PREFIX]
make install


# install boost (version 1.75 exact)
wget https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz
tar -xf boost_1_75_0.tar.gz
```

In the case that you install these libraries in your local, it is recommended to create a cmaplap/usr directory and specify cmaplap/usr in the PREFIX above.

<br>

# Compiling

The binaries will be created in `DeepBKZ/bin`.

```bash
mkdir build
cd build
cmake .. (options)
make
```

When you want to compile using MPI, set `CMAKE_CXX_COMPILER` to mpi compiler, e.g. `cmake .. -DCMAKE_CXX_COMPILER=mpicxx`.

## Created binaries

- `bin/fcmapdeepbkz`: shared memory version of CMAP-DeepBKZ
- `bin/paracmapdeepbkz`: distributed memory version of CMAP-DeepBKZ

## CMake Options

- `-DBOOST_DIR`: (required) directory of boost 1.75
- `-DCMAKE_BUILD_TYPE=Debug`: compile with debug mode
- `-DSHARED_MEMORY_ONLY=ON`: compile only shared-memory version
- `-DCMAKE_CXX_COMPILER=XXX`: use XXX mpi compiler (e.g. mpicxx)

**Note**
If you do not use mpi compiler (e.g. -DCMAKE_CXX_COMPILER=gcc), the build will fail under `SHARED_MEMORY_ONLY` option is OFF, so you have to use the option `-DSHARED_MEMORY_ONLY=ON`.


### Examples

- compile only shared memory version: `cmake .. -DBOOST_DIR=/xxx/boost_1_75_0 -DSHARED_MEMORY_ONLY=ON`
- compile both shared and distributed memory version: `cmake .. -DBOOST_DIR=/xxx/boost_1_75_0 -DCMAKE_CXX_COMPILER=mpicxx`
- compile with debug mode (add -g, and remove -O3): `cmake .. -DBOOST_DIR=/xxx/boost_1_75_0 -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Debug`


## Check

```sh
python test.py
```

<br>

# Usage

This solver is the parallel solver for lattice basis reduction.
All workers execute lattice basis reduction, and a supervisor shares a part of a lattice basis.

## Shared memory version

`./bin/fcmapdeepbkz settingfile matrixfile -sth yy` (e.g. `./bin/fcmapdeepbkz settings/default.set storage/sample_mats/dim80.txt -sth 3`)

matrixfile is the basis file in the same format as the SVP Challenge instance. 

**options**

- sth [Int] : the number of solver threads used

## Distributed memory version

`mpirun -np yy ./bin/paracmapdeepbkz settingfile matrixfile` (e.g. `mpirun -np 3 ./bin/paracmapdeepbkz settings/default.set storage/sample_mats/dim80.txt`)

**options**

- np [Int] : the number of solver process + 1


## Sequential version (test for DeepBKZ algorithms)

`./bin/seqcmaplap -i (instance file) -a exdeepbkz -b (blocksize)` (e.g. `./bin/seqcmaplap -i ./storage/sample_mats/dim80.txt -a exdeepbkz -b 30`)

<br>

# Bibliography

<a id="ref.cmaplap"></a>
\[Tat+21a\] Nariaki Tateiwa, Yuji Shinano, Keiichiro Yamamura, Akihiro Yoshida, Shizuo Kaji, Masaya Yasuda, and Katsuki Fujisawa. “CMAP-LAP: Configurable massively parallel solver for lattice problems”. In: 2021 IEEE 28th International Conference on High Performance Computing, Data, and Analytics (HiPC). IEEE. 2021, pp. 42–52.

<a id="ref.cmapdeepbkz"></a>
\[Tat+21b\] Nariaki Tateiwa, Yuji Shinano, Masaya Yasuda, Shizuo Kaji, Keiichiro Yamamura, and Katsuki Fujisawa. Massively parallel sharing lattice basis reduction. eng. Tech. rep. 21-38. Takustr. 7, 14195 Berlin: ZIB, 2021.
