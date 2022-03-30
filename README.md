# CMAP-LAP
CMAP-LAP is a C++ & MPI parallel framework for the lattice problems \[[Tat+21a](#ref.cmaplap)\].
This framework is based on the "UG framework" \[[Ug](#ref.Ug)\].

The CMAP-LAP framework has a supervisor-worker style, and it is designed to run in a massively parallel environment with more than a hundred thousand cores.
Load Coordinator, a supervisor process, distributes tasks, collects data from workers, and distributes data to workers.
Checkpoint and restart mechanisms are also implemented for long-term execution.

The lattice algorithms and data sharing strategies implemented in the CMAP-LAP worker process can be customized.

This repository includes a CMAP-TEST solver \[[Tat+20](#mapsvp), [Tat+21a](#ref.cmaplap)\] that is the test configurations of CMAP-LAP and a CMAP-DeepBKZ solver \[[Tat+21b](#ref.cmapdeepbkz)\] that is a configuration for parallelizing lattice-based reduction.

Framework and solvers provide a distributed memory parallel program with MPI and shared memory one with communication between threads.
The shared memory parallel programs are designed to mimic the behavior of the distributed memory one for debugging it.

These frameworks and solvers are distributed under the LGPL version 3 or later, in accordance with the license of the UG framework. Commercial licenses are available through <licenses@zib.de>.
In addition, some programs independent of the UG framework are distributed under the MIT License.

<br>

# Table of contents
- [CMAP-LAP](#cmap-lap)
- [Table of contents](#table-of-contents)
- [Setup](#setup)
  - [Setup on Docker](#setup-on-docker)
  - [Setup Manually](#setup-manually)
- [Compiling](#compiling)
  - [Created binaries](#created-binaries)
  - [CMake Options](#cmake-options)
  - [Check](#check)
- [Usage](#usage)
  - [CMAP-TEST](#cmap-test)
    - [Shared memory version](#shared-memory-version)
    - [Distributed memory version](#distributed-memory-version)
  - [CMAP-DeepBKZ](#cmap-deepbkz)
  - [Common Options (CMAP-TEST and CMAP-DeepBKZ)](#common-options-cmap-test-and-cmap-deepbkz)
  - [Sequential version (test for lattice algorithms)](#sequential-version-test-for-lattice-algorithms)
- [Contributors](#contributors)
- [Bibliography](#bibliography)

<br>

# Setup

We provide Dockerfile for building virtual environments and making binaries. If you do not use docker, you set up third-party libraries using `setup.sh` or install them manually.

## Setup on Docker

If you use docker to set up, you just have to run the following building command. The docker image will be created from ubuntu:20.04.

```
docker build -t cmaplap .
```

Then you run the docker image.


## Setup Manually

If you do not use docker, you have to install third-party libraries as follows.

- CMake (version >= 3.18)
- NTL
- Eigen

**Note** We have checked our framework and sovlers with CMake(v3.22.2), NTL(v11.5.1), Eigen(v3.4.0).
The commands to install these libraries are as follows.
```bash
# install cmake (version 3.22.3)
wget https://github.com/Kitware/CMake/releases/download/v3.22.3/cmake-3.22.3.tar.gz
tar -xf cmake-3.22.3.tar.gz
cd cmake-3.22.3
./bootstrap [--prefix=PREFIX] [--parallel=PARALLEL]
make
make install

# install NTL (version 11.5.1)
wget https://libntl.org/ntl-11.5.1.tar.gz
tar -xf ntl-11.5.1.tar.gz
cd ntl-11.5.1/src
./configure [PREFIX=PREFIX]
make
make install

# install Eigen (version 3.4.0)
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xf eigen-3.4.0.tar.gz
cd eigen-3.4.0
mkdir build
cd build
cmake .. [-DCMAKE_INSTALL_PREFIX=PREFIX]
make install
```

<br>

# Compiling

The binaries will be created in `DeepBKZ/bin`.

```bash
mkdir build
cd build
cmake .. (same options)
make
```

When you want to compile a parallel version, set `CMAKE_CXX_COMPILER` to mpi compiler, e.g. `cmake .. -DCMAKE_CXX_COMPILER=mpicxx`.

## Created binaries

- `bin/seqcmaplap` : sequential version (for testing lattice algoritms)
- `bin/fcmaptest`: shared memory version of CMAP-TEST
- `bin/paracmaptest`: distributed memory version of CMAP-TEST
- `bin/fcmapdeepbkz`: shared memory version of CMAP-DeepBKZ
- `bin/paracmapdeepbkz`: distributed memory version of CMAP-DeepBKZ


## CMake Options

- `-DCMAKE_BUILD_TYPE=Debug`: compile with -g option and debug flag for Eigen library
- `-DSHARED_MEMORY_ONLY=ON`: compile only shared-memory version
- `-DCMAKE_CXX_COMPILER=XXX`: use XXX mpi compiler (e.g. mpicxx)

**Note**
If you do not use mpi compiler (e.g. -DCMAKE_CXX_COMPILER=gcc), the build will fail under `SHARED_MEMORY_ONLY` option is OFF, so you have to use the option `-DSHARED_MEMORY_ONLY=ON`.

## Check

```sh
python test.py
```

<br>

# Usage

## CMAP-TEST

This is a test configuration solver of CMAP-LAP.
In this solver, workers can execute DeepBKZ, ENUM, or Sieve algorithms, and a supervisor collects and distributes short lattice vectors among solvers.

### Shared memory version

`./bin/fcmaptest settingfile matrixfile -sth yy`

- matrixfile is the basis file in the same format as the SVP Challenge instance.

**example**

`./bin/fcmaptest settings/default.set storage/sample_mats/dim80.txt -sth 3`

**options**

- `-sth [Int]` : the number of solver threads used

<br>

### Distributed memory version

`mpirun -np yy ./bin/paracmaptest settingfile matrixfile`

**example**

`mpirun -np 3 ./bin/paracmaptest settings/default.set storage/sample_mats/dim80.txt`

**options**

- `-np [Int]` : the number of solver process + 1


## CMAP-DeepBKZ

This solver is the parallel solver for lattice basis reduction.
All workers execute lattice basis reduction, and a supervisor shares a part of a lattice basis.

The commands for the shared memory version and distributed memory one are the same as CMAP-TEST, but binaries are `./bin/fcmapdeepbkz` (for shared memory) and `./bin/paracmapdeepbkz` (for distributed memory).


## Common Options (CMAP-TEST and CMAP-DeepBKZ)

- `-q` : suppress screen messages
- `-w <prefix_warm>` : warm start file prefix ( prefix_warm_nodes.gz and prefix_warm_solution.txt are read )
- `-fsol <solution file>` : specify output solution file
- `-isol <intial solution file>` : specify initial solution file


## Sequential version (test for lattice algorithms)

`./bin/seqcmaplap [options] `

**example**

`./bin/seqcmaplap -i ./storage/sample_mats/dim80.txt`

**options**

```
-i [basis file]
-f [param file]
-a [algorithm]
-b [beta (bkz;default 20)]
-c [csv log file path]
-p [prob (enum; default 1, auto -1)]
-s [seed (randomize; default 0)]
-t [time limit[s]]
-l [lower bound]
-o [output file path]
-q quiet output
-h show help
```

**algorithm**
```
deeplll, deepbkz(default)
enum, subenum
gausssieve
```

<br>

# Contributors

The following people have contributed to this project.

- Nariaki Tateiwa
- Yuji Shinano
- Masaya Yasuda
- Shizuo Kaji
- Keiichiro Yamamura
- Akihiro Yoshida
- Katsuki Fujisawa

<br>

# Bibliography

<a id="ref.mapsvp"></a>
\[Tat+20\] Nariaki Tateiwa, Yuji Shinano, Satoshi Nakamura, Akihiro Yoshida, Shizuo Kaji, Masaya Yasuda, and Katsuki Fujisawa. “Massive parallelization for finding shortest lattice vectors based on ubiquity generator framework”. In: SC20: International Conference for High Performance Computing, Networking, Storage and Analysis. IEEE. 2020, pp. 1–15.

<a id="ref.cmaplap"></a>
\[Tat+21a\] Nariaki Tateiwa, Yuji Shinano, Keiichiro Yamamura, Akihiro Yoshida, Shizuo Kaji, Masaya Yasuda, and Katsuki Fujisawa. “CMAP-LAP: Configurable massively parallel solver for lattice problems”. In: 2021 IEEE 28th International Conference on High Performance Computing, Data, and Analytics (HiPC). IEEE. 2021, pp. 42–52.

<a id="ref.cmapdeepbkz"></a>
\[Tat+21b\] Nariaki Tateiwa, Yuji Shinano, Masaya Yasuda, Shizuo Kaji, Keiichiro Yamamura, and Katsuki Fujisawa. Massively parallel sharing lattice basis reduction. eng. Tech. rep. 21-38. Takustr. 7, 14195 Berlin: ZIB, 2021.

<a id="ref.Ug"></a>
\[Ug\] UG: Ubiquity Generator framework. http://ug.zib.de/.
