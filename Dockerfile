# -------------------------------------------------------------
# $ docker build -t <image_name> . [--build-arg PARALLEL=<int>]
# -------------------------------------------------------------
FROM ubuntu:20.04
RUN mkdir -p usr/lib
RUN apt-get update -y
RUN apt-get install -y sudo
RUN apt-get install -y tzdata
RUN apt-get install -y build-essential
RUN apt-get install -y libgmp3-dev zlib1g-dev libmpfr-dev openmpi-bin \
    libopenmpi-dev libreadline-dev libboost-dev \
    libntl-dev libopenblas-dev libeigen3-dev \
    libgsl-dev libgtest-dev libssl-dev \
    libblas-dev liblapack-dev liblapacke-dev
RUN apt-get install -y vim git gdb valgrind bash-completion wget

ARG PARALLEL=1

ENV HOME /cmaplap
WORKDIR ${HOME}

# copy local files into image
COPY . .

# install cmake with 3.22.3 veresion
RUN mkdir -p Library usr
RUN cd Library && wget https://github.com/Kitware/CMake/releases/download/v3.22.3/cmake-3.22.3.tar.gz
RUN cd Library && tar -xf cmake-3.22.3.tar.gz
RUN cd Library/cmake-3.22.3 && \
    ./bootstrap --prefix=$HOME/usr --parallel=${PARALLEL} && \
    make -j ${PARALLEL} && \
    make install

RUN echo "export PATH=${HOME}/usr/bin:\$PATH" > .bashrc
RUN echo source /etc/bash_completion >> .bashrc
