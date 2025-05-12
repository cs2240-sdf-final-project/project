FROM debian:bookworm-slim
RUN apt update && apt install -y\
    sudo \
    build-essential\
    cmake\
    ninja-build\
    curl\
    wget\
    llvm-16-dev\
    clang-16\
    libomp-16-dev\
    libopenmpi-dev\
    openmpi-bin\
    git\
    libzstd-dev\
    imagemagick\
    ffmpeg
WORKDIR /usr/local/src
RUN git clone https://github.com/EnzymeAD/Enzyme
RUN cd /usr/local/src/Enzyme &&\
    git checkout v0.0.173 &&\
    mkdir -p /usr/local/src/Enzyme/enzyme/build &&\
    cd /usr/local/src/Enzyme/enzyme/build &&\
    cmake -G Ninja .. -DLLVM_DIR=/usr/lib/llvm-16/lib/cmake/llvm -DLLVM_EXTERNAL_LIT=/usr/lib/llvm-16/build/utils/lit/lit.py &&\
    ninja -j8
