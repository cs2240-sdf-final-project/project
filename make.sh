#!/bin/bash



docker build -t cs2240-final-build -f - . <<EOF || exit 1
FROM debian:bookworm-slim
RUN apt update && apt install -y neovim sudo build-essential curl wget llvm-16-dev clang-16 git cmake ninja-build libzstd-dev
WORKDIR /usr/local/src
RUN git clone https://github.com/EnzymeAD/Enzyme
WORKDIR /usr/local/src/Enzyme
RUN git checkout v0.0.173 && mkdir -p /usr/local/src/Enzyme/enzyme/build && cd /usr/local/src/Enzyme/enzyme/build && cmake -G Ninja .. -DLLVM_DIR=/usr/lib/llvm-16/lib/cmake/llvm -DLLVM_EXTERNAL_LIT=/usr/lib/llvm-16/build/utils/lit/lit.py && ninja -j1
RUN apt install -y imagemagick libomp-16-dev ffmpeg openmpi-bin libopenmpi-dev
EOF

cont=$(docker run -d --name cs2240-final-builder -v$(pwd):/usr/local/src/project cs2240-final-build sleep infinity)

docker exec --workdir=/usr/local/src/project $cont make $@

docker kill $cont
docker rm $cont
