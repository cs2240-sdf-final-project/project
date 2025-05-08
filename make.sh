#!/bin/bash

# docker build -t cs2240-final-build . || exit 1
cont=$(docker run -d --name cs2240-final-builder -v$(pwd):/usr/local/src/project cs2240-final-build sleep infinity)
docker exec --workdir=/usr/local/src/project $cont make $@
docker kill $cont
docker rm $cont
