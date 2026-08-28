#!/bin/sh
sudo docker run --rm -v $(pwd):/work -w /work docker.1ms.run/acculix/solvers:latest "$@" 
