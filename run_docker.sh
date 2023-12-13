#!/bin/bash

echo "Building and running the docker file."


docker run -it -v "$(pwd):/tp" $(docker build -q .)