#!/bin/bash
docker run -it \
        --platform linux/amd64 \
	--rm \
	-w /injection/ \
	-p 8989:8989 \
        -v $(pwd)/src/:/injection/src \
        -v $(pwd)/python/:/injection/python \
        -v $(pwd)/data/:/injection/data \
        -v $(pwd)/bin/:/injection/bin \
        -v $(pwd)/build/:/injection/build \
        -v $(pwd)/makefile:/injection/makefile \
        -v $(pwd)/jupyter_launch.sh:/injection/jupyter_launch.sh \
	gwat:dev
