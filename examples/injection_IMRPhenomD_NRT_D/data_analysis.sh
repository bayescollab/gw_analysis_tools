#!/bin/bash
docker run  \
        --platform linux/amd64 \
	-p 8990:8989\
	--rm \
	--name gw-injection-example-jupyter \
	-w /injection/ \
        -v $(pwd)/src:/injection/8169056/src \
        -v $(pwd)/data:/injection/8169056/data \
        -v $(pwd)/bin:/injection/8169056/bin \
        -v $(pwd)/build:/injection/8169056/build \
        -v $(pwd)/python:/injection/8169056/python \
        -v $(pwd)/makefile:/injection/8169056/makefile \
        -v $(pwd)/jupyter_launch.sh:/injection/8169056/jupyter_launch.sh \
	scottperkins/gwat.jupyter 
