FROM debian:10.9

WORKDIR /gwat

ENV CPATH $CPATH:/usr/local/include:/gwat/build_docker/install/include
ENV LIBRARY_PATH $LIBRARY_PATH:/usr/local/lib64:/usr/local/lib:/gwat/build_docker/install/lib
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/local/lib64:/usr/local/lib:/gwat/build_docker/install/lib

RUN apt update && apt install -y \
build-essential \
make \
cmake \
git  \
wget \
libfftw3-dev \
libgsl-dev \
libeigen3-dev \
cmake-curses-gui \
libhdf5-dev \
vim \
automake \
autoconf \
libtool \
devscripts \
python3 \
python3-pip \
swig \
doxygen


RUN pip3 install \
numpy \
matplotlib \
scipy \
pandas


## LALSUITE FRAME install ##
RUN mkdir /usr/local/LALFRAME
WORKDIR /usr/local/LALFRAME
RUN git clone https://git.ligo.org/virgo/virgoapp/Fr.git
WORKDIR /usr/local/LALFRAME/Fr
RUN cmake .
RUN make && make install 

## LALSUITE METAIO install ##
RUN mkdir /usr/local/LALMETAIO
WORKDIR /usr/local/LALMETAIO
RUN wget http://software.ligo.org/lscsoft/source/metaio-8.4.0.tar.gz
RUN tar -zxf metaio-8.4.0.tar.gz  
WORKDIR /usr/local/LALMETAIO/metaio-8.4.0
RUN ./configure --prefix=/usr/local
RUN make && make install 

## HDF5 install ##
#RUN mkdir /usr/local/hdf5 
#WORKDIR /usr/local/hdf5
#RUN wget https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_0/source/hdf5-1.12.0.tar.gz
#RUN tar -zxf hdf5-1.12.0.tar.gz  
#WORKDIR /usr/local/hdf5/hdf5-1.12.0/
#RUN ./configure --prefix=/usr/local && make check && make install && make check-install

## LALSUITE install ##
WORKDIR /usr/local
RUN git clone https://git.ligo.org/lscsoft/lalsuite.git
WORKDIR /usr/local/lalsuite
RUN ./00boot
RUN ./configure --prefix=/usr/local
RUN make -j && make install 


## ADOL-C install ##
RUN mkdir /usr/local/adolc
WORKDIR /usr/local/adolc
RUN git clone https://github.com/coin-or/ADOL-C.git
RUN mkdir ADOL-C/build
WORKDIR /usr/local/adolc/ADOL-C/build
RUN ../configure --prefix=/usr/local --with-openmp-flag=-fopenmp
RUN make && make install

WORKDIR /gwat
RUN mkdir build_docker
WORKDIR /gwat/build_docker
#RUN cmake ..  
