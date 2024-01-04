# Use a base image with the necessary dependencies
FROM ubuntu:latest

# Add the required paths to environment variables
ENV BLAS_PATH /etc/alternatives/blas.pc-x86_64-linux-gnu
ENV CBLAS_PATH /etc/alternatives/cblas.h-x86_64-linux-gnu
ENV LAPACK_PATH /etc/alternatives/lapack.pc-x86_64-linux-gnu

# Install build dependencies
RUN apt-get update && \
    apt-get install -y build-essential liblapacke-dev

# Set the working directory
WORKDIR /src

# Copy the source code and header files into the container
COPY src/tp_env.c .
COPY src/tp_poisson1D_direct.c .
#COPY src/tp_poisson1D_iter.c .
COPY src/lib_poisson1D.c .
#COPY src/lib_poisson1D_richardson.c .
COPY include ./include

# Compile the program
RUN gcc -o tp_env tp_env.c -I./include -llapacke -llapack -lblas -lm
RUN gcc -o tp_poisson1d_direct tp_poisson1D_direct.c lib_poisson1D.c -I./include -llapacke -llapack -lblas -lm
#RUN gcc -o tp_poisson1d_iter tp_poisson1D_iter.c lib_poisson1D_richardson.c lib_poisson1D.c -I./include -llapacke -llapack -lblas -lm

# Run the program