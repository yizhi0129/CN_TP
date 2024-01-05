CC=gcc -Wall -g
LIBSLOCAL=-L/usr/lib -llapacke -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include
OPTCLOCAL=-fPIC -march=native