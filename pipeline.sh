#!/bin/sh

echo -e "Preparing environment..."

mkdir ./data
mkdir ./plots

export OMP_NUM_THREADS=3

echo -e "Compiling MD program...\n"

g++ -fopenmp ./md.cc -lm &&

echo -e "Running MD program...\n"

time ./a.out && 

echo -e "Plotting data...\n"

python ./plot.py 

echo -e "Completion Successful: Please refer to ./plots for .png plots, and ./data for data outputs!"
