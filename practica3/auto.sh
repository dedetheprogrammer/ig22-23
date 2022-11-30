#!/bin/bash
# ulimit -c unlimited <- Core dump.
g++ main.cpp -g -I ../libs -O3 -Wall -std=c++11 -o ray_tracer && ./ray_tracer