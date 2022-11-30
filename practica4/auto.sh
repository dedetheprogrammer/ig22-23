#!/bin/bash
g++ main.cpp -g -I ../libs -O3 -pthread -Wall -std=c++11 -o path_tracer && ./path_tracer