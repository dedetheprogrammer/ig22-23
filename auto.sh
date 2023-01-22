#!/bin/bash
# usage ./auto.sh <clamp_val or 0> <equalize_val or 0> <gamma_val or 0> [render]
if [ "$4" == "render" ]; then
    g++ main.cpp -g -I./libs -O3 -pthread -std=c++17 -Wall -o path_tracer && ./path_tracer
    if [ "$?" != 0 ]; then exit 1; fi
fi
# if [ ! -e tone_mapper* ]; then
#     g++ ../practica2/main.cpp -g -I../libs -O3 -pthread -Wall -std=c++11 -o tone_mapper
# fi
if [ "$1" != 0 ]; then c="--clamp $1"; fi
if [ "$2" != 0 ]; then e="--equalize $2"; fi
if [ "$3" != 0 ]; then g="--gamma $3"; fi
./tone_mapper -i ./new_scene.ppm -o ./scene.ppm -c $c $e $g
