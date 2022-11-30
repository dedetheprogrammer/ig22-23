#!/bin/bash
# Cambiar tone mapper: gamma puede recibir dos parametros que sea el coeficiente de equalizacion
# y el de gamma.
./tone_mapper -i new_scene.ppm -o scene.ppm --equalize 1.0 --gamma 2.2