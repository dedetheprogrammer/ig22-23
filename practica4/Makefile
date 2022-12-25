#---------------------------------------------------------------
# File: Path tracer makefile.
# Author: Devid Dokash.
#---------------------------------------------------------------
# Macros:
CC=g++ # Compiler
#---------------------------------------------------------------
# Flags and options:
CXXFLAGS= -g -I. -I${LIBS} -O3 -pthread -std=c++17 -Wall 
#---------------------------------------------------------------
# Variables:
MAIN=main
EXEC=path_tracer
#---------------------------------------------------------------
# Libraries dependencies:
LIBS=../libs
OBJECTS=${LIBS}/objects
#---------------------------------------------------------------
# Linking:
${EXEC}: ${EXEC}.o ${OBJECTS}.o
	${CC} ${EXEC}.o ${OBJECTS}.o ${CXXFLAGS} -o ${EXEC}
#---------------------------------------------------------------
# Compiling:
${EXEC}.o: ${MAIN}.cpp
	${CC} -c ${MAIN}.cpp ${CXXFLAGS}

${OBJECTS}.o: ${OBJECTS}.hpp ${OBJECTS}.cpp
	${CC} -c ${OBJECTS}.cpp ${CXXFLAGS} -o ${OBJECTS}.o
	
