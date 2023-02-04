# First variables
PROJECT_NAME = als-math
BUILD_DIR = /tmp/build_tmp/${PROJECT_NAME}
INCLUDE_DIR = /usr/include/${PROJECT_NAME}
LIB_DIR = /usr/lib
BIN_DIR = /usr/bin
CXX = g++
CXXFLAGS = -Wall -Wextra -Wpedantic -fPIC -O3 -pthread
LIBRARY_DEPENDENCIES = -l als-basic-utilities

all: ${BUILD_DIR}/libals-math.so

${BUILD_DIR}/Instantiation/Instantiation_double.o: ./Instantiation/Instantiation_double.cpp\
		Vector.cpp\
		Vector.hpp\
		Matrix.cpp\
		Matrix.hpp\
		DualNumbers.cpp\
		DualNumbers.hpp
	mkdir -p ${BUILD_DIR}/Instantiation
	${CXX} ${CXXFLAGS} -o ${BUILD_DIR}/Instantiation/Instantiation_double.o -c ./Instantiation/Instantiation_double.cpp

${BUILD_DIR}/Instantiation/Instantiation_complex.o: ./Instantiation/Instantiation_complex.cpp\
		Vector.cpp\
		Vector.hpp
	mkdir -p ${BUILD_DIR}/Instantiation
	${CXX} ${CXXFLAGS} -o ${BUILD_DIR}/Instantiation/Instantiation_complex.o -c ./Instantiation/Instantiation_complex.cpp
	
${BUILD_DIR}/libals-math.so: ${BUILD_DIR}/Instantiation/Instantiation_double.o\
		${BUILD_DIR}/Instantiation/Instantiation_complex.o\
		${BUILD_DIR}/ODESolvers.o\
		${BUILD_DIR}/Integration.o\
		${BUILD_DIR}/AlgebraicSolvers.o\
		${BUILD_DIR}/Interpolation.o
	${CXX} -shared ${CXXFLAGS} ${LIBRARY_DEPENDENCIES} -o ${BUILD_DIR}/libals-math.so\
		${BUILD_DIR}/Instantiation/Instantiation_double.o\
		${BUILD_DIR}/Instantiation/Instantiation_complex.o\
		${BUILD_DIR}/ODESolvers.o\
		${BUILD_DIR}/Integration.o\
		${BUILD_DIR}/AlgebraicSolvers.o\
		${BUILD_DIR}/Interpolation.o

install: ${BUILD_DIR}/libals-math.so
	mkdir -p ${INCLUDE_DIR}
	install -T AlgebraicSolvers.hpp ${INCLUDE_DIR}/AlgebraicSolvers.hpp
	install -T DualNumbers.hpp ${INCLUDE_DIR}/DualNumbers.hpp
	install -T Integration.hpp ${INCLUDE_DIR}/Integration.hpp
	install -T Interpolation.hpp ${INCLUDE_DIR}/Interpolation.hpp
	install -T Matrix.hpp ${INCLUDE_DIR}/Matrix.hpp
	install -T ODESolvers.hpp ${INCLUDE_DIR}/ODESolvers.hpp
	install -T Conj.hpp ${INCLUDE_DIR}/Conj.hpp
	install -T Vector.hpp ${INCLUDE_DIR}/Vector.hpp
	install -T ${BUILD_DIR}/libals-math.so ${LIB_DIR}/libals-math.so
	rm -r ${BUILD_DIR}

${BUILD_DIR}/%.o: %.cpp %.hpp
	mkdir -p $(@D)
	${CXX} ${CXXFLAGS} -o $@ -c $<

    
