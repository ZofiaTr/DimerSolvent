#------------------
#   Compilateur
#------------------
CC=g++

#------------------
#   Librairies
#------------------
ITINCLUDE3=   -I/usr/include/gsl -I/usr/include  #s'il faut mettre lien vers des fichiers include/ : portable
#ITINCLUDE3= -I/usr/include/ -I/home/stoltz/local/include/gsl/  -I/home/stoltz/local/include/ #-- cluster --
ITINCLUDE = ${ITINCLUDE3}
ITLIB= -L/usr/lib/ # emplacement des librairies requises : portable
#ITLIB= -L/usr/lib/ -L/home/stoltz/local/lib/ #-- cluster --

#----------------------------
#   Options de compilation
#----------------------------
INCLUDE_PATH=-I./ ${ITINCLUDE}
LDFLAGS=${ITLIB}
LIBRARY=-lstdc++ -lc -lm #-lgsl -lgslcblas 
WARN = -Wno-deprecated 
CPPFLAGS= -O1 -ffast-math -funroll-loops ${WARN} ${INCLUDE_PATH}
DEFINES = -Wall

.cpp.o: ; ${CC} ${DEFINES} ${INCLUDE_PATH} ${WARN} -g -c $*.cpp

#----------------------------
#   Fonctionalités du code
#----------------------------

OBJ = particle.cpp SBDTypeRandom.cpp Simulation.cpp main.cpp NeighborSearchGrid.cpp NeighborSearchGridPBC.cpp

ham: ${OBJ}
	${CC} ${CPPFLAGS}  ${DEFINES}  ${OBJ} -o ham++ ${LDFLAGS} ${LIBRARY}

clean:
	/bin/rm -rf *.o *.exe *~ err *++ data/*

