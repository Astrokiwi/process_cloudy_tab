CC       =  gcc
CXX      =  g++

#MPICHLIB =  -L/usr/lib/openmpi/lib


#LIBS   = -lm $(MPICHLIB)
#INCL = $(HDF5INCL)

OPTIONS = -O3 -ldl -Wall -mcmodel=medium
#OPTIONS = -O0 -g3 -ldl -Wall -mcmodel=medium
#CFLAGS = $(LIBS) $(OPTIONS) ${INCL}

.SUFFIXES : .o .c .cpp .h

.cpp.o :
	${CXX} ${OPTIONS} ${INCL} -c  ${LIBS} $<

.c.o :
	${CC} ${OPTIONS} ${INCL} -c  ${LIBS} $<

interplib:
	swig -python -c++ -o _tab_interp_module.cc tab_interp.i
# 	python3 setup.py build_ext --inplace
	python setup.py build_ext --inplace

INTERP_OBJS=dust_temp_interp.o

test: ${INTERP_OBJS}
	${CXX} ${OPTIONS} ${INCL} ${INTERP_OBJS} ${LIBS} -o ../tab_interp_test

clean :
	rm *.o
