FOR = gfortran -IMODF -JMODF -O3 -ffast-math -funroll-loops
#FOR = gfortran -IMODF -JMODF -fbounds-check
#FOR = ifort -fast -IMODF -module MODF -heap-arrays
#FOR = ifort -check -traceback -IMODF -module MODF -heap-arrays

EXE = xlmg

LIB = ../../LIBRARY

OBJ = OF/global.o OF/utilities.o OF/problem.o OF/linearmultigrid.o OF/lmgdriver.o

$(EXE): $(OBJ)
	$(FOR) $(OBJ) -o $(EXE)

OF/global.o: $(LIB)/global.f90
	$(FOR) -c $(LIB)/global.f90 -o OF/global.o

OF/utilities.o: $(LIB)/utilities.f90
	$(FOR) -c $(LIB)/utilities.f90 -o OF/utilities.o

OF/problem.o: problem.f90
	$(FOR) -c problem.f90 -o OF/problem.o

OF/linearmultigrid.o: $(LIB)/linearmultigrid.f90
	$(FOR) -c $(LIB)/linearmultigrid.f90 -o OF/linearmultigrid.o

OF/lmgdriver.o: ./lmgdriver.f90
	$(FOR) -c ./lmgdriver.f90 -o OF/lmgdriver.o

reset:
	rm $(EXE) MODF/* OF/*

remove:
	rm OUT/*.dat ./output.dat
