FOR = gfortran -IMODF -JMODF -O3 -ffast-math -funroll-loops
#FOR = gfortran -IMODF -JMODF -fbounds-check
#FOR = ifort -fast -IMODF -module MODF -heap-arrays
#FOR = ifort -CB -traceback -IMODF -module MODF -heap-arrays

EXE = xnlmg

LIB = ../../LIBRARY

OBJ = OF/global.o OF/utilities.o OF/problem.o OF/nonlinearmultigrid.o OF/nlmgdriver.o

$(EXE): $(OBJ)
	$(FOR) $(OBJ) -o $(EXE)

OF/global.o: $(LIB)/global.f90
	$(FOR) -c $(LIB)/global.f90 -o OF/global.o

OF/utilities.o: $(LIB)/utilities.f90
	$(FOR) -c $(LIB)/utilities.f90 -o OF/utilities.o

OF/problem.o: problem.f90
	$(FOR) -c problem.f90 -o OF/problem.o

OF/nonlinearmultigrid.o: $(LIB)/nonlinearmultigrid.f90
	$(FOR) -c $(LIB)/nonlinearmultigrid.f90 -o OF/nonlinearmultigrid.o

OF/nlmgdriver.o: ./nlmgdriver.f90
	$(FOR) -c ./nlmgdriver.f90 -o OF/nlmgdriver.o

reset:
	rm $(EXE) MODF/* OF/*

remove:
	rm OUT/*.dat ./output.dat
