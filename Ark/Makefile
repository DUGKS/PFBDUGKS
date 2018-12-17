vpath %.h Header:Lib:Constant
vpath %.c SourceC:Lib

Target = PFDUGKS_Basilisk
QModel = D2Q9AC

CC = qcc

CFLAGS = -fopenmp 

MPFLAG = -lm

#DIR_SRC = $(shell find ${ROOT_SRC} -maxdepth 1 -type d)

DIR_SRC = ./

SRC = $(wildcard *.c)

DIR_OBJ = ./Obj

DIR_OUT = ../FlowField/global ../FlowField/Convergence

subObjs = $(patsubst %.c,${DIR_OBJ}/%.o,$(notdir ${SRC}))

objects = $(subObjs) ${DIR_OBJ}/${QModel}.o

# ALL:

# 	@echo "SRC:" $(notdir ${SRC})

# $@ : target file
# $^ : all of the dependent file

# $(Target):$(objects)

# 	$(CC) $(objects) -o $@ ${MPFLAG}

# ${DIR_OBJ}/%.o:%.c

# 	$(CC) $(CFLAGS) -c $< -o $@ ${MPFLAG}

# ${DIR_OBJ}/${QModel}.o:${QModel}.c

# 	$(CC) $(CFLAGS) -c $< -o $@ ${MPFLAG}

$(Target) : main.c
	$(CC) $(CFLAGS) main.c -o $(Target) ${MPFLAG}

.PHONY:clean
clean:
	find -name "${Target}" -exec mv {} ~/.trash/ \;
	find ${DIR_OUT} -name "*.dat" -exec rm {} \; 

# headers = ZeroDeclaration.h Mesh_2D.h ZeroConstant.h $(QModel).h ZeroReference.h ZeroFlip.h

# main : $(objects)
# 	$(CC) -o main $(objects) -fopenmp

# main.o : main.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 main.cpp

# MeshConstruct.o : MeshConstruct.cpp MeshConstructFunc.h $(headers)
# 	$(CC) -c -O3 -std=c++11 MeshConstruct.cpp

# Preprocess.o : Preprocess.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 Preprocess.cpp

# DUGKS2DSolver.o : DUGKS2DSolver.cpp NewtonCotes.h $(headers)
# 	$(CC) -c -O3 -std=c++11 DUGKS2DSolver.cpp -fopenmp

# $(QModel).o : $(QModel).cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 $(QModel).cpp

# Mesh_2D.o : Mesh_2D.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 Mesh_2D.cpp

# Output.o : Output.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 Output.cpp

# .PHONY : clean
# clean :
# 	rm pexe $(objects)
