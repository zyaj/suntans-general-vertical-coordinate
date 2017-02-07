SUNTANSHOME=../../main
include $(SUNTANSHOME)/Makefile.in

ifneq ($(MPIHOME),)
  CC = $(MPIHOME)/bin/mpicc
  MPIDEF = 
  MPIINC = -I$(MPIHOME)
else
  CC = gcc
  MPIDEF = -DNOMPI
  MPIINC = 
endif

ifneq ($(PARMETISHOME),)
  PARMETISINC = -I$(PARMETISHOME)/ParMETISLib
endif

LD = $(CC) 
CFLAGS = 
MATHLIB = -lm

EXEC = overbump
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
NUMPROCS = 2

all:	data 

test:	data
	sh $(EXEC).sh $(NUMPROCS)

test0.05:	r0.05 data
	sh $(EXEC).sh $(NUMPROCS)

test0.1:	r0.1 data
	sh $(EXEC).sh $(NUMPROCS)

test0.2:	r0.2 data
	sh $(EXEC).sh $(NUMPROCS)

test0.4:	r0.4 data
	sh $(EXEC).sh $(NUMPROCS)

data:	$(SUN)

r0.05:
	cp rundata/r0.05/cells.dat rundata/cells.dat
	cp rundata/r0.05/edges.dat rundata/edges.dat
	cp rundata/r0.05/points.dat rundata/points.dat

r0.1:
	cp rundata/r0.1/cells.dat rundata/cells.dat
	cp rundata/r0.1/edges.dat rundata/edges.dat
	cp rundata/r0.1/points.dat rundata/points.dat

r0.2:
	cp rundata/r0.2/cells.dat rundata/cells.dat
	cp rundata/r0.2/edges.dat rundata/edges.dat
	cp rundata/r0.2/points.dat rundata/points.dat

r0.4:
	cp rundata/r0.4/cells.dat rundata/cells.dat
	cp rundata/r0.4/edges.dat rundata/edges.dat
	cp rundata/r0.4/points.dat rundata/points.dat		

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o state.o
	cp initialization.o boundaries.o state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

clean:
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
