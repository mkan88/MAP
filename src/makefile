IDIR =../include
CC=mpicc
CFLAGS=-Wall -std=c99 -fopenmp -O3 -I$(IDIR) -IC:/GSL/include


ODIR=.
LDIR =../lib

LIBS=-lm -lgsl -lgslcblas -L$(LDIR) -LC:/GSL/lib

_DEPS = diffusionmatrix.h maths_functions.h particles.h stochastic_force.h moving_on.h forces.h initial_finalisation.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o diffusionmatrix.o maths_functions.o particles.o stochastic_force.o moving_on.o forces.o initial_finalisation.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

args =

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) ${args}

../bin/executable: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) ${args}

.PHONY: clean run

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

run:
	mpiexec -n $(n) ../bin/executable.exe ${args}
