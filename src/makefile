IDIR =../include
CC=mpicc
CFLAGS=-Wall -std=c99 -fopenmp -O3 -I$(IDIR) -IC:/GSL/include
NVCC=nvcc -ccbin g++
NVCC_CFLAGS=-I$(IDIR) -IC:/GSL/include -I/cm/shared/apps/cuda-7.5.18/NVIDIA_CUDA-7.5_Samples/common/inc/ -m64


SMS ?= 20 30 35 37 50 52

ifeq ($(SMS),)
$(info >>> WARNING - no SM architectures have been specified - waiving sample <<<)
SAMPLE_ENABLED := 0
endif

ifeq ($(GENCODE_FLAGS),)
# Generate SASS code for each SM architecture listed in $(SMS)
$(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))

# Generate PTX code from the highest SM architecture in $(SMS) to guarantee forward-compatibility
HIGHEST_SM := $(lastword $(sort $(SMS)))
ifneq ($(HIGHEST_SM),)
GENCODE_FLAGS += -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)
endif
endif


ODIR=.
LDIR =../lib

LIBS=-lm -lgsl -lgslcblas -L$(LDIR) -LC:/GSL/lib
NVCC_LIBS=-lcusolver -lcublas -lcusparse
CC_CUDA_LIBS=-lcudart

_DEPS = diffusionmatrix.h maths_functions.h particles.h moving_on.h forces.h initial_finalisation.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
_NVCC_DEPS = stochastic_force.h
NVCC_DEPS = $(patsubst %,$(IDIR)/%,$(_NVCC_DEPS))

_OBJ = main.o diffusionmatrix.o maths_functions.o particles.o moving_on.o forces.o initial_finalisation.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
_NVCC_OBJ = stochastic_force.o
NVCC_OBJ = $(patsubst %,$(ODIR)/%,$(_NVCC_OBJ))

args =

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) ${args}

$(ODIR)/%.o: %.cpp $(NVCC_DEPS)
	$(NVCC) $(NVCC_CFLAGS) $(GENCODE_FLAGS) -c -o $@ $<

../bin/executable: $(OBJ) $(NVCC_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(NVCC_FLAGS) $(LIBS) $(NVCC_LIBS) ${CC_CUDA_LIBS} ${args}

.PHONY: clean run

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

run:
	mpiexec -n $(n) ../bin/executable.exe ${args}
