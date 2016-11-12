ifeq ($(strip $(PETSC_DIR)),)
	SRCSF = src/petsc_solver_stub.F90
	PETSC_INC=
else
	SRCSF = src/petsc_solver.F90
	PETSC_INC=-I $(PETSC_DIR)/include
endif

BUILDDIR=build
COREDIR=../../model_core/build
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) $(PETSC_INC) $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
