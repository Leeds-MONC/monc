SRCSF = src/fftw.F90 src/pencilfft.F90 src/fftsolver.F90

BUILDDIR=build
COREDIR=../../model_core/build
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) -I $(FFTW_DIR)/include $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
