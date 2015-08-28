SRCSF = src/iobridge.F90

BUILDDIR=build
COREDIR=../../model_core/build
IODIR=../../io/build
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) -I $(IODIR) $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
