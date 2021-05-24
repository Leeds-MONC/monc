SRCSF = src/profile_diagnostics.F90

BUILDDIR=build
COREDIR=../../model_core/build
TVDDIR=../tvdadvection/build

FFLAGS=-I $(BUILDDIR) -I $(COREDIR) -I $(TVDDIR) $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
