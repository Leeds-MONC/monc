SRCSF = src/dephy_forcings.F90

BUILDDIR=build
COREDIR=../../model_core/build
LOWERBCDIR=../lowerbc/build
SETFLUXLOOKDIR=../setfluxlook/build
GRIDMANAGERDIR=../gridmanager/build
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) -I $(LOWERBCDIR) -I $(SETFLUXLOOKDIR) -I $(GRIDMANAGERDIR) $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
