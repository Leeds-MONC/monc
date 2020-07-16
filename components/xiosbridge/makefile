SRCSF = src/xiosbridge.F90

BUILDDIR=build
COREDIR=../../core/build
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) -I $(XIOS_DIR)/inc $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
