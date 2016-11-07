SRCSF = src/pressuresource.F90

BUILDDIR=build
COREDIR=../../model_core/build
FFLAGS_SPECIAL=$(filter-out -O3,$(COMPILERFFLAGS))
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) -O1 $(FFLAGS_SPECIAL)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
