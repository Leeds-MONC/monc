ifeq ($(strip $(CASIM_DIR)),)
	SRCSF = src/casim_stub.F90
else
	SRCSF = src/casim.F90
endif

BUILDDIR=build
COREDIR=../../model_core/build
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) -I $(CASIM_DIR)/build $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
