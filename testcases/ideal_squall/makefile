SRCSF = src/ideal_squall.F90

BUILDDIR=build
NETCDF_DIR?=/usr
COREDIR=../../model_core/build
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) -I $(NETCDF_DIR)/include $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
