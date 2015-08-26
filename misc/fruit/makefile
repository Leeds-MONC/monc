BUILDDIR=build
FFLAGS=-I build
FRUIT_SRC= src/fruit_util.f90 src/fruit.f90
FRUIT_OBJS=$(patsubst %.f90,$(BUILDDIR)/%.o,$(FRUIT_SRC))

GNU: FTN=gfortran
GNU: FFLAGS+=-J $(BUILDDIR)
GNU: create-build-dirs $(FRUIT_OBJS) 

Cray: FTN=ftn
Cray: create-build-dirs $(FRUIT_OBJS)

Intel: FTN=ifort
Intel: FFLAGS+= -module $(BUILDDIR)
Intel: create-build-dirs $(FRUIT_OBJS)

$(FRUIT_OBJS): $(BUILDDIR)/%.o : %.f90
	$(FTN) $(FFLAGS) -c $< -o $(BUILDDIR)/$(notdir $@)

create-build-dirs:
	rm -Rf $(BUILDDIR) ; mkdir -p $(BUILDDIR)
