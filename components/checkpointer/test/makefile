BUILDDIR=build
FRUIT_DIR=../../../misc/fruit/build
TEST_SRC=
MONC_SRC=
FFLAGS=-I $(FRUIT_DIR) -I $(BUILDDIR)
LFLAGS=
SRCSF = $(MONC_SRC) $(TEST_SRC)
TARGET=
FTN=ftn
OBJSF = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

SRCS=$(SRCSF)
OBJS=$(OBJSF)

GNU: FFLAGS += -cpp -DTEST_MODE -J $(BUILDDIR) -c
GNU: compile-test

Cray: FFLAGS +=  -e Z f03 -e m -J $(BUILDDIR) -DTEST_MODE -c
Cray: compile-test
	rm *.i *.mod

Intel: FFLAGS += -fpp -free -DTEST_MODE -c -std03 -module $(BUILDDIR) -Tf
Intel: compile-test

compile-test: create-build-dirs test-driver

create-build-dirs:
	rm -rf $(BUILDDIR); mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)

test-driver: $(OBJS)
	$(FTN) $(LFLAGS) -o $(TARGET) $(FRUIT_DIR)/*.o $(addprefix $(BUILDDIR)/, $(notdir $(OBJS)))
