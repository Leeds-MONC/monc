BUILDDIR=build
MONC_CORE_BUILDDIR=../core/build
CHECKPOINT_COMPONENT_BUILDDIR=../components/checkpointer/build

DEFINES = -D_XOPEN_SOURCE=700
PLATFORM = -D__DARWIN

CFLAGS= -fPIC $(PLATFORM)
FFLAGS= -fPIC -I $(BUILDDIR) -I $(MONC_CORE_BUILDDIR) -I $(CHECKPOINT_COMPONENT_BUILDDIR) -I $(NETCDF_DIR)/include $(COMPILERRECURSIVE) $(COMPILERFFLAGS) -DENFORCE_THREAD_SAFETY $(PLATFORM)
LFLAGS=

CSOURCES := $(wildcard src/forthread/*.c)
FORTHREAD_F= src/forthread/forthread_types.F90 src/forthread/forthread_ciface.F90 src/forthread/forthread_data.F90 src/forthread/forthread.F90
SRCSF= $(FORTHREAD_F) src/ioclient.F90 src/saxparser.F90  src/configurationparser.F90 src/threadpool.F90  src/datautils.F90 src/mpicommunication.F90 src/diagnostics/inter-io/inter-io-specifics.F90 src/diagnostics/inter-io/reduction-inter-io.F90 src/diagnostics/inter-io/broadcast-inter-io.F90 src/diagnostics/inter-io/allreduction-inter-io.F90 src/diagnostics/inter-io/global-callback.F90 src/writers/writer_types.F90 src/writers/file_types/netcdf_filetype.F90 src/writers/time_manipulation/instantaneous_manipulation.F90  src/writers/time_manipulation/timeaveraged_manipulation.F90 src/writers/writer_federator.F90 src/writers/writer_field_manager.F90  src/diagnostics/operators/arithmetic-operator.F90 src/diagnostics/operators/reductionlocation-operator.F90 src/diagnostics/operators/operator.F90 src/diagnostics/diagnostic_federator.F90 src/ioserver.F90

OBJSF = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))
OBJSEXEC = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSEXEC))
COBJECTS := $(patsubst %.c,$(BUILDDIR)/%.o,$(CSOURCES))

SRCS=$(SRCSF)
OBJS=$(OBJSF)

all: create-build-dirs create-forthread $(OBJS) create-shared create-static
	

create-build-dirs:
	mkdir -p $(BUILDDIR)

create-forthread: $(COBJECTS)
#	$(MAKE) -C src/forthread
#	mv src/forthread/*.o src/forthread/*.mod build/.

$(COBJECTS):  $(BUILDDIR)/%.o : %.c
	$(CC) $(DEFINES) $(CFLAGS)  -c $< -o $(BUILDDIR)/$(notdir $@)

create-shared:
	$(FTN) -shared -Wl,-soname,libmonc_io.so -o libmonc_io.so $(BUILDDIR)/*.o

create-static:
	ar rcs libmonc_io.a $(BUILDDIR)/*.o

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $(DBG) $< -o $(BUILDDIR)/$(notdir $@)

$(OBJSEXEC) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) -cpp -J $(BUILDDIR) -c $(DBG) $< -o $(BUILDDIR)/$(notdir $@)

clean:
	rm -Rf $(BUILDDIR)/*

clean-build:
	rm -Rf $(BUILDDIR)
