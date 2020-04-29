SRCSF = src/socrates_couple.F90

BUILDDIR=build
COREDIR=../../core/build
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

ESRAD_RADIANCE_CORE=/data/local/fra23/MONC/es_radiation/src/radiance_core/
ESRAD_MODULE_CORE=src/modules_core/
ESRAD_CORE_DIR_LOCAL=src/um_edward_slingo_core_src
COPY_ESRAD_RADCORE=$(wildcard $(ESRAD_RADIANCE_CORE)/*.F90)
COPY_ESRAD_MODCORE=$(wildcard $(ESRAD_MODULE_CORE)/*.F90)
ESRAD_OBJS= $(patsubst %.F90, $(BUILDDIR)/%.o, $(notdir $(wildcard $(ESRAD_CORE_DIR_LOCAL)/*.F90)))

all: create-esrad-dirs copy_radcore copy_modcore create-build-dirs $(ESRAD_OBJS) $(OBJS)

copy_radcore:
	cp -f ${COPY_ESRAD_RADCORE} $(ESRAD_DIR_LOCAL)

copy_modcore:
	cp -f ${COPY_ESRAD_MODCORE} $(ESRAD_DIR_LOCAL)

create-esrad-dirs:
	mkdir -p $(ESRAD_DIR_LOCAL)	

create-build-dirs:
	mkdir -p $(BUILDDIR)

include $(ESRAD_DIR)/Mkdepend

$(ESRAD_OBJS) : $(BUILDDIR)/%.o : $(ESRAD_DIR_LOCAL)/%.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
