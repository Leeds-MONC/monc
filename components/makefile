# List of components to compile is the directories in this top level with makefiles in them
COMPONENTS=buoyancy casim casim_profile_dgs cfltest checkpointer clearsourceterms conditional_diagnostics_column conditional_diagnostics_whole coriolis damping debugger decomposition diagnostics_3d diffusion diverr fftsolver flux_budget forcing gridmanager haloswapper iobridge iterativesolver iterativesolver_single_prec kidreader lateral_bcs lowerbc lwrad_exponential meanprofiles modelsynopsis pdf_analysis petsc_solver pressuresource tvdadvection profile_diagnostics pstep pwadvection randomnoise scalar_diagnostics set_consistent_lowbc setfluxlook simplecloud simplesetup smagorinsky socrates_couple stepfields steppingdirection subgrid_profile_diagnostics swapsmooth terminationcheck thadvection tracers viscosity xiosbridge dephy_forcings

COMPONENT_HEADER_FILE=componentheaders.autogen
COMPONENT_REGISTRATIONS_FILE=componentregistrations.autogen
BUILDDIR=build

.PHONY: all $(COMPONENTS)
     
all: create-build-dirs cleancomponents $(COMPONENTS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

clean: cleancomponents	

cleancomponents:
	rm -Rf */$(BUILDDIR)/*
	rm -Rf $(BUILDDIR)/*
	rm -f $(COMPONENT_HEADER_FILE) $(COMPONENT_REGISTRATIONS_FILE)
	touch $(COMPONENT_HEADER_FILE) $(COMPONENT_REGISTRATIONS_FILE)

clean-build:
	rm -Rf */$(BUILDDIR)  
	rm -Rf $(BUILDDIR)

$(COMPONENTS):
	cd $@; $(MAKE) ; cp $(BUILDDIR)/* ../$(BUILDDIR)/.
	@echo "use $@_mod, only : $@_get_descriptor" >> $(COMPONENT_HEADER_FILE)
	@echo "call add_component(component_descriptions, $@_get_descriptor())" >> $(COMPONENT_REGISTRATIONS_FILE)
