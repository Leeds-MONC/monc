# List of components to compile is the directories in this top level with makefiles in them
COMPONENTS=$(subst /,,$(sort $(dir $(wildcard */makefile))))

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
