# List of test cases to compile is the directories in this top level with makefiles in them
TESTCASES = $(subst /,,$(sort $(dir $(wildcard */makefile))))

TESTCASE_HEADER_FILE=testcaseheaders.autogen
TESTCASE_REGISTRATIONS_FILE=testcaseregistrations.autogen
BUILDDIR=build

.PHONY: all $(TESTCASES)
     
all: create-build-dirs cleantestcases $(TESTCASES)

create-build-dirs:
	mkdir -p $(BUILDDIR)

clean: cleantestcases

cleantestcases:
	rm -Rf */$(BUILDDIR)/*
	rm -Rf $(BUILDDIR)/*
	rm -f $(TESTCASE_HEADER_FILE) $(TESTCASE_REGISTRATIONS_FILE)
	touch $(TESTCASE_HEADER_FILE) $(TESTCASE_REGISTRATIONS_FILE)

clean-build:
	rm -Rf */$(BUILDDIR) 
	rm -Rf $(BUILDDIR)   

$(TESTCASES):
	cd $@; $(MAKE) ; cp $(BUILDDIR)/* ../$(BUILDDIR)/.
	@echo "use $@_mod, only : $@_get_descriptor" >> $(TESTCASE_HEADER_FILE)
	@echo "call add_component(component_descriptions, $@_get_descriptor())" >> $(TESTCASE_REGISTRATIONS_FILE)
