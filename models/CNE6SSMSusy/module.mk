DIR          := models/CNE6SSMSusy
MODNAME      := CNE6SSMSusy
SARAH_MODEL  := NE6SSM

CNE6SSMSusy_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

CNE6SSMSusy_MK     := \
		$(DIR)/module.mk

CNE6SSMSusy_SUSY_BETAS_MK := \
		$(DIR)/susy.mk

CNE6SSMSusy_BETAS_MK := \
		$(CNE6SSMSusy_SUSY_BETAS_MK)

CNE6SSMSusy_SLHA_INPUT := \
		$(DIR)/LesHouches.in.CNE6SSM

CNE6SSMSusy_GNUPLOT := \
		$(DIR)/CNE6SSMSusy_plot_rgflow.gnuplot

CNE6SSMSusy_TARBALL := \
		$(MODNAME).tar.gz

LIBCNE6SSMSusy_SRC := \
		$(DIR)/CNE6SSMSusy_info.cpp \
		$(DIR)/CNE6SSMSusy_susy_parameters.cpp

EXECNE6SSMSusy_SRC :=

LIBCNE6SSMSusy_HDR := \
		$(DIR)/CNE6SSMSusy_constraint_handler.hpp \
		$(DIR)/CNE6SSMSusy_convergence_tester.hpp \
		$(DIR)/CNE6SSMSusy_high_scale_constraint.hpp \
		$(DIR)/CNE6SSMSusy_info.hpp \
		$(DIR)/CNE6SSMSusy_initial_guesser.hpp \
		$(DIR)/CNE6SSMSusy_input_parameters.hpp \
		$(DIR)/CNE6SSMSusy_low_scale_constraint.hpp \
		$(DIR)/CNE6SSMSusy_model.hpp \
		$(DIR)/CNE6SSMSusy_susy_parameters.hpp \
		$(DIR)/CNE6SSMSusy_utilities.hpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBCNE6SSMSusy_SRC += \
		$(DIR)/CNE6SSMSusy_two_scale_constraint_handler.cpp \
		$(DIR)/CNE6SSMSusy_two_scale_convergence_tester.cpp \
		$(DIR)/CNE6SSMSusy_two_scale_high_scale_constraint.cpp \
		$(DIR)/CNE6SSMSusy_two_scale_initial_guesser.cpp \
		$(DIR)/CNE6SSMSusy_two_scale_input_parameters.cpp \
		$(DIR)/CNE6SSMSusy_two_scale_low_scale_constraint.cpp \
		$(DIR)/CNE6SSMSusy_two_scale_model.cpp

EXECNE6SSMSusy_SRC +=

LIBCNE6SSMSusy_HDR += \
		$(DIR)/CNE6SSMSusy_spectrum_generator.hpp \
		$(DIR)/CNE6SSMSusy_two_scale_constraint_handler.hpp \
		$(DIR)/CNE6SSMSusy_two_scale_convergence_tester.hpp \
		$(DIR)/CNE6SSMSusy_two_scale_high_scale_constraint.hpp \
		$(DIR)/CNE6SSMSusy_two_scale_initial_guesser.hpp \
		$(DIR)/CNE6SSMSusy_two_scale_input_parameters.hpp \
		$(DIR)/CNE6SSMSusy_two_scale_low_scale_constraint.hpp \
		$(DIR)/CNE6SSMSusy_two_scale_model.hpp
endif

ifneq ($(findstring semianalytic,$(ALGORITHMS)),)
LIBCNE6SSMSusy_SRC +=

EXECNE6SSMSusy_SRC +=

LIBCNE6SSMSusy_HDR +=

endif

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(CNE6SSMSusy_SUSY_BETAS_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CNE6SSMSusy_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

# remove duplicates in case all algorithms are used
LIBCNE6SSMSusy_SRC := $(sort $(LIBCNE6SSMSusy_SRC))
EXECNE6SSMSusy_SRC := $(sort $(EXECNE6SSMSusy_SRC))

LIBCNE6SSMSusy_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCNE6SSMSusy_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCNE6SSMSusy_SRC)))

EXECNE6SSMSusy_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECNE6SSMSusy_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECNE6SSMSusy_SRC)))

LIBCNE6SSMSusy_DEP := \
		$(LIBCNE6SSMSusy_OBJ:.o=.d)

EXECNE6SSMSusy_DEP := \
		$(EXECNE6SSMSusy_OBJ:.o=.d)

LIBCNE6SSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

METACODE_STAMP_CNE6SSMSusy := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_CNE6SSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBCNE6SSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(CNE6SSMSusy_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNE6SSMSusy_SRC) $(CNE6SSMSusy_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNE6SSMSusy_HDR) $(CNE6SSMSusy_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXECNE6SSMSusy_SRC) $(CNE6SSMSusy_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(CNE6SSMSusy_MK) $(CNE6SSMSusy_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(CNE6SSMSusy_BETAS_MK) $(CNE6SSMSusy_INSTALL_DIR)
ifneq ($(CNE6SSMSusy_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(CNE6SSMSusy_SLHA_INPUT) $(CNE6SSMSusy_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(CNE6SSMSusy_GNUPLOT) $(CNE6SSMSusy_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBCNE6SSMSusy_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCNE6SSMSusy_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBCNE6SSM)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(CNE6SSMSusy_TARBALL) \
		$(LIBCNE6SSMSusy_SRC) $(LIBCNE6SSMSusy_HDR) \
		$(EXECNE6SSMSusy_SRC) \
		$(CNE6SSMSusy_MK) $(CNE6SSMSusy_BETAS_MK) \
		$(CNE6SSMSusy_SLHA_INPUT) $(CNE6SSMSusy_GNUPLOT)

$(LIBCNE6SSMSusy_SRC) $(LIBCNE6SSMSusy_HDR) $(EXECNE6SSMSusy_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CNE6SSMSusy)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CNE6SSMSusy): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CNE6SSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_CNE6SSM)"
		@echo "Note: to regenerate CNE6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CNE6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CNE6SSMSusy):
		@true
endif

$(LIBCNE6SSMSusy_DEP) $(EXECNE6SSMSusy_DEP) $(LIBCNE6SSMSusy_OBJ) $(EXECNE6SSMSusy_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCNE6SSMSusy_DEP) $(EXECNE6SSMSusy_DEP) $(LIBCNE6SSMSusy_OBJ) $(EXECNE6SSMSusy_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBCNE6SSM): $(LIBCNE6SSMSusy_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBCNE6SSMSusy_DEP) $(EXECNE6SSMSusy_DEP)
ALLSRC += $(LIBCNE6SSMSusy_SRC) $(EXECNE6SSMSusy_SRC)
ALLLIB += $(LIBCNE6SSM)
ALLEXE +=
