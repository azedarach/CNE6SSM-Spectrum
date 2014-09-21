DIR          := models/lowNE6SSM
MODNAME      := lowNE6SSM
SARAH_MODEL  := NE6SSM

lowNE6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

lowNE6SSM_MK     := \
		$(DIR)/module.mk

lowNE6SSM_TWO_SCALE_MK := \
		$(DIR)/two_scale_susy.mk \
		$(DIR)/two_scale_soft.mk

lowNE6SSM_SLHA_INPUT := \


lowNE6SSM_GNUPLOT := \
		$(DIR)/lowNE6SSM_plot_rgflow.gnuplot \
		$(DIR)/lowNE6SSM_plot_spectrum.gnuplot

LIBlowNE6SSM_SRC :=
EXElowNE6SSM_SRC :=

LIBlowNE6SSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBlowNE6SSM_SRC += \
		$(DIR)/lowNE6SSM_info.cpp \
		$(DIR)/lowNE6SSM_slha_io.cpp \
		$(DIR)/lowNE6SSM_physical.cpp \
		$(DIR)/lowNE6SSM_utilities.cpp \
		$(DIR)/lowNE6SSM_two_scale_convergence_tester.cpp \
		$(DIR)/lowNE6SSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/lowNE6SSM_two_scale_initial_guesser.cpp \
		$(DIR)/lowNE6SSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/lowNE6SSM_two_scale_model.cpp \
		$(DIR)/lowNE6SSM_two_scale_susy_parameters.cpp \
		$(DIR)/lowNE6SSM_two_scale_soft_parameters.cpp \
		$(DIR)/lowNE6SSM_two_scale_susy_scale_constraint.cpp
EXElowNE6SSM_SRC += \
		$(DIR)/run_lowNE6SSM.cpp \
		$(DIR)/scan_lowNE6SSM.cpp
LIBlowNE6SSM_HDR += \
		$(DIR)/lowNE6SSM_convergence_tester.hpp \
		$(DIR)/lowNE6SSM_high_scale_constraint.hpp \
		$(DIR)/lowNE6SSM_info.hpp \
		$(DIR)/lowNE6SSM_initial_guesser.hpp \
		$(DIR)/lowNE6SSM_input_parameters.hpp \
		$(DIR)/lowNE6SSM_low_scale_constraint.hpp \
		$(DIR)/lowNE6SSM_model.hpp \
		$(DIR)/lowNE6SSM_physical.hpp \
		$(DIR)/lowNE6SSM_slha_io.hpp \
		$(DIR)/lowNE6SSM_spectrum_generator.hpp \
		$(DIR)/lowNE6SSM_susy_scale_constraint.hpp \
		$(DIR)/lowNE6SSM_utilities.hpp \
		$(DIR)/lowNE6SSM_two_scale_convergence_tester.hpp \
		$(DIR)/lowNE6SSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/lowNE6SSM_two_scale_initial_guesser.hpp \
		$(DIR)/lowNE6SSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/lowNE6SSM_two_scale_model.hpp \
		$(DIR)/lowNE6SSM_two_scale_soft_parameters.hpp \
		$(DIR)/lowNE6SSM_two_scale_susy_parameters.hpp \
		$(DIR)/lowNE6SSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(DIR)/two_scale_susy.mk
-include $(DIR)/two_scale_soft.mk
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(DIR)/two_scale_susy.mk: run-metacode-$(MODNAME)
		@true
$(DIR)/two_scale_soft.mk: run-metacode-$(MODNAME)
		@true
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
LIBlowNE6SSM_SRC := $(sort $(LIBlowNE6SSM_SRC))
EXElowNE6SSM_SRC := $(sort $(EXElowNE6SSM_SRC))

LIBlowNE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBlowNE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBlowNE6SSM_SRC)))

EXElowNE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXElowNE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXElowNE6SSM_SRC)))

LIBlowNE6SSM_DEP := \
		$(LIBlowNE6SSM_OBJ:.o=.d)

EXElowNE6SSM_DEP := \
		$(EXElowNE6SSM_OBJ:.o=.d)

LIBlowNE6SSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_lowNE6SSM_OBJ := $(DIR)/run_lowNE6SSM.o
RUN_lowNE6SSM_EXE := $(DIR)/run_lowNE6SSM.x

SCAN_lowNE6SSM_OBJ := $(DIR)/scan_lowNE6SSM.o
SCAN_lowNE6SSM_EXE := $(DIR)/scan_lowNE6SSM.x

METACODE_STAMP_lowNE6SSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

SARAH_MODEL_FILES_lowNE6SSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		distclean-$(MODNAME) run-metacode-$(MODNAME)

all-$(MODNAME): $(LIBlowNE6SSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(lowNE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBlowNE6SSM_SRC) $(lowNE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBlowNE6SSM_HDR) $(lowNE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXElowNE6SSM_SRC) $(lowNE6SSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(lowNE6SSM_MK) $(lowNE6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(lowNE6SSM_TWO_SCALE_MK) $(lowNE6SSM_INSTALL_DIR)
ifneq ($(lowNE6SSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(lowNE6SSM_SLHA_INPUT) $(lowNE6SSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(lowNE6SSM_GNUPLOT) $(lowNE6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBlowNE6SSM_DEP)
		-rm -f $(EXElowNE6SSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBlowNE6SSM_OBJ)
		-rm -f $(EXElowNE6SSM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBlowNE6SSM)
		-rm -f $(RUN_lowNE6SSM_EXE)
		-rm -f $(SCAN_lowNE6SSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBlowNE6SSM_SRC) $(LIBlowNE6SSM_HDR) $(EXElowNE6SSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_lowNE6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_lowNE6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_lowNE6SSM)
		$(MATH) -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_lowNE6SSM)"
		@echo "Note: to regenerate lowNE6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_lowNE6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_lowNE6SSM):
		@true
endif

$(LIBlowNE6SSM_DEP) $(EXElowNE6SSM_DEP) $(LIBlowNE6SSM_OBJ) $(EXElowNE6SSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBlowNE6SSM_DEP) $(EXElowNE6SSM_DEP) $(LIBlowNE6SSM_OBJ) $(EXElowNE6SSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBlowNE6SSM): $(LIBlowNE6SSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_lowNE6SSM_EXE): $(RUN_lowNE6SSM_OBJ) $(LIBlowNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS) 

$(SCAN_lowNE6SSM_EXE): $(SCAN_lowNE6SSM_OBJ) $(LIBlowNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(FLIBS) $(THREADLIBS) 

ALLDEP += $(LIBlowNE6SSM_DEP) $(EXElowNE6SSM_DEP)
ALLSRC += $(LIBlowNE6SSM_SRC) $(EXElowNE6SSM_SRC)
ALLLIB += $(LIBlowNE6SSM)
ALLEXE += $(RUN_lowNE6SSM_EXE) $(SCAN_lowNE6SSM_EXE)
