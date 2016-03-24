DIR          := models/CNE6SSM
MODNAME      := CNE6SSM
SARAH_MODEL  := NE6SSM

CNE6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

CNE6SSM_MK     := \
		$(DIR)/module.mk

CNE6SSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

CNE6SSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

CNE6SSM_BETAS_MK := \
		$(CNE6SSM_SUSY_BETAS_MK) \
		$(CNE6SSM_SOFT_BETAS_MK)

CNE6SSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.CNE6SSM

CNE6SSM_GNUPLOT := \
		$(DIR)/CNE6SSM_plot_rgflow.gnuplot \
		$(DIR)/CNE6SSM_plot_spectrum.gnuplot

CNE6SSM_TARBALL := \
		$(MODNAME).tar.gz

LIBCNE6SSM_SRC := \
		$(DIR)/CNE6SSM_higgs_upper_bound.cpp \
		$(DIR)/CNE6SSM_info.cpp \
		$(DIR)/CNE6SSM_mass_eigenstates.cpp \
		$(DIR)/CNE6SSM_physical.cpp \
		$(DIR)/CNE6SSM_slha_io.cpp \
		$(DIR)/CNE6SSM_soft_parameters.cpp \
		$(DIR)/CNE6SSM_susy_parameters.cpp \
		$(DIR)/CNE6SSM_utilities.cpp \
		$(DIR)/CNE6SSM_scan_utilities.cpp

EXECNE6SSM_SRC :=

LIBCNE6SSM_HDR := \
		$(DIR)/CNE6SSM_constraint_handler.hpp \
		$(DIR)/CNE6SSM_convergence_tester.hpp \
		$(DIR)/CNE6SSM_higgs_upper_bound.hpp \
		$(DIR)/CNE6SSM_high_scale_constraint.hpp \
		$(DIR)/CNE6SSM_info.hpp \
		$(DIR)/CNE6SSM_initial_guesser.hpp \
		$(DIR)/CNE6SSM_input_parameters.hpp \
		$(DIR)/CNE6SSM_low_scale_constraint.hpp \
		$(DIR)/CNE6SSM_mass_eigenstates.hpp \
		$(DIR)/CNE6SSM_model.hpp \
		$(DIR)/CNE6SSM_model_slha.hpp \
		$(DIR)/CNE6SSM_physical.hpp \
		$(DIR)/CNE6SSM_slha_io.hpp \
		$(DIR)/CNE6SSM_soft_parameters.hpp \
		$(DIR)/CNE6SSM_susy_parameters.hpp \
		$(DIR)/CNE6SSM_susy_scale_constraint.hpp \
		$(DIR)/CNE6SSM_utilities.hpp \
		$(DIR)/CNE6SSM_scan_utilities.hpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBCNE6SSM_SRC += \
		$(DIR)/CNE6SSM_scan_parameters.cpp \
		$(DIR)/CNE6SSM_two_scale_constraint_handler.cpp \
		$(DIR)/CNE6SSM_two_scale_convergence_tester.cpp \
		$(DIR)/CNE6SSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/CNE6SSM_two_scale_initial_guesser.cpp \
		$(DIR)/CNE6SSM_two_scale_input_parameters.cpp \
		$(DIR)/CNE6SSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/CNE6SSM_two_scale_model.cpp \
		$(DIR)/CNE6SSM_two_scale_model_slha.cpp \
		$(DIR)/CNE6SSM_two_scale_susy_scale_constraint.cpp
EXECNE6SSM_SRC += \
		$(DIR)/gridscan_CNE6SSM.cpp \
		$(DIR)/gridscan_rge_coeffs_CNE6SSM.cpp \
		$(DIR)/rge_coefficients_CNE6SSM.cpp \
		$(DIR)/run_CNE6SSM.cpp \
		$(DIR)/run_cmd_line_CNE6SSM.cpp \
		$(DIR)/scan_CNE6SSM.cpp
LIBCNE6SSM_HDR += \
		$(DIR)/CNE6SSM_scan_parameters.hpp \
		$(DIR)/CNE6SSM_spectrum_generator.hpp \
		$(DIR)/CNE6SSM_two_scale_constraint_handler.hpp \
		$(DIR)/CNE6SSM_two_scale_convergence_tester.hpp \
		$(DIR)/CNE6SSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/CNE6SSM_two_scale_initial_guesser.hpp \
		$(DIR)/CNE6SSM_two_scale_input_parameters.hpp \
		$(DIR)/CNE6SSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/CNE6SSM_two_scale_model.hpp \
		$(DIR)/CNE6SSM_two_scale_model_slha.hpp \
		$(DIR)/CNE6SSM_two_scale_susy_scale_constraint.hpp
endif

ifneq ($(findstring semianalytic,$(ALGORITHMS)),)
LIBCNE6SSM_SRC += \
		$(DIR)/CNE6SSM_semi_two_scale_constraint_handler.cpp \
		$(DIR)/CNE6SSM_semi_two_scale_convergence_tester.cpp \
		$(DIR)/CNE6SSM_semi_two_scale_high_scale_constraint.cpp \
		$(DIR)/CNE6SSM_semi_two_scale_initial_guesser.cpp \
		$(DIR)/CNE6SSM_semi_two_scale_input_parameters.cpp \
		$(DIR)/CNE6SSM_semi_two_scale_low_scale_constraint.cpp \
		$(DIR)/CNE6SSM_semi_two_scale_model.cpp \
		$(DIR)/CNE6SSM_semi_two_scale_model_slha.cpp \
		$(DIR)/CNE6SSM_semi_two_scale_susy_scale_constraint.cpp \
		$(DIR)/CNE6SSM_susy_two_scale_convergence_tester.cpp

EXECNE6SSM_SRC += \
		$(DIR)/run_semianalytic_CNE6SSM.cpp \
		$(DIR)/gridscan_semianalytic_CNE6SSM.cpp \
		$(DIR)/scan_semianalytic_CNE6SSM.cpp \
		$(DIR)/generate_semianalytic_points.cpp

ifneq ($(findstring addons/susyhd_call,$(ADDONS)),)
EXECNE6SSM_SRC += \
		$(DIR)/run_susyhd_CNE6SSM.cpp \
		$(DIR)/get_susyhd_higgs_mass.cpp
endif

LIBCNE6SSM_HDR += \
		$(DIR)/CNE6SSM_semi_constraint_handler.hpp \
		$(DIR)/CNE6SSM_semi_convergence_tester.hpp \
		$(DIR)/CNE6SSM_semi_high_scale_constraint.hpp \
		$(DIR)/CNE6SSM_semi_initial_guesser.hpp \
		$(DIR)/CNE6SSM_semi_input_parameters.hpp \
		$(DIR)/CNE6SSM_semi_low_scale_constraint.hpp \
		$(DIR)/CNE6SSM_semi_model.hpp \
		$(DIR)/CNE6SSM_semi_model_slha.hpp \
		$(DIR)/CNE6SSM_semi_susy_scale_constraint.hpp \
		$(DIR)/CNE6SSM_semianalytic_spectrum_generator.hpp \
		$(DIR)/CNE6SSM_semi_two_scale_constraint_handler.hpp \
		$(DIR)/CNE6SSM_semi_two_scale_convergence_tester.hpp \
		$(DIR)/CNE6SSM_semi_two_scale_high_scale_constraint.hpp \
		$(DIR)/CNE6SSM_semi_two_scale_initial_guesser.hpp \
		$(DIR)/CNE6SSM_semi_two_scale_input_parameters.hpp \
		$(DIR)/CNE6SSM_semi_two_scale_low_scale_constraint.hpp \
		$(DIR)/CNE6SSM_semi_two_scale_model.hpp \
		$(DIR)/CNE6SSM_semi_two_scale_model_slha.hpp \
		$(DIR)/CNE6SSM_semi_two_scale_susy_scale_constraint.hpp \
		$(DIR)/CNE6SSM_susy_convergence_tester.hpp \
		$(DIR)/CNE6SSM_susy_two_scale_convergence_tester.hpp
endif

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(CNE6SSM_SUSY_BETAS_MK)
-include $(CNE6SSM_SOFT_BETAS_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CNE6SSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CNE6SSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
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
LIBCNE6SSM_SRC := $(sort $(LIBCNE6SSM_SRC))
EXECNE6SSM_SRC := $(sort $(EXECNE6SSM_SRC))

LIBCNE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCNE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCNE6SSM_SRC)))

EXECNE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECNE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECNE6SSM_SRC)))

LIBCNE6SSM_DEP := \
		$(LIBCNE6SSM_OBJ:.o=.d)

EXECNE6SSM_DEP := \
		$(EXECNE6SSM_OBJ:.o=.d)

LIBCNE6SSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

GRIDSCAN_CNE6SSM_OBJ := $(DIR)/gridscan_CNE6SSM.o
GRIDSCAN_CNE6SSM_EXE := $(DIR)/gridscan_CNE6SSM.x

GRIDSCAN_RGE_CNE6SSM_OBJ := $(DIR)/gridscan_rge_coeffs_CNE6SSM.o
GRIDSCAN_RGE_CNE6SSM_EXE := $(DIR)/gridscan_rge_coeffs_CNE6SSM.x

RGE_COEFF_CNE6SSM_OBJ := $(DIR)/rge_coefficients_CNE6SSM.o
RGE_COEFF_CNE6SSM_EXE := $(DIR)/rge_coefficients_CNE6SSM.x

RUN_CNE6SSM_OBJ := $(DIR)/run_CNE6SSM.o
RUN_CNE6SSM_EXE := $(DIR)/run_CNE6SSM.x

RUN_CMD_LINE_CNE6SSM_OBJ := $(DIR)/run_cmd_line_CNE6SSM.o
RUN_CMD_LINE_CNE6SSM_EXE := $(DIR)/run_cmd_line_CNE6SSM.x

SCAN_CNE6SSM_OBJ := $(DIR)/scan_CNE6SSM.o
SCAN_CNE6SSM_EXE := $(DIR)/scan_CNE6SSM.x

RUN_SEMI_CNE6SSM_OBJ := $(DIR)/run_semianalytic_CNE6SSM.o
RUN_SEMI_CNE6SSM_EXE := $(DIR)/run_semianalytic_CNE6SSM.x

RUN_SUSYHD_CNE6SSM_OBJ := $(DIR)/run_susyhd_CNE6SSM.o
RUN_SUSYHD_CNE6SSM_EXE := $(DIR)/run_susyhd_CNE6SSM.x

GRIDSCAN_SEMI_CNE6SSM_OBJ := $(DIR)/gridscan_semianalytic_CNE6SSM.o
GRIDSCAN_SEMI_CNE6SSM_EXE := $(DIR)/gridscan_semianalytic_CNE6SSM.x

SCAN_SEMI_CNE6SSM_OBJ := $(DIR)/scan_semianalytic_CNE6SSM.o
SCAN_SEMI_CNE6SSM_EXE := $(DIR)/scan_semianalytic_CNE6SSM.x

GET_HIGGS_MASS_OBJ := $(DIR)/get_susyhd_higgs_mass.o
GET_HIGGS_MASS_EXE := $(DIR)/get_susyhd_higgs_mass.x

GENERATE_POINTS_OBJ := $(DIR)/generate_semianalytic_points.o
GENERATE_POINTS_EXE := $(DIR)/generate_semianalytic_points.x

METACODE_STAMP_CNE6SSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

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
		install -d $(CNE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNE6SSM_SRC) $(CNE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNE6SSM_HDR) $(CNE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXECNE6SSM_SRC) $(CNE6SSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(CNE6SSM_MK) $(CNE6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(CNE6SSM_BETAS_MK) $(CNE6SSM_INSTALL_DIR)
ifneq ($(CNE6SSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(CNE6SSM_SLHA_INPUT) $(CNE6SSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(CNE6SSM_GNUPLOT) $(CNE6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBCNE6SSM_DEP)
		-rm -f $(EXECNE6SSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCNE6SSM_OBJ)
		-rm -f $(EXECNE6SSM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBCNE6SSM)
		-rm -f $(RUN_CNE6SSM_EXE)
		-rm -f $(RUN_CMD_LINE_CNE6SSM_EXE)
		-rm -f $(SCAN_CNE6SSM_EXE)
		-rm -f $(GRIDSCAN_CNE6SSM_EXE)
		-rm -f $(GRIDSCAN_RGE_CNE6SSM_EXE)
		-rm -f $(RGE_COEFF_CNE6SSM_EXE)
		-rm -f $(RUN_SEMI_CNE6SSM_EXE)
		-rm -f $(RUN_SUSYHD_CNE6SSM_EXE)
		-rm -f $(GRIDSCAN_SEMI_CNE6SSM_EXE)
		-rm -f $(SCAN_SEMI_CNE6SSM_EXE)
		-rm -f $(GET_HIGGS_MASS_EXE)
		-rm -f $(GENERATE_POINTS_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(CNE6SSM_TARBALL) \
		$(LIBCNE6SSM_SRC) $(LIBCNE6SSM_HDR) \
		$(EXECNE6SSM_SRC) \
		$(CNE6SSM_MK) $(CNE6SSM_BETAS_MK) \
		$(CNE6SSM_SLHA_INPUT) $(CNE6SSM_GNUPLOT)

$(LIBCNE6SSM_SRC) $(LIBCNE6SSM_HDR) $(EXECNE6SSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CNE6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CNE6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CNE6SSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_CNE6SSM)"
		@echo "Note: to regenerate CNE6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CNE6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CNE6SSM):
		@true
endif

$(LIBCNE6SSM_DEP) $(EXECNE6SSM_DEP) $(LIBCNE6SSM_OBJ) $(EXECNE6SSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(MLINKFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCNE6SSM_DEP) $(EXECNE6SSM_DEP) $(LIBCNE6SSM_OBJ) $(EXECNE6SSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_STATIC_LIBS),yes)
$(LIBCNE6SSM): $(LIBCNE6SSM_OBJ)
		$(MAKELIB) $@ $^
else
$(LIBCNE6SSM): $(LIBCNE6SSM_OBJ)
		$(MAKELIB) $@ $^ $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)
endif

$(GRIDSCAN_CNE6SSM_EXE): $(GRIDSCAN_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(GRIDSCAN_RGE_CNE6SSM_EXE): $(GRIDSCAN_RGE_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RGE_COEFF_CNE6SSM_EXE): $(RGE_COEFF_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_CNE6SSM_EXE): $(RUN_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_CMD_LINE_CNE6SSM_EXE): $(RUN_CMD_LINE_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(SCAN_CNE6SSM_EXE): $(SCAN_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_SEMI_CNE6SSM_EXE): $(RUN_SEMI_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(RUN_SUSYHD_CNE6SSM_EXE): $(RUN_SUSYHD_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBMATHLINK) $(LIBSUSYHDLINK)
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(MLINKLIBS) $(EXTRA_MLINK_LIBS)

$(GRIDSCAN_SEMI_CNE6SSM_EXE): $(GRIDSCAN_SEMI_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(SCAN_SEMI_CNE6SSM_EXE): $(SCAN_SEMI_CNE6SSM_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(GET_HIGGS_MASS_EXE): $(GET_HIGGS_MASS_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS)) $(LIBMATHLINK) $(LIBSUSYHDLINK)
		$(CXX) -Wl,-no-as-needed -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(MLINKLIBS) $(EXTRA_MLINK_LIBS)

$(GENERATE_POINTS_EXE): $(GENERATE_POINTS_OBJ) $(LIBCNE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

ALLDEP += $(LIBCNE6SSM_DEP) $(EXECNE6SSM_DEP)
ALLSRC += $(LIBCNE6SSM_SRC) $(EXECNE6SSM_SRC)
ALLLIB += $(LIBCNE6SSM)
ifneq ($(findstring two_scale,$(ALGORITHMS)),)
ALLEXE += $(GRIDSCAN_CNE6SSM_EXE) $(GRIDSCAN_RGE_CNE6SSM_EXE) $(RGE_COEFF_CNE6SSM_EXE) $(RUN_CNE6SSM_EXE) $(RUN_CMD_LINE_CNE6SSM_EXE) $(SCAN_CNE6SSM_EXE)
endif
ifneq ($(findstring semianalytic,$(ALGORITHMS)),)
ALLEXE += $(RUN_SEMI_CNE6SSM_EXE) $(GRIDSCAN_SEMI_CNE6SSM_EXE) $(SCAN_SEMI_CNE6SSM_EXE) $(GENERATE_POINTS_EXE)
ifneq ($(findstring addons/susyhd_call,$(ADDONS)),)
ALLEXE += $(RUN_SUSYHD_CNE6SSM_EXE) $(GET_HIGGS_MASS_EXE)
endif
endif
