DIR		:= config
MODNAME		:= config

CONFIG_HDR := \
		$(DIR)/config.h

CONFIG_MK    := \
		$(DIR)/module.mk

CONFIG_TMPL  := \
		$(DIR)/abspathx.mk \
		$(DIR)/config.h.in \
		$(DIR)/ne6ssm-config.in \

CONFIG_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MAKEFILE_IN  := \
		$(DIR)/Makefile.in

REMOVE_EXPORT_MARKERS := \
		$(DIR)/remove_export_markers.sh

.PHONY:		all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(CONFIG_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(CONFIG_TMPL) $(CONFIG_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(CONFIG_MK) $(CONFIG_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MAKEFILE_IN) $(CONFIG_INSTALL_DIR) -m u=rw,g=r,o=r
		$(INSTALL_STRIPPED) $(REMOVE_EXPORT_MARKERS) $(CONFIG_INSTALL_DIR) -m u=rwx,g=r,o=r
		$(INSTALL_STRIPPED) $(INSTALL_STRIPPED) $(CONFIG_INSTALL_DIR) -m u=rwx,g=r,o=r
endif

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)
		-rm -f $(CONFIG_HDR)

clean::		clean-$(MODNAME)

distclean::	distclean-$(MODNAME)
