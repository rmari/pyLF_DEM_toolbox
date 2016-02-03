# load Makefile inputs
makeconfig = config/Makefile_config.mk
gotconfig = 0
ifneq ("$(wildcard $(makeconfig))","")
        include $(makeconfig)
        gotconfig = 1
endif

install:
	cp *.py $(INSTALL_DIR)
