# load Makefile inputs

PYTHON_DIR = $(shell env python -c "import site;print(site.getsitepackages()[0])")

install:
	cp *.py $(INSTALL)
	cp LFDEM_confgen.py $(PYTHON_DIR)
