THISDIR = $(CURDIR)
THISOS = $(shell uname -s)

init.sh : jobTree exe
	echo \
	export PATH=$(THISDIR)/bin:$(THISDIR)/utilities:\$${PATH} > init.sh
	echo \
	if [ -n "\$${PYTHONPATH+x}" ] >> init.sh
	echo \
	then >> init.sh
	echo \
	  export PYTHONPATH=$(THISDIR):$(THISDIR)/bin:\$${PYTHONPATH} >> init.sh
	echo \
	else >> init.sh
	echo \
	  export PYTHONPATH=$(THISDIR):$(THISDIR)/bin >> init.sh
	echo \
	fi >> init.sh
	echo \
	setenv PATH $(THISDIR)/bin:$(THISDIR)/utilities:\$${PATH} > init.csh
	echo \
	if \$$?PYTHONPATH then >> init.csh
	echo \
	  setenv PYTHONPATH $(THISDIR):$(THISDIR)/bin:\$${PYTHONPATH} >> init.csh
	echo \
	else >> init.csh
	echo \
	  setenv PYTHONPATH $(THISDIR):$(THISDIR)/bin >> init.csh
	echo \
	endif >> init.csh

jobTree : sonLib
	git clone git://github.com/benedictpaten/jobTree.git
	cd jobTree; make

sonLib :
	git clone git://github.com/benedictpaten/sonLib.git

exe :
	mkdir exe
	if (test -d /hive); then \
	cd exe; cp /hive/groups/cancerGB/paradigm/exe/newEmSpec/paradigm /hive/groups/cancerGB/paradigm/current/Paradigm/ParadigmDistribution/collectParameters .; \
	fi
	if (test -d /inside); then \
	cd exe; cp /inside/grotto/users/sng/bin/Paradigm/paradigm /inside/grotto/users/sng/bin/Paradigm/collectParameters .; \
	fi
	if (! test -e exe/paradigm); then \
	if [ $(THISOS) == Darwin ]; then \
	cd exe; cp ../public/exe/collectParameters ../public/exe/MACOSX/paradigm .; \
	elif [ $(THISOS) == Linux ]; then \
	cd exe; cp ../public/exe/collectParameters ../public/exe/LINUX/paradigm .; \
	else \
	echo "paradigm not compiled for os $(THISOS)"; \
	fi \
	fi
	ln -s $(THISDIR)/exe/paradigm bin/
	ln -s $(THISDIR)/exe/collectParameters bin/

pathmark-scripts :
	cd ..; git clone git://github.com/sng87/pathmark-scripts.git
	ln -s ../pathmark-scripts pathmark-scripts
	cd pathmark-scripts; make

galaxy : pathmark-scripts
	mkdir -p paradigm_module
	cp bin/* paradigm_module/
	cp -r pathmark-scripts/bin/* paradigm_module/

clean :
	rm -rf bin/paradigm bin/collectParameters pathmark-scripts jobTree sonLib exe init.sh init.csh
	cd test; make clean
