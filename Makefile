THISDIR = $(CURDIR)

init.sh : jobTree exe
	echo \
	export PATH=$(THISDIR)/bin:$(THISDIR)/exe:\$${PATH} > init.sh
	echo \
	export PYTHONPATH=$(THISDIR):$(THISDIR)/bin:\$${PYTHONPATH} >> init.sh
	echo \
	setenv PATH $(THISDIR)/bin:$(THISDIR)/exe:\$${PATH} > init.csh
	echo \
	setenv PYTHONPATH $(THISDIR):$(THISDIR)/bin:\$${PYTHONPATH} >> init.csh

jobTree : sonLib
	git clone git://github.com/benedictpaten/jobTree.git
	cd jobTree; make

sonLib :
	git clone git://github.com/benedictpaten/sonLib.git

exe :
	mkdir exe
	if (test -d /hive); \
	then \
	cd exe; cp /hive/groups/cancerGB/paradigm/exe/newEmSpec/paradigm /hive/groups/cancerGB/paradigm/current/Paradigm/ParadigmDistribution/collectParameters .; \
	fi
	if (test -d /inside); \
	then \
	cd exe; cp /inside/grotto/users/sng/bin/Paradigm/paradigm /inside/grotto/users/sng/bin/Paradigm/collectParameters .; \
	fi
	if (! test -e exe/paradigm); \
	then \
	cd exe; cp ../public/exe/paradigm .; \
	fi

pathmark-scripts :
	git clone git://github.com/sng87/pathmark-scripts.git; \
	cd pathmark-scripts; make

clean :
	rm -rf pathmark-scripts jobTree sonLib exe init.sh init.csh
	cd test; make clean
