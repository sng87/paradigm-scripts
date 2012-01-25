THISDIR = $(CURDIR)
USERNAME ?= $(LOGNAME)
PUBLICHOST ?= cancer2.cse.ucsc.edu
FIREWALLHOST ?= tcga1

init.sh : jobForest exe
	echo \
	export PATH=$(THISDIR)/bin:$(THISDIR)/exe:\$${PATH} > init.sh
	echo \
	export PYTHONPATH=$(THISDIR)/jobForest:$(THISDIR)/bin:\$${PYTHONPATH} >> init.sh
	echo \
	setenv PATH $(THISDIR)/bin:$(THISDIR)/exe:\$${PATH} > init.csh
	echo \
	setenv PYTHONPATH $(THISDIR)/jobForest:$(THISDIR)/bin:\$${PYTHONPATH} >> init.csh

jobForest :
	git clone git://github.com/kellrott/jobForest.git

exe :
	mkdir exe
	if (! test -d /inside); \
	then cd exe; scp $(USERNAME)@$(PUBLICHOST):{/hive/groups/cancerGB/paradigm/exe/newEmSpec/paradigm,/hive/groups/cancerGB/paradigm/current/Paradigm/ParadigmDistribution/collectParameters} .; \
	else cd exe; scp $(USERNAME)@$(FIREWALLHOST):{/inside/grotto/users/sng/bin/Paradigm/paradigm,/inside/grotto/users/sng/bin/Paradigm/collectParameters} .; \
	fi

clean :
	rm -rf jobForest exe init.sh init.csh
