#!/usr/bin/env python
"""wrapParadigm.py: handles setup and running of paradigm on multiple cohorts and/or pathways

Usage:
  wrapParadigm.py [options] attachment file:path [attachment file:path ...]

Options:
   -d dir               dogma directory (take top of config from here)
   -p dir               the directory containing pathway files 
   -b flt;flt[,flt;flt] boundaries for discretization, use comma to specify different
                        boundaries per evidence (default 0.333;0.667)
"""
## Written by: Sam Ng
import getopt, os, os.path, re, sys
from optparse import OptionParser
from jtParadigm import *

from jobTree.src.bioio import logger
from jobTree.src.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

binDir = os.path.realpath(os.path.dirname(sys.argv[0]))
dogmaDir = "%s/dogmas" % (binDir.rstrip("/bin"))
pathwayDir = "%s/pathways" % (binDir.rstrip("/bin"))
exeDir = "%s/exe" % (binDir.rstrip("/bin"))

dogmaDefault = "standard"
pathwayDefault = "global_five3_v2"

paradigmExec = "%s/paradigm" % (exeDir)
inferSpec = "method=BP,updates=SEQFIX,tol=1e-9,maxiter=10000,logdomain=0"

class prepareParadigm(Target):
    def __init__(self, evidSpec, disc, paramFile, nullBatches, paradigmExec, inferSpec, dogmaLib, pathwayLib, em, directory):
        Target.__init__(self, time=10000)
        self.evidSpec = evidSpec
        self.disc = disc
        self.paramFile = paramFile
        self.nullBatches = nullBatches
        self.paradigmExec = paradigmExec
        self.inferSpec = inferSpec
        self.dogmaLib = dogmaLib
        self.pathwayLib = pathwayLib
        self.em = em
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        if self.paramFile is not None:
            cmd = "prepareParadigm.py -b \"%s\" -t %s -s same -n %s -i %s -e %s -d %s -p %s %s >& jt.err" % (self.disc, self.paramFile, self.nullBatches, self.inferSpec, self.paradigmExec, self.dogmaLib, self.pathwayLib, self.evidSpec)
        else:
            cmd = "prepareParadigm.py -b \"%s\" -s same -n %s -i %s -e %s -d %s -p %s %s >& jt.err" % (self.disc, self.nullBatches, self.inferSpec, self.paradigmExec, self.dogmaLib, self.pathwayLib, self.evidSpec)
        system(cmd)
        self.setFollowOnTarget(jtParadigm(self.em, self.directory))

class jtParadigm(Target):
    def __init__(self, em, directory):
        Target.__init__(self, time=10000)
        self.em = em
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        if self.em:
            self.addChildTarget(ExpectationIteration(0, 0.001, "%s" % (self.directory)))
        else:
            self.addChildTarget(FinalRun(0, "%s" % (self.directory)))

def wrapParadigm():
    ## parse arguments
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help = "Add as a child of jobFile rather " +
                      "than making a new jobTree")
    parser.add_option("-d", "--dogma", dest="dogmaPath", default="")
    parser.add_option("-p", "--pathway", dest="pathwayPath", default="")
    parser.add_option("-b", "--boundaries", dest="discBound", default="")
    parser.add_option("-n", "--nulls", dest="nullBatches", default="")
    parser.add_option("-t", "--storedparam", dest="paramFile", default="")
    parser.add_option("-s", "--skipem", action="store_false", dest="runEM", default=True)
    options, args = parser.parse_args()
    print "Using Batch System '" + options.batchSystem + "'"
   
    evidList = args 
    if (len(evidList) % 2 == 1) | (len(evidList) == 0):
        sys.stderr.write("ERROR: incorrect number of arguments\n")
        sys.exit(1)
    
    if len(options.discBound) == 0:
        disc = "0.3333;0.6667"
    else:
        disc = options.discBound
    if len(options.dogmaPath) == 0:
        dogma = "%s/%s" % (dogmaDir, dogmaDefault)
    else:
        dogma = options.dogmaPath
    if len(options.pathwayPath) == 0:
        pathway = "%s/%s" % (pathwayDir, pathwayDefault)
    else:
        pathway = options.pathwayPath
    if len(options.nullBatches) == 0:
        nullBatches = 0
    else:
        nullBatches = int(options.nullBatches)
    if len(options.paramFile) == 0:
        paramFile = None
    else:
        paramFile = options.paramFile
    runEM = options.runEM
    logger.info("options: " + str(options))
    
    ## run
    logger.info("starting prepare")
    s = Stack(prepareParadigm(" ".join(evidList), disc, paramFile, nullBatches, paradigmExec, inferSpec, dogma, pathway, runEM, os.getcwd()))
    if options.jobFile:
        s.addToJobFile(options.jobFile)
    else:
        if options.jobTree == None:
            options.jobTree = "./.jobTree"
        
        failed = s.startJobTree(options)
        if failed:
            print ("%d jobs failed" % failed)
        else:
            logger.info("Run complete!")
            system("rm -rf .lastjobTree")
            system("mv .jobTree .lastjobTree")

if __name__ == "__main__":
    from wrapParadigm import *
    if os.path.exists(".jobTree"):
        print "WARNING: .jobTree directory already exists"
    wrapParadigm()
