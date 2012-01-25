#!/usr/bin/env python
"""wrapParadigm.py: 

Usage:
  wrapParadigm.py [options] attachment file:path [attachment file:path ...]
"""
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
    def __init__(self, evidSpec, paradigmExec, dogmaLib, pathwayLib, directory):
        Target.__init__(self, time=10000)
        self.evidSpec = evidSpec
        self.paradigmExec = paradigmExec
        self.dogmaLib = dogmaLib
        self.pathwayLib = pathwayLib
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        cmd = "prepareParadigm.py -b \"0.3333;0.6667\" -s same -n 0 -i %s -e %s -d %s -p %s %s" % (inferSpec, self.paradigmExec, self.dogmaLib, self.pathwayLib, self.evidSpec)
        system(cmd)
        self.setFollowOnTarget(jtParadigm(self.directory))

class jtParadigm(Target):
    def __init__(self, directory):
        Target.__init__(self, time=10000)
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        self.addChildTarget(ExpectationIteration(0, 0.001, "%s" % (self.directory)))

def wrapParadigm():
    ## parse arguments
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help = "Add as a child of jobFile rather " +
                      "than making a new jobTree")
    parser.add_option("-d", "--dogma", dest="dogmaPath", default="")
    parser.add_option("-p", "--pathway", dest="pathwayPath", default="")
    options, args = parser.parse_args()
    print "Using Batch System '" + options.batchSystem + "'"
   
    evidList = args 
    if (len(evidList) % 2 == 1) | (len(evidList) == 0):
        sys.stderr.write("ERROR: incorrect number of arguments\n")
        sys.exit(1)
    
    if len(options.dogmaPath) == 0:
        dogma = "%s/%s" % (dogmaDir, dogmaDefault)
    else:
        dogma = options.dogmaPath
    if len(options.pathwayPath) == 0:
        pathway = "%s/%s" % (pathwayDir, pathwayDefault)
    else:
        pathway = options.pathwayPath
    
    logger.info("options: " + str(options))
    
    ## run
    logger.info("starting first iteration")
    s = Stack(prepareParadigm(" ".join(evidList), paradigmExec, dogma, pathway, os.getcwd()))
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
