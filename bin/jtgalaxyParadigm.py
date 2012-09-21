#!/usr/bin/env python
"""
jtgalaxyParadigm.py: handles setup and running of paradigm on multiple cohorts and/or pathways
"""
## Written by: Sam Ng
import getopt, os, os.path, re, sys
from optparse import OptionParser
from jtParadigm import *

from jobTree.src.bioio import logger
from jobTree.src.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

basedir = os.path.dirname(os.path.abspath(__file__))

basedogma = os.path.join(basedir, "d_standard.zip")
basepathway = os.path.join(basedir, "p_global_five3_v2.zip")

paradigmExec = os.path.join(basedir, "paradigm")
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
    parser = OptionParser(usage = "%prog [options] attachment file:path [attachment file:path ...]")
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help = "Add as a child of jobFile rather " +
                      "than making a new jobTree")
    parser.add_option("-w", "--workdir", dest="workdir", help="Common Work directory", default="./")
    parser.add_option("-n", "--nulls", dest="nulls", help="Number of Null Samples", default="5")
    parser.add_option("-d", "--dogma", dest="dogmazip", help="Path to PARADIGM Dogma Specification", default=basedogma)
    parser.add_option("-p", "--pathway", dest="pathwayzip", help="Path to PARADIGM Pathway Specification", default=basepathway)
    parser.add_option("-b", "--boundaries", dest="disc", help="Data Discretization Bounds", default="0.33;0.67")
    parser.add_option("-t", "--storedparam", dest="param", help="Initial Parameter Starting Point", default=None)
    parser.add_option("-s", "--skipem", action="store_false", dest="em", help="Skip Running EM", default=True)
    
    parser.add_option("--fr", "--filter-real", dest="filtered_real", help="Filtered Output", default=None)
    parser.add_option("--fa", "--filter-all", dest="filtered_all", help="Filtered Output", default=None)
    parser.add_option("--ur", "--unfilter-real", dest="unfiltered_real", help="Filtered Output", default=None)
    parser.add_option("--ua", "--unfilter-all", dest="unfiltered_all", help="Filtered Output", default=None)
    
    options, args = parser.parse_args()
    logger.info("options: " + str(options))
    print "Using Batch System '" + options.batchSystem + "'"
    
    evidList = args 
    if (len(evidList) % 2 == 1) | (len(evidList) == 0):
        sys.stderr.write("ERROR: incorrect number of arguments\n")
        sys.exit(1)
    
    workdir = os.path.abspath(options.workdir)
    nullBatches = int(options.nullBatches)
    dogmaZip=os.path.abspath(options.dogmazip)
    pathwayZip=os.path.abspath(options.pathwayzip)
    disc=options.disc
    paramFile=os.path.abspath(options.param) if args.param is not None else None
    runEM = options.runEM
    
    dogmaLib = os.path.join(workdir, "dogma")
    pathwayLib = os.path.join(workdir, "pathway")
    system("unzip %s -d %s" % (self.dogmaZip, dogmaLib))
    system("unzip %s -d %s" % (self.pathwayZip, pathwayLib))

    ## run
    logger.info("starting prepare")
    s = Stack(prepareParadigm(" ".join(evidList), disc, paramFile, nullBatches, paradigmExec, inferSpec, dogmaLib, pathwayLib, runEM, workdir))
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
    from jtgalaxyParadigm import *
    if os.path.exists(".jobTree"):
        print "WARNING: .jobTree directory already exists"
    wrapParadigm()
