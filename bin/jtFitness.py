#!/usr/bin/env python
"""jtFitness.py:

Usage:
  jtFitness.py [options] attachment file:path [attachment file:path ...]

Options:
   -d dir               dogma directory (take top of config from here)
   -p dir               the directory containing pathway files 
   -b flt;flt[,flt;flt] boundaries for discretization, use comma to specify different
                        boundaries per evidence (default 0.333;0.667)
"""
## Written By: Sam Ng
## Last Updated: ##/##/####
import math, os, os.path, sys, random, re, types

from wrapParadigm import prepareParadigm
from wrapParadigm import jtParadigm

from optparse import OptionParser
from copy import deepcopy

from jobTree.src.bioio import logger
from jobTree.src.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

## default variables
nDecoys = 1

## executables and directories
binDir = os.path.realpath(os.path.dirname(sys.argv[0]))
dogmaDir = "%s/dogmas" % (binDir.rstrip("/bin"))
pathwayDir = "%s/pathways" % (binDir.rstrip("/bin"))
exeDir = "%s/exe" % (binDir.rstrip("/bin"))

dogmaDefault = "standard"
pathwayDefault = "global_five3_v2"

paradigmExec = "%s/paradigm" % (exeDir)
inferSpec = "method=BP,updates=SEQFIX,tol=1e-9,maxiter=10000,logdomain=0"

def openAnyFile(inf):
    """performs an open() on a file, url, or stream"""
    if type(inf) == types.StringType:
        if inf.startswith("http"):
            import urllib2
            stream = urllib2.urlopen(inf)
        elif os.path.exists(inf):
            stream = open(inf, 'r')
        else:
            log("Could not open %s no file, url, or stream" % (inf), die = True)
    else:
        stream = inf
    return stream

def retColumns(inf, delim = "\t"):
    """returns the columns of a .tsv"""
    f = openAnyFile(inf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank header\n", die = True)
    line = line.rstrip("\r\n")
    return(re.split(delim, line)[1:])

def rwCRSData(outf, inf, delim = "\t", null = "NA", useCols = None, useRows = None, colMap = {}, rowMap = {}, enumerateRows = False):
    """reads and writes .tsv"""
    colFeatures = []
    ## read header
    f = openAnyFile(inf)
    o = open(outf, "w")
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    lineLength = len(pline)
    colIndex = {}
    for i, col in enumerate(pline[1:]):
        colIndex[col] = i
        if useCols != None:
            if col not in useCols:
                continue
        colFeatures.append(col)
    o.write("id%s\n" % (delim+delim.join(colFeatures)))
    ## read and write data
    rowCount = 0
    if enumerateRows:
        m = open("%s.idmap" % (outf))
    for line in f:
        if line.isspace():
            continue
        rowCount += 1
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        if pline[0] in rowMap:
            mrow = rowMap[pline[0]]
        else:
            mrow = pline[0]
        if useRows != None:
            if mrow not in useRows:
                continue
        if len(pline) != lineLength:
            log("ERROR: length of line does not match the rest of the file\n", die = True)
        if enumerateRows:
            o.write("r%s" % (rowCount))
            m.write("r%s\t%s\n" % (rowCount, mrow))
        else:
            o.write("%s" % (mrow))
        for col in colFeatures:
            if col in colMap:
                mcol = colMap[col]
            else:
                mcol = col
            if pline[colIndex[mcol]+1] == "":
                o.write("%s" % (delim+null))
            else:            
                o.write("%s" % (delim+pline[colIndex[mcol]+1]))
        o.write("\n")
    f.close()
    o.close()
    if enumerateRows:
        m.close()

class jtCmd(Target):
    def __init__(self, cmd, directory):
        Target.__init__(self, time=1000)
        self.cmd = cmd
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        system(self.cmd)

class branchDecoys(Target):
    def __init__(self, evidSpec, disc, paradigmExec, inferSpec, dogmaLib, pathwayLib, nDecoys, directory):
        Target.__init__(self, time=10000)
        self.evidSpec = evidSpec
        self.disc = disc
        self.paradigmExec = paradigmExec
        self.inferSpec = inferSpec
        self.dogmaLib = dogmaLib
        self.pathwayLib = pathwayLib
        self.nDecoys = nDecoys
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        evidList = zip(re.split("\s", self.evidSpec)[0::2], re.split("\s", self.evidSpec)[1::2])
        
        ## assert files are in data/
        for i in evidList:
            assert(re.split(":", i[1])[1].startswith("data"))
        ## copy files over to internal directory
        system("mkdir paradigm-real")
        system("mkdir paradigm-real/data")
        for i in evidList:
            system("cp %s paradigm-real/data" % (re.split(":", i[1])[1]))
        ## kick off real run
        # self.addChildTarget(prepareParadigm(self.evidSpec, self.disc, 0, self.paradigmExec, self.inferSpec, self.dogmaLib, self.pathwayLib, "%s/paradigm-real" % (self.directory)))
        
        ## iterate through decoys
        for decoy in range(1, self.nDecoys+1):
            system("mkdir paradigm-decoy%s" % (decoy))
            system("mkdir paradigm-decoy%s/data" % (decoy))
            ## decoys will consist of shuffled methylation data
            for i in evidList:
                if i[0] != "promoter":
                    system("cp %s paradigm-decoy%s/data" % (re.split(":", i[1])[1], decoy))
                else:
                    if i[1].startswith("rawFile"):
                        ## permute cols to directory
                        colFeatures = retColumns(re.split(":", i[1])[1])
                        permuteFeatures = random.sample(colFeatures, len(colFeatures))
                        colMap = {}
                        for j in range(0, len(colFeatures)):
                            colMap[colFeatures[j]] = permuteFeatures[j]
                        rwCRSData("paradigm-decoy%s/%s" % (decoy, re.split(":", i[1])[1]), re.split(":", i[1])[1], colMap = colMap)
                    else:
                        ## permute rows to directory
                        rowFeatures = retRows(re.split(":", i[1])[1])
                        permuteFeatures = random.sample(rowFeatures, len(rowFeatures))
                        rowMap = {}
                        for j in range(0, len(rowFeatures)):
                            rowMap[rowFeatures[j]] = permuteFeatures[j]
                        rwCRSData("paradigm-decoy%s/%s" % (decoy, re.split(":", i[1])[1]), re.split(":", i[1])[1], rowMap = rowMap)
            ## kick off decoy run
            self.addChildTarget(prepareParadigm(self.evidSpec, self.disc, 0, self.paradigmExec, self.inferSpec, self.dogmaLib, self.pathwayLib, "%s/paradigm-decoy%s" % (self.directory, decoy)))
            
        ## no repeating or folds yet

        
class branchRepeats(Target):
    def __init__(self, nRepeats, directory):
        Target.__init__(self, time=10000)
        self.nRepeats = nRepeats
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        for i in range(1, self.nRepeats+1):
            ## create partitioning
            self.addChildTarget(branchFolds())

class branchFolds(Target):
    def __init__(self, directory):
        Target.__init__(self, time=10000)
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        ## prepareParadigm
        ## addChildTarget(jtTrain())
        ## addFollowOnTarget(jtTest())


def jtFitness():
    ## parse arguments
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help = "Add as a child of jobFile rather " +
                      "than making a new jobTree")
    parser.add_option("-d", "--dogma", dest="dogmaPath", default="")
    parser.add_option("-p", "--pathway", dest="pathwayPath", default="")
    parser.add_option("-b", "--boundaries", dest="discBound", default="")
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
    
    ## clean
    if len(args) == 1:
        if args[0] == "clean":
            print "rm -rf .jobTree paradigm-*"
            os.system("rm -rf .jobTree paradigm-*")
            sys.exit(0)
    
    ## run
    logger.info("options: " + str(options))
    s = Stack(branchDecoys(" ".join(evidList), disc, paradigmExec, inferSpec, dogma, pathway, nDecoys, os.getcwd()))
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
    from jtFitness import *
    if os.path.exists(".jobTree"):
        print "WARNING: .jobTree directory already exists"
    jtFitness()
