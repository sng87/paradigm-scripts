#!/usr/bin/env python
"""jtFitness.py: handles setup and running paradigm in a cross validation framework for calculation of fitness

Usage:
  jtFitness.py [options] attachment file:path [attachment file:path ...]

Options:
   -d dir               dogma directory (take top of config from here)
   -p dir               the directory containing pathway files 
   -b flt;flt[,flt;flt] boundaries for discretization, use comma to specify different
                        boundaries per evidence (default 0.333;0.667)
"""
## Written By: Sam Ng
import math, os, os.path, sys, random, re, types

from wrapParadigm import prepareParadigm
from wrapParadigm import jtParadigm
from jtParadigm import *

from optparse import OptionParser
from copy import deepcopy

from jobTree.src.bioio import logger
from jobTree.src.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

## default variables
mFolds = 5

## executables and directories
binDir = os.path.realpath(os.path.dirname(sys.argv[0]))
dogmaDir = "%s/dogmas" % (binDir.rstrip("/bin"))
pathwayDir = "%s/pathways" % (binDir.rstrip("/bin"))
exeDir = "%s/exe" % (binDir.rstrip("/bin"))

dogmaDefault = "standard"
pathwayDefault = "global_five3_v2"

paradigmExec = "%s/paradigm" % (exeDir)
inferSpec = "method=BP,updates=SEQFIX,tol=1e-9,maxiter=10000,logdomain=0"

def log(msg, die = False):
    """logger function"""
    sys.stderr.write(msg)
    if die:
        sys.exit(1)

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

def retRows(inf, delim = "\t"):
    """returns the rows of a .tsv"""
    rows = []
    f = openAnyFile(inf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank header\n", die = True)
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        rows.append(re.split(delim, line)[0])
    return(rows)

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

class branchFolds(Target):
    def __init__(self, evidSpec, disc, paramFile, paradigmExec, inferSpec, dogmaLib, pathwayLib, shuffleNode, nShuffle, mFolds, directory):
        Target.__init__(self, time=10000)
        self.evidSpec = evidSpec
        self.disc = disc
        self.paramFile = paramFile
        self.paradigmExec = paradigmExec
        self.inferSpec = inferSpec
        self.dogmaLib = dogmaLib
        self.pathwayLib = pathwayLib
        self.shuffleNode = shuffleNode
        self.nShuffle = nShuffle
        self.mFolds = mFolds
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        evidList = zip(re.split("\s", self.evidSpec)[0::2], re.split("\s", self.evidSpec)[1::2])
        
        ## assert files are in data/
        for i in evidList:
            assert(re.split(":", i[1])[1].startswith("data"))
        
        ## check if new run
        if not os.path.exists("fold1"):            
            ## find sample overlap
            dataSamples = None
            for i in evidList:
                if i[1].startswith("rawFile"):
                    if dataSamples is None:
                        dataSamples = set(retRows(re.split(":", i[1])[1]))
                    else:
                        dataSamples = dataSamples & set(retRows(re.split(":", i[1])[1]))
                else:
                    if dataSamples is None:
                        dataSamples = set(retColumns(re.split(":", i[1])[1]))
                    else:
                        dataSamples = dataSamples & set(retColumns(re.split(":", i[1])[1]))
            dataSamples = list(dataSamples)
            
            ## pick samples
            foldSamples = {}
            for f in range(1, self.mFolds+1):
                foldSamples[f] = []
            selectSamples = deepcopy(dataSamples)
            while len(selectSamples) > 0:
                for f in range(1, mFolds+1):
                    if len(selectSamples) > 0:
                        foldSamples[f].append(selectSamples.pop(random.randint(0,len(selectSamples)-1)))
        
            ## create directories and data
            for f in range(1, self.mFolds+1):
                system("mkdir fold%s" % (f))
                system("mkdir fold%s/train" % (f))
                system("mkdir fold%s/train/data" % (f))
                system("mkdir fold%s/test" % (f))
                system("mkdir fold%s/test/data" % (f))
                trainSamples = list(set(dataSamples) - set(foldSamples[f]))
                testSamples = foldSamples[f]
                for i in evidList:
                    if i[1].startswith("rawFile"):
                        rwCRSData("fold%s/train/%s" % (f, re.split(":", i[1])[1]), re.split(":", i[1])[1], useRows = trainSamples)
                        rwCRSData("fold%s/test/%s" % (f, re.split(":", i[1])[1]), re.split(":", i[1])[1], useRows = testSamples)
                    else:
                        rwCRSData("fold%s/train/%s" % (f, re.split(":", i[1])[1]), re.split(":", i[1])[1], useCols = trainSamples)
                        rwCRSData("fold%s/test/%s" % (f, re.split(":", i[1])[1]), re.split(":", i[1])[1], useCols = testSamples)
        
        ## kick off runs
        for f in range(1, self.mFolds+1):
            self.addChildTarget(branchTrain(self.evidSpec, self.disc, self.paramFile, self.paradigmExec, self.inferSpec, self.dogmaLib, self.pathwayLib, self.shuffleNode, self.nShuffle, "%s/fold%s" % (self.directory, f)))

class branchTrain(Target):
    def __init__(self, evidSpec, disc, paramFile, paradigmExec, inferSpec, dogmaLib, pathwayLib, shuffleNode, nShuffle, directory):
        Target.__init__(self, time=10000)
        self.evidSpec = evidSpec
        self.disc = disc
        self.paramFile = paramFile
        self.paradigmExec = paradigmExec
        self.inferSpec = inferSpec
        self.dogmaLib = dogmaLib
        self.pathwayLib = pathwayLib
        self.shuffleNode = shuffleNode
        self.nShuffle = nShuffle
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        evidList = zip(re.split("\s", self.evidSpec)[0::2], re.split("\s", self.evidSpec)[1::2])
        
        ## check if new run
        if not os.path.exists("train/real"):
            ## copy files over to internal directory
            system("mkdir train/real")
            system("mkdir train/real/data")
            system("mkdir test/real")
            system("mkdir test/real/data")
            for i in evidList:
                system("cp train/%s train/real/data" % (re.split(":", i[1])[1]))
                system("cp test/%s test/real/data" % (re.split(":", i[1])[1]))
        
        ## iterate over shuffles
        for shuffle in range(1, self.nShuffle+1):
            ## check if new run
            if not os.path.exists("train/shuffle%s" % (shuffle)):
                ## create shuffle files
                system("mkdir train/shuffle%s" % (shuffle))
                system("mkdir train/shuffle%s/data" % (shuffle))
                system("mkdir test/shuffle%s" % (shuffle))
                system("mkdir test/shuffle%s/data" % (shuffle))
                ## shuffles will consist of shuffled methylation data
                for i in evidList:
                    if i[0] != self.shuffleNode:
                        system("cp train/%s train/shuffle%s/data" % (re.split(":", i[1])[1], shuffle))
                        system("cp test/%s test/shuffle%s/data" % (re.split(":", i[1])[1], shuffle))
                    else:
                        if i[1].startswith("rawFile"):
                            ## permute cols to directory
                            colFeatures = retColumns("train/%s" % (re.split(":", i[1])[1]))
                            permuteFeatures = random.sample(colFeatures, len(colFeatures))
                            colMap = {}
                            for j in range(0, len(colFeatures)):
                                colMap[colFeatures[j]] = permuteFeatures[j]
                            rwCRSData("train/shuffle%s/%s" % (shuffle, re.split(":", i[1])[1]), "train/%s" % (re.split(":", i[1])[1]), colMap = colMap)
                            rwCRSData("test/shuffle%s/%s" % (shuffle, re.split(":", i[1])[1]), "test/%s" % (re.split(":", i[1])[1]), colMap = colMap)
                        else:
                            ## permute rows to directory
                            rowFeatures = retRows("train/%s" % (re.split(":", i[1])[1]))
                            permuteFeatures = random.sample(rowFeatures, len(rowFeatures))
                            rowMap = {}
                            for j in range(0, len(rowFeatures)):
                                rowMap[rowFeatures[j]] = permuteFeatures[j]
                            rwCRSData("train/shuffle%s/%s" % (shuffle, re.split(":", i[1])[1]), "train/%s" % (re.split(":", i[1])[1]), rowMap = rowMap)
                            rwCRSData("test/shuffle%s/%s" % (shuffle, re.split(":", i[1])[1]), "test/%s" % (re.split(":", i[1])[1]), rowMap = rowMap)
        
        ## kick off train.real run
        if not os.path.exists("train/real/outputFiles"):
            if os.path.exists("%s/real" % (self.pathwayLib)):
                self.addChildTarget(prepareParadigm(self.evidSpec, self.disc, self.paramFile, 0, self.paradigmExec, self.inferSpec, self.dogmaLib, "%s/real" % (self.pathwayLib), "%s/train/real" % (self.directory)))
            else:
                self.addChildTarget(prepareParadigm(self.evidSpec, self.disc, self.paramFile, 0, self.paradigmExec, self.inferSpec, self.dogmaLib, self.pathwayLib, "%s/train/real" % (self.directory)))
        
        ## kick off train.shuffle runs
        for shuffle in range(1, self.nShuffle+1):
            if not os.path.exists("train/shuffle%s/outputFiles" % (shuffle)):
                if os.path.exists("%s/shuffle%s" % (self.pathwayLib, shuffle)):
                    self.addChildTarget(prepareParadigm(self.evidSpec, self.disc, self.paramFile, 0, self.paradigmExec, self.inferSpec, self.dogmaLib, "%s/shuffle%s" % (self.pathwayLib, shuffle), "%s/train/shuffle%s" % (self.directory, shuffle)))
                else:
                    self.addChildTarget(prepareParadigm(self.evidSpec, self.disc, self.paramFile, 0, self.paradigmExec, self.inferSpec, self.dogmaLib, self.pathwayLib, "%s/train/shuffle%s" % (self.directory, shuffle)))
        
        ## kick off test.real and test.shuffle runs
        self.setFollowOnTarget(branchTest(self.evidSpec, self.disc, self.paradigmExec, self.inferSpec, self.dogmaLib, self.pathwayLib, self.nShuffle, self.directory))
        
class branchTest(Target):
    def __init__(self, evidSpec, disc, paradigmExec, inferSpec, dogmaLib, pathwayLib, nShuffle, directory):
        Target.__init__(self, time=10000)
        self.evidSpec = evidSpec
        self.disc = disc
        self.paradigmExec = paradigmExec
        self.inferSpec = inferSpec
        self.dogmaLib = dogmaLib
        self.pathwayLib = pathwayLib
        self.nShuffle = nShuffle
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        evidList = zip(re.split("\s", self.evidSpec)[0::2], re.split("\s", self.evidSpec)[1::2])
        
        ## kick off test.real prepareParadigm
        if os.path.exists("%s/real" % (self.pathwayLib)):
            cmd = "prepareParadigm.py -b \"%s\" -s same -n %s -i %s -e %s -d %s -p %s %s" % (self.disc, 0, self.inferSpec, self.paradigmExec, self.dogmaLib, "%s/real" % (self.pathwayLib), self.evidSpec)
        else:
            cmd = "prepareParadigm.py -b \"%s\" -s same -n %s -i %s -e %s -d %s -p %s %s" % (self.disc, 0, self.inferSpec, self.paradigmExec, self.dogmaLib, self.pathwayLib, self.evidSpec)
        self.addChildTarget(jtCmd(cmd, "%s/test/real" % (self.directory)))
        
        ## kick off test.real prepareParadigm
        for shuffle in range(1, self.nShuffle+1):
            if os.path.exists("%s/shuffle%s" % (self.pathwayLib, shuffle)):
                cmd = "prepareParadigm.py -b \"%s\" -s same -n %s -i %s -e %s -d %s -p %s %s" % (self.disc, 0, self.inferSpec, self.paradigmExec, self.dogmaLib, "%s/shuffle%s" % (self.pathwayLib, shuffle), self.evidSpec)
            else:
                cmd = "prepareParadigm.py -b \"%s\" -s same -n %s -i %s -e %s -d %s -p %s %s" % (self.disc, 0, self.inferSpec, self.paradigmExec, self.dogmaLib, self.pathwayLib, self.evidSpec)
            self.addChildTarget(jtCmd(cmd, "%s/test/shuffle%s" % (self.directory, shuffle)))
        
        ## copy params
        self.setFollowOnTarget(runTest(self.nShuffle, self.directory))
        

class runTest(Target):
    def __init__(self, nShuffle, directory):
        Target.__init__(self, time=10000)
        self.nShuffle = nShuffle
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## kick off test.real run
        trainDir = "%s/train/real" % (self.directory)
        testDir = "%s/test/real" % (self.directory)
        system("rm -f %s/params*" % (testDir))
        index = 0
        while os.path.exists("%s/params%s.txt" % (trainDir, index)):
            system("cp %s/params%s.txt %s" % (trainDir, index, testDir))
            index += 1
        self.addChildTarget(FinalRun(index-1, testDir))
        
        ## kick off test.shuffle runs
        for shuffle in range(1, self.nShuffle+1):
            trainDir = "%s/train/shuffle%s" % (self.directory, shuffle)
            testDir = "%s/test/shuffle%s" % (self.directory, shuffle)
            index = 0
            while os.path.exists("%s/params%s.txt" % (trainDir, index)):
                system("cp %s/params%s.txt %s" % (trainDir, index, testDir))
                index += 1
            self.addChildTarget(FinalRun(index-1, testDir))

def jtFitness():
    ## parse arguments
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help = "Add as a child of jobFile rather " +
                      "than making a new jobTree")
    parser.add_option("-d", "--dogma", dest="dogmaPath", default="")
    parser.add_option("-p", "--pathway", dest="pathwayPath", default="")
    parser.add_option("-b", "--boundaries", dest="discBound", default="")
    parser.add_option("-s", "--shuffle", dest="shuffleNode", default="")
    parser.add_option("-n", "--nulls", dest="nNulls", default="10")
    parser.add_option("-t", "--storedparam", dest="paramFile", default="")
    
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
        if not dogma.startswith("/"):
            dogma = "%s/%s" % (os.getcwd(), dogma)        
    if len(options.pathwayPath) == 0:
        pathway = "%s/%s" % (pathwayDir, pathwayDefault)
    else:
        pathway = options.pathwayPath
        if not pathway.startswith("/"):
            pathway = "%s/%s" % (os.getcwd(), pathway)
    if len(options.shuffleNode) == 0:
        shuffleNode = "genome"
    else:
        shuffleNode = options.shuffleNode
    nShuffle = int(options.nNulls)
    if len(options.paramFile) == 0:
        paramFile = None
    else:
        paramFile = options.paramFile

    ## clean
    if len(args) == 1:
        if args[0] == "clean":
            print "rm -rf .jobTree fold*"
            os.system("rm -rf .jobTree fold*")
            sys.exit(0)
    
    ## run
    logger.info("options: " + str(options))
    s = Stack(branchFolds(" ".join(evidList), disc, paramFile, paradigmExec, inferSpec, dogma, pathway, shuffleNode, nShuffle, mFolds, os.getcwd()))
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
