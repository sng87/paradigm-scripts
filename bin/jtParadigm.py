#!/usr/bin/env python
## Written by: Charles Vaske
## Modified by: Sam Ng
import sys
import os
import os.path
import re
import resource
import glob

from optparse import OptionParser

from jobTree.src.bioio import logger
from jobTree.src.bioio import system

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

collectParamsExec = "collectParameters"
mergeSwarm = "mergeSwarmFiles.py"
mergeMerge = "merge_merged.py"

class ParadigmCmd(Target):
    def __init__(self, command, cwd):
        Target.__init__(self, time=1000)
        self.cmd = command
        self.cwd = cwd

    def run(self):
        os.chdir(self.cwd)
        resource.setrlimit(resource.RLIMIT_CORE, (0,0))
        system(self.cmd)

class MaximizationIteration(Target):
    def __init__(self, iteration, tolerance, cwd):
        Target.__init__(self, time=10000)
        self.iteration = iteration
        self.tolerance = tolerance
        self.cwd = cwd
    def readLL(self, filename):
        f = open(filename, "r")
        topline = f.readline().rstrip()
        f.close()
        m = re.search("logZ=([0-9.e+-]*)", topline)
        return float(m.group(1))

    def emHasTerminated(self):
        if self.iteration < 2:
            return False
        prevLL = self.readLL("params%i.txt" % (self.iteration))
        currLL = self.readLL("params%i.txt" % (self.iteration + 1))
        decrease = ((prevLL - currLL) / currLL)
        logger.info("LL: %5g, Decrease: %3g" % (currLL, 100*decrease))
        return decrease < self.tolerance

    def run(self):
        os.chdir(self.cwd)
        cmd = "%s -p outputFilesEM/*learn* " % collectParamsExec
        if (os.path.exists("mask.expectations")):
            cmd += " mask.expectations "
        cmd += "| %s -o params%i.txt /dev/stdin " \
                       % (collectParamsExec, self.iteration + 1)
        if (os.path.exists("mask.params")):
            cmd += " mask.params "
        system(cmd)
        if self.emHasTerminated():
            self.setFollowOnTarget(FinalRun(self.iteration + 1, self.cwd))
        else:
            self.setFollowOnTarget(ExpectationIteration(self.iteration + 1, 
                                                        self.tolerance, self.cwd))
        
        

class ExpectationIteration(Target):
    def __init__(self, iteration, tolerance, cwd):
        Target.__init__(self, time=1000)
        self.iteration = iteration
        self.tolerance = tolerance
        self.cwd = cwd

    def run(self):
        os.chdir(self.cwd)
        system("rm -f params.txt")
        system("ln -s params%i.txt params.txt" % self.iteration)

        system("mkdir -p outputFilesEM%i" % self.iteration)
        system("rm -f outputFilesEM")
        system("ln -s outputFilesEM%i outputFilesEM" % self.iteration)

        sys.stderr.write("Current directory: " + os.getcwd() + "\n")
        jfile = open("jobsEM.list", "r")
        for job in jfile:
            self.addChildTarget(ParadigmCmd(job, self.cwd))
        jfile.close()
        self.setFollowOnTarget(MaximizationIteration(self.iteration, 
                                                     self.tolerance, self.cwd))

class FinalRun(Target):
    def __init__(self, iteration, cwd):
        Target.__init__(self, time=10000)
        self.iteration = iteration
        self.cwd = cwd

    def run(self):
        os.chdir(self.cwd)
        system("rm -f params.txt")
        system("ln -s params%i.txt params.txt" % self.iteration)
        system("mkdir -p outputFiles")

        jfile = open("jobs.list", "r")
        for job in jfile:
            self.addChildTarget(ParadigmCmd(job, self.cwd))
        jfile.close()
        self.setFollowOnTarget(Merge(self.cwd))

class Merge(Target):
    def __init__(self, cwd):
        Target.__init__(self, time=10000)
        self.cwd = cwd
    def run(self):
        os.chdir(self.cwd)
        system("mkdir -p mergeFiles")
        system("%s outputFiles mergeFiles" % mergeSwarm)
        mergeFiles = glob.glob("mergeFiles/*transpose*")
        system("%s bioInt mergeFiles/ merge_merged.tab" % mergeMerge)
        #if len(mergeFiles) == 1: # a global pathway
        #    system("mv merge_merged.tab merge_merged_unfiltered.tab")
        #    system("filterFeatures.py -n merge_merged_unfiltered.tab 1,0.01" +
        #           "> merge_merged.tab")

def commandAvailable(executable):
    return 0 == os.system("which %s > /dev/null 2> /dev/null" % executable)

def main():
    ## Make sure we're in the right type of directory
    assert os.path.exists("jobs.list")
    assert os.path.exists("jobsEM.list")
    assert os.path.exists("config.txt")
    assert os.path.exists("configEM.txt")
    assert os.path.exists("params0.txt")

    assert commandAvailable(collectParamsExec)
    assert commandAvailable(mergeSwarm)
    assert commandAvailable(mergeMerge)

    ##
    ## Parse options
    ##
    parser = OptionParser()
    Stack.addJobTreeOptions(parser) # so that the stack will work
    parser.add_option("--jobFile", help="Add as a child of jobFile rather " +
                      "than making a new jobTree")
    options, args = parser.parse_args()
    print "Using Batch System '" + options.batchSystem + "'"
    assert len(args) == 0 or len(args) == 1

    tolerance = 0.001
    if len(args) == 1:
        tolerance = float(args[0])

    logger.info("options: " + str(options))

    ##
    ## Run
    ##
    logger.info("starting first EM iteration")
    s = Stack(ExpectationIteration(0, tolerance, os.getcwd()))
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

if __name__ == "__main__":
    from jtParadigm import *
    main()
