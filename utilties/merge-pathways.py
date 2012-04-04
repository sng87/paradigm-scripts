#!/usr/bin/env python
"""merge-pathways.py: 

Usage:
  merge-pathways.py [options] output input [input, ...]

Options:
  -q            run quietly
"""
## Written By: Sam Ng
## Last Updated: ##/##/####
import os, sys, getopt, re
import mPathway

verbose = True

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg, die = False):
    if verbose:
        sys.stderr.write(msg)
    if die:
        sys.exit(1)

def syscmd(cmd):
    log("running:\n\t"+cmd+"\n")
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % exitstatus
        sys.exit(10)
    log("... done\n")

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) < 2:
        log("ERROR: incorrect number of arguments", die = True)
    
    prefix = args[0]
    inputArguments = args[1:]
    
    global verbose
    for o, a in opts:
        if o == "-q":
            verbose = False
    
    ## execute
    inputPathways = []
    for element in inputArguments:
        if os.path.isdir(element):
            for file in os.listdir(element):
                if file.endswith("pathway.tab"):
                    inputPathways.append(file)
        elif element.endswith("pathway.tab"):
            inputPathways.append(element)
    
    ## append pathways
    outPathway = mPathway.Pathway({}, {})
    for file in inputPathways:
        (nodes, interactions) = mPathway.rPathway(file)
        appendPathway = mPathway.Pathway(nodes, interactions)
        outPathway = mPathway.combinePathways(outPathway, appendPathway)
        
    ## write pathways
    mPathway.wPathway(prefix, outPathway.nodes, outPathway.interactions)

if __name__ == "__main__":
    main(sys.argv[1:])
