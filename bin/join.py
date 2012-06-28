#!/usr/bin/env python
"""join.py: 

Usage:
  join.py [options] file1 file2 [file3 ...]

Options:
  -h            header
  -q            run quietly
"""
import os, os.path, sys, getopt, re

delim = "\t"
verbose = True

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg, die = False):
    if (verbose):
        sys.stderr.write(msg)
    if die:
        sys.exit(1)

def readFile(inFile, header = True):
    dataWidth = None
    dataMap = {}
    f = open(inFile, "r")
    if header:
        line = f.readline()
        if line.isspace():
            log("ERROR: missing header\n", die = True)
        pline = re.split(delim, line.rstrip("\n\r"))
        dataMap["HEADER"] = pline
        dataWidth = len(pline[1:])
    for line in f:
        if line.isspace():
            continue
        pline = re.split(delim, line.rstrip("\n\r"))
        if dataWidth is None:
            dataWidth = len(pline[1:])
        assert(len(pline[1:]) == dataWidth)
        dataMap[pline[0]] = pline[1:]
    f.close()
    return (dataMap, dataWidth)

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "hq")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) > 0:
        files = args
    else:
        files = []
        for i in sys.stdin:
           files.append(i.rstrip("\n\r"))
    
    if len(files) < 2:
        print "incorrect number of arguments"
        usage(1)

    header = False
    global verbose
    for o, a in opts:
        if o == "-h":
            header = True
        elif o == "-q":
            verbose = False
    
    ## read files
    fileData = {}
    fileWidth = {}
    for file in files:
        (fileData[file], fileWidth[file]) = readFile(file, header = header)
    features = list(set(fileData[files[0]].keys()) - set(["HEADER"]))
    features.sort()
    
    ## output
    if header:
        lineElements = [fileData[files[0]]["HEADER"][0]]
        for file in files:
            lineElements += fileData[file]["HEADER"][1:]
        print "%s" % (delim.join(lineElements))
    for feature in features:
        lineElements = []
        for file in files:
            if feature in fileData[file]:
                lineElements += fileData[file][feature]
            else:
                lineElements += ["" for i in range(fileWidth[file])]
        print "%s" % (feature + delim + delim.join(lineElements))

if __name__ == "__main__":
    main(sys.argv[1:])
