#!/usr/bin/env python
import os
import sys
import getopt
import time
import random

verbose = True
minSampleFrac = 0.2
minGeneFrac = 0.2
outputPrefix = "na"
samplePrefix = "na_iter"
numberingOffset = 0
batches = 1
sameSample = False
trueFileDir = ''
geneIntersection = True

def usage(code=0):
    print "createNullFiles.py: create null data files from tuples of data"
    print "    This will create a new output file for each input file"
    print "    where matrix entry i, j in each output file comes from a"
    print "    corresponding i', j' from the input files."
    print "IMPORTANT: to create permutations of genes within the same "
    print "    sample, pass the string \"same\" as the first argument"
    print "Usage:"
    print "   createNullFiles.py [options] nullSamples file file2 file3..."
    print ""
    print "Options: "
    print "   -p string  prefix for output files (defaults to 'na')"
    print "   -b int     number of batches (appended to prefix) in filenames"
    print "   -o int     offset for sample numbering (defaults to 0)"
    print "   -s float   minimum fraction of samples present in all files"
    print "   -g float   minimum fraction of genes present in all files"
    print "   -t dir     write restricted true files here"
    print "   -u         use union of gene list rather than intersection"
    print "   -q         don't output logging information"
    if code != None:
        sys.exit(code)

class NamedMatrix(dict):
    @classmethod
    def fromFile(cls, filename, sep="\t"):
        fh = open(filename, "r")
        header = fh.readline().rstrip("\r\n").split(sep)
        self = cls(header[0], header[1:])
        for l in fh:
            vals = l.rstrip("\r\n").split(sep)
            self.addSample(vals[0], vals[1:])
        fh.close()
        return self
    def __init__(self, corner, colnames):
        self._corner = corner
        self.__setColNames(colnames)
    def __setColNames(self, colnames):
        self._colnames = colnames
        self._nameToCol = dict(zip(self._colnames, range(len(self._colnames))))
    def addSample(self, sampleName, vals):
        # print "adding sample %s with %i vals" % (sampleName, len(vals))
        self[sampleName] = vals
    def restrictColumns(self, columnList):
        l = len(self._colnames)
        colOrder = [self._nameToCol.get(name, l) for name in columnList]
        for row in self.keys():
            data = self[row]
            data.append("NA")
            self[row] = [data[col] for col in colOrder]
        self.__setColNames(columnList)
    def describe(self):
        print "%i rows, %i columns" % (len(self.keys()), len(self._colnames))
    def writeLineToFile(self, fh, label, vals, sep="\t"):
        fh.write(label)
        for v in vals:
            fh.write(sep)
            fh.write(v)
        fh.write("\n")
    def writeToFile(self, filename, sep="\t"):
        fh = open(filename, "w")
        self.writeLineToFile(fh, self._corner, self._colnames, sep)
        for row in self.keys():
            # print "writing row %s" % row
            self.writeLineToFile(fh, row, self[row], sep)
        fh.close()

def outputFileName(outputPrefix, fn):
    path = outputPrefix + os.path.basename(fn)
    return path

def log(msg):
    if (verbose):
        sys.stderr.write(msg)

def createNullFiles(numSamples, files):
    matrices = [NamedMatrix.fromFile(filename) for filename in files]
    samples = set(matrices[0].keys())
    if len(matrices) > 1:
        samples = samples.intersection(*[m.keys() for m in matrices[1:]])
    samples = list(samples)
    sampleFraction = [float(len(samples)) / len(m.keys()) for m in matrices]
    sFrac = [str(f) for f in sampleFraction]
    log("%i samples, fraction of each matrix: %s\n" % (len(samples), 
                                                       " ".join(sFrac)))
    if sameSample:
        numSamples = len(samples)
    
    genes = set(matrices[0]._colnames)
    if len(matrices) > 1:
        if geneIntersection:
            genes = genes.intersection(*[m._colnames for m in matrices[1:]])
        else:
            genes = genes.union(*[m._colnames for m in matrices[1:]])
    genes = list(genes)
    geneFraction = [float(len(genes)) / len(m._colnames) for m in matrices]
    gFrac = [str(f) for f in geneFraction]
    log("%i genes, fraction of each matrix: %s\n" % (len(genes), 
                                                     " ".join(gFrac)))
    
    if min(sampleFraction) < minSampleFrac:
        print "Need at least %f overlap in samples" % minSampleFrac
        sys.exit(3)
    if min(geneFraction) < minGeneFrac:
        print "Need at least %f overlap in genes" % minGeneFrac
        sys.exit(4)

    print "restricting down to %d genes" % len(genes)
    for m in matrices:
        m.restrictColumns(genes)

    for b in range(batches):
        print "writing batch %i" % (b+1)
        nullM = [NamedMatrix(m._corner, genes) for m in matrices]
        prefix = outputPrefix
        numOffset = 1+ numberingOffset + b*numSamples
        if batches > 1:
            prefix += "_" + str(b + 1) +  "_"
        for i in range(numSamples):
            geneIndices = range(len(genes))
            if sameSample:
                random.shuffle(geneIndices)
                randomG = geneIndices
                randomS = [samples[i] for s in geneIndices]
                sample = samplePrefix + "_" + str(b+1) + "_" + samples[i]
            else:
                randomG = [random.choice(geneIndices) for gi in geneIndices]
                randomS = [random.choice(samples) for s in geneIndices]
                sample = samplePrefix + str(i + numOffset)
            for null, m in zip(nullM, matrices):
                data = [m[s][g] for s, g in zip(randomS, randomG)]
                null.addSample(sample, data)
                
        for fn, m in zip(files, nullM):
            m.describe()
            m.writeToFile(outputFileName(prefix, fn))

    if trueFileDir != '':
        for fn, m in zip(files, matrices):
            restrictM = NamedMatrix(m._corner, m._colnames)
            for s in samples: restrictM.addSample(s, m[s]) 
            outFileName = os.path.join(trueFileDir, os.path.basename(fn))
            restrictM.writeToFile(outFileName)

def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "s:b:g:qp:o:t:u")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)

    global geneIntersection
    global minSampleFrac, minGeneFrac, outputPrefix, numberingOffset, verbose
    global batches, trueFileDir, sameSample
    for o, a in opts:
        if o == "-s":
            minSampleFrac = float(a)
        if o == "-b":
            batches = int(a)
        if o == "-g":
            minGeneFrac = float(a)
        if o == "-p":
            outputPrefix = a
        if o == "-o":
            numberingOffset = int(a)
        if o == "-t":
            trueFileDir = a
        if o == "-q":
            verbose = False
        if o == "-u":
            geneIntersection = False
            minGeneFrac = 0.0

    if (len(args) < 2):
        print "Not enough arguments: specify number of samples and >=1 file"
        usage(1)

    if args[0] == "same":
        sameSample = True
        numSamples = 0
    else:
        try:
            numSamples = int(args[0])
        except ValueError:
            print "First argument must be an integer"
            usage(1)

    createNullFiles(numSamples, args[1:])

if __name__ == "__main__":
    main(sys.argv)
