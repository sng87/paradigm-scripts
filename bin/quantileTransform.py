#!/usr/bin/env python
import sys
import getopt
import time
import array
import csv
import math
import ctypes
import ctypes.util

libc = ctypes.cdll.LoadLibrary(ctypes.util.find_library("c"))

verbose = False;

def usage(code=0):
    print "quantileTransform.py: rank transform an evidence file for Paradigm"
    print ""
    print "Usage:"
    print "   transformEvidence.py [options] filename"
    print "Options:"
    print "  -r int   number of header rows"
    print "  -c int   number of header columns "
    print "  -v       print progress reports to stderr"
    if (code != None):
        sys.exit(code)

def rankTransform(rank, total):
    return rank / float(total)

def log(msg):
    if (verbose):
        sys.stderr.write(msg)

def printElapsed(start):
    ms = int((time.time() - start) * 1000)
    log(" %ims elapsed\n" % ms)


def py_cmp_float(a, b):
    return (a > b) - (a < b)

CMPFUNC = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float))

cmp_float = CMPFUNC(py_cmp_float)

def csort(buf):
    """
    This is an inplace sort of an array.array float class using the C qsort function
    """
    addr, count = buf.buffer_info()
    libc.qsort( ctypes.cast(addr, ctypes.POINTER(ctypes.c_float)), count, ctypes.sizeof(ctypes.c_float), cmp_float)     


def transformFile(fh, sep="\t"):
    startTime = time.time()
    log("reading file...")
    
    floatValues = None
    cols = None
    rows = None
    totalValues = 0
    reader = csv.reader(fh, delimiter=sep)
    for row in reader:
        if cols is None:
            cols = row[1:]
            numCols = len(cols)
            rows = []
            floatValues = array.array('f')
        else:
            rows.append(row[0])
            assert(len(row)-1 == numCols)
            for val in row[1:]:
                try:
                    floatValues.append(float(val))
                    totalValues += 1
                except ValueError:
                    floatValues.append(float('nan'))

    numRows = len(rows)
    if (numRows == 0):
        print "Empty input"
        exit(10)

    log("read %i rows %i columns" % (numRows, numCols))
    printElapsed(startTime)

    if totalValues == 0:
        assert False, "did not read any values"
    log("%f%% missing data\n" % (100 - 100*float(totalValues)
                      / ((numRows) * (numCols))))
    printElapsed(startTime)
    
    rankDict = dict()
    log("sorting float values...")
    i = 0
    sortedValues = array.array('f', floatValues) 
    
    csort(sortedValues)
    
    for val in sortedValues:
        if not math.isnan(val):
            rankDict[val] = rankTransform(i, totalValues)
            i += 1

    log("rank transformed all data\n")
    printElapsed(startTime)

    def rowString(j):
        def matrixVal(i):
            val = floatValues[i*numCols + j]
            if val in rankDict: return ("%5g" % rankDict[val])
            else: return "NA"
        return "\t".join(map(matrixVal, range(numRows)))

    print "samples\t%s" % ("\t".join(rows))
    for j in range(numCols):
        print "%s\t%s" % (cols[j], rowString(j))
    log("wrote transposed matrix")
    printElapsed(startTime)

def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "v")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)

    if (len(args) == 0):
        fh = sys.stdin
    elif (len(args) == 1):
        fh = open(args[0], 'r')
    else:
        usage(1)

    global verbose
    for o, a in opts:
        if o == "-v":
            verbose = True;
        else:
            assert False, "unhandled option"
    
    transformFile(fh)

if __name__ == "__main__":
    main(sys.argv)
