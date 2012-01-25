#!/usr/bin/python
import sys
import getopt
import time

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

def transformFile(fh, startRow, startCol, sep="\t"):
    startTime = time.time()
    log("reading file...")

    stringValues = [l.rstrip().split(sep) for l in fh] # everything in file

    numRows = len(stringValues)
    if (numRows == 0):
        print "Empty input"
        exit(10)
    numCols = len(stringValues[0])
    log("read %i rows %i columns" % (numRows, numCols))
    printElapsed(startTime)

    floatValues = [] # everything that's not missing data

    log("parsing float values...")
    for i, vals in enumerate(stringValues):
        if (i < startRow): continue
        if i > 0: assert len(vals) == len(stringValues[i-1]), \
                ("lines %i and %i have uneven columns" % 
                           (startRow + i + 1, startRow + i))
        if startCol >= len(vals):
            assert False, ("not enough columns in line %i" % i + startRow + 1)
        for j, val in enumerate(vals):
            if (j < startCol): continue
            try:
                fval = float(val)
                if fval != fval:
                    raise ValueError
                floatValues.append(val)
            except ValueError:
                pass

    totalValues = len(floatValues)
    if totalValues == 0:
        assert False, "did not read any values"
    log("%f%% missing data\n" % (100 - 100*float(totalValues)
                      / ((numRows - startRow) * (numCols - startCol))))
    printElapsed(startTime)

    floatValues.sort(key=lambda x: float(x))
    log("sorted all values\n")
    printElapsed(startTime)
    rankDict = dict()
    for rank, v in enumerate(floatValues):
        rankDict[v] = rankTransform(rank, totalValues)

    log("rank transformed all data\n")
    printElapsed(startTime)

    # for i, vals in enumerate(stringValues):
    #     for j, val in enumerate(vals):
    #         if j > 0: print sep,
    #         if (j < startCol or parseFloat(val) == None):
    #             print val,
    #         else:
    #             print rankDict[(i,j)],
    #     print ""


    def rowString(j):
        def matrixVal(i):
            val = stringValues[i][j]
            if val in rankDict: return ("%5g" % rankDict[val])
            else: return val
        return "\t".join(map(matrixVal, range(numRows)))

    for j in range(numCols):
        print rowString(j)
    log("wrote transposed matrix")
    printElapsed(startTime)

def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "r:c:v")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)

    if (len(args) == 0):
        fh = sys.stdin
    elif (len(args) == 1):
        fh = open(args[0], 'r')
    else:
        usage(1)

    rows = 1
    columns = 1
    global verbose
    for o, a in opts:
        if o == "-r":
            rows = int(a)
        elif o == "-c":
            columns = int(a)
        elif o == "-v":
            verbose = True;
        else:
            assert False, "unhandled option"
    
    transformFile(fh, rows, columns)

if __name__ == "__main__":
    main(sys.argv)
