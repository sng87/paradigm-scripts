#!/usr/bin/env python
"""prepareParadigm.py: gather data files and create job lists

Creates in the current directory:
  clusterFiles/  the specified database tables or tab files
                 will be rank transformed and also pathway 
                 files are placed in this directory

  configEM.txt  the configuration file for EM runs
  jobsEM.list   the list of parasol jobs for EM

  config.txt     the configuration for the final run
  jobs.list      the list of parasol jobs fro the final run

Usage:
  prepareParadigm.py [options] attach1 evid1 [attach2 evid2 ...]

Evidence is specified as:
    table:<hg18_tablename>   where <hg18_tablename> is a table in hg18
    tableNoNorm:<hg18        like table:, but without median centering genes
    bioInt:<tablename>       where <tablename> is a table in bioInt
    http:<url>               read file from URL
    http:<url>               read file from URL, no quantile transform
    file:<tabfile_name>      where <tabfile_name> is the path to a file
    rawFile:<tabfile>        where <tabfile> is a file that won't be run
                             through quantile transformation

Options:
   -n int               the number of null batches (500 each) to create (default: 2)
   -s int               size of null batches (or 'same', see createNullSamples.py)
   -e path              path to paradigm exe (default: /hive/users/$USER/bin)
   -d dir               dogma directory (take top of config from here)
   -p dir               the directory containing pathway files 
                        (default is /hive/groups/cancerGB/paradigm/pathwayfiles/v1)
   -t str               use initial parameters stored from another params.txt
   -i string            inference parameters 
                        (default is method=JTREE,updates=HUGIN,verbose=1)
   -c options           options to pass to createNullFiles.py (quote them all)
   -b flt;flt[,flt;flt] boundaries for discretization, use comma to specify different
                        boundaries per evidence (default 0.333;0.667)
   -q                   run quietly, don't output status
"""
## Written by: Charles Vaske
## Modified by: Sam Ng
import os, sys, glob, getopt, re, subprocess, math, json

###
### Experiments to find the path of the currently executing script
###
# print "os.path.curdir:", os.path.curdir
# print "os.path.realpath(os.path.curdir):", os.path.realpath(os.path.curdir)
# print "sys.argv[0]:", sys.argv[0]
# print "os.path.dirname(sys.argv[0]):", os.path.dirname(sys.argv[0])
# print "os.path.realpath(os.path.dirname(sys.argv[0])):", os.path.realpath(os.path.dirname(sys.argv[0]))
scriptDirectory = os.path.realpath(os.path.dirname(sys.argv[0]))
evidenceTypes = {
    "table": (("(cd %s/../extractData; " + 
               "./extractData -median -norm=median %%s)" +
               "| cut -f 1,3- | %s/quantileTransform.py /dev/stdin")
              % (scriptDirectory, scriptDirectory)),
    "tableNoNorm": (("(cd %s/../extractData; ./extractData -median %%s)" +
                     "| cut -f 1,3- | %s/quantileTransform.py /dev/stdin")
                    % (scriptDirectory, scriptDirectory)),
    "bioInt": (("%s/grabData.sh %%s" +
                " | %s/quantileTransform.py /dev/stdin")
               % (scriptDirectory, scriptDirectory)),
    "file":  (("cat %%s" +
               " | %s/quantileTransform.py /dev/stdin") % scriptDirectory),
    "http": "curl http:%s | quantileTransform.py /dev/stdin",
    "rawhttp": "curl http:%s | transpose.py /dev/stdin /dev/stdout",
    "rawFile":  "cat %s"
    }

dataDir = "clusterFiles"
#outputEmDir = "outputFilesEM"
#outputDir = "outputFiles"

verbose = True
standardAttach = ["genome", "mRNA", "protein", "active"]
paradigmExec = "/hive/users/" + os.getenv("USER") + "/bin/paradigm"

newSpecStyle = True
dryrun = False

nullOptions = ""
nullBatches = 2
nullBatchSize = 500

targetJobLength = 45 # seconds

disc = "0.333;0667"
paramFile = ""

### If dogmaDir is defined, files in that directory are copied to the 
### destination directory.  Additionally, if the files 'configTop' or
### 'configTopEM' are in that directory, their contents replace
### the following two variables (don't forget the %s for inference!)
dogmaDir = ''
configTop = """pathway [max_in_degree=5,param_file=params.txt]
inference [%s]\n"""
configTopEM = """pathway [max_in_degree=5,param_file=params.txt]
inference [%s]
em_step ["""
inference="method=JTREE,updates=HUGIN,verbose=1"
configEmLine = "em [max_iters=1,log_z_tol=1e-10]\n"
def configELine(evidStub):
    if "disc" in evidStub:
        b = evidStub["disc"]
    else:
        b = disc
    return "evidence [suffix=%s,node=%s,disc=%s,epsilon=0.01,epsilon0=0.2]\n" \
        % (evidStub["suffix"], evidStub["attachment"], b)

paramHeader="> parameters em_iters=0 logZ=-1e300\n"

mutationParams="""> shared CondProbEstimation [pseudo_count=1,target_dim=3,total_dim=9] %s=codeMut
0.0
0.2
0.8
0.9998
0.0001
0.0001
0.0
0.2
0.8
"""

mutationParamsMask="""> mask shared CondProbEstimation [pseudo_count=1,target_dim=3,total_dim=9] %s=codeMut
0
nan
nan
nan
nan
nan
0
nan
nan
"""

def norm(array):
    s = sum(array)
    return [v/float(s) for v in array]

def initParams(numBins, reverse=False):
    base = [v for v in norm(range(1, 1+numBins))]
    up   = [str(v) for v in base]
    if (reverse):
        up.reverse()
    zero = norm([base[min(i, len(base) - i - 1)] for i in range(len(base))])
    zero = [str(v) for v in zero]
    down = [v for v in up]
    down.reverse()
    return "\n".join(down + zero + up) + "\n"

def readParams(paramFile):
    storedParams = {}
    f = open(paramFile, "r")
    f.readline()
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            attachment = re.split("=", line)[-1]
            storedParams[attachment] = ""
        else:
            storedParams[attachment] += "%s\n" % (line)
    f.close()
    return storedParams

def writeBaseParamsFile(pfilename, evidence, storedParams = {}):
    writeHeader = not os.path.exists(pfilename)
    pfile = open(pfilename, "a")
    if writeHeader:
        pfile.write(paramHeader)
    for e in evidence:
        bins = len(e["disc"].split(";")) + 1
        global newSpecStyle
        if newSpecStyle:
            spec = e["attachment"]
        else:
            spec = "-obs>"
        if e["attachment"] != "codeMut":
            pfile.write("> shared CondProbEstimation [pseudo_count=1,target_dim=%i,total_dim=%i] %s=%s\n" % (bins, 3*bins, e["suffix"], spec))
            if e["attachment"] in storedParams:
                pfile.write(storedParams[e["attachment"]])
            else:
                pfile.write(initParams(bins, reverse=("reversed" in e)))
        else:
            pfile.write(mutationParams % e["suffix"])
            if os.path.exists("mask.params"):
                mfile = open("mask.params", mode="a")
                mfile.write(mutationParamsMask % e["suffix"])
                mfile.close()
    pfile.close()

def readPathwayTiming(directory):
    tfile = open(directory + "/timings.tab", "r")
    samplesline = tfile.readline().rstrip();
    m = re.search('^#\s*samples\s*(\d+)\s*', samplesline)
    if not m:
        print "missing samples line on pathway timings"
        sys.exit(1)
    samples = int(m.group(1))
    result = {}
    for line in tfile:
        timestring, pathway = line.rstrip().split("\t")
        result[pathway] = float(timestring) / samples
    return result

def numBuckets(pathway, samples, timings, targetLength):
    if pathway not in timings:
        return 1
    length = timings[pathway] * samples
    return min(int(math.ceil(length / targetLength)), samples)

def usage(code=0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg):
    if (verbose):
        sys.stderr.write(msg)

def evidenceStreamCommand(evidenceSpec):
    (type, sep, name) = evidenceSpec.partition(":")
    if type not in evidenceTypes:
        print "Unrecognized evidence spec \"%s:\"" % evidenceSpec 
        usage(1)
    return (evidenceTypes[type] % name)

def evidenceStub(attachment, evidspec, index):
    (type, sep, name)= evidspec.partition(":")
    bname = os.path.basename(name)
    (where, sep, options) = attachment.partition(":")
    if (len(re.split(",", disc)) > 1):
        eviddisc = re.split(",", disc)[index]
    else:
        eviddisc = disc
    stub = {"attachment" : where, 
            "spec" : evidspec, 
            "suffix" : bname, 
            "outputFile" : dataDir + "/" + bname,
            "disc" : eviddisc}
    if options != "":
        stub.update(json.loads(options))
    return stub

def readFileLineNumber(filename):
    wcArgs = ["sh", "-c", "cat %s | wc -l" % filename]
    return int(subprocess.Popen(wcArgs, 
                                stdout=subprocess.PIPE).communicate()[0])


def syscmd(cmd):
    log("running:\n    " + cmd + "\n")
    if dryrun:
        exitstatus = 0
    else:
        exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % exitstatus
        sys.exit(10)
    log(" ... done\n")

def mkdir(dirname):
    try:
        os.mkdir(dirname)
    except OSError, err:
        if err.strerror == "File exists":
            print "WARNING: directory %s already exists" % dirname
        else:
            print "couldn't make directory %s: %s" % (dirname, str(err))
            sys.exit(1)

def prepareParadigm(args):
    pathwayDir = '/hive/groups/cancerGB/paradigm/pathwayfiles/v1'
    try:
        opts, args = getopt.getopt(args, "p:n:e:qc:b:s:t:i:d:")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)

    if len(args) < 2 or len(args) % 2 != 0:
        print "need an even number of arguments"
        usage(1)

    global paradigmExec, dryrun, nullOptions, disc
    global nullBatches, nullBatchSize, paramFile, inference, dogmaDir
    global configTop, configTopEM
    for o, a in opts:
        if o == "-p":
            pathwayDir = a
        elif o == "-d":
            dogmaDir = a
            fn = os.path.join(dogmaDir, "configTop")
            if os.path.exists(fn):
                cfile = open(fn)
                configTop = cfile.read()
                cfile.close()
            fn = os.path.join(dogmaDir, "configTopEM")
            if os.path.exists(fn):
                cfile = open(fn)
                configTopEM = cfile.read().rstrip("\n")
                cfile.close()
        elif o == "-n":
            nullBatches = int(a)
        elif o == "-s":
            if a == "same":
                nullBatchSize = a
            else:
                nullBatchSize = int(a)
        elif o == "-e":
            paradigmExec = a
        elif o == "-q":
            verbose = False
        elif o == "-c":
            nullOptions = a
        elif o == "-b":
            disc = a
        elif o == "-t":
            paramFile = a
        elif o == "-i":
            inference = a
        
    log("Making sub-directories\n")
    mkdir(dataDir)

    evidence = [evidenceStub(a,e,i) for a, e, i in zip(args[0::2], args[1::2], range(len(args[0::2])))]
    for e in evidence:
        log("Evidence:\n")
        for k in e.keys():
            log("    %s\t%s\n" % (k,e[k]))
    
    for e in evidence:
        if (e["attachment"] not in standardAttach):
            print "WARNING: %s is non-standard: " % e["attachment"]
            print "         standard attachments are: " + str(standardAttach)
        cmd = evidenceStreamCommand(e["spec"]) + " > " + e["suffix"]
        syscmd(cmd)

    cmd = "%s/createNullFiles.py %s -t %s -p %s/na_batch -b %i %s " % \
        (scriptDirectory, nullOptions, dataDir, dataDir, 
         nullBatches, str(nullBatchSize)) \
        + " ".join([e["suffix"] for e in evidence])
    syscmd(cmd)

    # minus 1 for header
    samples = readFileLineNumber(evidence[0]["outputFile"]) - 1

    log("Writing config file for EM\n")
    confFile = open("configEM.txt", "w")
    confFile.write("# " + " ".join(sys.argv) + "\n")
    confFile.write(configTopEM % inference)
    if newSpecStyle:
        confFile.write(",".join([e["suffix"] + "=" + e["attachment"]
                                 for e in evidence]))
    else:
        confFile.write(",".join([e["suffix"] + "=-obs>" for e in evidence]))

    confFile.write("]\n")
    confFile.write(configEmLine)
    [confFile.write(configELine(e)) for e in evidence]
    confFile.close()

    log("Writing config file for final run\n")
    confFile = open("config.txt", "w")
    confFile.write(configTop % inference)
    [confFile.write(configELine(e)) for e in evidence]
    confFile.close()

    if dogmaDir:
        log("Copying dogma files\n")
        syscmd("cp %s/* ." % dogmaDir)

    log("Copying pathway files\n")
    syscmd("cp %s/*_pathway.tab %s" % (pathwayDir, dataDir))
    timings = readPathwayTiming(pathwayDir)
    
    pathFiles = glob.glob(dataDir + "/*_pathway.tab")

    log("writing EM jobs list\n")
    jfile = open("jobsEM.list", "w")
    for p in pathFiles:
        pathway = os.path.basename(p)
        buckets = numBuckets(pathway, samples, timings, targetJobLength)
        pid = pathway[0:-len("_pathway.tab")]
        if buckets == 1:
            emOut = "outputFilesEM/" + pid + "_learned_parameters.fa"
            jfile.write("%s -p %s -c configEM.txt -b %s/ -e %s\n" % 
                        (paradigmExec, p, dataDir, emOut))
        else:
            for b in range(buckets):
                bpid = "%s_b%i_%i" % (pid, b, buckets)
                emOut = "outputFilesEM/" + bpid + "_learned_parameters.fa"
                out = "outputFilesEM/" + bpid + "_output.fa"
                c = "%s -p %s -c configEM.txt -b%s/ -e %s -s %i,%i\n" % \
                    (paradigmExec, p, dataDir, emOut, b, buckets)
                jfile.write(c)

    jfile.close()

    log("writing jobs list\n")
    jfile = open("jobs.list", "w")
    for p in pathFiles:
        pathway = os.path.basename(p)
        buckets = numBuckets(pathway, samples, timings, targetJobLength)
        pid = pathway[0:-len("_pathway.tab")]
        if (buckets == 1):
            out = "outputFiles/" + pid + "_output.fa"
            jfile.write("%s -p %s -c config.txt -b %s/ -o %s\n" % 
                        (paradigmExec, p, dataDir, out))
        else:
            for b in range(buckets):
                bpid = "%s_b%i_%i" % (pid, b, buckets)
                out = "outputFiles/" + bpid + "_output.fa"
                c = "%s -p %s -c config.txt -b %s/ -o %s -s %i,%i\n" % \
                    (paradigmExec, p, dataDir, out, b, buckets)
                jfile.write(c)
        for b in range(1, 1 + nullBatches):
            if nullBatchSize == "same":
                numNullSamples = samples
            else:
                numNullSamples = nullBatchSize
            buckets = numBuckets(pathway, numNullSamples, 
                                 timings, targetJobLength)
            if buckets == 1:
                out = "outputFiles/" + pid + "_batch_" + str(b) + "_output.fa"
                c = "%s -p %s -c config.txt -b %s/na_batch_%i_ -o %s\n" % \
                    (paradigmExec, p, dataDir, b, out)
                jfile.write(c)
            else:
                for k in range(buckets):
                    bpid = "%s_b%i_%i" % (pid, k, buckets)
                    out = "outputFiles/%s_batch_%s_output.fa" % (bpid, str(b))
                    batch = "%s/na_batch_%i_" % (dataDir, b)
                    c = "%s -p %s -c config.txt -b %s -o %s -s %i,%i\n" % \
                        (paradigmExec, p, batch, out, k, buckets)
                    jfile.write(c)
    jfile.close()
    
    if len(paramFile) > 0:
        writeBaseParamsFile("params0.txt", evidence, storedParams = readParams(paramFile))
    else:
        writeBaseParamsFile("params0.txt", evidence)
    
    syscmd("ln -s params0.txt params.txt")

if __name__ == "__main__":
    prepareParadigm(sys.argv[1:])


