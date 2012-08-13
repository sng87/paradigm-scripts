#!/usr/bin/env python
"""circlePlot.py: 

Usage:
  circlePlot.py [options] outputDir inputFile [inputFile ...]

Options:
  -s str        list file containing samples to include
  -f str        list file containing features to include
  -o str        feature;file[,file ...] or feature
  -c str        file to use as center colors
  -l            print the feature identifier in the circle or not (default: FALSE)
  -q            run quietly
"""
## Written By: Steve Benz and Zack Sanborn
## Modified By: Sam Ng
## Last Updated: 10/10/2011
import getopt, math, os, sys, re
from matplotlib import *
use('Agg')
from pylab import *
from random import random
import mData

verbose = True
tstep = 0.01

class rgb:
    def __init__(self,r,g,b):
        self.r = int(round(r))
        self.g = int(round(g))
        self.b = int(round(b))
        
        if self.r > 255:
            self.r = 255
        elif self.r < 0:
            self.r = 0
        if self.g > 255:
            self.g = 255
        elif self.g < 0:
            self.g = 0
        if self.b > 255:
            self.b = 255
        elif self.b < 0:
            self.b = 0
        
    def tohex(self):
        r = self.r
        g = self.g
        b = self.b
        hexchars = "0123456789ABCDEF"
        return "#" + hexchars[r / 16] + hexchars[r % 16] + hexchars[g / 16] + hexchars[g % 16] + hexchars[b / 16] + hexchars[b % 16]

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

def scmp(a, b, feature, dataList):
    dataFeature = feature
    if (a not in dataList[0]) & (b in dataList[0]):
        return(1)
    elif (a in dataList[0]) & (b not in dataList[0]):
        return(-1)
    elif (b not in dataList[0]) & (a not in dataList[0]):
        return(0)
    if dataFeature not in dataList[0][a]:
        if "*" in dataList[0][a]:
            dataFeature = "*"
        else:
            return(0)
    val = cmp(dataList[0][a][dataFeature], dataList[0][b][dataFeature])
    if val == 0:
        if len(dataList) > 1:
            val = scmp(a, b, feature, dataList[1:])
        else:
            return(0)
    return(val)

def polar(r, val):
    theta = -2.0 * math.pi * val + math.pi/2.0
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    return x, y

def getColor(val, minVal, maxVal, minColor = rgb(0, 0, 255), zeroColor = rgb(255, 255, 255), maxColor = rgb(255, 0, 0)):
    try:
        fval = float(val)
        if fval != fval:
            raise ValueError
    except ValueError:
        col = rgb(200,200,200)
        return col.tohex()
    if fval < 0.0:
        if fval < minVal:
            fval = -1.0
        else:
            fval = fval / minVal
        col = minColor
    else:
        if fval > maxVal:
            fval = 1.0
        else:
            fval = fval/maxVal
        col = maxColor
    r = fval * float(col.r - zeroColor.r) + zeroColor.r
    g = fval * float(col.g - zeroColor.g) + zeroColor.g
    b = fval * float(col.b - zeroColor.b) + zeroColor.b
    try:
        col = rgb(r,g,b)
    except ValueError:
        col = rgb(200,200,200)
    return col.tohex()

def plotScale(imgFile, minVal, maxVal):
    imgSize = (2, 4)
    fig = plt.figure(figsize=imgSize, dpi=100, frameon=True, facecolor='w')
    for i in range(10):
        val = minVal+i*(maxVal-minVal)/10
        col = getColor(val, minVal, maxVal)
        X = [float(i)/10, float(i+1)/10, float(i+1)/ 10, float(i)/10, float(i)/10]
        Y = [1, 1, 0, 0, 1]
        fill(X, Y, col, lw = 1, ec = col)
    savefig(imgFile)
    close()

def plotCircle(imgFile, label = "", centerCol = rgb(255, 255, 255).tohex(), circleCols = [[rgb(200, 200, 200).tohex()]], innerRadTotal=0.2, outerRadTotal=0.5, width = 5):
    ## image settings
    imgSize = (width, width)
    fig = plt.figure(figsize=imgSize, dpi=100, frameon=True, facecolor='w')
    axes([0, 0, 1, 1], frameon=True, axisbg='w')
    axis('off')
    circleWid = (outerRadTotal-innerRadTotal)/float(len(circleCols))
    
    ## color center
    outerRad = innerRadTotal
    outerRad -= .01
    X = []
    Y = []
    x, y = polar(outerRad, 0)
    X.append(x)
    Y.append(y)
    ti = 0
    while ti < 1:
        x, y = polar(outerRad, ti)
        X.append(x)
        Y.append(y)
        ti += tstep
        if ti > 1:
            break
    x, y = polar(outerRad, 1)
    X.append(x)
    Y.append(y)
    fill(X, Y, centerCol, lw = 1, ec = centerCol)
    
    ## color rings
    for i in range(len(circleCols)):
        innerRad = (i*circleWid)+innerRadTotal
        outerRad = ((i+1)*circleWid)+innerRadTotal-.01
        for j in range(len(circleCols[i])):
            t0 = float(j)/len(circleCols[i])
            t1 = float(j+1)/len(circleCols[i])
            X = []
            Y = []
            x, y = polar(innerRad, t0)
            X.append(x)
            Y.append(y)
            ti = t0
            while ti < t1:
                x, y = polar(outerRad, ti)
                X.append(x)
                Y.append(y)
                ti += tstep
                if ti > t1:
                    break
            x, y = polar(outerRad, t1)
            X.append(x)
            Y.append(y)
            ti = t1
            while ti > t0:
                x, y = polar(innerRad, ti)
                X.append(x)
                Y.append(y)
                ti -= tstep
                if ti < t0:
                    break
            x, y = polar(innerRad, t0)
            X.append(x)
            Y.append(y)
            fill(X, Y, circleCols[i][j], lw = 1, ec = circleCols[i][j])
    
    ## save image
    text(0, 0, label, ha='center', va='center')
    xlim(-0.5, 0.5)
    ylim(-0.5, 0.5)
    savefig(imgFile)
    close()

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "s:f:o:c:lq")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    if len(args) < 2:
        usage(2)
    
    outputDir = args[0].rstrip("/")
    circleFiles = args[1:]
    
    sampleFile = None
    featureFile = None
    orderFeature = None
    centerFile = None
    printLabel = False
    global verbose
    for o, a in opts:
        if o == "-s":
            sampleFile = a
        elif o == "-f":
            featureFile = a
        elif o == "-o":
            sa = re.split(";", a)
            if len(sa) == 1:
                orderFeature = sa[0]
                orderFiles = []
            else:
                orderFeature = sa[0]
                orderFiles = re.split(",", sa[1])
        elif o == "-c":
            centerFile = a
        elif o == "-l":
            printLabel = True
        elif o == "-q":
            verbose = False
    
    ## execute
    samples = []
    features = []
    if sampleFile != None:
        samples = mData.rList(sampleFile)
    if featureFile != None:
        features = mData.rList(featureFile)
    
    ## read circleFiles
    circleData = []
    circleColors = []
    for i in range(len(circleFiles)):
        (data, cols, rows) = mData.rCRSData(circleFiles[i], retFeatures = True)
        circleData.append(data)
        minCol = rgb(0, 0, 255)
        zerCol = rgb(255, 255, 255)
        maxCol = rgb(255, 0, 0)
        if circleFiles[i].endswith("meth"):
            maxCol = rgb(0, 0, 255)
            minCol = rgb(255, 0, 0)
            log("Color: meth\n")
        elif circleFiles[i].startswith("mut."):
            maxCol = rgb(0, 0, 0)
            minCol = rgb(255, 255, 255)
            log("Color: mut\n")
        circleColors.append( (minCol, zerCol, maxCol) )
        if sampleFile == None:
            samples = list(set(cols) | set(samples))
        if featureFile == None:
            features = list(set(rows) | set(features))
    
    ## read centerFile
    centerData = None
    if centerFile != None:
        centerData = mData.r2Col(centerFile, header = True)
        
    ## sort
    if orderFeature != None:
        if len(orderFiles) > 0:
            orderData = []
            orderColors = []
            for i in range(len(orderFiles)):
                orderData.append(mData.rCRSData(orderFiles[i]))
                minCol = rgb(255, 255, 255)
                zerCol = rgb(255, 255, 255)
                maxCol = rgb(0, 0, 0)
                orderColors.append( (minCol, zerCol, maxCol) )
        else:
            orderData = circleData
        samples.sort(lambda x, y: scmp(x, y, orderFeature, orderData))
        
        ## cohort png
        if len(orderFiles) > 0:
            imgFile = "%s/Cohort.png" % (outputDir)
            label = "Cohort"
            centerCol = rgb(255, 255, 255).tohex()
            circleCols = []
            for i in range(len(orderData)):
                ringCols = []
                ringVals = []
                for sample in samples:
                    if sample in orderData[i]:
                        if orderFeature in orderData[i][sample]:
                            ringVals.append(orderData[i][sample][orderFeature])
                        elif "*" in orderData[i][sample]:
                            ringVals.append(orderData[i][sample]["*"])
                minVal = min([-0.01]+mData.floatList(ringVals))
                maxVal = max([0.01]+mData.floatList(ringVals))
                for sample in samples:
                    if sample in orderData[i]:
                        if orderFeature in orderData[i][sample]:
                            ringCols.append(getColor(orderData[i][sample][orderFeature], minVal, maxVal, minColor = orderColors[i][0], zeroColor = orderColors[i][1], maxColor = orderColors[i][2]))
                        elif "*" in orderData[i][sample]:
                            ringCols.append(getColor(orderData[i][sample]["*"], minVal, maxVal, minColor = orderColors[i][0], zeroColor = orderColors[i][1], maxColor = orderColors[i][2]))
                        else:
                            ringCols.append(rgb(200, 200, 200).tohex())
                    else:
                        ringCols.append(rgb(200, 200, 200).tohex())
                circleCols.append(ringCols)
            plotCircle(imgFile, label = label, centerCol = centerCol, circleCols = circleCols, innerRadTotal=0.2, outerRadTotal=0.5, width = 5)
        
    ## plot images
    for feature in features:
        log("Drawing %s\n" % (feature))
        imgName = re.sub("[/:]", "_", feature)
        if len(imgName) > 100:
            imgName = imgName[:100]
        imgFile = "%s/%s.png" % (outputDir, imgName)
        label = ""
        if printLabel:
            label = feature
        centerCol = rgb(255, 255, 255).tohex()
        if centerData != None:
            if feature in centerData:
                minVal = min([-0.01]+mData.floatList(centerData.values()))
                maxVal = max([0.01]+mData.floatList(centerData.values()))
                centerCol = getColor(centerData[feature], minVal, maxVal)
                log("\t%s,%s,%s,%s\n" % (centerData[feature],minVal,maxVal,centerCol))
        circleCols = []
        for i in range(len(circleData)):
            ringCols = []
            ringVals = []
            for sample in samples:
                if sample in circleData[i]:
                    if feature in circleData[i][sample]:
                        ringVals.append(circleData[i][sample][feature])
                    elif "*" in circleData[i][sample]:
                        ringVals.append(circleData[i][sample]["*"])
            minVal = min([-0.01]+mData.floatList(ringVals))
            maxVal = max([0.01]+mData.floatList(ringVals))
            for sample in samples:
                if sample in circleData[i]:
                    if feature in circleData[i][sample]:
                        ringCols.append(getColor(circleData[i][sample][feature], minVal, maxVal, minColor = circleColors[i][0], zeroColor = circleColors[i][1], maxColor = circleColors[i][2]))
                    elif "*" in circleData[i][sample]:
                        ringCols.append(getColor(circleData[i][sample]["*"], minVal, maxVal, minColor = circleColors[i][0], zeroColor = circleColors[i][1], maxColor = circleColors[i][2]))
                    else:
                        ringCols.append(rgb(200, 200, 200).tohex())
                else:
                    ringCols.append(rgb(200, 200, 200).tohex())
            circleCols.append(ringCols)
        plotCircle(imgFile, label = label, centerCol = centerCol, circleCols = circleCols, innerRadTotal=0.2, outerRadTotal=0.5, width = 5)

if __name__ == "__main__":
    main(sys.argv[1:])
