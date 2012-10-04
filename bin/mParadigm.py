"""
mPathway.py: Python module for handling pathways
Written By: Sam Ng
"""
import re, sys, os
from copy import deepcopy

class Pathway:
    def __init__(self, nodes, interactions, pid = None):
        self.nodes = nodes
        self.interactions = interactions
        self.pid = pid
    def removeNode(self, node):
        (self.nodes, self.interactions) = removeNode(node, self.nodes, self.interactions)
    def selfTest(self):
        ## check for unannotated features
        for source in self.interactions.keys():
            if source not in self.nodes:
                if source.endswith("(complex)"):
                    log("WARNING: %s not defined, default to complex\n" % (source))
                    self.nodes[source] = "complex"
                elif source.endswith("(family)"):
                    log("WARNING: %s not defined, default to family\n" % (source))
                    self.nodes[source] = "family"
                elif source.endswith("(abstract)"):
                    log("WARNING: %s not defined, default to abstract\n" % (source))
                    self.nodes[source] = "abstract"
                else:
                    log("WARNING: %s not defined, default to protein\n" % (source))
                    self.nodes[source] = "protein"
            for target in self.interactions[source].keys():
                if target not in self.nodes:
                    if target.endswith("(complex)"):
                        log("WARNING: %s not defined, default to complex\n" % (target))
                        self.nodes[target] = "complex"
                    elif target.endswith("(family)"):
                        log("WARNING: %s not defined, default to family\n" % (target))
                        self.nodes[target] = "family"
                    elif target.endswith("(abstract)"):
                        log("WARNING: %s not defined, default to abstract\n" % (target))
                        self.nodes[target] = "abstract"
                    else:
                        log("WARNING: %s not defined, default to protein\n" % (target))
                        self.nodes[target] = "protein"
        ## check for invalid links
        legalLinks = ["-t>", "-t|", "-a>", "-a|", "-ap>", "-ap|", "component>", "member>"]
        for source in self.interactions.keys():
            for target in self.interactions[source].keys():
                for link in re.split(";", self.interactions[source][target]):
                    if link not in legalLinks:
                        log("WARNING: illegal link %s\t%s\t%s found\n" % (source, target, link))
        ## report number of separate components
        components = sortConnected(self)
        log("found %s subpathway components\n" % (len(components)))
        (self.nodes, self.interactions) = constructInteractions(components[0], self.nodes, self.interactions)
    def wAliasMap(self, hugof, outf):
        hugoList = rList(hugof)
        if self.pid != None:
            prefix = "%s_" % (self.pid)
        else:
            prefix = ""
        f = open(outf, "w")
        for node in (set(hugoList) & set(self.nodes.keys())):
            f.write("%s\t%s\n" % (prefix+node, node))
        f.close()

def log(msg, die = False):
    """logger function"""
    sys.stderr.write(msg)
    if die:
        sys.exit(1)

def rPathway(inf, reverse = False, retProteins = False, delim = "\t"):
    """read UCSC pathway format"""
    proteins = set()
    readPathway = Pathway(dict(), dict())
    f = open(inf, "r")
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        if len(pline) == 2:
            readPathway.nodes[pline[1]] = pline[0]
            if pline[0] == "protein":
                proteins.update([pline[1]])
        elif len(pline) == 3:
            if reverse:
                if pline[1] not in readPathway.interactions:
                    readPathway.interactions[pline[1]] = dict()
                if pline[0] not in readPathway.interactions[pline[1]]:
                    readPathway.interactions[pline[1]][pline[0]] = pline[2]
                else:
                    readPathway.interactions[pline[1]][pline[0]] += ";"+pline[2]
            else:
                if pline[0] not in readPathway.interactions:
                    readPathway.interactions[pline[0]] = dict()
                if pline[1] not in readPathway.interactions[pline[0]]:
                    readPathway.interactions[pline[0]][pline[1]] = pline[2]
                else:
                    readPathway.interactions[pline[0]][pline[1]] += ";"+pline[2]
        else:
            print >> sys.stderr, "ERROR: line length not 2 or 3: \"%s\"" % (line)
            sys.exit(1)
    f.close()
    if retProteins:
        return(readPathway.nodes, readPathway.interactions, proteins)
    else:
        return(readPathway.nodes, readPathway.interactions)

def wPathway(outf, outNodes, outInteractions, useNodes = None):
    """write UCSC pathway format"""
    f = open(outf, "w")
    if useNodes == None:
        useNodes = outNodes.keys()
    for i in useNodes:
        if i not in outNodes:
            continue
        f.write("%s\t%s\n" % (outNodes[i], i))
    for i in useNodes:
        if i not in outInteractions:
            continue
        for j in outInteractions[i].keys():
            if j not in useNodes:
                continue
            for k in re.split(";", outInteractions[i][j]):
                f.write("%s\t%s\t%s\n" % (i, j, k))
    f.close()

def rSIF(inf, typef = None, reverse = False):
    """read .sif"""
    readPathway = Pathway(dict(), dict())
    inNodes = dict()                            #Dictionary with (A : type)
    inInteractions = dict()                     #Dictionary with (A : (B : interaction))
    nodeMap = dict()
    if typef != None:
        nodeMap = r2Col(typef, delim = " = ", header = True)
    f = open(inf, "r")
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split("\s*\t\s*", line)
        if pline[0] not in inNodes:
            if pline[0] in nodeMap:
                inNodes[pline[0]] = nodeMap[pline[0]]
            else:
                inNodes[pline[0]] = "concept"
        if pline[2] not in inNodes:
            if pline[2] in nodeMap:
                inNodes[pline[2]] = nodeMap[pline[2]]
            else:
                inNodes[pline[2]] = "concept"
        if reverse:
            if pline[2] not in inInteractions:
                inInteractions[pline[2]] = dict()
            if pline[0] not in inInteractions[pline[2]]: 
                inInteractions[pline[2]][pline[0]] = pline[1]
            else:
                inInteractions[pline[2]][pline[0]] += ";"+pline[1]
        else:
            if pline[0] not in inInteractions:
                inInteractions[pline[0]] = dict()
            if pline[2] not in inInteractions[pline[0]]:
                inInteractions[pline[0]][pline[2]] = pline[1]
            else:
                inInteractions[pline[0]][pline[2]] += ";"+pline[1]
    f.close()
    return(inNodes, inInteractions)

def wSIF(writeFile, writeInteractions, useNodes = None):
    """write .sif"""
    f = open(writeFile, "w")
    if useNodes == None:
        for i in writeInteractions.keys():
            for j in writeInteractions[i].keys():
                for k in re.split(";", writeInteractions[i][j]):
                    f.write("%s\t%s\t%s\n" % (i, k, j))
    else:
        for i in useNodes:
            if i not in writeInteractions:
                continue
            for j in writeInteractions[i].keys():
                if j not in useNodes:
                    continue
                for k in re.split(";", writeInteractions[i][j]):
                    f.write("%s\t%s\t%s\n" % (i, k, j))
    f.close()

def reverseInteractions(interactions):
    """reverse interaction mapping"""
    rinteractions = {}
    for source in interactions.keys():
        for target in interactions[source].keys():
            if target not in rinteractions:
                rinteractions[target] = {}
            rinteractions[target][source] = interactions[source][target]
    return(rinteractions)

def getComponentMap(pNodes, pInteractions):
    """create the dictionary componentMap from interaction map"""
    rpInteractions = reverseInteractions(pInteractions)
    componentMap = dict()
    for i in pNodes.keys():
        if pNodes[i] != "complex":
            continue
        componentMap[i] = []
        if i not in rpInteractions:
            continue
        for j in rpInteractions[i]:
            if rpInteractions[i][j] == "component>":
                componentMap[i].append(j)
    return(componentMap)

def constructInteractions(nodeList, refNodes, refInteractions):
    """select concepts from list and construct Pathway"""
    outPathway = Pathway({}, {})
    for i in nodeList:
        outPathway.nodes[i] = refNodes[i]
        if i in refInteractions:
            for j in refInteractions[i].keys():
                if j in nodeList:
                    if i not in outPathway.interactions:
                        outPathway.interactions[i] = dict()
                    outPathway.interactions[i][j] = refInteractions[i][j]
    return(outPathway.nodes, outPathway.interactions)

def combinePathways(currentPathway, appendPathway, exclude = []):
    """combine Pathways"""
    for source in appendPathway.interactions.keys():
        if source in exclude:
            continue
        if source not in currentPathway.nodes:
            currentPathway.nodes[source] = appendPathway.nodes[source]
        for target in appendPathway.interactions[source].keys():
            if target in exclude:
                continue
            if target not in currentPathway.nodes:
                currentPathway.nodes[target] = appendPathway.nodes[target]
            if source not in currentPathway.interactions:
                currentPathway.interactions[source] = dict()
            if target not in currentPathway.interactions[source]:
                currentPathway.interactions[source][target] = appendPathway.interactions[source][target]
    return(currentPathway)

def sortConnected(inPathway):
    rinteractions = reverseInteractions(inPathway.interactions)
    index = 1
    mapNets = dict()
    sortedNets = []
    seenNodes = set()
    for i in inPathway.nodes.keys():
        if i in seenNodes:
            continue
        borderNodes = [i]
        currentNet = [i]
        while len(borderNodes) > 0:
            if borderNodes[0] in rinteractions:
                for j in rinteractions[borderNodes[0]].keys():
                    if j not in seenNodes:
                        seenNodes.update([j])
                        borderNodes.append(j)
                        currentNet.append(j)
            if borderNodes[0] in inPathway.interactions:
                for j in inPathway.interactions[borderNodes[0]].keys():
                    if j not in seenNodes:
                        seenNodes.update([j])
                        borderNodes.append(j)
                        currentNet.append(j)
            borderNodes.pop(0)
        if ("__DISCONNECTED__" not in currentNet):
            mapNets[index] = deepcopy(currentNet)
            index += 1
    indexList = mapNets.keys()
    netScore = dict()
    for i in indexList:
        netScore[i] = len(mapNets[i])
    indexList.sort(lambda x, y: cmp(netScore[y], netScore[x]))
    for i in indexList:
        sortedNets.append(mapNets[i])
    return(sortedNets)

def flattenPathway(inPathway):
    """expands complexes into their respective gene components"""
    allowedNodes = ["abstract", "family", "miRNA", "protein", "rna"]
    outPathway = Pathway({}, {})
    ## read and search componentMap for protein components
    componentMap = getComponentMap(inPathway.nodes, inPathway.interactions)
    for entity in componentMap.keys():
        seenNodes = set()
        elements = []
        expand = deepcopy(componentMap[entity])
        while len(expand) > 0:
            if expand[0] in seenNodes:
                expand.pop(0)
                continue
            seenNodes.update([expand[0]])
            if inPathway.nodes[expand[0]] == "protein":
                elements.append(expand[0])
            elif expand[0] in componentMap:
                expand += deepcopy(componentMap[expand[0]])
            expand.pop(0)
        componentMap[entity] = elements
    ## iterate over all interactions
    for source in inPathway.interactions.keys(): 
        for target in inPathway.interactions[source].keys():
            ## update interactions map
            if inPathway.nodes[source] in allowedNodes:
                if inPathway.nodes[target] in allowedNodes:
                    if inPathway.interactions[source][target] not in ["component>", "member>"]:
                        if source not in outPathway.nodes:
                            outPathway.nodes[source] = inPathway.nodes[source]
                        if target not in outPathway.nodes:
                            outPathway.nodes[target] = inPathway.nodes[target]
                        if source not in outPathway.interactions:
                            outPathway.interactions[source] = {}
                        outPathway.interactions[source][target] = inPathway.interactions[source][target]
                elif target in componentMap:
                    for element in componentMap[target]:
                        if source != element:
                            if source not in outPathway.nodes:
                                outPathway.nodes[source] = inPathway.nodes[source]
                            if element not in outPathway.nodes:
                                outPathway.nodes[element] = inPathway.nodes[element]
                            if source not in outPathway.interactions:
                                outPathway.interactions[source] = {}
                            if inPathway.interactions[source][target] == "component>":
                                outPathway.interactions[source][element] = "-a>"
                            else:
                                outPathway.interactions[source][element] = inPathway.interactions[source][target]
            elif source in componentMap:
                if inPathway.nodes[target] in allowedNodes:
                    for element in componentMap[source]:
                        if element not in outPathway.nodes:
                            outPathway.nodes[element] = inPathway.nodes[element]
                        if target not in outPathway.nodes:
                            outPathway.nodes[target] = inPathway.nodes[target]
                        if element not in outPathway.interactions:
                            outPathway.interactions[element] = {}
                        outPathway.interactions[element][target] = inPathway.interactions[source][target]
                elif target in componentMap:
                    continue
    return(outPathway)

def getDownstream(node, distance, interactions):
    """returns downstream neighbors of distance"""
    seenNodes = set([node])
    borderNodes = [node]
    frontierNodes = []
    for dist in range(distance):
        while len(borderNodes) > 0:
            currNode = borderNodes.pop()
            if currNode in interactions:
                for i in interactions[currNode].keys():
                    if i not in seenNodes:
                        seenNodes.update([i])
                        frontierNodes.append(i)
        borderNodes = deepcopy(frontierNodes)
        frontierNodes = list()
    return(list(seenNodes))

def getUpstream(node, distance, interactions):
    """returns downstream neighbors of distance"""
    rinteractions = reverseInteractions(interactions)
    seenNodes = set([node])
    borderNodes = [node]
    frontierNodes = []
    for dist in range(distance):
        while len(borderNodes) > 0:
            currNode = borderNodes.pop()
            if currNode in rinteractions:
                for i in rinteractions[currNode].keys():
                    if i not in seenNodes:
                        seenNodes.update([i])
                        frontierNodes.append(i)
        borderNodes = deepcopy(frontierNodes)
        frontierNodes = []
    return(list(seenNodes))

def getNeighbors(node, distance, interactions):
    """returns upstream and downstream neighbors of distance"""
    rinteractions = reverseInteractions(interactions)
    seenNodes = set([node])
    borderNodes = [node]
    frontierNodes = []
    for d in range(distance):
        while len(borderNodes) > 0:
            currNode = borderNodes.pop()
            if currNode in interactions:
                for i in interactions[currNode].keys():
                    if i not in seenNodes:
                        seenNodes.update([i])
                        frontierNodes.append(i)
            if currNode in rinteractions:
                for i in rinteractions[currNode].keys():
                    if i not in seenNodes:
                        seenNodes.update([i])
                        frontierNodes.append(i)
        borderNodes = deepcopy(frontierNodes)
        frontierNodes = []
    return(seenNodes)
    
def wNodeAttributes(pNodes, scoreMap = None, directory = "."):
    """write cytoscape node attribute files"""
    ## TYPE.NA
    typef = open("%s/TYPE.NA" % (directory), "w")
    typef.write("TYPE (class=java.lang.String)\n")
    for node in pNodes.keys():
        typef.write("%s = %s\n" % (node, pNodes[node]))
    typef.close()
    ## LABEL.NA
    labelf = open("%s/LABEL.NA" % (directory), "w")
    labelf.write("LABEL (class=java.lang.String)\n")
    for node in pNodes.keys():
        if pNodes[node] == "protein":
            labelf.write("%s = %s\n" % (node, node))
        else:
            labelf.write("%s = %s\n" % (node, ""))
    labelf.close()
    ## *_SCORE.NA
    if scoreMap != None:
        for element in scoreMap.keys():
            scoref = open("%s/%s_SCORE.NA" % (directory, element), "w")
            scoref.write("SCORE (class=java.lang.Double)\n")
            for node in scoreMap[element].keys():
                scoref.write("%s = %s\n" % (node, scoreMap[element][node]))
            scoref.close()
   
def removeNode(node, pNodes, pInteractions):
    """remove a node and its interactions from a Pathway"""
    rpInteractions = reverseInteractions(pInteractions)
    del pNodes[node]
    if node in pInteractions:
        for element in pInteractions[node].keys():
            del pInteractions[node][element]
            if len(pInteractions[node].keys()) == 0:
                del pInteractions[node]
            del rpInteractions[element][node]
            if len(rpInteractions[element].keys()) == 0:
                del rpInteractions[element]
    if node in rpInteractions:
        for element in rpInteractions[node].keys():
            del pInteractions[element][node]
            if len(pInteractions[element].keys()) == 0:
                del pInteractions[element]
            del rpInteractions[node][element]
            if len(rpInteractions[node].keys()) == 0:
                del rpInteractions[node]
    return(pNodes, pInteractions)

def isTerminal(node, interactions):
    rinteractions = reverseInteractions(interactions)
    nLinks = 0
    if node in interactions:
        nLinks += len(interactions[node].keys())
    if node in rinteractions:
        nLinks += len(rinteractions[node].keys())
    if nLinks > 1:
        return False
    else:
        return True

def numChildren(node, interactions):
    if node in interactions:
        return(len(interactions[node].keys()))
    else:
        return(0)

def isTranscriptional(node, interactions):
    rinteractions = reverseInteractions(interactions)
    if node in rinteractions:
        for target in rinteractions[node].keys():
            for type in re.split(";", rinteractions[node][target]):
                if type.startswith("-t"):
                    return(True)
    return(False)
    
def hasTranscriptional(node, interactions, retTargets = False):
    transcriptionalTargets = set()
    seenNodes = set([node])
    unresolvedNodes = [node]
    while len(unresolvedNodes) > 0:
        currNode = unresolvedNodes.pop(0)
        if currNode in interactions:
            for target in interactions[currNode].keys():
                if target in seenNodes:
                    continue
                seenNodes.update([target])
                unresolvedNodes.append(target)
                for type in re.split(";", interactions[currNode][target]):
                    if type.startswith("-t"):
                        transcriptionalTargets.update([target])
    transcriptionalTargets = list(transcriptionalTargets)
    if len(transcriptionalTargets) > 0:
        if retTargets:
            return(transcriptionalTargets)
        else:
            return(True)
    else:
        if retTargets:
            return([])
        else:
            return(False)

def getListIndices(inItem, inList):
    """returns indices of the occurence of inItem in inList"""
    indices = []
    for i, item in enumerate(inList):
        if item == inItem:
            indices.append(i)
    return(indices)

def filterComplexesByGeneSupport(allNodes, forInteractions, revInteractions, typeMap, componentMap, threshold = 0.5):
    """remove complexes by percent support"""
    ## mark complexes in allNodes as unvisitedComplexes
    unvisitedComplexes = set()
    for i in allNodes:
        if allNodes[i] == "complex":
            unvisitedComplexes.update([i])
    def keepMajority(complex, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions, typeMap, componentMap, threshold = 0.5):
        unvisitedComplexes.discard(complex)
        complexSubunits = []
        otherSubunits = []
        for i in componentMap[complex]:
            if typeMap[i] == "complex":
                complexSubunits.append(i)
            else:
                otherSubunits.append(i)
        n = len(complexSubunits)+len(otherSubunits)
        complexSubunitsInNet = []
        otherSubunitsInNet = []
        for i in complexSubunits:
            if i in allNodes:
                if i in recursedComplexes:
                    indices = getListIndices(i, recursedComplexes)
                    log("WARNING: complex loop %s\n" % (recursedComplexes[indices[-1]:]))
                    continue
                recursedComplexes.append(i)
                (logical, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions) = keepMajority(i, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions, typeMap, componentMap, threshold = threshold)
                if logical:
                    complexSubunitsInNet.append(i)
        for i in otherSubunits:
            if i in allNodes:
                otherSubunitsInNet.append(i)
        m = len(complexSubunitsInNet)+len(otherSubunitsInNet)
        if m <= threshold*n:
            (allNodes, forInteractions, revInteractions) = removeNode(complex, allNodes, forInteractions)
            return(False, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions)
        else:
            return(True, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions)

    ## visit complexes while there are still unvisitedComplexes
    while len(unvisitedComplexes) > 0:
        complex = list(unvisitedComplexes)[0]
        recursedComplexes = [complex]
        (logical, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions) = keepMajority(complex, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions, typeMap, componentMap, threshold = threshold)
    return(allNodes, forInteractions)
    
def shortestPath(source, target, interactions):
    shortestPaths = []
    allPaths = [[source]]
    while len(shortestPaths) == 0:
        nextPaths = []
        while len(allPaths) > 0:
            currentPath = allPaths.pop()
            if currentPath[-1] not in interactions:
                continue
            for intermediate in interactions[currentPath[-1]]:
                if intermediate in currentPath:
                    continue
                if intermediate == target:
                    shortestPaths.append(currentPath + [intermediate])
                    nextPaths.append(currentPath + [intermediate])
                else:
                    nextPaths.append(currentPath + [intermediate])
        assert(len(nextPaths) > 0)
        allPaths = deepcopy(nextPaths)
    return (shortestPaths)

def rList(inf, header = False):
    """read 1 column list"""
    inList = []
    f = open(inf, "r")
    if header:
        f.readline()
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\t\r\n")
        inList.append(line)
    f.close()
    return(inList)

def r2Col(inf, appendData = {}, delim = "\t", null = "NA", header = False):
    """read 2 column data"""
    inData = deepcopy(appendData)
    f = open(inf, "r")
    if header:
        line = f.readline()
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        if len(pline[1]) == 0:
            pline[1] = null
        if len(pline) != 2:
            log("ERROR: Length of data line is not 2\n", die = True)
        inData[pline[0]] = pline[1]
    f.close()
    return(inData)
