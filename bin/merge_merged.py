#!/usr/bin/env python

import os
import sys
import fnmatch
#from hgSQL import hgSQL

resultMatrix = {} # first key is entity, then sample
sampleList = []

def usage():
	print "Usage: "+sys.argv[0]+" db merged_dir outputfile"
	sys.exit(0)

def getFilesMatching(baseDir, patterns):
	list = []
	
	for root, dirs, files in os.walk(baseDir):
		for file in files:
			ptr = os.path.join(root, file)
			for pattern in patterns:
				if fnmatch.fnmatch(ptr, pattern):
					list.append(ptr)
	return list

def addFileToResult(file):
	pid = file[:-4].split("_").pop()
	fh = open(file,"r")
	header = fh.readline().strip("\n").split("\t")
	sampleOrder = header[1:]
	#print sampleOrder
	for sample in sampleOrder:
		if sample.startswith("na_") == False  and sample.split(" ").pop(0) not in sampleList:
			sampleList.append(sample.split(" ").pop(0))
	for line in fh:
		dataA = line.strip("\n").split("\t")
		entity = dataA.pop(0)
		resultMatrix[pid+"_"+entity] = {}
		for i in range(len(dataA)):
			sample = sampleOrder[i] ## sng ##
			if sample.startswith("na_") == False:
				resultMatrix[pid+"_"+entity][sampleOrder[i].split(" ").pop(0)] = dataA[i]
	fh.close()

def main(db,directory,outputFile):
	files = getFilesMatching(directory, ["*_transpose_*"])

	print "Loading data",
	sys.stdout.flush()
	for f in files:
		pid = f[:-4].split("_").pop()
		if pid == "example":
			continue
		addFileToResult(f)
		print ".",
		sys.stdout.flush()
		
	#print "Done reading data, converting sample ids...",
	sys.stdout.flush()	
	# open the sql connection
	#sql = hgSQL()
	#connection = sql.connect("localDb","bioIntTCGAOV")
	#connection = sql.connect("localDb",db)
	#transSampleList = []
	#for i in range(len(sampleList)):
	#	(samplePrefix,sampleID) = sampleList[i].split("_")
	#	if samplePrefix == "sample":
	#		query = "select name from samples where id = "+sampleID
	#	elif samplePrefix =="patient":
	#		query = "select patient_name as name from samples where patient_id = "+sampleID
	#	else:
	#		print "Error detecting sample type!  Make sure the header starts with sample_ or patient_"
	#		sys.exit(1)
	#		
	#	result = connection.execute(query)
	#	row = result.fetchone()
	#	transSampleList.append(row['name'])
	## close the sql connections cause we're done now
	#connection.close()
	#print "done."
	print "Printing Results..."
	resultFile = open(outputFile,"w")
	resultFile.write("pid_entity\t")
	resultFile.write("\t".join(sampleList)+"\n")
	for entity in resultMatrix:
		resultFile.write(entity)
		for sample in sampleList:
			resultFile.write("\t")
			#print sample
			#print resultMatrix[entity]
			#sys.exit(0)
			if sample in resultMatrix[entity]:
				resultFile.write(resultMatrix[entity][sample])
		resultFile.write("\n")
				
	resultFile.close()
	

if __name__ == "__main__":
	if len(sys.argv) != 4:
		usage()

	db = sys.argv[1]
	directory = sys.argv[2]
	outputFile = sys.argv[3]
	main(db,directory,outputFile)

