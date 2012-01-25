#!/usr/bin/env python

import sys, os

args = list(sys.argv)
args.pop(0) # remove the script name

if args[0] == "-n":
	filterNA = True
	args.pop(0)
else:
	filterNA = False

filter = args[1]
(count,cutoff) = filter.split(",")
count = int(count)
cutoff = float(cutoff)

file = open(args[0])
header = file.readline()
sys.stdout.write(header)

headerA = header.split("\t")
headerA.pop(0)

for line in file:
	lineA = line.strip("\n").split("\t")
	lineA.pop(0)
	currCount = 0
	#for val in lineA:
	for i in range(len(lineA)):
		val = lineA[i]
		if filterNA and headerA[i][:3] == "na_":
			continue
		if abs(float(val)) >= cutoff:
			currCount += 1
		if currCount >= count:
			sys.stdout.write(line)
			sys.stdout.flush()
			break
file.close()


