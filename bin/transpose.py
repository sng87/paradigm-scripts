#!/usr/bin/env python
import string,sys

if (len(sys.argv))!=3:
    print "python transpose.py extractDataIn transposeOut-Paradigm"
    print "       the output file skips the probe column"
    sys.exit()

fin= open(sys.argv[1],'r')
fout= open(sys.argv[2],'w')

matrix=[]
for line in fin.readlines():
    data = string.split(line.strip(),'\t')
    matrix.append(data)

row=len(matrix)
col= len(matrix[0])

#header
for i in range (0, row-1):
    fout.write(matrix[i][0]+'\t')
fout.write(matrix[i+1][0]+'\n')

#body
for j in range (1,col):
    for i in range(0, row-1):
        fout.write(matrix[i][j]+'\t')
    fout.write(matrix[i+1][j]+'\n')

fin.close()
fout.close()
