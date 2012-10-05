#!/usr/bin/env python

import sys
import array
from math import isnan, sqrt


def median(inList):
    """calculates median"""
    cList = copy(inList)
    cList.sort()
    if len(cList) == 0:
        median = float("nan")
    else:
        if len(cList)%2 == 1:
            median = cList[len(cList)/2]
        else:
            median = (cList[len(cList)/2]+cList[(len(cList)/2)-1])/2.0
    return (median)

class FloatMatrix:
    def __init__(self):
        self.corner_name = "probe"
        self.data = None
        self.nrows = None
        self.ncols = None
        self.rowmap = None
        self.colmap = None

    def read(self, handle):
        header = None
        for line in handle:
            row = line.rstrip().split("\t")
            if header is None:
                header = row
                self.data = array.array("f")
                self.colmap = {}
                self.rowmap = {}
                self.ncols = len(row) - 1
                self.nrows = 0
                for i, c in enumerate(row[1:]):
                    self.colmap[c] = i
            else:
                if len(row) - 1 != self.ncols:
                    raise DataException("Misformed matrix")
                self.rowmap[row[0]] = len(self.rowmap)
                a = []
                for v in row[1:]:
                    try:
                        a.append(float(v))
                    except ValueError:
                        a.append(float('Nan'))
                self.data.extend(a)
                self.nrows += 1

    def init_blank(self, rows, cols):
        self.data = array.array("f")
        self.colmap = {}
        for i,c in enumerate(cols):
            self.colmap[c] = i
        self.rowmap = {}
        for i,r in enumerate(rows):
            self.rowmap[r] = i
        self.ncols = len(cols)
        self.nrows = len(rows)
        for i in range(self.nrows):
            self.data.extend([float('nan')] * self.ncols)

    def get_value(self, row_name, col_name):
        return self.data[ self.rowmap[row_name] * self.ncols + self.colmap[col_name] ]

    def set_value(self, row_name, col_name, value):
        self.data[ self.rowmap[row_name] * self.ncols + self.colmap[col_name] ] = value
    
    def get_row(self, row_name):
        return self.data[ self.rowmap[row_name] * self.ncols :  (self.rowmap[row_name]+1) * self.ncols ]

    def get_cols(self):
        out = self.colmap.keys()
        return sorted(out, key=self.colmap.get)
    
    def has_row(self, row):
        return row in self.rowmap 

    def has_col(self, col):
        return col in self.colmap 

    def get_rows(self):
        out = self.rowmap.keys()
        return sorted(out, key=self.rowmap.get)
    
    def write(self, handle, missing='NA'):
        write = csv.writer(handle, delimiter="\t", lineterminator='\n')
        col_list = self.get_cols()
        
        write.writerow([self.corner_name] + col_list)
        for rowName in self.rowmap:
            out = [rowName]
            row = self.get_row(rowName)
            for col in col_list:
                val = row[self.colmap[col]]
                if val is None or math.isnan(val):
                    val = missing
                else:
                    val = "%.5f" % (val)
                out.append(val)
            write.writerow(out)        
            
    def toRmatrix(self, r):
        out = r.matrix(self.data, ncol=self.ncols, dimnames=[ self.get_rows(), self.get_cols() ], byrow=True)
        return out

    def row_median_shift(self, normal_list=None):
        if normal_list is None:
            """
            Shift all values in matrix by median value of row (probe)
            """
            out = FloatMatrix()
            out.init_blank(cols=self.get_cols(), rows=self.get_rows())
            
            for row in self.get_rows():
                medianVal = median(self.get_row(row))
                if not isnan(medianVal):
                    for col in self.get_cols():
                        val = self.get_value(col_name=col, row_name=row)
                        if not isnan(val):  
                            out.set_value(col_name=col, row_name=row, value=val-medianVal)
            return out
        else:            
            """
            Shift all non-normal values by median value for  in matrix by median value of row (probe)
            """
            out = FloatMatrix()
            out.init_blank(cols=self.get_cols(), rows=self.get_rows())

            col_list = self.get_cols()
            normal_cols = []    
            for norm in normal_list:
                if norm in col_list:
                    normal_cols.append(col_list.index(norm))
            
            if len(normal_cols) == 0:
                raise TypeError("No normals found in matrix")
            
            for row in self.get_rows():
                r = self.get_row(row)
                normal_row = []
                for i in normal_cols:
                    normal_row.append(r[i])
                medianVal = median(normal_row)
                if medianVal == "NA":
                    medianVal = 0
                
                for col in self.get_cols():
                    val = self.get_value(col_name=col, row_name=row)
                    if not isnan(val):  
                        out.set_value(col_name=col, row_name=row, value=val-medianVal)
            return out

    def row_average_shift(self):
        out = FloatMatrix()
        out.init_blank(cols=self.get_cols(), rows=self.get_rows())
        
        for row in self.get_rows():
            s = 0.0
            i = 0.0
            for v in self.get_row(row):
                if not isnan(v):
                    s += v
                    i += 1.0
            
            if i == 0:
                aveVal = 0
            else:
                aveVal = s / i
            
            for col in self.get_cols():
                val = self.get_value(col_name=col, row_name=row)
                if not isnan(val):  
                    out.set_value(col_name=col, row_name=row, value=val-aveVal)
        return out
    
    def row_std_normalize(self):
        out = FloatMatrix()
        out.init_blank(cols=self.get_cols(), rows=self.get_rows())
        
        for row in self.get_rows():
            s = 0.0
            s2 = 0.0
            i = 0.0
            for v in self.get_row(row):
                if not isnan(v):
                    s += v
                    s2 += v * v
                    i += 1.0
            
            if i == 0:
                std_val = 1.0
            else:
                ave_val = s / i
                std_val = sqrt((s2 / i) - (ave_val * ave_val))
                if std_val == 0.0:
                    std_val = 1.0
            
            for col in self.get_cols():
                val = self.get_value(col_name=col, row_name=row)
                if not isnan(val):  
                    out.set_value(col_name=col, row_name=row, value=val / std_val )
        return out
            

def isect(a, b):
    out = []
    for v in a:
        if v in b:
            out.append(v)
    return out

def uniq(e):
    out = {}
    for v in e:
        out[v] = True
    return out.keys()

class MatrixTupleIter:
    
    def __init__(self, matrix1, matrix2, emit_names=False):
        self.matrix1 = matrix1
        self.matrix2 = matrix2
        self.emit_names = emit_names
    
    def __iter__(self):
        for row in self.matrix1.get_rows():
            if self.matrix2.has_row(row):
                for col in self.matrix1.get_cols():
                    if self.matrix2.has_col(col):
                        if self.emit_names:
                            yield (row, col, self.matrix1.get_value(row_name=row, col_name=col), self.matrix2.get_value(row_name=row, col_name=col))
                        else:
                            yield (self.matrix1.get_value(row_name=row, col_name=col), self.matrix2.get_value(row_name=row, col_name=col))


def pearsonr(vals):
    
    n = 0.0
    sum_x = 0.0
    sum_y = 0.0
    sum_x_sq = 0.0
    sum_y_sq = 0.0
    psum = 0.0
    
    for x,y in vals:
        if not isnan(x) and not isnan(y):
            n += 1
            sum_x += x
            sum_y += y
            sum_x_sq += x * x
            sum_y_sq += y * y
            psum += x * y
    num = psum - (sum_x * sum_y / n)
    den = pow((sum_x_sq - pow(sum_x,2) / n) * (sum_y_sq - pow(sum_y,2) / n), 0.5)
    if den == 0: return 0
    return num / den


def maxdiffs(vals, top_count=10):
    out = []
    for row, col, x, y in vals:
        diff = abs(x-y)
        if len(out) < top_count:
            out.append((row, col, diff))
        elif diff > out[-1][2]:
            out[-1] = (row, col, diff)
            out.sort(key=lambda x:x[2], reverse=True )
    return out

if __name__ == "__main__":
    handle = open(sys.argv[1])
    matrix_1 = FloatMatrix()
    matrix_1.read(handle)
    handle.close()
    
    handle = open(sys.argv[2])
    matrix_2 = FloatMatrix()
    matrix_2.read(handle)
    handle.close()
    
    col_isect_size = len(isect(matrix_1.get_cols(), matrix_2.get_cols()))
    row_isect_size = len(isect(matrix_1.get_rows(), matrix_2.get_rows()))
    print "Sample Intersection size", float(col_isect_size) / float(len(matrix_1.get_cols())), float(col_isect_size) / float(len(matrix_2.get_cols()))
    print "Probe Intersection size", float(row_isect_size) / float(len(matrix_1.get_rows())), float(row_isect_size) / float(len(matrix_2.get_rows()))
    
    count = 0.0
    x_nan_count = 0.0
    y_nan_count = 0.0
    for x, y in MatrixTupleIter(matrix_1, matrix_2):
        if isnan(x):
            x_nan_count += 1.0
        if isnan(y):
            y_nan_count += 1.0
        count += 1.0
    print "Missing data: %f %f" % ( x_nan_count / count, y_nan_count / count )
        
    print "Pearson Corrleation Coeff", pearsonr( MatrixTupleIter(matrix_1, matrix_2))
    
    ave_matrix_1 = matrix_1.row_average_shift()
    ave_matrix_2 = matrix_2.row_average_shift()
    
    print "Average Shift Pearson Corrleation Coeff", pearsonr( MatrixTupleIter(ave_matrix_1, ave_matrix_2))
    
    print "Zscore Shift Pearson Corrleation Coeff", pearsonr( MatrixTupleIter(ave_matrix_1.row_std_normalize(), ave_matrix_2.row_std_normalize()))
    
    print "Max Differences:"
    for v in maxdiffs( MatrixTupleIter(matrix_1, matrix_2, True)):
        if v[2] > 0.0:
            print "", v[0], v[1], v[2] 
    
    
    
