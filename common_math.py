#!/usr/bin/env python
#home-made math tool
#accept [ [x1,y11,y12..].[x2,y21,y22..]..] data or [y1,y2,y3...]
#some function only accept two-column data
import copy,math
from math import *
import list_utils as lu

def f_Data_Sum(data,i1,filter=None):
    '''
    Get sum of one column
    :param i1: the column index ( first is 0). If data is just a value then use -1.
    :param filter: ignore all data which filter(row) == True:
    '''
    if ( not isinstance(data,list)):
        raise TypeError,"Data must be list type."

    result = 0
    if ( i1 !=-1):
        for row in data:
            if ( filter != None ):
                if ( not filter(row) ):
                    continue
            result += row[i1]
    else:
        for row in data:
            if ( filter != None ):
                if ( not filter(row) ):
                    continue
            result += row

    return result
   

def f_Data_SumProduct(data,i1,i2,filter=None):
    '''
    Get sum product of two column data
    :param filter: ignore all data which filter(row) == True:
    '''
    if ( not isinstance(data,list)):
        raise TypeError,"Data must be list type."

    result = 0
    for row in data:
        if ( filter != None ):
            if ( not filter(row) ):
                continue
        result += row[i1]*row[i2]
    return result

def f_Data_AvgAndMeanSquareError(data,i1=-1,filter=None):
    '''
    Get mean square error of one column data
    :param i1: the column index of each row; if -1, then assume the row is directly a data ( not data list of [x,y1,y2..]
    :param filter: ignore all data which filter(row) == True:
    '''
    if ( not isinstance(data,list)):
        raise TypeError,"Data must be list type"
    if ( i1 != -1 and not isinstance(data[0],list) ):
        raise TypeError,"Column index must be used with data row but not single data value"


    nCount = 0
    dSum = 0
    dSumSquare = 0
    
    if ( i1 == -1):
        for row in data:
            if ( filter != None ):
                if ( not filter(row) ):
                    continue
            nCount += 1
            dSum += row
            dSumSquare += row**2
    else:
        for row in data:
            if ( filter != None ):
                if ( not filter(row) ):
                    continue
            nCount += 1
            dSum += row[i1]
            dSumSquare += row[i1]**2

    return dSum*1.0/nCount,(dSumSquare-dSum**2)*1.0/nCount

def f_Data_Interpolate(data,nTime,method=2):
    '''
    Interpolate data to specific times as original.
    Note, every point will be extended including first but except last
    Note, error of rectangle method is much larger than trapezodial method. Especially, if interpolate n points ( symmetric as a_x = a_{n-x} ) to 2n points ( or 4n, 6n...), the result will never be symmetric, but to 3n ( 5n, 7n ..) is exactly symmetric. For trapezodial method there is no such a problem.
    :param data: must be a list of two-element list
    :param method: 0 means midpoint-rectangle method, 1 means trapezodial method ( or linear spline ), 2means ^2 spline, 3 means ^3 spline ( cubic spline , not implented yet)
    '''
    if ( nTime == 1):
        return copy.deepcopy(data)
    result = []
    fInterval = data[1][0] - data[0][0]
    if ( method == 0): # expand [-x,a1] [0,a2] [x,a3]  to ... a1 [ 0-x/2,a2] ... [0+x/2,a2] a3... ( midpoint )
        for i in range(0,len(data)-1):
            d1 = data[i]
            d2 = data[i+1]
            for j in range(0,nTime / 2 + 1 ):
                result.append([d1[0] + fInterval / nTime * j,d1[1]])
            for j in range(nTime / 2 + 1, nTime):
                result.append([d1[0] + fInterval / nTime * j,d2[1]])
    elif ( method == 1): #expand linear ( trapezodial )
        i = 0
#        for i in range(0,len(data)-1):
        while i < len(data)-1:
            d1 = data[i]
            d2 = data[i+1]
            for j in range(0,nTime):
                result.append([d1[0] + fInterval / nTime * j,d2[1]*j/nTime + d1[1]*(nTime-j)/nTime])
            i += 1
        result.append( [ data[-1][0],data[-1][1] ] )
    elif ( method == 2): # expand with ^2
        #1st derivative is key parameter
        det0 = lu.f_Matrix_det([ [data[0][0]**2,data[0][0],1],[data[1][0]**2,data[1][0],1],[data[2][0]**2,data[2][0],1]])
        z2 = 2 * data[0][0] * lu.f_Matrix_det([ [data[0][1],data[0][0],1],[data[1][1],data[1][0],1],[data[2][1],data[2][0],1] ]) / det0 + lu.f_Matrix_det( [ [data[0][0]**2,data[0][1],1],[data[0][0]**2,data[0][1],1],[data[0][0]**2,data[0][1],1]]) / det0

        i = 0
        while ( i < len(data)-1 ):
        #for i in range(0,len(data)-1):
            d1 = data[i]
            d2 = data[i+1]
            z1 = z2
            z2 = -z1 + 2 * (d2[1]-d1[1])/(d2[0]-d1[0])

            for j in range(0,nTime):
                dx = fInterval / nTime *j
                result.append([d1[0] + dx, d1[1] + dx * z1 + (z2-z1)/2/(d2[0]-d1[0])*dx**2 ])
            i += 1
        result.append([data[-1][0],data[-1][1]])

    return result

def f_Data_Save(data,stFileName):
    '''
    Write a multi-column data to a file
    WARNING: NO overwrite prompt!
    '''
    f = open(stFileName,'w')
    for row in data:
        f.write(str(row)[1:-1].replace(" ",' '))
        f.write('\n')
    f.close()


class NumericFunction():
    '''
    Store a list in format [x,y] with asceding order of x. Distance between x[n] and x[n+1] for any x should be same to achieve best performance.
    Simplest piecewise function method: only the nearest point after the specific x is used. When the number of points is large enough it give nearly the same result to complex piecewise function method, and very easy to implented
    Trapz method: use line interpolation of x[n] and x[n+1] 
    '''
    fErr = 0.000001
    def __init__(self,data):
        if ( isinstance(data,str)):
            fIn = open(data,'r')
            data = [ [float(x) for x in y.split('\t')] for y in fIn.readlines()]
            fIn.close()
        self.data = copy.deepcopy(data)
        self.len = len(self.data)
        self.xmin = self.data[0][0]
        self.xmax = self.data[-1][0]
        self.step = self.data[1][0] - self.data[0][0] # only meaningful when data is uniform
#add additional line for convinience
        self.data.append(copy.deepcopy(self.data[-1]))
        self.data[-1][0] = 2*self.data[-2][0] - self.data[-3][0]

    def __Index__(self,x1,bLeft=True):
        '''
        use two-way search to find the max index with data[index].x <= x or data[index].x > x
        :param bLeft: if true ,return data <= x ; else return data > x
        '''
        #global to
        data = self.data
        #print(x1)
        #avoid numerical error at the border ( not used )
        #if ( x1 < self.xmin):
        #    if ( self.xmin-x1 < fErr):
        #        x1 = self.xmin
        #elif ( x1 > self.xmax):
        #    if ( x1-self.xmax < fErr):
        #        x1 = self.xmax

        if ( bLeft):
            nIndex = int(math.floor( (x1 - self.xmin) / self.step ))
        else:
            nIndex = int(math.ceil( (x1 - self.xmin) / self.step ))

        if ( nIndex > self.len -1  or nIndex < 0 ): #overflow
            #try to avoid numerical error caused exception ( not used ) 
            if ( abs(x1-self.xmin) < NumericFunction.fErr):
                nIndex = 0
            elif ( abs(x1-self.xmax) < NumericFunction.fErr):
                nIndex =  self.len-1
            else:
                #raise ValueError,"x ( %f ) is not inclulded in the data %f~%f" % (x1,self.xmin,self.xmax)
                raise ValueError,"x ( %-20.14f ) at index %d is not inclulded in the data %-20.14f~%-20.14f index 0~%d" % (x1,nIndex,self.xmin,self.xmax,self.len-1)
        return nIndex

    def Get(self,x1):
        '''
        Get the value at x
        '''
        data = self.data
        result = 0.0
        nIndex = self.__Index__(x1)
        #rectangle method
        #return data[nIndex][1]
        #Trapezodial rule
        return (data[nIndex][1] * ( data[nIndex+1][0]-x1) + data[nIndex+1][1]* ( x1 - data[nIndex][0]) ) / (data[nIndex+1][0] - data[nIndex][0])

    def Set(self,x,y):
        '''
        Set the value at x to y
        '''
        data = self.data
        result = 0.0
        data[self.__Index__(x)][1] = y

    def Change(self,x,y):
        '''
        Change the value at x to f(x) + y
        '''
        data = self.data
        result = 0.0
        data[self.__Index__(x)][1] += y

    def ChangeRange(self,x1,x2,func):
        '''
        Change all value in (x1,x2) by function
        :param func: a function accept x and return y
        '''
        data = self.data
        i1 = self.__Index__(x1,False) #ensure  xmin <= i1 <= i2 < xmax
        i2 = self.__Index__(x2,True) + 1 
        #print("Input range:",x1,x2 ,"Change range:",i1,i2-1,self.data[i1][0],self.data[i2-1][0])
        #print(x1,x2,i1,i2,data[i2-1][0]-12.2703222374)
#        for j in range(i,len(data)-1):
#            if ( x2 < data[i][0]):
#                break
        for j in range(i1,i2):
            #print(data[j][0])
            #if ( data[j][0] == 26.57297726 ):
            #    print("Del ",func(data[j][0]))
            data[j][1] += func(data[j][0])

    def Integrate(self,x1,x2,func):
        '''
        Integrate func(x) from x1 to x2 on a numerical defined function data.
        Assumed the data is thick enough. the last is directly discarded.
        Note: the border is not guaranteed to work!
        :param func: a function accept [x,y] parameter and return a value
        '''
        #global to
        #print("to:%f" %to)
        data = self.data
        result = 0.0
        #print(time.time()-to)
        i = self.__Index__(x1)
        #print(time.time()-to)
        i2 = self.__Index__(x2)+1
        #print(time.time()-to)
        #print(i,i2)
        #for j in range(i,len(data)-1):
        #    if ( x2 < data[i][0]):
        #        break
        for j in range(i+2,i2-1):
        #left rectangle method
            result += func(data[j]) * ( data[j+1][0] -data[j][0])
        #Add tail
        result += func(data[i+1]) * ( data[i+2][0] / 2 + data[i+1][0]/2 -x1)
        result += func(data[i2-1]) * ( x2 - data[i2-1][0] / 2 - data[i2-2][0] /2 )
        #right rectangle method
            #result += func(data[j+1]) * ( data[j+1][0] -data[j][0])
        #trapezodial method
            #result += (func(data[j]) + func(data[j+1])) * ( data[j+1][0] -data[j][0])/2
        #print(time.time()-to)

        return result

    def Write(self,stFileName):
        '''
        Save data in tab-seperated text
        '''
        fOut = open(stFileName,'w')
        for i in range(0,self.len):
            fOut.write("%16.10f\t%16.10f\n" % ( self.data[i][0],self.data[i][1]))
        fOut.close()


