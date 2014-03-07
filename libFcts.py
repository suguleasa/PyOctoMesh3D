
import SimpleITK as sitk
from globalVars import *
from math import sqrt, floor, copysign
from numpy import *
import scipy
import numpy

class Coordinate(object):
    def __init__(self,x=-1,y=-1,z=-1):
        self.x = x
        self.y = y
        self.z = z
        
def find_mid_point(p1, p2):
    xMid = int( (p1.x + p2.x)/2.0 )
    yMid = int( (p1.y + p2.y)/2.0 )
    zMid = int( (p1.z + p2.z)/2.0 )
    return Coordinate(xMid,yMid,zMid)


## Find distance between two points
def find_distance(p1, p2):
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    dz = p2.z - p1.z
    return ( sqrt(dx**2 + dy**2 + dz**2) )


def bresenham_line3d(p1,p2):
    # adapted from the Matlab code by Jimmy Shen
    # http://www.mathworks.com/matlabcentral/fileexchange/21057-3d-bresenhams-line-generation
    nx = abs(p2.x-p1.x)+1
    ny = abs(p2.y-p1.y)+1
    nz = abs(p2.z-p1.z)+1
    d = max(nx,ny,nz)
    
    X = scipy.zeros(d)
    Y = scipy.zeros(d)
    Z = scipy.zeros(d)

    x1 = p1.x
    y1 = p1.y
    z1 = p1.z

    x2 = p2.x
    y2 = p2.y
    z2 = p2.z

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1

    ax = abs(dx)*2
    ay = abs(dy)*2
    az = abs(dz)*2

    sx = copysign(1,dx)
    sy = copysign(1,dy)
    sz = copysign(1,dz)

    x = x1;
    y = y1;
    z = z1;
    idx = 0;

    if(ax>=max(ay,az)):           # x dominant
        yd = ay - ax/2;
        zd = az - ax/2;

        while(1):
            X[idx] = x;
            Y[idx] = y;
            Z[idx] = z;
            idx = idx + 1;

            if(x == x2)  :      # end
                break;

            if(yd >= 0) :       # move along y
                y = y + sy;
                yd = yd - ax;

            if(zd >= 0) :       # move along z
                z = z + sz;
                zd = zd - ax;

            x  = x  + sx;        # move along x
            yd = yd + ay;
            zd = zd + az;
        
    elif(ay>=max(ax,az)):        # y dominant
        xd = ax - ay/2
        zd = az - ay/2

        while(1):
            X[idx] = x
            Y[idx] = y
            Z[idx] = z
            idx = idx + 1

            if(y == y2) :       # end
                break;

            if(xd >= 0) :      # move along x
                x = x + sx
                xd = xd - ay

            if(zd >= 0) :       # move along z
                z = z + sz;
                zd = zd - ay

            y  = y  + sy        # move along y
            xd = xd + ax
            zd = zd + az
        
    elif(az>=max(ax,ay))  :     # z dominant
        xd = ax - az/2
        yd = ay - az/2

        while(1):
            X[idx] = x
            Y[idx] = y
            Z[idx] = z
            idx = idx + 1

            if(z == z2) :       # end
                break;

            if(xd >= 0)  :      # move along x
                x = x + sx;
                xd = xd - az

            if(yd >= 0)  :      # move along y
                y = y + sy
                yd = yd - az

            z  = z  + sz       # move along z
            xd = xd + ax
            yd = yd + ay
   
    return [X,Y,Z]

## This function draws a line between two points in space: pinit and pend
def draw_line(image,pinit,pend):

    [X,Y,Z] = bresenham_line3d(pinit,pend)
    for i in range(0,len(X)):
        image.SetPixel( int(X[i]), int(Y[i]), int(Z[i]), VAL)
 
## This function returns 1 if the pixel value is within a given
## bin with boundaries: [lowLim,highLim], and 0 otherwise
def in_bin_i(pxVal, lowLim, highLim):
    return (lowLim <= pxVal and pxVal < highLim)

## This function checks if two given pixels belong in the same bin
def is_in_same_bin(pxVal1, pxVal2):
    isHomogeneous = False

    for i in range(1,len(binBnd)):
        lowLim = binBnd[i-1]
        highLim = binBnd[i]
        isHomogeneous = isHomogeneous or ( in_bin_i(pxVal1, lowLim, highLim) and in_bin_i(pxVal2, lowLim, highLim))
    return isHomogeneous

        
## This function computes the probability that there is an inclusion
## of a different material inside thie element.
## Given that all 4 corners are found to be homogeneous, or to belong in to
## the same bin, then is there any inclusion inside this element that is
## larger than a certain area? If so, return true, otherwise return false
def has_inclusions(image,p1,p2,p3,p4,p5):

    pxVal1 = image.GetPixel(p1.x,p1.y,p1.z)
    
#    print ' p1,p2,p3,p4,p5'
#    print p1.x, p1.y, p1.z
#    print p2.x, p2.y, p2.z
#    print p3.x, p3.y, p3.z
#    print p4.x, p4.y, p4.z
#    print p5.x, p5.y, p5.z
    
    xHigh = p2.x 
    xLow = p1.x    
    yHigh = p4.y
    yLow = p1.y
    zHigh = max(p1.z, p5.z)
    zLow = min(p1.z, p5.z)
    
    areaElem = abs( (p4.y - p1.y) * (p2.x - p1.x) * (p1.z - p5.z))
    nr_samples = int( log(PROB) / log (abs(areaElem - AREA_INCLUSION)/areaElem) )
    
    for i in range (1,nr_samples):
        rx = random.randint(xLow,xHigh)
        ry = random.randint(yLow,yHigh)
        rz = random.randint(zLow,zHigh)
        samplePixel = image.GetPixel(rx,ry,rz)
        if (0 == is_in_same_bin(pxVal1,samplePixel)):
            return True

    return False

## This function finds the midpoint between two points in space: p1 and p2
def find_mid_point(p1, p2):
    xMid = int( (p1.x + p2.x)/2.0 )
    yMid = int( (p1.y + p2.y)/2.0 )
    zMid = int( (p1.z + p2.z)/2.0 )
    return Coordinate(xMid,yMid,zMid)


## This function checks if the four corners belong to the same
## pixel value (same bin that has a range of pixel values)
## check if the 4 corners are all in the same bin of colors or not
def eight_corners_test(pxVal1, pxVal2, pxVal3, pxVal4, pxVal5, pxVal6, pxVal7, pxVal8):

    isHomogeneous = False
    for i in range(len(binBnd)-1):
        lowLim = binBnd[i-1]
        highLim = binBnd[i]
        isHomogeneous = ( isHomogeneous or 
                         (
                          in_bin_i(pxVal1,lowLim,highLim) and in_bin_i(pxVal2,lowLim,highLim) 
                          and in_bin_i(pxVal3,lowLim,highLim) and in_bin_i(pxVal4,lowLim,highLim)
                          and in_bin_i(pxVal5,lowLim,highLim) and in_bin_i(pxVal6,lowLim,highLim) 
                          and in_bin_i(pxVal7,lowLim,highLim) and in_bin_i(pxVal8,lowLim,highLim)
                          )
                         )

    return isHomogeneous

## Check if the ends of a line are from the same bin
def ends_in_same_bin(image, p1, p2):
    
#    print p2.x, p2.y, p2.z
#    pxVal1 = image.GetPixel(int(p1.x), int(p1.y), int(p1.z));
#    pxVal2 = image.GetPixel(int(p2.x), int(p2.y), int(p2.z));
    pxVal1 = image.GetPixel(p1.x, p1.y, p1.z)
    pxVal2 = image.GetPixel(p2.x, p2.y, p2.z)   
        
    val1 = pxVal1 == pxVal2
    val2 = is_in_same_bin(pxVal1, pxVal2)
    
    return (val1 or val2)

def linear_search(image,bbegin,eend):
        
        list_nodes = []
        begin = Coordinate(bbegin.x,bbegin.y, bbegin.z)
        end = Coordinate(eend.x,eend.y, eend.z)
        
        old = Coordinate(bbegin.x,bbegin.y, bbegin.z)
        dist = find_distance(begin,end)

        if bbegin.x == eend.x and bbegin.y == eend.y and dist>2: # vertical line: z changes, x & y stay the same
            
            next = Coordinate(begin.x,begin.y, begin.z+1)
            while next.z <= end.z:
                if not(ends_in_same_bin(image,next,old)):
                    
                    list_nodes = list_nodes + [Coordinate(next.x,next.y, next.z)]
                old = Coordinate(next.x,next.y, next.z)
                next = Coordinate(next.x,next.y, next.z+1)
        
        else:
            if bbegin.y == eend.y and bbegin.z == eend.z and dist>2: # horizontal line: x changes, y & z stay the same
        
                next = Coordinate(begin.x+1,begin.y, begin.z)            
                
                while next.x <= end.x :
    
                    if not(ends_in_same_bin(image, next,old)):
                        list_nodes = list_nodes + [Coordinate(next.x,next.y, next.z)]
                        
                    old = Coordinate(next.x,next.y, next.z)
                    next = Coordinate(next.x+1,next.y, next.z)        
        
            else:
                if bbegin.x == eend.x and bbegin.z == eend.z and dist>2: # 90 degrees line: y changes, x & z stay the same
        
                    next = Coordinate(begin.x,begin.y+1, begin.z)            
                    
                    while next.y <= end.y :
        
                        if not(ends_in_same_bin(image, next,old)):
                            list_nodes = list_nodes + [Coordinate(next.x,next.y, next.z)]
                            
                        old = Coordinate(next.x,next.y, next.z)
                        next = Coordinate(next.x,next.y+1, next.z)   
#                else:
#                    print 'bbegin', bbegin.x, bbegin.y, bbegin.z
#                    print 'eend', eend.x, eend.y, eend.z     
                    
        return list(list_nodes)


#def non_linear_search(image,bbegin,eend):
#    
#    [X,Y,Z] = bresenham_line3d(bbegin,eend)
#    
#    list_nodes = []
#    dist = find_distance(bbegin, eend)
#    
#    if dist > 2: # if the end and beginning are more than 2 pixels apart
#        print 'len - ', len(X)
#        for i in range(0,len(X)-2):
#        
#            old = Coordinate( int(X[i]), int(Y[i]), int(Z[i]))
#            next = Coordinate( int(X[i])+1, int(Y[i])+1, int(Z[i])+1)
##            print old.x, old.y, old.z
##            print next.x, next.y, next.z
#            if not(ends_in_same_bin(image,next,old)) :                    
#
#                list_nodes = list_nodes + [Coordinate(next.x,next.y, next.z)]
#  
#    return list_nodes

def in_child_k(cubes,L):
    p1,p2,p3,p4,p5,p6,p7,p8 = cubes
    if p1.x <= L.x and L.x <= p2.x and p1.y <= L.y and L.y <= p4.y and p1.z <= L.z and L.z<p5.z:
        return True
    
    return False

def calc_plane_res(x, y, z):
    # 1 = a*x + b*y + c*z
    a = numpy.column_stack((x, y, z))
    (coeffs, resid,rank, sing_vals) = numpy.linalg.lstsq(a, numpy.ones_like(x))
    return resid

def draw_plane_connections(image, l1,l2,l3,l4, L1,L2,L3,L4):
#    print l1, l2, l3 ,l4
#    print L1, L2, L3, L4
    
    if l1 == 0 and l2 == 0:
        draw_line(image,L1,L2)
    if l1 == 0 and l3 == 0:
        draw_line(image,L1,L3)
    if l1 == 0 and l4 == 0:
#        print L1, len(L1)
#        print L4, len(L4)
#        print L1.x, L1.y, L1.z
#        print L4.x, L4.y, L4.z
        draw_line(image,L1,L4)
    if l2 == 0 and l3 == 0:
        draw_line(image,L2,L3)
    if l2 == 0 and l4 == 0:
#        print L2, len(L2)
#        print L4
#        print L2.x, L2.y, L2.z
#        print L4.x, L4.y, L4.z       
        draw_line(image,L2,L4)
    if l3 == 0 and l4 == 0:
        draw_line(image,L3,L4)
    
    
    