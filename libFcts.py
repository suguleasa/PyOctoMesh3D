# Developer: Elena Caraba, copyright 2014

import SimpleITK as sitk
from globalVars import *
from math import sqrt, floor, copysign
from numpy import *
import scipy
import numpy
from bresenham import *
import types


D = {}
D['0'] = {
          'D': {'Quadrant':'2', 'Direction':'H'},
          'U': {'Quadrant':'2', 'Direction':'U'},          
          'R': {'Quadrant':'1', 'Direction':'H'},
          'L': {'Quadrant':'1', 'Direction':'L'},
          
          'B': {'Quadrant':'4', 'Direction':'H'},
          'F': {'Quadrant':'4', 'Direction':'F'}
          
          }

D['1'] = {
          'D': {'Quadrant':'3', 'Direction':'H'},
          'U': {'Quadrant':'3', 'Direction':'U'},          
          'R': {'Quadrant':'0', 'Direction':'R'},
          'L': {'Quadrant':'0', 'Direction':'H'},
                    
          'B': {'Quadrant':'5', 'Direction':'H'},
          'F': {'Quadrant':'5', 'Direction':'F'}

          }

D['2'] = {
          'D': {'Quadrant':'0', 'Direction':'D'},
          'U': {'Quadrant':'0', 'Direction':'H'},
          'R': {'Quadrant':'3', 'Direction':'H'},
          'L': {'Quadrant':'3', 'Direction':'L'},
          
          'B': {'Quadrant':'6', 'Direction':'H'},
          'F': {'Quadrant':'6', 'Direction':'F'}                    

          }

D['3'] = {
          'D': {'Quadrant':'1', 'Direction':'D'},
          'U': {'Quadrant':'1', 'Direction':'H'},
          'R': {'Quadrant':'2', 'Direction':'R'},
          'L': {'Quadrant':'2', 'Direction':'H'},

          'B': {'Quadrant':'7', 'Direction':'H'},
          'F': {'Quadrant':'7', 'Direction':'F'}        
                    
          }

D['4'] = {
          'D': {'Quadrant':'6', 'Direction':'H'},
          'U': {'Quadrant':'6', 'Direction':'U'},
          'R': {'Quadrant':'5', 'Direction':'H'},
          'L': {'Quadrant':'5', 'Direction':'L'},
                    
          'B': {'Quadrant':'0', 'Direction':'B'},
          'F': {'Quadrant':'0', 'Direction':'H'}        
                    
          }

D['5'] = {
          'D': {'Quadrant':'7', 'Direction':'H'},
          'U': {'Quadrant':'7', 'Direction':'U'},
          'R': {'Quadrant':'4', 'Direction':'R'},
          'L': {'Quadrant':'4', 'Direction':'H'},
                    
          'B': {'Quadrant':'1', 'Direction':'B'},
          'F': {'Quadrant':'1', 'Direction':'H'}        
                    
          }

D['6'] = {
          'D': {'Quadrant':'4', 'Direction':'D'},
          'U': {'Quadrant':'4', 'Direction':'H'},
          'R': {'Quadrant':'7', 'Direction':'H'},
          'L': {'Quadrant':'7', 'Direction':'L'},
          
          'B': {'Quadrant':'2', 'Direction':'B'},
          'F': {'Quadrant':'2', 'Direction':'H'}        
                    
          }

D['7'] = {
          'D': {'Quadrant':'5', 'Direction':'D'},
          'U': {'Quadrant':'5', 'Direction':'H'},
          'R': {'Quadrant':'6', 'Direction':'R'},
          'L': {'Quadrant':'6', 'Direction':'H'},
          
          'B': {'Quadrant':'3', 'Direction':'B'},
          'F': {'Quadrant':'3', 'Direction':'H'}        
                    
          }

class Coordinate(object):
    # 3D coordinate class
    def __init__(self,x=-1,y=-1,z=-1):
        self.x = x
        self.y = y
        self.z = z
        self.all = [x,y,z]

class Coordinate2D(object):
    # 2D coordinate class
    def __init__(self,x=-1,y=-1):
        self.x = x
        self.y = y
        self.all = [x,y]
        
            
def search_in(my_list,pi,pj,inImage):
    # search in a given list for coordinates
    Lk_list1 = [[x[0],x[1]] in [[pi,pj],] for x in my_list]  
    Lk_list2 = [[x[0],x[1]] in [[pj,pi],] for x in my_list] 
    if True in Lk_list1: #if we found it in the list:
        Lk_ind = Lk_list1.index(True)
        Lk = my_list[Lk_ind][2]
    else: 
        if True in Lk_list2:
            Lk_ind = Lk_list2.index(True)
            Lk = my_list[Lk_ind][2]
        else:
            Lk = linear_search(inImage,pi,pj)
            my_list.append([pi,pj,Lk])
    return Lk
                            
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

#    [X,Y,Z] = bresenham_line3d(pinit,pend)
    list_bnd = bND([pinit.x,pinit.y, pinit.z], [pend.x, pend.y, pend.z])
    [X,Y,Z] = map(list,map(None,*list_bnd))
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
## of a different material inside this element.
## Given that all 4 corners are found to be homogeneous, or to belong in to
## the same bin, then is there any inclusion inside this element that is
## larger than a certain area? If so, return true, otherwise return false
def has_inclusions(image,p1,p2,p3,p4,p5):

    pxVal1 = image.GetPixel(p1.x,p1.y,p1.z)
    
    xHigh = p2.x 
    xLow = p1.x    
    yHigh = p4.y
    yLow = p1.y
    zHigh = max(p1.z, p5.z)
    zLow = min(p1.z, p5.z)
    
    areaElem = abs( (p4.y - p1.y) * (p2.x - p1.x) * (p1.z - p5.z))
    
    nr_samples = int( log(PROB) / log (abs(areaElem - VOL_INCLUSION)/areaElem) )
        
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

## Check if the ends of a line are from the same bin
def ends_in_same_bin_2d(image, p1, p2, same_x, same_y, same_z, L):
    if same_x == True:
        pxVal1 = image.GetPixel(L.x, int(p1.x), int(p1.y))
        pxVal2 = image.GetPixel(L.x, int(p2.x), int(p2.y))
        
    if same_y == True:
        pxVal1 = image.GetPixel(int(p1.x), L.y, int(p1.y))
        pxVal2 = image.GetPixel(int(p2.x), L.y, int(p2.y))
        
    if same_z == True:
        pxVal1 = image.GetPixel(int(p1.x), int(p1.y), L.z)
        pxVal2 = image.GetPixel(int(p2.x), int(p2.y), L.z)
        
    val1 = pxVal1 == pxVal2
    val2 = is_in_same_bin(pxVal1, pxVal2)
    return (val1 or val2)

def linear_search(image,bbegin,eend):
        # linear search for interface along a line defined by bbegin and eend coordinates
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
                    
        return list(list_nodes)



def in_child_k(cubes,L):
    # check  if coordinate L is inside element cubes
    p1,p2,p3,p4,p5,p6,p7,p8 = cubes
    if p1.x <= L.x and L.x <= p2.x and p1.y <= L.y and L.y <= p4.y and p1.z <= L.z and L.z<=p5.z:
        return True
    
    return False

def calc_plane_residual(x, y, z):
    #calculate residual and coefficients of planar least square approximation
    # 1 = a*x + b*y + c*z
    a = numpy.column_stack((x, y, z))
    (coeffs, resid,rank, sing_vals) = numpy.linalg.lstsq(a, numpy.ones_like(x))
    return [resid,coeffs]

def draw_plane_connections(image, l1,l2,l3,l4, L1,L2,L3,L4):
    # draw the planar least square approx of interface
    if l1 == 0 and l2 == 0 and l3 == 1 and l4 == 1:
        if coords_not_equal(L1[0],L2[0]):
            draw_line(image,L1[0],L2[0])
        
    if l1 == 0 and l3 == 0:
        if coords_not_equal(L1[0],L3[0]):
            draw_line(image,L1[0],L3[0])
    
    if l1 == 0 and l4 == 0 and l2 == 1 and l3 == 1:
        if coords_not_equal(L1[0],L4[0]):
            draw_line(image,L1[0],L4[0])
    
    if l2 == 0 and l3 == 0 and l1 == 1 and l4 == 1:
        if coords_not_equal(L2[0],L3[0]):
            draw_line(image,L2[0],L3[0])
            
    if l2 == 0 and l4 == 0 and l1 == 1 and l3 == 1:
        if coords_not_equal(L2[0],L4[0]):
            draw_line(image,L2[0],L4[0])
            
    if l3 == 0 and l4 == 0 and l1 == 1 and l2 == 1:
        if coords_not_equal(L3[0],L4[0]):
            draw_line(image,L3[0],L4[0])

    
def set_interval(imSize,level):
    # create uniform grid based on level of recursion
    my_arr = [0,imSize-1]
    
    for i in range(0,level):
        mdpt = numpy.zeros((len(my_arr)-1,1))
        new_arr = numpy.zeros(( len(my_arr) + len(mdpt), 1  ))
        for j in range(0,len(my_arr)-1):
            mdpt[j] = (my_arr[j] + my_arr[j+1] ) // 2

        for k in range(0,len(new_arr)):
            
            if  k % 2 == 0:
                new_arr[k] = my_arr[k/2]
            else:
            
                new_arr[k] = mdpt[k/2]
        my_arr = new_arr
        
    return my_arr 

  
def convert_to_base_8(n):
    # converting number n in base 10 to base 8

    result = []
    (n,remainder) = divmod(n,8)
    result.append(str(remainder))
    while n:
        (n,remainder) = divmod(n,8)
        result.append(str(remainder))
        
    result.reverse()
    
    return ''.join(result)

def tomorton(x,y,z):
# adapted to include the 3D case, from:
#  http://www.thejach.com/view/id/207
  x = bin(x)[2:]
  lx = len(x)
  
  y = bin(y)[2:]
  ly = len(y)
  
  z = bin(z)[2:]
  lz = len(z)
  
  L = max(lx, ly, lz)
  m = 0
  
  for j in xrange(1, L+1):
    # note: ith bit of x requires x[lx - i] since our bin numbers are big endian
    xi = int(x[lx-j]) if j-1 < lx else 0
    yi = int(y[ly-j]) if j-1 < ly else 0
    zi = int(z[lz-j]) if j-1 < lz else 0
  
    m += 2**(3*j)*xi + 2**(3*j+1)*yi + 2**(3*j+2)*zi
  
  return m/8

def find_neighbor(index,direction):
    # finds neighbor given one direction letter
    loc = str(index)
    llist_str = list(loc)
    for i in range(len(loc)-1,-1,-1):
        
        new_quadrant =  D[str(loc[i])][direction]['Quadrant']
        new_direction = D[str(loc[i])][direction]['Direction']
        if new_direction != 'H':
            direction = new_direction
            llist_str[i] = str(new_quadrant)
        else:
            llist_str[i] = str(new_quadrant)
            return str("".join(llist_str))

    return str("".join(llist_str))
  
def find_neighbor_of(index, direction, masterNode):
    # find neighbor of element given index and direction
    root = get_node_by_id(masterNode, index)
    p1,p2,p3,p4,p5,p6,p7,p8 = root.cube
    
    if len(direction) == 1:
        # we only do one lookup pass through the table
        
        # check we are not going outside the image
        if direction == 'L' and p1.x == 0:
            return ""
        if direction == 'R' and p2.x == root.imsize[0]:
            return ""
        if direction == 'F' and p4.y == root.imsize[1]:
            return ""
        if direction == 'B' and p1.y == 0:
            return ""
        if direction == 'D' and p5.z == root.imsize[2]:
            return ""
        if direction == 'U' and p1.z == 0:
            return ""
                
        return find_neighbor(index,direction)
    if len(direction) == 2:
        # we do two passes through the table
        
        # check we are not going outside the image
        if direction == 'LF' and p1.x == 0 and  p4.y == root.imsize[1]:
            return ""
        if direction == 'RF' and p2.x == root.imsize[0] and p4.y == root.imsize[1]:
            return ""
        if direction == 'LB' and p1.x == 0 and p1.y == 0:
            return ""
        if direction == 'RB' and p2.x == root.imsize[0] and p1.y == 0:
            return ""
        if direction == 'UF' and p1.z == 0 and p4.y == root.imsize[1]:
            return ""
        if direction == 'DF' and p5.z == root.imsize[2] and p4.y == root.imsize[1]:
            return ""
        if direction == 'UB' and p1.z == 0 and p1.y == 0:
            return ""
        if direction == 'DB' and p5.z == root.imsize[2] and p1.y == 0:
            return ""
        if direction == 'RU' and p2.x == root.imsize[0] and p1.z == 0:
            return ""
        if direction == 'RD' and p2.x == root.imsize[0] and p5.z == root.imsize[2]:
            return ""
        if direction == 'LU' and p1.x == 0 and p1.z == 0:
            return ""
        if direction == 'LD' and p1.x == 0 and p5.z == root.imsize[2]:
            return ""
        
        pass1 = find_neighbor(index,direction[0])
        pass2 = find_neighbor(pass1,direction[1])
        return pass2
    
    if len(direction) == 3:
        # we do three passes through the table
        
        # check we are not going outside the image
        if direction == 'LUB' and p1.x == 0 and p1.z == 0 and p1.y == 0:
            return ""
        if direction == 'RUB' and p2.x == root.imsize[0] and p1.z == 0 and p1.y == 0:
            return ""
        if direction == 'LUF' and p1.x == 0 and p1.z == 0 and p4.y == root.imsize[1]:
            return ""
        if direction == 'RUF' and p2.x == root.imsize[0] and p1.z == 0 and p4.y == root.imsize[1]:
            return ""
        if direction == 'LDF' and p1.x == 0 and p5.z == root.imsize[2] and p4.y == root.imsize[1]:
            return ""
        if direction == 'RDF' and p2.x == root.imsize[0] and p5.z == root.imsize[2] and p4.y == root.imsize[1]:
            return ""
        if direction == 'LDB' and p1.x == 0 and p5.z == root.imsize[2] and p1.y == 0:
            return ""
        if direction == 'RDB' and p2.x == root.imsize[0] and p5.z == root.imsize[2] and p1.y == 0:
            return ""
        pass1 = find_neighbor(index,direction[0])
        pass2 = find_neighbor(pass1,direction[1])
        pass3 = find_neighbor(pass2,direction[2])
        return pass3
    
          
def get_node_by_id(node,id):
        # returns node, given index 
        # index could be ['00101'], thus it's a list
        
        index = id[0]
        ll = len(index)

        p = node
        for i in range(0,ll):
            p = p.children[int(index[i])]
            
        return p  
    
def get_list_of_nodes(tree, root, masterNode,llist):

        if root.has_children == False:
            llist.append([root.index]) 
            
        if root.children[0] != None:
            get_list_of_nodes(tree,root.children[0],masterNode,llist)
        if root.children[1] != None:
            get_list_of_nodes(tree,root.children[1],masterNode,llist)
        if root.children[2] != None:
            get_list_of_nodes(tree,root.children[2],masterNode,llist)
        if root.children[3] != None:
            get_list_of_nodes(tree,root.children[3],masterNode,llist)
        if root.children[4] != None:
            get_list_of_nodes(tree,root.children[4],masterNode,llist)
        if root.children[5] != None:
            get_list_of_nodes(tree,root.children[5],masterNode,llist)
        if root.children[6] != None:
            get_list_of_nodes(tree,root.children[6],masterNode,llist)
        if root.children[7] != None:
            get_list_of_nodes(tree,root.children[7],masterNode,llist)
            
        return llist
    
    
    
def it_exists(index,masterNode):   
    # check if node of given index exists in the tree masterNode    
    llen = len(index)
    child = masterNode
    for i in range(0,llen):
        if child.children[int(index[i])].has_children == False:
            return False
        child = child.children[int(index[i])]
    return True            

def get_node_of_neighbor(root,my_ind,neigh_ind):
    
    p = root
    for i in range(0,len(my_ind)):
        p = p.parent
        
    r = p
    for j in range(0,len(neigh_ind)):
        r = r.children[int(neigh_ind[j])]
    return r
            
def neigh_has_children(root, masterNode, direction): 
    
    neigh_index = str(find_neighbor_of(root.index,direction,masterNode)) 
    if len(neigh_index)<1:
        return False   
    # checking to see if the neighbor exists or is a ghost
    if it_exists(neigh_index, masterNode):
        node_neighbor = get_node_of_neighbor(root, root.index, neigh_index)
        if node_neighbor.has_children == True:
            return True
    return False
           
def ghost_nodes_enrichment_nodes(tree, root, masterNode):
# this function triggers refinement when both an interface node and a hanging node
# are present on an edge or face of an element
        p1,p2,p3,p4,p5,p6,p7,p8 = root.cube
            
        up_has_children = False
        down_has_children = False
        left_has_children = False
        right_has_children = False
        back_has_children = False
        front_has_children = False
        
        rb_has_children = False
        rf_has_children = False
        lb_has_children = False
        lf_has_children = False
        
        ub_has_children = False
        uf_has_children = False
        db_has_children = False
        df_has_children = False
        
        ld_has_children = False
        lu_has_children = False
        ru_has_children = False
        rd_has_children = False
        
        # if root has no children look at his neighbors
        # if all of them do have children, 
        # root needs to be subdivided
        if root.has_children == False:
            
            up_has_children = neigh_has_children(root,masterNode,'U')
            down_has_children = neigh_has_children(root,masterNode,'D')
            left_has_children = neigh_has_children(root,masterNode,'L')
            right_has_children = neigh_has_children(root,masterNode,'R')
            back_has_children = neigh_has_children(root,masterNode,'B')
            front_has_children = neigh_has_children(root,masterNode,'F')
            
            rb_has_children = neigh_has_children(root,masterNode,'RB')
            rf_has_children = neigh_has_children(root,masterNode,'RF')
            lb_has_children = neigh_has_children(root,masterNode,'LB')
            lf_has_children = neigh_has_children(root,masterNode,'LF')
            
            ub_has_children = neigh_has_children(root,masterNode,'UB')
            uf_has_children = neigh_has_children(root,masterNode,'UF')
            db_has_children = neigh_has_children(root,masterNode,'DB')
            df_has_children = neigh_has_children(root,masterNode,'DF')
            
            ld_has_children = neigh_has_children(root,masterNode,'LD')
            lu_has_children = neigh_has_children(root,masterNode,'LU')
            ru_has_children = neigh_has_children(root,masterNode,'RU')
            rd_has_children = neigh_has_children(root,masterNode,'RD')
                    
            # if the node has interface nodes and its neighbors have children, 
            # we need to subdivide
            if (len(root.enrichNodes) > 0 and 
                (up_has_children == True or down_has_children == True or
                 left_has_children == True or right_has_children == True or
                 back_has_children == True or down_has_children == True or
                 
                 rb_has_children == True or rf_has_children == True or
                 lb_has_children == True or lf_has_children == True or
                 
                 ub_has_children == True or uf_has_children == True or
                 db_has_children == True or df_has_children == True or
                 
                 ld_has_children == True or lu_has_children == True or
                 ru_has_children == True or rd_has_children == True
                 
                 )):

                root.divideOnce()      

        if root.children[0] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[0],masterNode)
        if root.children[1] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[1],masterNode)
        if root.children[2] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[2],masterNode)
        if root.children[3] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[3],masterNode)            
        if root.children[4] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[4],masterNode)
        if root.children[5] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[5],masterNode)
        if root.children[6] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[6],masterNode)
        if root.children[7] != None:
            ghost_nodes_enrichment_nodes(tree,root.children[7],masterNode)
            
def coords_not_equal(p1,p2):
    if p1.x == p2.x and p1.y == p2.y and p1.z == p2.z:
        return 0
    return 1

def neigh_has_grandchildren(root, masterNode, direction, whichChildren): 
    
    neigh_index = str(find_neighbor_of(root.index,direction,masterNode))    
    if len(neigh_index) <1:
        return False
    if it_exists(neigh_index, masterNode):
        node_neighbor = get_node_of_neighbor(root, root.index, neigh_index)
        if node_neighbor.has_children == True:
            has_kids = False
            for i in range(0, len(whichChildren)):
                if node_neighbor.children[whichChildren[i]].has_children:
                    i_has_kids = True
                else:
                    i_has_kids = False
                has_kids = has_kids or i_has_kids
            return has_kids
    return False
            
def tree_balance(tree, root,masterNode):
   
        p1,p2,p3,p4,p5,p6,p7,p8 = root.cube
            
        # if the root node has children, let's look if the neighbors do, 
        # and if they do, do their children have children?
        # 1 - irregularity rule            
        if root.has_children == False:
            if neigh_has_grandchildren(root, masterNode, 'U', [4,5,6,7]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'D', [0,1,2,3]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'F', [0,1,4,5]) == True:
                    root.divideOnce()                                            
            if neigh_has_grandchildren(root, masterNode, 'B', [2,3,6,7]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'L', [1,3,5,7]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'R', [0,2,4,6]) == True:
                    root.divideOnce()
                    
            if neigh_has_grandchildren(root, masterNode, 'UB', [6,7]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'UF', [4,5]) == True:
                    root.divideOnce()                                                            
            if neigh_has_grandchildren(root, masterNode, 'LU', [5,7]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'RU', [4,6]) == True:
                    root.divideOnce()
                    
            if neigh_has_grandchildren(root, masterNode, 'DB', [2,3]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'DF', [0,1]) == True:
                    root.divideOnce()                                                            
            if neigh_has_grandchildren(root, masterNode, 'LD', [1,3]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'RD', [0,2]) == True:
                    root.divideOnce()
                    
                    
            if neigh_has_grandchildren(root, masterNode, 'LF', [1,5]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'RF', [0,4]) == True:
                    root.divideOnce()                                                            
            if neigh_has_grandchildren(root, masterNode, 'LB', [3,7]) == True:
                    root.divideOnce()
            if neigh_has_grandchildren(root, masterNode, 'RB', [2,6]) == True:
                    root.divideOnce()
                                                                                                                        
        if root.children[0] != None:
            tree_balance(tree,root.children[0],masterNode)
        if root.children[1] != None:
            tree_balance(tree,root.children[1],masterNode)
        if root.children[2] != None:
            tree_balance(tree,root.children[2],masterNode)
        if root.children[3] != None:
            tree_balance(tree,root.children[3],masterNode)
        if root.children[4] != None:
            tree_balance(tree,root.children[4],masterNode)
        if root.children[5] != None:
            tree_balance(tree,root.children[5],masterNode)
        if root.children[6] != None:
            tree_balance(tree,root.children[6],masterNode)
        if root.children[7] != None:
            tree_balance(tree,root.children[7],masterNode)
            
def k_neighbor_rule(tree, root,masterNode):
# checks if k neighobors have children, if so, then subdivide
# part 1: no more than k1 faces can have children (not grandchildren)
# part 2: no more than k2: 2 faces and one edge neighbor associated with edge e can have children 
        p1,p2,p3,p4,p5,p6,p7,p8 = root.cube
        
        up_has_children = False
        down_has_children = False
        left_has_children = False
        right_has_children = False
        back_has_children = False
        front_has_children = False
        
        rb_has_children = False
        rf_has_children = False
        lb_has_children = False
        lf_has_children = False
        
        ub_has_children = False
        uf_has_children = False
        db_has_children = False
        df_has_children = False
        
        ld_has_children = False
        lu_has_children = False
        ru_has_children = False
        rd_has_children = False
        

        k1_counter = 0
        k2_counter = 0
        # if root has no children look at his neighbors
        # if all of them do have children, 
        # root needs to be subdivided
        if root.has_children == False:
            
            up_has_children = neigh_has_children(root,masterNode,'U')
            if up_has_children:
                k1_counter += 1
            
            down_has_children = neigh_has_children(root,masterNode,'D')
            if down_has_children:
                k1_counter += 1
                
            left_has_children = neigh_has_children(root,masterNode,'L')
            if left_has_children:
                k1_counter += 1
                
            right_has_children = neigh_has_children(root,masterNode,'R')
            if right_has_children:
                k1_counter += 1
                
            back_has_children = neigh_has_children(root,masterNode,'B')
            if back_has_children:
                k1_counter += 1
                
            front_has_children = neigh_has_children(root,masterNode,'F')
            if front_has_children:
                k1_counter += 1
                
            # PART 1 constraint   
            if k1_counter >= k1_CONST:
                root.divideOnce()
            
            # PART 2 constraint
            rb_has_children = neigh_has_children(root,masterNode,'RB')
            if ( (rb_has_children and right_has_children and back_has_children)
                ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
            
            rf_has_children = neigh_has_children(root,masterNode,'RF')
            if ( (rf_has_children and right_has_children and front_has_children)
                ):
                # edge is nearly regular
                k2_counter += 1
                
            lb_has_children = neigh_has_children(root,masterNode,'LB')
            if ( (lb_has_children and left_has_children and back_has_children)
                 ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
            
            
            lf_has_children = neigh_has_children(root,masterNode,'LF')
            if ( (lf_has_children and left_has_children and front_has_children)
                 ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
            
                
                
            ub_has_children = neigh_has_children(root,masterNode,'UB')
            if ( (ub_has_children and up_has_children and back_has_children)
                  ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
            
                
            uf_has_children = neigh_has_children(root,masterNode,'UF')
            if ( (uf_has_children and up_has_children and  front_has_children)
                  ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
                            
            db_has_children = neigh_has_children(root,masterNode,'DB')
            if ( (db_has_children and down_has_children and back_has_children)
                 ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
                            
            df_has_children = neigh_has_children(root,masterNode,'DF')
            if ( (df_has_children and down_has_children and front_has_children)
                  ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
                        
            
            ld_has_children = neigh_has_children(root,masterNode,'LD')
            if ( (ld_has_children and left_has_children and down_has_children)
                 ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
                            
            lu_has_children = neigh_has_children(root,masterNode,'LU')
            if ( (lu_has_children and left_has_children and up_has_children)
                 ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
                            
            ru_has_children = neigh_has_children(root,masterNode,'RU')
            if ( (ru_has_children and right_has_children and up_has_children)
                  ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
                            
            rd_has_children = neigh_has_children(root,masterNode,'RD')
            if ( (rd_has_children and right_has_children and down_has_children)
                ):
                # edge is nearly regular
                k2_counter += 1
            if k2_counter >= k2_CONST:
                root.divideOnce()                
                        
                                                                                                                        
        if root.children[0] != None:
            k_neighbor_rule(tree,root.children[0],masterNode)
        if root.children[1] != None:
            k_neighbor_rule(tree,root.children[1],masterNode)
        if root.children[2] != None:
            k_neighbor_rule(tree,root.children[2],masterNode)
        if root.children[3] != None:
            k_neighbor_rule(tree,root.children[3],masterNode)
        if root.children[4] != None:
            k_neighbor_rule(tree,root.children[4],masterNode)
        if root.children[5] != None:
            k_neighbor_rule(tree,root.children[5],masterNode)
        if root.children[6] != None:
            k_neighbor_rule(tree,root.children[6],masterNode)
        if root.children[7] != None:
            k_neighbor_rule(tree,root.children[7],masterNode)

def opposite_direction(direc):
# swap direction - look opposite
    N_len = len(direc)    
    direction = list(direc)    
    new_dir = list(' '* N_len)
    
    for i in range(0, N_len):
        if direction[i] == 'U':
            new_dir[i] = 'D'
        if direction[i] == 'D':
            new_dir[i] = 'U'
            
        if direction[i] == 'F':
            new_dir[i] = 'B'
        if direction[i] == 'B':
            new_dir[i] = 'F'
            
        if direction[i] == 'L':
            new_dir[i] = 'R'
        if direction[i] == 'R':
            new_dir[i] = 'L'
            
    return ''.join(new_dir)
    
#=======================================================================
#http://stackoverflow.com/questions/5666222/3d-line-plane-intersection
#=======================================================================
# generic math functions
def add_v3v3(a, b):
    return [a[0] + b[0],
            a[1] + b[1],
            a[2] + b[2]]


def sub_v3v3(a, b):
    return [a[0] - b[0],
            a[1] - b[1],
            a[2] - b[2]]


def dot_v3v3(a, b):
    return (a[0] * b[0] +
            a[1] * b[1] +
            a[2] * b[2])


def mul_v3_fl(a, f):
    a[0] *= f
    a[1] *= f
    a[2] *= f


# intersection function
def isect_line_plane_v3(p0, p1, p_co, p_no, epsilon=0.000001):
    """
    p0, p1: define the line
    p_co, p_no: define the plane:
        p_co is a point on the plane (plane coordinate).
        p_no is a normal vector defining the plane direction; does not need to be normalized.

    return a Vector or None (when the intersection can't be found).
    """

    u = sub_v3v3(p1, p0)
    w = sub_v3v3(p0, p_co)
    dot = dot_v3v3(p_no, u)

    if abs(dot) > epsilon:
        # the factor of the point between p0 -> p1 (0 - 1)
        # if 'fac' is between (0 - 1) the point intersects with the segment.
        # otherwise:
        #  < 0.0: behind p0.
        #  > 1.0: infront of p1.
        fac = -dot_v3v3(p_no, w) / dot
        mul_v3_fl(u, fac)
        return add_v3v3(p0, u)
    else:
        # The segment is parallel to plane
        return None
            

def calc_internal_pts(x,y,z,root,Npts):
    # function to determin internal points to be used with Coons NURBS patch
    # x,y,z: 3 vector lists of coordinates of intersection points
    # len(x) = len(y) = len(z)  is between 3 and 6
    p1,p2,p3,p4,p5,p6,p7,p8 = root.cube
    
    # 1 = a*x + b*y + c*z
    A = numpy.column_stack((x, y, z))
    (coeffs, resid,rank, sing_vals) = numpy.linalg.lstsq(A, numpy.ones_like(x))
    a = coeffs[0]
    b = coeffs[1]
    c = coeffs[2]
    Normal = [a,b,c]
    Centroid = [ (p1.x+p2.x)/2.0, (p1.y + p4.y)/2.0, (p1.z + p5.z)/2.0]
    
    #intersection points of the plane with the box
    int_pts = []
    p12_Int = isect_line_plane_v3([p1.x,p1.y,p1.z], [p2.x, p2.y, p2.z], Centroid, Normal)
    if p12_Int != None:
        if in_child_k(root.cube, Coordinate(p12_Int[0], p12_Int[1], p12_Int[2])):
            int_pts.append(p12_Int)
            
    p15_Int = isect_line_plane_v3([p1.x,p1.y,p1.z], [p5.x, p5.y, p5.z], Centroid, Normal)
    if p15_Int != None:
        if in_child_k(root.cube, Coordinate(p15_Int[0], p15_Int[1], p15_Int[2])):
            int_pts.append(p15_Int)
    
    p14_Int = isect_line_plane_v3([p1.x,p1.y,p1.z], [p4.x, p4.y, p4.z], Centroid, Normal)
    if p14_Int != None:
        if in_child_k(root.cube, Coordinate(p14_Int[0], p14_Int[1], p14_Int[2])):
            int_pts.append(p14_Int)
            
    
    p32_Int = isect_line_plane_v3([p3.x,p3.y,p3.z], [p2.x, p2.y, p2.z], Centroid, Normal)
    if p32_Int != None:
        if in_child_k(root.cube, Coordinate(p32_Int[0], p32_Int[1], p32_Int[2])):
            int_pts.append(p32_Int)
            
    p65_Int = isect_line_plane_v3([p6.x,p6.y,p6.z], [p5.x, p5.y, p5.z], Centroid, Normal)
    if p65_Int != None:
        if in_child_k(root.cube, Coordinate(p65_Int[0], p65_Int[1], p65_Int[2])):
            int_pts.append(p65_Int)
    
    p85_Int = isect_line_plane_v3([p8.x,p8.y,p8.z], [p5.x, p5.y, p5.z], Centroid, Normal)
    if p85_Int != None:
        if in_child_k(root.cube, Coordinate(p85_Int[0], p85_Int[1], p85_Int[2])):
            int_pts.append(p85_Int)
            
    p34_Int = isect_line_plane_v3([p3.x,p3.y,p3.z], [p4.x, p4.y, p4.z], Centroid, Normal)
    if p34_Int != None:
        if in_child_k(root.cube, Coordinate(p34_Int[0], p34_Int[1], p34_Int[2])):
            int_pts.append(p34_Int)
            
    
    p37_Int = isect_line_plane_v3([p3.x,p3.y,p3.z], [p7.x, p7.y, p7.z], Centroid, Normal)
    if p37_Int != None:
        if in_child_k(root.cube, Coordinate(p37_Int[0], p37_Int[1], p37_Int[2])):
            int_pts.append(p37_Int)
    
    p67_Int = isect_line_plane_v3([p6.x,p6.y,p6.z], [p7.x, p7.y, p7.z], Centroid, Normal)
    if p67_Int != None:
        if in_child_k(root.cube, Coordinate(p67_Int[0], p67_Int[1], p67_Int[2])):
            int_pts.append(p67_Int)
            
    p62_Int = isect_line_plane_v3([p6.x,p6.y,p6.z], [p2.x, p2.y, p2.z], Centroid, Normal)
    if p62_Int != None:
        if in_child_k(root.cube, Coordinate(p62_Int[0], p62_Int[1], p62_Int[2])):
            int_pts.append(p62_Int)
            
    p87_Int = isect_line_plane_v3([p8.x,p8.y,p8.z], [p7.x, p7.y, p7.z], Centroid, Normal)
    if p87_Int != None:
        if in_child_k(root.cube, Coordinate(p87_Int[0], p87_Int[1], p87_Int[2])):
            int_pts.append(p87_Int)
            
    p84_Int = isect_line_plane_v3([p8.x,p8.y,p8.z], [p4.x, p4.y, p4.z], Centroid, Normal)
    if p84_Int != None:
        if in_child_k(root.cube, Coordinate(p84_Int[0], p84_Int[1], p84_Int[2])):
            int_pts.append(p84_Int)
                    
    
    if len(int_pts) == 4: # if there are 4 intersection points
        ptA = int_pts[0]
        ptB = int_pts[1]
        ptC = int_pts[2]
        ptD = int_pts[3]
        
        xmin = min( ptA[0], ptB[0], ptC[0], ptD[0])
        xmax = max( ptA[0], ptB[0], ptC[0], ptD[0])
        
        ymin = min( ptA[1], ptB[1], ptC[1], ptD[1])
        ymax = max( ptA[1], ptB[1], ptC[1], ptD[1])

        zmin = min( ptA[2], ptB[2], ptC[2], ptD[2])
        zmax = max( ptA[2], ptB[2], ptC[2], ptD[2])
        
        
        # create grid of points of size Npts x Npts
        XX = numpy.linspace(xmin, xmax, num=Npts)
        XX = XX[1:-1]
        
        YY = numpy.linspace(ymin, ymax, num=Npts)
        YY = YY[1:-1]
        
        ZZ = numpy.linspace(zmin, zmax, num=Npts)
        ZZ = ZZ[1:-1]
        

        # to be implemented: NORMALS through these points
        
def line_face_intersection(centroid, N, pCorner1, pCorner2, pCorner3):
# find the intersection of the ray with a plane (face of the box)
#     begin_line, end_line define the ray through the box
#     pCorner is one of the coordinates of the face plane
    
    # face Back
    begin_line = [centroid.x,centroid.y,centroid.z]
    end_line= [centroid.x+N[0], centroid.y+N[1], centroid.z+N[2]]
    
    # one point on the plane - a corner of the face
    p_co = [pCorner1.x, pCorner1.y, pCorner1.z]
    # compute the vectors defining the plane
    pc1_pc2 = [pCorner2.x - pCorner1.x, pCorner2.y - pCorner1.y, pCorner2.z - pCorner1.z]
    pc1_pc3 = [pCorner3.x - pCorner1.x, pCorner3.y - pCorner1.y, pCorner3.z - pCorner1.z]
    # compute the normal
    p_no = list(numpy.cross(pc1_pc2, pc1_pc3))
    
    return isect_line_plane_v3(begin_line, end_line, p_co, p_no, epsilon=0.000001)

def cmp_floats(a,b):
    abs_diff = abs(a-b)
    if abs_diff < 1e-12:
        return True
    else:
        return False

def intervalcheck(a,b,c):
    """
    Returns whether a <= b <= c is True or False
    """
    if cmp_floats(a,b) == True or cmp_floats(b,c) == True:
        return True
    if a<b and b<c:
        return True
    else:
        return False
    
def which_face(L, root):
            
        p1,p2,p3,p4,p5,p6,p7,p8 = root.cube        
        x1 = p1.x
        x2 = p2.x
        y1 = p1.y
        y2 = p4.y
        z1 = p1.z
        z2 = p5.z
        surface_name = []

        L.x = int(L.x)
        L.y = int(L.y)
        L.z = int(L.z)
        if L.x == x1:
            surface_name.append('L')
        if L.x == x2:
            surface_name.append('R')
        if L.y == y1:
            surface_name.append('B')
        if L.y == y2:
            surface_name.append('F')
        if L.z == z1:
            surface_name.append('U')
        if L.z == z2:
            surface_name.append('D')

        
        return_id = ''
        
        for j in range(len(surface_name)):
            return_id = return_id + surface_name[j] + ''
            
        return return_id
    
            
def line_box_intersection(root, centroid, N):
    inters = []
    p1,p2,p3,p4,p5,p6,p7,p8 = root.cube
    
    #back face
    inters_B = line_face_intersection(centroid, N, p1, p2, p5)
    if inters_B != None:
        # check to see if the intersection point is indeed in the plane
        if in_child_k(root.cube, Coordinate(inters_B[0], inters_B[1], inters_B[2])):
            inters = inters + [inters_B]
        
    #front face
    inters_F = line_face_intersection(centroid, N, p4, p3, p8)
    if inters_F != None:
        # check to see if the intersection point is indeed in the plane
        if in_child_k(root.cube, Coordinate(inters_F[0], inters_F[1], inters_F[2])):
            inters = inters + [inters_F]

    #left face
    inters_L = line_face_intersection(centroid, N, p1, p4, p5)
    if inters_L != None:
        # check to see if the intersection point is indeed in the plane
        if in_child_k(root.cube, Coordinate(inters_L[0], inters_L[1], inters_L[2])):
            inters = inters + [inters_L]
        
    #right face
    inters_R = line_face_intersection(centroid, N, p2, p3, p6)
    if inters_R != None:
        # check to see if the intersection point is indeed in the plane
        if in_child_k(root.cube, Coordinate(inters_R[0], inters_R[1], inters_R[2])):
            inters = inters + [inters_R]
        
    #up face
    inters_U = line_face_intersection(centroid, N, p1, p2, p4)
    if inters_U != None:
        # check to see if the intersection point is indeed in the plane
        if in_child_k(root.cube, Coordinate(inters_U[0], inters_U[1], inters_U[2])):
            inters = inters + [inters_U]

    #down face
    inters_D = line_face_intersection(centroid, N, p5, p6, p8)
    if inters_D != None:
        # check to see if the intersection point is indeed in the plane
        if in_child_k(root.cube, Coordinate(inters_D[0], inters_D[1], inters_D[2])):
            inters = inters + [inters_D]                 
        
    return inters       
                                
def stress_concentration_constraint(tree_list, masterNode, image):

    n = len(tree_list)

    full_list = []
    
    # for each node in the tree:
    for i in range(0,n):
        root_i = get_node_by_id(masterNode,tree_list[i])   
        p1,p2,p3,p4,p5,p6,p7,p8 = root_i.cube
        # for each non-hom node in the tree
        if len(root_i.enrichNodes)>1:

            x_list_c = []
            y_list_c = []
            z_list_c = []
            for j in range(0, len(root_i.enrichNodes)):
                L = root_i.enrichNodes[j] 
                if isinstance(L,list):
                    L = L[0]
                x_list_c.append(L.x)
                y_list_c.append(L.y)
                z_list_c.append(L.z)
            
            if len(x_list_c)>0:
                x_sum = 0.0
                y_sum = 0.0
                z_sum = 0.0
                len_list = len(x_list_c)
                for j in range(0, len_list):
                    x_sum += x_list_c[j]
                    y_sum += y_list_c[j]
                    z_sum += z_list_c[j]
                    
                x_centroid = x_sum / len_list
                y_centroid = y_sum / len_list
                z_centroid = z_sum / len_list
                centroid = Coordinate(int(x_centroid), int(y_centroid), int(z_centroid))
                
                [res,coeffs] = calc_plane_residual(x_list_c, y_list_c, z_list_c)
                N = coeffs
                dx = abs(p1.x - p2.x)
                dy = abs(p1.y - p4.y)
                dz = abs(p1.z - p5.z)
                

                inters = line_box_intersection(root_i, centroid, N)
                if len(inters)<2:
                    break
                currentIndex1 = root_i.index
                currentIndex2 = root_i.index

                one_way = inters[0]
                which_dir1 =  which_face(Coordinate(one_way[0],one_way[1],one_way[2]), root_i)
                counter1 = 0
                list1 = []
                list1.append(root_i.index)
                

                    
                other_way = inters[1]
                which_dir2 = which_face(Coordinate(other_way[0],other_way[1],other_way[2]), root_i)
                counter2 = 0
                list2 = []
                list2.append(root_i.index)
  
  
                neigh_index1 = find_neighbor_of(currentIndex1, which_dir1, masterNode)
                neigh_index2 = find_neighbor_of(currentIndex2, which_dir2, masterNode)
                
                new_inters1 = inters[0] 
                while counter1 <= N_ELEMS:
                    
                    # find the index of the neighbor in the direction given by which_dir1
                    neigh_index1 = find_neighbor_of(currentIndex1, which_dir1, masterNode)
                    if len(neigh_index1) < 1:
                        # we are going outside the image
                        break
                    
                    # if the neighbor exists: either my size or 1/2 smaller
                    if [str(neigh_index1)] in tree_list:
                        # get the neighbor node
                        neigh1 = get_node_by_id(masterNode,[str(neigh_index1)])

                        # if it is 1/2 smaller
                        if neigh1.has_children:
                            
                            # find out in which box is the intersection point:
                            for ii in range(0,8):
                                child_box_index1 = str(str(neigh_index1) + str(ii))
                                child_box1 = get_node_by_id( masterNode, [str(child_box_index1)] )
                                if in_child_k(child_box1.cube, Coordinate(new_inters1[0], new_inters1[1], new_inters1[2]) ):
                                    new_box1 = child_box1
                                    box_index1 = child_box_index1
                                
                        else: # same size with me
                            box_index1 = str(neigh_index1)
                            new_box1 = neigh1
                    else: # the neighbor index doesn't exist, must be 2 times bigger
                        box_index1 = str(neigh_index1[:-1])
                        new_box1 = get_node_by_id(masterNode,[str(box_index1)])
                     
                    if not(new_box1.ishomog):
                        #neighbor is non-homogeneous, add box_index to list and quit            
                        list1.append(box_index1)
                        break  
                   
                    else:
                        # Neighbor is homogeneous 
                                            
                        # get the intersection points of the ray/line/normal with the new box
                        inters_pts = line_box_intersection(new_box1, centroid, N)
                        # if at the end of the image, then no intersection points:
                        if len(inters_pts) < 2:
                            list1.append(box_index1)
                            break 
                        
                        inters1 = inters_pts[0]
                        inters2 = inters_pts[1]
                        # find out on which faces these intersection points lie
                        id1_1 = which_face(Coordinate(inters1[0], inters1[1], inters1[2]), new_box1)
                        id1_2 = which_face(Coordinate(inters2[0], inters2[1], inters2[2]), new_box1)
                        
                        # choose the direction not used for the new direction
                        if opposite_direction(which_dir1) != id1_1:
                            which_dir1 = id1_1
                            new_inters1 = inters1
                        else:
                            which_dir1 = id1_2
                            new_inters1 = inters2
                        
                        # currentIndex becomes that of the new box / neighbor
                        currentIndex1 = str(box_index1)
                        # add the new index to the list, since it is homogeneous
                        list1.append(box_index1)
                        counter1 += 1
                 
                full_list.append(list1)   
                
                new_inters2 = inters[1]
                while counter2 <= N_ELEMS:
                    
                    # find the index of the neighbor in the direction given by which_dir1
                    neigh_index2 = find_neighbor_of(currentIndex2, which_dir2, masterNode)
                    
                    if len(neigh_index2) < 1:
                        # we are going outside the image
                        break
                    
                    # if the neighbor exists: either my size or 1/2 smaller
                    if [str(neigh_index2)] in tree_list:
                        # get the neighbor node
                        neigh2 = get_node_by_id(masterNode,[str(neigh_index2)])

                        # if it is 1/2 smaller
                        if neigh2.has_children:
                            
                            # find out in which box is the intersection point:
                            for ii in range(0,8):
                                child_box_index2 = str(str(neigh_index2) + str(ii))
                                child_box2 = get_node_by_id( masterNode, [str(child_box_index2)] )
                                if in_child_k(child_box2.cube, Coordinate(new_inters2[0], new_inters2[1], new_inters2[2]) ):
                                    new_box2 = get_node_by_id( masterNode, [str(child_box_index2)] )
                                    box_index2 = child_box_index2
                        else: # same size with me
                            box_index2 = str(neigh_index2)
                            new_box2 = get_node_by_id(masterNode,[str(neigh_index2)])
                    else: # the neighbor index doesn't exist, must be 2 times bigger
                        box_index2 = str(neigh_index2[:-1])
                        new_box2 = get_node_by_id(masterNode,[str(box_index2)])
                    
                      
                    if not(new_box2.ishomog):
                        #neighbor is non-homogeneous, add box_index to list and quit            
                        list2.append(box_index2)
                        break  
                   
                    else:
                        # Neighbor is homogeneous 
                        # get the intersection points of the ray/line/normal with the new box
                        inters_pts = line_box_intersection(new_box2, centroid, N)
                        # if at the end of the image, then no intersection points:
                        if len(inters_pts) < 2:
                            list2.append(box_index2)
                            break  
                        inters1 = inters_pts[0]
                        inters2 = inters_pts[1]
                        
                        # find out on which faces these intersection points lie
                        id2_1 = which_face(Coordinate(inters1[0], inters1[1], inters1[2]), new_box2)
                        id2_2 = which_face(Coordinate(inters2[0], inters2[1], inters2[2]), new_box2)
                        
                        # choose the direction not used for the new direction
                        if opposite_direction(which_dir2) != id2_1:
                            which_dir2 = id2_1
                            new_inters2 = inters1
                        else:
                            which_dir2 = id2_2
                            new_inters2 = inters2
                        # currentIndex becomes that of the new box / neighbor
                        currentIndex2 = str(box_index2)
                        # add the new index to the list, since it is homogeneous
                        list2.append(box_index2)
                                    
                    
                    counter2 += 1

                full_list.append(list2)
      
    return full_list
                

def divide_high_stress_elements(full_list, masterNode, tree_list_of_nodes):
    # processing the list with elements in between interfaces
    
    extra_list = []
    for k in range(0, len(full_list)):
        llist = full_list[k]
        if len(llist) < STRESS_MIN + 2:
            last_index = str(llist[-1])

            nodee = get_node_by_id(masterNode, ['0160'])

            last_node = get_node_by_id(masterNode, [str(last_index)])

            if [str(last_index)] in tree_list_of_nodes and len(last_node.enrichNodes)>0: 
                if len(llist) == 2:
                    node1 = get_node_by_id(masterNode, [str(llist[0])])
                    node1.divideOnce()
                    node2 = get_node_by_id(masterNode, [str(llist[1])])
                    node2.divideOnce()
                if len(llist) == 2 + 1:
                    node1 = get_node_by_id(masterNode, [str(llist[1])])
                    node1.divideOnce()
                    extra_list.append(node1)
                    
                        
                if len(llist) == 2 + 2:
                    node1 = get_node_by_id(masterNode, [str(llist[1])])
                    node1.divideOnce()
                    node2 =  get_node_by_id(masterNode, [str(llist[2])])
                    node2.divideOnce()
                    
                if len(llist) == 2 + 3:
                    node1 = get_node_by_id(masterNode, [str(llist[1])])
                    node2 = get_node_by_id(masterNode, [str(llist[2])])
                    node3 = get_node_by_id(masterNode, [str(llist[3])])
                    
                    node1.divideOnce()
                    node2.divideOnce()
                    node3.divideOnce()

    for j in range(0, len(extra_list)):
        node1 = extra_list[j]
        node1.children[0].divideOnce()
        node1.children[1].divideOnce()
        node1.children[2].divideOnce()
        node1.children[3].divideOnce()
        node1.children[4].divideOnce()
        node1.children[5].divideOnce()
        node1.children[6].divideOnce()
        node1.children[7].divideOnce()
                    

def set_homog(masterNode,llist):
    n = len(llist)
    # for each element 
    for i in range(0,n):
        root = get_node_by_id(masterNode,llist[i])


def N_12(t):
    if t<1.0/2.0:
        N_12 = (1 - 2 * t) * (1 - 2 * t)
    else:
        N_12 =  0.0
    return N_12

def N_22(t):        
    if t<1.0/2.0:
        N_22 = 2 * t * (2 - 3 * t)
    else:
        N_22 = 2 * (1 - t) * (1 - t) 
    return N_22

def N_32(t):
    if t<1.0/2.0:
        N_32 =  2 * t * t
    else:
        N_32 = (8*t - 6*t*t - 2)
    return N_32

def N_42(t):
    if t<1.0/2.0:
        N_42 = 0.0
    else:
        N_42 = (2*t - 1) * (2*t -1) 
    return  N_42
                                            
def Nurbs_control_points(Q):
    
    Q = numpy.matrix( [ [Q[0][0], Q[0][1]], [Q[1][0], Q[1][1]], [Q[2][0], Q[2][1]], [Q[3][0], Q[3][1]] ])
    Q = Q[numpy.argsort(Q.A[:,0])]
    R = numpy.matrix([[1.0, 0.0, 0.0, 0.0],[ 1.0/9.0, 2.0/3.0, 2.0/9.0, 0.0 ], [  0.0, 2.0/9.0, 2.0/3.0, 1.0/9.0],[ 0.0, 0.0, 0.0, 1.0] ] )
    P = numpy.linalg.solve(R,Q)
    
    return P

def Nurbs_basis_fcts(t,P):
    
    if t < 1.0/2.0:
        N = lambda t: P[0] * (1 - 2 * t) * (1 - 2 * t) + P[1] * 2 * t * (2 - 3 * t) + P[2] * 2 * t * t
    else:
        N = lambda t: P[1] * 2 * (1 - t) * (1 - t) + P[2] * (8*t - 6*t*t - 2) + P[3] * (2*t - 1) * (2*t -1) 
    return N

def Nurbs_NW_case(image,p1,p2,p3,p4,L1,L4, same_x, same_y, same_z, L):
    
    x_is_F_of_y = False
    
    if (abs(L1.x-p1.x) < abs(L4.y - p1.y)):
        x_is_F_of_y = True
    else:
        x_is_F_of_y = False
    
    if (x_is_F_of_y == False):
        
        # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third12 = Coordinate2D( double ((L1.x - p1.x) / 3.0) + p1.x, p1.y);
        p1third34 = Coordinate2D( p1third12.x, p4.y);
        E = log_search_2d(image,p1third12,p1third34, same_x, same_y, same_z, L);
        
        p2thirds12 = Coordinate2D( double (2.0 * (L1.x - p1.x) / 3.0) + p1.x, p1.y);
        p2thirds34 = Coordinate2D(p2thirds12.x, p4.y);
        F = log_search_2d(image,p2thirds12,p2thirds34, same_x, same_y, same_z, L);
        
        # Step 2. Build the sample points vector   
        Q = [ [L4.x, L4.y], [E.x, E.y], [F.x, F.y], [L1.x, L1.y]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(L1.x-p1.x)
        
        # Step 4. Search for the interface along a vertical line
        p12 = Coordinate2D( (p1.x + L1.x)/2.0, p1.y)
        p34 = Coordinate2D( p12.x, p4.y)
        C = log_search_2d(image,p12,p34, same_x, same_y, same_z, L)
         
        Nx = Nurbs_basis_fcts(0.5,P[:,0])
        Ny = Nurbs_basis_fcts(0.5,P[:,1])
        NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )
               
    else:
        
        # Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
        p1third14 = Coordinate2D( p1.x, float ((L4.y - p1.y) / 3.0) + p1.y);
        p1third23 = Coordinate2D( p2.x, p1third14.y);
        E = log_search_2d(image,p1third14,p1third23, same_x, same_y, same_z, L);
        
        p2thirds14 = Coordinate2D( p1.x, float (2.0 * (L4.y - p1.y) / 3.0) + p1.y);
        p2thirds23 = Coordinate2D( p2.x, p2thirds14.y);
        F = log_search_2d(image,p2thirds14,p2thirds23, same_x, same_y, same_z, L);

        # Step 2. Build the sample points vector   
        Q = [ [L4.y, L4.x], [F.y, F.x], [E.y, E.x], [L1.y, L1.x]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)

        t_step = 1.0 / abs(L4.y - p1.y)
        
        # Step 4. Search for the interface along a vertical line
        p14 = Coordinate2D( p1.x, (p1.y + L4.y) / 2.0 )
        p23 = Coordinate2D( p2.x, p14.y )
        C = log_search_2d(image,p14,p23, same_x, same_y, same_z, L)
        
        Nx = Nurbs_basis_fcts(0.5,P[:,1])
        Ny = Nurbs_basis_fcts(0.5,P[:,0])
        NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )
        
    t_step = T_STEP
    t = arange(0, 1+t_step, t_step)
    
    if find_distance2D(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]


def Nurbs_NE_case(image,p1,p2,p3,p4,L1,L2, same_x, same_y, same_z, L):
    
    x_is_F_of_y = False
        
    if (abs(p2.x-L1.x) < abs(L2.y - p2.y)):
        x_is_F_of_y = True;
    else:
        x_is_F_of_y = False;

    if( x_is_F_of_y == False):
        # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third12 = Coordinate2D( float ((p2.x - L1.x) / 3.0) + L1.x, p1.y);
        p1third34 = Coordinate2D( p1third12.x, p4.y );
        E = log_search_2d(image,p1third12,p1third34, same_x, same_y, same_z, L);

        p2thirds12  = Coordinate2D( float (2.0 * (p2.x - L1.x) / 3.0) + L1.x, p1.y);
        p2thirds34 = Coordinate2D( p2thirds12.x, p4.y);
        F = log_search_2d(image,p2thirds12,p2thirds34, same_x, same_y, same_z, L);
        
        # Step 2. Build the sample points vector   
        Q = [ [L1.x, L1.y], [E.x, E.y], [F.x, F.y], [L2.x, L2.y]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(p2.x - L1.x)
        
        p12 = Coordinate2D( (L1.x + p2.x)/2.0, p2.y )
        p34 = Coordinate2D( p12.x, p3.y )
        C = log_search_2d(image, p12, p34, same_x, same_y, same_z, L)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,0])
        Ny = Nurbs_basis_fcts(0.5,P[:,1])
        NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )
        
        
    else:
        # Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
        p1third14 = Coordinate2D( p1.x, float((L2.y - p2.y) / 3.0) + p2.y );
        p1third23 = Coordinate2D( p2.x, p1third14.y);
        E = log_search_2d(image,p1third14,p1third23, same_x, same_y, same_z, L);

        p2thirds14 = Coordinate2D( p1.x, float (2.0 * (L2.y - p2.y) / 3.0) + p2.y );
        p2thirds23 = Coordinate2D( p2.x, p2thirds14.y );
        F = log_search_2d(image,p2thirds14,p2thirds23, same_x, same_y, same_z, L);
        
        # Step 2. Build the sample points vector   
        Q = [ [L1.y, L1.x], [E.y, E.x], [F.y, F.x], [L2.y, L2.x]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(L2.y - p2.y)
        
        p14 = Coordinate2D( p1.x, (p2.y + L2.y) / 2.0 )
        p23 = Coordinate2D( p2.x, p14.y )
        C = log_search_2d(image, p14, p23, same_x, same_y, same_z, L)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,1])
        Ny = Nurbs_basis_fcts(0.5,P[:,0])
        NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )
        
    t_step = T_STEP
    t = arange(0, 1+t_step, t_step)
    
    if find_distance2D(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]
    

def Nurbs_SE_case(image,p1,p2,p3,p4,L2,L3, same_x, same_y, same_z, L):
    
    x_is_F_of_y = False
    if( abs(L3.x - p3.x) < abs(p3.y - L2.y) ):
        x_is_F_of_y = True;

    else:
        x_is_F_of_y = False;

    if (x_is_F_of_y == False):
        # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third12 = Coordinate2D( float ((p3.x - L3.x) / 3.0) + L3.x, p1.y);
        p1third34 = Coordinate2D( p1third12.x, p4.y);
        E = log_search_2d(image,p1third12,p1third34, same_x, same_y, same_z, L);

        p2thirds12  = Coordinate2D( float (2.0 * (p3.x - L3.x) / 3.0) + L3.x, p1.y );
        p2thirds34 = Coordinate2D( p2thirds12.x, p4.y );
        F = log_search_2d(image,p2thirds12,p2thirds34, same_x, same_y, same_z, L);

        # Step 2. Build the sample points vector   
        Q = [ [L3.x, L3.y], [E.x, E.y], [F.x, F.y], [L2.x, L2.y]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(p3.x - L3.x)
        
        p12 = Coordinate2D( (L3.x + p3.x)/2.0, p2.y)
        p34 = Coordinate2D( p12.x, p3.y)
        C = log_search_2d(image,p12,p34, same_x, same_y, same_z, L)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,0])
        Ny = Nurbs_basis_fcts(0.5,P[:,1])
        NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )
        
    else:
        # Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
        p1third14 = Coordinate2D( p1.x, float ((p3.y - L2.y) / 3.0) + L2.y );
        p1third23 = Coordinate2D( p2.x, p1third14.y);
        E = log_search_2d(image,p1third14,p1third23, same_x, same_y, same_z, L);

        p2thirds14 = Coordinate2D( p1.x, float (2.0 * (p3.y - L2.y) / 3.0) + L2.y );
        p2thirds23 = Coordinate2D( p2.x, p2thirds14.y );
        F = log_search_2d(image,p2thirds14,p2thirds23, same_x, same_y, same_z, L);

        # Step 2. Build the sample points vector   
        Q = [ [L3.y, L3.x], [F.y, F.x], [E.y, E.x], [L2.y, L2.x]]
        
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
    
        t_step = 1.0 / abs(p3.y - L2.y)
        
        p14 = Coordinate2D( p4.x, (L2.y + p3.y) / 2.0 )
        p23 = Coordinate2D( p3.x, p14.y)    
        C = log_search_2d(image,p14,p23, same_x, same_y, same_z, L)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,1])
        Ny = Nurbs_basis_fcts(0.5,P[:,0])
        NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )
    
    t_step = T_STEP
    t = arange(0, 1+t_step, t_step)
    
    if find_distance2D(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]
                    


def Nurbs_SW_case(image,p1,p2,p3,p4,L3,L4, same_x, same_y, same_z, L):
    
    x_is_F_of_y = False
    
    # Step 2. Find intersection of line between L4,L3 and the 45 degree line
    if ( abs(L3.x - p4.x) < abs(p4.y - L4.y) ):
        x_is_F_of_y = True;
    else: 
        x_is_F_of_y = False;


    if (x_is_F_of_y == False):
        # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
        p1third12 = Coordinate2D( float ((p4.x - L3.x) / 3.0) + L3.x, p1.y);
        p1third34 = Coordinate2D( p1third12.x, p4.y);
        E = log_search_2d(image,p1third12,p1third34, same_x, same_y, same_z, L);
        
        p2thirds12 = Coordinate2D(float (2.0 * (p4.x - L3.x) / 3.0) + L3.x, p1.y);
        p2thirds34 = Coordinate2D( p2thirds12.x, p4.y);
        F = log_search_2d(image, p2thirds12, p2thirds34, same_x, same_y, same_z, L);

        # Step 2. Build the sample points vector   
        Q = [ [L4.x, L4.y], [F.x, F.y], [E.x, E.y], [L3.x, L3.y]]

        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
        
        t_step = 1.0 / abs(p4.x - L3.x)
        
        p12 = Coordinate2D( (p4.x + L3.x)/2.0, p1.y)
        p34 = Coordinate2D( p12.x, p4.y )
        C = log_search_2d(image,p12,p34, same_x, same_y, same_z, L)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,0])
        Ny = Nurbs_basis_fcts(0.5,P[:,1])
        NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )
        
    else:
        # Step 1. Search for the interface along horizontal lines at 1/3 and 2/3
        p1third14 = Coordinate2D( p1.x, float ((p4.y - L4.y) / 3.0) + L4.y);
        p1third23 = Coordinate2D( p2.x, p1third14.y);
        E = log_search_2d(image,p1third14,p1third23, same_x, same_y, same_z, L);

        p2thirds14 = Coordinate2D(p1.x, float (2.0 * (p4.y - L4.y) / 3.0) + L4.y );
        p2thirds23 = Coordinate2D(p2.x, p2thirds14.y);
        F = log_search_2d(image, p2thirds14,p2thirds23, same_x, same_y, same_z, L);

        # Step 2. Build the sample points vector   
        Q = [ [L3.y, L3.x], [F.y, F.x], [E.y, E.x], [L4.y, L4.x]]
    
        # Step 3. Determine the control points
        P = Nurbs_control_points(Q)
    
        t_step = 1.0 / abs(p4.y - L4.y)
        
        p14 = Coordinate2D( p1.x, (L4.y + p4.y) / 2.0)
        p23 = Coordinate2D(p3.x, p14.y)
        C = log_search_2d(image,p14,p23, same_x, same_y, same_z, L)
            
        Nx = Nurbs_basis_fcts(0.5,P[:,1])
        Ny = Nurbs_basis_fcts(0.5,P[:,0])
        NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )

    t_step = T_STEP
    t = arange(0, 1+t_step, t_step)
    
    if find_distance2D(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]
    


def Nurbs_vertical_case(image,p1,p2,p3,p4,L1,L3, same_x, same_y, same_z, L):
    
    x_is_F_of_y = True
    
    # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
    p1third14 = Coordinate2D( p1.x, float ((p4.y - p1.y) / 3.0) + p1.y);
    p1third23 = Coordinate2D( p2.x, float ((p3.y - p2.y) / 3.0) + p2.y);
    E = log_search_2d(image,p1third14,p1third23, same_x, same_y, same_z, L);

    p2thirds14 = Coordinate2D(p1.x, float (2.0 * (p4.y - p1.y) / 3.0) + p1.y );
    p2thirds23 = Coordinate2D(p2.x, float (2.0 * (p3.y - p2.y) / 3.0) + p2.y );
    F = log_search_2d(image,p2thirds14,p2thirds23, same_x, same_y, same_z, L);

    # Step 2. Build the sample points vector   
    Q = [ [L3.y, L3.x], [E.y, E.x], [F.y, F.x], [L1.y, L1.x]]
    # Step 3. Determine the control points
    P = Nurbs_control_points(Q)

    t_step = 1.0 / abs(p1.y - p4.y)
    
    p14 = Coordinate2D(p1.x, floor( (p1.y + p4.y)/2.0) )
    p23 = Coordinate2D(p2.x, floor( (p2.y + p3.y)/2.0) )
    C = log_search_2d(image, p14, p23, same_x, same_y, same_z, L)
        
    Nx = Nurbs_basis_fcts(0.5,P[:,1])
    Ny = Nurbs_basis_fcts(0.5,P[:,0])
    NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )
       
    t_step = T_STEP
    t = arange(0, 1+t_step, t_step)
    
    if find_distance2D(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]
    
def Nurbs_horizontal_case(image,p1,p2,p3,p4,L2,L4, same_x, same_y, same_z, L):
    
    x_is_F_of_y = False    

    # Step 1. Search for the interface along vertical lines at 1/3 and 2/3
    p1third12 = Coordinate2D( float ((p2.x - p1.x) / 3.0) + p1.x, p1.y);
    p1third34 = Coordinate2D( float ((p3.x - p4.x) / 3.0) + p4.x, p4.y);
    E = log_search_2d(image,p1third12,p1third34, same_x, same_y, same_z, L);

    p2thirds12 = Coordinate2D( float (2.0 * (p2.x - p1.x) / 3.0) + p1.x, p1.y);
    p2thirds34 = Coordinate2D( float (2.0 * (p3.x - p4.x) / 3.0) + p4.x, p4.y);
    F = log_search_2d(image,p2thirds12,p2thirds34, same_x, same_y, same_z, L);

    # Step 2. Build the cubic polynomial
    xi = [L4.x, E.x, F.x, L2.x];
    yi = [L4.y, E.y, F.y, L2.y];
        
            
    # Step 2. Build the sample points vector   
    Q = [ [L4.x, L4.y], [E.x, E.y], [F.x, F.y], [L2.x, L2.y]]
    
    # Step 3. Determine the control points
    P = Nurbs_control_points(Q)

    t_step = T_STEP
    t_step = 1.0 / abs(p1.x - p2.x)
    
    p12 = Coordinate2D( (p1.x + p2.x)/2.0, p1.y)
    p34 = Coordinate2D( (p3.x + p4.x)/2.0, p4.y)
    C = log_search_2d(image,p12,p34, same_x, same_y, same_z, L)
        
    Nx = Nurbs_basis_fcts(0.5,P[:,0])
    Ny = Nurbs_basis_fcts(0.5,P[:,1])
    NC = Coordinate2D( int(Nx(0.5)) , int(Ny(0.5)) )
            
    t = arange(0, 1+t_step, t_step)
    
    if find_distance2D(C,NC) <= TOL_NURBS:
        return [t,P,x_is_F_of_y, True]
    else:
        return [t,P,x_is_F_of_y, False]

    
def check_nurbs_on_face(inImage, l1,l2,l3,l4, L1,L2,L3,L4, p1,p2,p3,p4):
    # check how well the NURBS is approximated
    same_x = False
    same_y = False
    same_z = False
    
    if p1.x == p2.x and p2.x == p3.x and p3.x == p4.x:
        np1 = Coordinate2D(p1.y, p1.z)
        np2 = Coordinate2D(p2.y, p2.z)
        np3 = Coordinate2D(p3.y, p3.z)
        np4 = Coordinate2D(p4.y, p4.z)
        same_x = True
        
    if p1.y == p2.y and p2.y == p3.y and p3.y == p4.y:
        np1 = Coordinate2D(p1.x, p1.z)
        np2 = Coordinate2D(p2.x, p2.z)
        np3 = Coordinate2D(p3.x, p3.z)
        np4 = Coordinate2D(p4.x, p4.z)
        same_y = True

    if p1.z == p2.z and p2.z == p3.z and p3.z == p4.z:
        np1 = Coordinate2D(p1.x, p1.y)
        np2 = Coordinate2D(p2.x, p2.y)
        np3 = Coordinate2D(p3.x, p3.y)
        np4 = Coordinate2D(p4.x, p4.y)
        same_z = True
        
    
    # NW case
    if (l1==0 and l2==1 and l3==1 and l4==0) and L1.x != -100 and L4.x != -100:
        if coords_not_equal(L1,L4):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L1, L4)  

            [t,P,x_is_F_of_y,test_approx] = Nurbs_NW_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L1)
            if test_approx == False:
                return False
            
    # NE case
    if (l1==0 and l2==0 and l3==1 and l4==1) and L1.x != -100 and L2.x != -100:
        if coords_not_equal(L1,L2):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L1, L2)
            
            [t,P,x_is_F_of_y,test_approx] = Nurbs_NE_case(inImage,np1,np2,np3,np4,nP1,nP2,same_x, same_y, same_z, L1)
            if test_approx == False:
                return False
            
     # SE case                        
    if(l1==1 and l2==0 and l3==0 and l4==1) and L2.x != -100 and L3.x != -100:
        
        if coords_not_equal(L2,L3):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L3, L2)

            [t,P,x_is_F_of_y,test_approx] = Nurbs_SE_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L2)
            if test_approx == False:
                return False
            
     # SW case
    if (l1==1 and l2==1 and l3==0 and l4==0) and L3.x != -100 and L4.x != -100:
        if coords_not_equal(L3,L4):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L4, L3)
            
            [t,P,x_is_F_of_y,test_approx] = Nurbs_SW_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L3)
            if test_approx == False:
                return False
            
     # vertical case
    if (l1==0 and l2==1 and l3==0 and l4==1 ) and L1.x != -100 and L3.x != -100:
        if coords_not_equal(L1,L3):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L1, L3)
            
            [t,P,x_is_F_of_y,test_approx] = Nurbs_vertical_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L1)
            if test_approx == False:
                return False
            
   # horizontal case 
    if (l1==True and l2==False and l3==True and l4==False) and L2.x != -100 and L4.x != -100:  
        if coords_not_equal(L2,L4):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L4, L2)
             
            [t,P,x_is_F_of_y,test_approx] = Nurbs_horizontal_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L2)
            if test_approx == False:
                return False
    
    return True
            
def draw_nurbs_on_face(image, inImage, l1,l2,l3,l4, L1,L2,L3,L4, p1,p2,p3,p4):
    
    same_x = False
    same_y = False
    same_z = False
    
    if p1.x == p2.x and p2.x == p3.x and p3.x == p4.x:
        np1 = Coordinate2D(p1.y, p1.z)
        np2 = Coordinate2D(p2.y, p2.z)
        np3 = Coordinate2D(p3.y, p3.z)
        np4 = Coordinate2D(p4.y, p4.z)
        same_x = True
        
    if p1.y == p2.y and p2.y == p3.y and p3.y == p4.y:
        np1 = Coordinate2D(p1.x, p1.z)
        np2 = Coordinate2D(p2.x, p2.z)
        np3 = Coordinate2D(p3.x, p3.z)
        np4 = Coordinate2D(p4.x, p4.z)
        same_y = True

    if p1.z == p2.z and p2.z == p3.z and p3.z == p4.z:
        np1 = Coordinate2D(p1.x, p1.y)
        np2 = Coordinate2D(p2.x, p2.y)
        np3 = Coordinate2D(p3.x, p3.y)
        np4 = Coordinate2D(p4.x, p4.y)
        same_z = True
        
    
    # NW case
    if (l1==0 and l2==1 and l3==1 and l4==0) and len(L1)>0 and len(L4)>0:
        L1 = L1[0]
        L4 = L4[0]
        if coords_not_equal(L1,L4):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L1, L4)  

            [t,P,x_is_F_of_y,test_approx] = Nurbs_NW_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L1)
            draw_nurbs(image,t,P,x_is_F_of_y,np1,np2,np4, same_x, same_y, same_z, L1)
                             
    # NE case
    if (l1==0 and l2==0 and l3==1 and l4==1) and len(L1) > 0 and len(L2) > 0:
        L1 = L1[0]
        L2 = L2[0]
        if coords_not_equal(L1,L2):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L1, L2)
            
            [t,P,x_is_F_of_y,test_approx] = Nurbs_NE_case(inImage,np1,np2,np3,np4,nP1,nP2,same_x, same_y, same_z, L1)
            draw_nurbs(image,t,P,x_is_F_of_y,np1,np2,np4, same_x, same_y, same_z, L1)
        
     # SE case                        
    if(l1==1 and l2==0 and l3==0 and l4==1) and len(L3) > 0 and len(L2) > 0:
        L2 = L2[0]
        L3 = L3[0]
        if coords_not_equal(L2,L3):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L3, L2)

            [t,P,x_is_F_of_y,test_approx] = Nurbs_SE_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L2)
            draw_nurbs(image,t,P,x_is_F_of_y,np1,np2,np4, same_x, same_y, same_z, L2)
    
     # SW case
    if (l1==1 and l2==1 and l3==0 and l4==0) and len(L3) > 0 and len(L4) > 0:
        L3 = L3[0]
        L4 = L4[0]
        if coords_not_equal(L3,L4):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L4, L3)
            
            [t,P,x_is_F_of_y,test_approx] = Nurbs_SW_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L3)
            draw_nurbs(image,t,P,x_is_F_of_y,np1,np2,np4, same_x, same_y, same_z, L3)
     
     # vertical case
    if (l1==0 and l2==1 and l3==0 and l4==1 ) and len(L1) > 0 and len(L3) > 0:
        L1 = L1[0]
        L3 = L3[0]
        if coords_not_equal(L1,L3):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L1, L3)
            
            [t,P,x_is_F_of_y,test_approx] = Nurbs_vertical_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L1)
            draw_nurbs(image,t,P,x_is_F_of_y,np1,np2,np4, same_x, same_y, same_z, L1)
        
   # horizontal case 
    if (l1==True and l2==False and l3==True and l4==False) and len(L4) > 0 and len(L2) > 0:  
        L2 = L2[0]
        L4 = L4[0]
        if coords_not_equal(L2,L4):
            [nP1, nP2] = sort_nPs(same_x, same_y, same_z, L4, L2)
             
            [t,P,x_is_F_of_y,test_approx] = Nurbs_horizontal_case(inImage,np1,np2,np3,np4,nP2,nP1,same_x, same_y, same_z, L2)
            draw_nurbs(image,t,P,x_is_F_of_y,np1,np2,np4, same_x, same_y, same_z, L2)  
        
def sort_nPs(same_x, same_y, same_z, L1, L4):
    if same_x == True:
        if L1.y <= L4.y:
            nP1 = Coordinate2D(L1.y, L1.z)
            nP2 = Coordinate2D(L4.y, L4.z)
        else:
            nP2 = Coordinate2D(L1.y, L1.z)
            nP1 = Coordinate2D(L4.y, L4.z)
            
    if same_y == True:
        if L1.x <= L4.x:
            nP1 = Coordinate2D(L1.x, L1.z)
            nP2 = Coordinate2D(L4.x, L4.z)
        else:
            nP2 = Coordinate2D(L1.x, L1.z)
            nP1 = Coordinate2D(L4.x, L4.z)
        
    if same_z == True:
        if L1.x <= L4.x:
            nP1 = Coordinate2D(L1.x, L1.y)
            nP2 = Coordinate2D(L4.x, L4.y)
        else:
            nP2 = Coordinate2D(L1.x, L1.y)
            nP1 = Coordinate2D(L4.x, L4.y)
    
    return [nP1, nP2]   
                    
def draw_nurbs(image,t,P,x_is_F_of_y,p1,p2,p4, same_x, same_y, same_z, L):
    px_col = 0
    
    imageSize = image.GetSize()
    
    Px = P[:,0]
    Py = P[:,1]
    
    for i in range(0,len(t)):
        if x_is_F_of_y == False:
            Nx = Nurbs_basis_fcts(t[i],Px)
            Ny = Nurbs_basis_fcts(t[i],Py)
        else:
            Nx = Nurbs_basis_fcts(t[i],Py)
            Ny = Nurbs_basis_fcts(t[i],Px)

        xloc = int(Nx(t[i]))
        yloc = int(Ny(t[i]))
        
        if ( (p1.x <= xloc and xloc <= p2.x and p1.y <= yloc and yloc <= p4.y) or 
         (p1.x <= yloc and yloc <= p2.x and p1.y <= xloc and xloc <= p4.y) ):
                if same_z == True:
                    image.SetPixel(xloc,yloc,L.z,px_col)
                if same_x == True:
                    image.SetPixel(L.x,xloc,yloc,px_col)
                if same_y == True:
                    image.SetPixel(xloc,L.y,yloc,px_col)
        
                  
def log_search_2d(image,bbegin,eend, same_x, same_y, same_z, L):
        
        begin = Coordinate2D(bbegin.x,bbegin.y);
        end = Coordinate2D(eend.x,eend.y);
        mid = Coordinate2D();

        mid.x = (begin.x + end.x)/2.0;
        mid.y = (begin.y + end.y)/2.0;
        dist = find_distance2D(begin,end)
        while dist>2 and not(ends_in_same_bin_2d(image,begin,end, same_x, same_y, same_z, L)):
            if ends_in_same_bin_2d(image,begin, mid, same_x, same_y, same_z, L):
                begin.x = mid.x
                begin.y = mid.y;
                mid.x = (begin.x + end.x)/2.0;
                mid.y = (begin.y + end.y)/2.0;
            else:
                end.x = mid.x;
                end.y = mid.y;
                mid.x = (begin.x + end.x)/2.0;
                mid.y = (begin.y + end.y)/2.0;
        
            dist = find_distance2D(begin,end)

        mid.x = int(mid.x)
        mid.y = int(mid.y)
        return mid
                    

## Find distance between two points
def find_distance2D(p1, p2):
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    return ( sqrt(dx**2 + dy**2) )                                           
