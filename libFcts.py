
import SimpleITK as sitk
from globalVars import *
from math import sqrt, floor, copysign
from numpy import *
import scipy
import numpy
from bresenham import *
from pvtrace import Geometry


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
    def __init__(self,x=-1,y=-1,z=-1):
        self.x = x
        self.y = y
        self.z = z
        self.all = [x,y,z]

            
def search_in(my_list,pi,pj,inImage):
    Lk_list1 = [[x[0],x[1]] in [[pi,pj],] for x in my_list] #[True, False, False, True, etc]
    Lk_list2 = [[x[0],x[1]] in [[pj,pi],] for x in my_list] #[True, False, False, True, etc]
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

def calc_plane_residual(x, y, z):
    # 1 = a*x + b*y + c*z
    a = numpy.column_stack((x, y, z))
    (coeffs, resid,rank, sing_vals) = numpy.linalg.lstsq(a, numpy.ones_like(x))
    return [resid,coeffs]

def draw_plane_connections(image, l1,l2,l3,l4, L1,L2,L3,L4):
#    print l1, l2, l3 ,l4
#    print L1, L2, L3, L4
    
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
  
def find_neighbor_of(index, direction):
    if len(direction) == 1:
        # we only do one lookup pass through the table
        return find_neighbor(index,direction)
    if len(direction) == 2:
        # we do two passes through the table
        pass1 = find_neighbor(index,direction[0])
        pass2 = find_neighbor(pass1,direction[1])
        return pass2
    if len(direction) == 3:
        # we do three passes through the table
        pass1 = find_neighbor(index,direction[0])
        pass2 = find_neighbor(pass1,direction[1])
        pass3 = find_neighbor(pass2,direction[2])
        return pass3
    
          
def get_node_by_id(node,id):
        # returns node, given index 
        # index could be ['00101'], thus it'a list
        
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
    
    neigh_index = str(find_neighbor_of(root.index,direction))    
    # checking to see if the west neighbor exists or is a ghost
    if it_exists(neigh_index, masterNode):
        node_neighbor = get_node_of_neighbor(root, root.index, neigh_index)
        if node_neighbor.has_children == True:
            return True
    return False
           
def ghost_nodes_enrichment_nodes(tree, root, masterNode):
# this function triggers refinement when both an interface node and a hanging node
# are present on an edge or face of an element
        p1,p2,p3,p4,p5,p6,p7,p8 = root.cube
            
#        if p1.x>=192 and p2.x<=254 and p1.y >=384 and p4.y<=446 and p1.z >=0 and p5.z <=62:
#            print root.index
            
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
                    
#            if p1.x>=192 and p2.x<=254 and p1.y >=384 and p4.y<=446 and p1.z >=0 and p5.z <=62:
#                print root.index
#                print up_has_children,down_has_children ,left_has_children ,right_has_children ,back_has_children , front_has_children
#                print rb_has_children, rf_has_children,lb_has_children,lf_has_children 
#                print ub_has_children,uf_has_children ,db_has_children,df_has_children
#                print  ld_has_children, lu_has_children, ru_has_children,rd_has_children 
                
#            print 'ever here?'        
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
#                print 'what about here?'

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
    
    neigh_index = str(find_neighbor_of(root.index,direction))    
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
            if ( (rb_has_children and righ_has_children and back_has_children)
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
            if ( (uf_has_children and up_has_children and  front_has_chidren)
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
# swap direction - look oppposite
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
    
    
def stress_concentration_constraint(tree_list, masterNode, image):

    n = len(tree_list)

    full_list = []
#     LIST = []
    
    # for each node in the tree:
    for i in range(0,n):
        root_i = get_node_by_id(masterNode,tree_list[i])   
        p1,p2,p3,p4,p5,p6,p7,p8 = root_i.cube
        # for each non-hom node in the tree
        if len(root_i.enrichNodes)>1:
#            print len(root_i.enrichNodes)
#            print root_i.enrichNodes[0], root_i.enrichNodes[1]

            x_list_c = []
            y_list_c = []
            z_list_c = []
            for j in range(0, len(root_i.enrichNodes)):
                L = root_i.enrichNodes[j] 
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
                centroid = Coordinate(x_centroid, y_centroid, z_centroid)
                
                [res,coeffs] = calc_plane_residual(x_list_c, y_list_c, z_list_c)
                N = coeffs
            
                dx = abs(p1.x - p2.x)
                dy = abs(p1.y - p5.y)
                dz = abs(p1.z - p4.z)
                box = Geometry.Box([p1.x,p1.y,p1.z], [dx,dy,dz])
                ray = Geometry.Ray(position=(centroid.x, centroid.y, centroid.z), direction=(N[0],N[1],N[2]))
                inters = box.intersection(ray)
                print inters
                one_way = inters[0]
                id1 = box.surface_identifier(one_way)
                counter1 = 0
                list1 = []
                
                other_way = inters[1]
                id2 = box.surface_identifier(other_way)
                counter2 = 0
                list2 = []
                                  
def stress_concentration_constraint2(tree_list, masterNode, image):

    n = len(tree_list)

    full_list = []
    
    # for each node in the tree:
    for i in range(0,n):
        root_i = get_node_by_id(masterNode,tree_list[i])    
         
        if len(root_i.enrichNodes) > 1: # for each non-hom node
        # root_i.index=='312':
            [side1,side2,side3,side4,neigh_list,dir_list,whichSide] = element_normal_intersection(root_i.enrichNodes[0], root_i.enrichNodes[1], root_i, image)
#             print '-----------------'
#             print root_i.index
#             root_i.printRect()
#             print root_i.enrichNodes[0].x, root_i.enrichNodes[0].y
#             print root_i.enrichNodes[1].x, root_i.enrichNodes[1].y
#             print find_neighbor_index_of(root_i.index,dir_list[0], masterNode, tree_list)
#             print find_neighbor_index_of(root_i.index,dir_list[1], masterNode, tree_list)
#             print dir_list
#             print '++++++++++++++++'

#             print neigh_list
            counter1 = 0
            counter2 = 0
            
            list1 = []
            list2 = []
            
            whichEdge1 = dir_list[0]
            whichEdge2 = dir_list[1]
            
            currentIndex1 = root_i.index
            currentIndex2 = root_i.index
            list1.append(currentIndex1)
            list2.append(currentIndex2)

#             print 'original ',whichSide[0].x, whichSide[0].y, whichSide[1].x, whichSide[1].y
            whichSidePrev1 = copy_list_of_sides(whichSide)

            while counter1 <= 4:
#                 print 'counter 1'
                neighs = find_neighbor_index_of(currentIndex1,whichEdge1, masterNode, tree_list)
                
                if len(neighs) == 2: # there are 2 neighbors sharing an edge with me
                    neigh1 = get_node_by_id(masterNode,[str(neighs[0])]) 
                    neigh2 = get_node_by_id(masterNode,[str(neighs[1])])

                    p1n1,p2n1,p3n1,p4n1 = neigh1.rect
                    
                    if whichEdge1 == 'U':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side1.x >= p1n1.x and side1.x <= p2n1.x:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])

                    if whichEdge1 == 'D':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side3.x >= p1n1.x and side3.x <= p2n1.x:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])

                    if whichEdge1 == 'L':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side4.y >= p1n1.y and side4.y <= p4n1.y:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])
                            
                    if whichEdge1 == 'R':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side2.y >= p1n1.y and side2.y <= p4n1.y:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])
                    
                else:
                    if len(neighs) == 0: # there are no neighbors of mine, perhaps we are at the margins of the image
                        break
                    else:
                        neighIndex = str(neighs[0])

                neigh_node = get_node_by_id(masterNode,[str(neighIndex)])
                if not(neigh_node.ishomog):
                    #neighbor is non-homogeneous:                    
                    list1.append(neighIndex)
                    break
                else:
                    counter1 += 1
                    [side1N,side2N,side3N,side4N,neigh_listN,dir_listN,whichSideN] = element_normal_intersection(root_i.enrichNodes[0], root_i.enrichNodes[1], neigh_node, image)
                    currentIndex1 = neighIndex
#                     whichEdge1 = swap_edges(whichEdge1, dir_listN)
#                     print 'counter 1, whichEdge=', whichEdge1, 'neigh index = ', neighIndex
#                     print 'dirlistN=', dir_listN, 'whichSideN=', whichSideN
#                     print side1N.x, side1N.y, side2N.x,side2N.y, side3N.x, side3N.y, side4N.x,side4N.y

#                     root_i.printRect()
#                     neigh_node.printRect()
                    
                    whichEdge1 = swap_directions(dir_listN,whichSidePrev1,whichSideN)
                    whichSidePrev1 = copy_list_of_sides(whichSideN)
                    list1.append(neighIndex)
                    side1.x = side1N.x
                    side1.y = side1N.y
                    side2.x = side2N.x
                    side2.y = side2N.y
                    side3.x = side3N.x
                    side3.y = side3N.y
                    side4.x = side4N.x
                    side4.y = side4N.y
                    
                 
#             print root_i.index
#             print list1   
            #===================================================================
#             divide_high_stress_elements(list1,masterNode)
            #===================================================================
            full_list.append(list1)
            
            [side1,side2,side3,side4,neigh_list,dir_list,whichSide] = element_normal_intersection(root_i.enrichNodes[0], root_i.enrichNodes[1], root_i,image)

            whichSidePrev2 = copy_list_of_sides(whichSide)
            dir_listN = dir_list
            while counter2 <= 4:
                
#                 newpts = []
#                 if side1.x != NoneINT and side1.y != NoneINT:
#                     newpts.append(side1)
#                 if side2.x != NoneINT and side2.y != NoneINT:
#                     newpts.append(side2)
#                 if side3.x != NoneINT and side3.y != NoneINT:
#                     newpts.append(side3)
#                 if side4.x != NoneINT and side4.y != NoneINT:
#                     newpts.append(side4)
#                     
#                 draw_line(newpts[0], newpts[1])

#                 
#                 print 'counter 2'
#                 print side1.x, side1.y, side2.x,side2.y, side3.x, side3.y, side4.x,side4.y
                
                neighs = find_neighbor_index_of(currentIndex2,whichEdge2, masterNode, tree_list)
                
#                 print neighs, currentIndex2, whichEdge2
                
                if len(neighs) == 2: # there are 2 neighbors sharing an edge with me
                    neigh1 = get_node_by_id(masterNode,[str(neighs[0])]) 
                    neigh2 = get_node_by_id(masterNode,[str(neighs[1])])

#                     print 'blah, blah',str(neighs[0]), str(neighs[1])
#                     print neigh2.index
                    
#                     print neigh1.index, neigh1.printRect()
#                     print whichEdge2
                    
                    p1n1,p2n1,p3n1,p4n1 = neigh1.rect
                    
                    if whichEdge2 == 'U':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side1.x >= p1n1.x and side1.x <= p2n1.x:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])

                    if whichEdge2 == 'D':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side3.x >= p1n1.x and side3.x <= p2n1.x:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])

                    if whichEdge2 == 'L':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side4.y >= p1n1.y and side4.y <= p4n1.y:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])
                            
                    if whichEdge2 == 'R':# or whichEdge1 == 'LU' or whichEdge1 ==  'RU':
                        if side2.y >= p1n1.y and side2.y <= p4n1.y:
                            neighIndex = str(neighs[0])
                        else:
                            neighIndex = str(neighs[1])
                    
                else:
                    if len(neighs) == 0: # there are no neighbors of mine, perhaps we are at the margins of the image
                        break
                    else:
                        neighIndex = str(neighs[0])

                neigh_node = get_node_by_id(masterNode,[str(neighIndex)])
                if not(neigh_node.ishomog):
                    #neighbor is non-homogeneous:                    
                    list2.append(neighIndex)
                    break
                else:
#                     print neighIndex, currentIndex2, neigh_node.ishomog,
                    counter2 += 1
                    [side1N,side2N,side3N,side4N,neigh_listN,dir_listN,whichSideN] = element_normal_intersection(root_i.enrichNodes[0], root_i.enrichNodes[1], neigh_node, image)
                    
#                     if len(dir_listN) < 2:
#                         list2.append(neighIndex)
#                         break
#                     else:
                    currentIndex2 = neighIndex
                    
                    neigh_node.printRect()
#                     print 'dirlist', dir_listN, neighIndex
#                     print side1N.x, side1N.y, side2N.x,side2N.y, side3N.x, side3N.y, side4N.x,side4N.y
#                     print 'old edge:', whichEdge2
#                     whichEdge2 = swap_edges(whichEdge2, dir_listN)
                    whichEdge2 = swap_directions(dir_listN,whichSidePrev2,whichSideN)
                    whichSidePrev2 = copy_list_of_sides(whichSideN)
#                     print 'new edge: ', whichEdge2
                    list2.append(neighIndex)
                    side1.x = side1N.x
                    side1.y = side1N.y
                    side2.x = side2N.x
                    side2.y = side2N.y
                    side3.x = side3N.x
                    side3.y = side3N.y
                    side4.x = side4N.x
                    side4.y = side4N.y
                    
#             print list2
            #===================================================================
#             masterNode = divide_high_stress_elements(list2,masterNode)
            #===================================================================
            full_list.append(list2)

    return full_list

def divide_high_stress_elements(full_list, masterNode,image):
    # processing the list with elements in between interfaces
    
    extra_list = []
    for k in range(0, len(full_list)):
        llist = full_list[k]
        
        if len(llist) < STRESS_MIN + 2:
            last_index = str(llist[-1])
            last_node = get_node_by_id(masterNode, [str(last_index)])
            if not(last_node.ishomog): 
            
                if len(llist) == 2:
                    node1 = get_node_by_id(masterNode, [str(llist[0])])
                    node1.divideOnce()
                    node2 = get_node_by_id(masterNode, [str(llist[1])])
                    node2.divideOnce()
                if len(llist) == 2 + 1:
    #                 print 'only one homog elem in between'
                    node1 = get_node_by_id(masterNode, [str(llist[1])])
                    node1.divideOnce()
                    extra_list.append(node1)
                    
                        
                if len(llist) == 2 + 2:
    #                 print 'only 2 homog elems in between'
                    node1 = get_node_by_id(masterNode, [str(llist[1])])
                    node1.divideOnce()
                    node2 =  get_node_by_id(masterNode, [str(llist[2])])
                    node2.divideOnce()
                    
                if len(llist) == 2 + 3:
    #                 print 'only 3 homog elems in between'
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
                    
def element_normal_intersection(pt1,pt2,node,image, centroid, coeffs):
    
    NoneINT = -9999
    
    ptN1 = Coordinate( centroid.x + coeffs[0], centroid.y + coeffs[1], centroid.z + coeffs[2])
    ptN2 = Coordinate( centroid.x - coeffs[0], centroid.y - coeffs[1], centroid.z - coeffs[2])

    
    if 0 <= ptN1.x and ptN1.x <= node.imsize[0] and 0 <= ptN1.y and ptN1.y <= node.imsize[1] and 0 <= ptN1.z and ptN1.z <= node.imsize[2]:
        ptN = Coordinate(ptN1.x, ptN1.y, ptN1.z)
    elif 0 <= ptN2.x and ptN2.x <= node.imsize[0] and 0 <= ptN2.y and ptN2.y <= node.imsize[1] and 0 <= ptN2.z and ptN2.z <= node.imsize[2]:
        ptN = Coordinate(ptN2.x, ptN2.y, ptN2.z)
        
    dx_m = ptN.x - centroid.x
    dy_m = ptN.y - centroid.y
    dz_m = ptN.z - centroid.z

    
    p1,p2,p3,p4,p5,p6,p7,p8 = node.cube

    
    #Compute the intersection of the normal with the 4 sides of an element
    if dx_m == 0: #vertical
        side1 = Coordinate(midpt.x, p1.y)
        side2 = Coordinate(NoneINT,NoneINT) # normal runs parallel with the edge
        side3 = Coordinate(midpt.x, p4.y)
        side4 = Coordinate(NoneINT,NoneINT) # normal runs parallel with the edge
#         draw_line(image, side1,side3)
    elif dy_m == 0: #horizontal
        side1 = Coordinate(NoneINT, NoneINT) # normal runs parallel with the edge
        side2 = Coordinate(p2.x, midpt.y)
        side3 = Coordinate(NoneINT,NoneINT) # normal runs parallel with the edge
        side4 = Coordinate(p1.x, midpt.y)
#         draw_line(image, side2,side4)
    else:
        m_slope = float(dy_m) / dx_m
        b = midpt.y - m_slope * midpt.x
#         print 'suntem aici!!!!', m_slope, b
#         node.printRect()
#         print 'side1', float(p1.y - b) / m_slope
#         print 'side4', float(p4.y - b) / m_slope 
        # SIDE 1
        if p1.x <= float(p1.y - b) / m_slope and float(p1.y - b) / m_slope <= p2.x:
            side1 = Coordinate( float(p1.y - b) / m_slope, p1.y)
            side1.x = int(side1.x)
            side1.y = int(side1.y)
#             draw_line(image, midpt, side1)
        else: # intersection happens outisde the element
            side1 = Coordinate(NoneINT, NoneINT)
        
        # SIDE 2    
        if p1.y <= m_slope * p2.x + b and  m_slope * p2.x + b <= p4.y:
            side2 = Coordinate( p2.x, m_slope * p2.x + b)
            side2.x = int(side2.x)
            side2.y = int(side2.y)
#             draw_line(image, midpt, (side2))
        else: # intersection happens outisde the element
            side2 = Coordinate(NoneINT, NoneINT)
        
        # SIDE 3
        if p1.x <= float(p4.y - b) / m_slope and float(p4.y - b) / m_slope <= p2.x:
            side3 = Coordinate( float(p4.y - b) / m_slope, p4.y)
            side3.x = int(side3.x)
            side3.y = int(side3.y)
#             draw_line(image, midpt, side3)
        else: # intersection happens outisde the element
            side3 = Coordinate(NoneINT, NoneINT)
        # SIDE 4    
        if p1.y <= m_slope * p1.x + b and m_slope * p1.x + b <= p4.y:
            side4 = Coordinate( p1.x, m_slope * p1.x + b)
            side4.x = int(side4.x)
            side4.y = int(side4.y)
#             draw_line(image, midpt, side4)
        else: # intersection happens outisde the element
            side4 = Coordinate(NoneINT, NoneINT)

    side1.x = int(side1.x)
    side1.y = int(side1.y)
    side2.x = int(side2.x)
    side2.y = int(side2.y)
    side3.x = int(side3.x)
    side3.y = int(side3.y)
    side4.x = int(side4.x)
    side4.y = int(side4.y)
    
    N_edge = 0
    S_edge = 0
    W_edge = 0
    E_edge = 0
    NW_edge = 0
    NE_edge = 0
    SE_edge = 0
    SW_edge = 0
    direction_list = []
    whichSide = []
    
    if (side1.x != NoneINT and side1.x != p1.x and side1.x != p2.x and
        side1.y != NoneINT and (side1.y == p1.y or side1.y == p2.y) ):
#         print 'normal going through North edge' 
        N_edge = 1
        direction_list.append('U')
        whichSide.append(side1)
        
    if (side2.x != NoneINT and (side2.x == p3.x or side2.x == p2.x) and
        side2.y != NoneINT and side2.y != p3.y and side2.y != p2.y):
#         print 'normal going through East edge' 
#         print p2.x, p3.x, p3.y, p2.y, side2.x, side2.y
        E_edge = 1
        direction_list.append('R')
        whichSide.append(side2)
        
    if (side3.x != NoneINT and side3.x != p3.x and side3.x != p4.x and
        side3.y != NoneINT and (side3.y == p3.y or side3.y == p4.y) ):
#         print 'normal going through South edge' 
        S_edge = 1
        direction_list.append('D')
        whichSide.append(side3)
        
    if (side4.x != NoneINT and (side4.x == p1.x or side4.x == p4.x) and
        side4.y != NoneINT and side4.y != p1.y and side4.y != p4.y):
#         print 'normal going through West edge' 
        W_edge = 1
        direction_list.append('L')
        whichSide.append(side4)
        
    if (side1.x == p2.x and side1.y == p2.y) or (side2.x == p2.x and side2.y == p2.y):
#         print 'normal going through North East corner' 
        NE_edge = 1
        direction_list.append('RU')
        if (side1.x == p2.x and side1.y == p2.y) :
            whichSide.append(side1)
        else:
            whichSide.append(side2)
        
    if (side1.x == p1.x and side1.y == p1.y) or (side4.x == p1.x and side4.y == p1.y):
#         print 'normal going through North West corner'
        NW_edge = 1
        direction_list.append('LU')
        if (side1.x == p1.x and side1.y == p1.y):
            whichSide.append(side1)
        else:
            whichSide.append(side4)
    
    if (side2.x == p3.x and side2.y == p3.y) or (side3.x == p3.x and side3.y == p3.y):
#         print 'normal going through South East corner' 
        SE_edge = 1
        direction_list.append('RD')
        if (side2.x == p3.x and side2.y == p3.y):
            whichSide.append(side2)
        else:
            whichSide.append(side3)
        
    if (side3.x == p4.x and side3.y == p4.y) or (side4.x == p4.x and side4.y == p4.y):
#         print 'normal going through South West edge' 
        SW_edge = 1
        direction_list.append('LD')
        if (side3.x == p4.x and side3.y == p4.y) :
            whichSide.append(side3)
        else:
            whichSide.append(side4)
        
    neigh_list = [NW_edge,N_edge,NE_edge,E_edge,SE_edge,S_edge,SW_edge,W_edge]
    return [side1,side2,side3,side4,neigh_list, direction_list,whichSide]

        
def set_homog(masterNode,llist):
    n = len(llist)
    # for each element 
    for i in range(0,n):
        root = get_node_by_id(masterNode,llist[i])

        if len(root.enrichNodes)>1:
            print root.index
                                            