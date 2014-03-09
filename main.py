import numpy
import SimpleITK as sitk
from libFcts import *
from globalVars import *
import dicom
import scipy
        
class Node():
    ROOT = 0
    BRANCH = 1
    LEAF = 2
    MAX_DEPTH = 0
    
    def __init__(self, parent, cube, inImage,outImage,imageSize):

        self.imsize = imageSize
        self.parent = parent
        self.children = [None, None, None, None,None,None,None,None]
        self.has_children = False

#        self.maxdepth = 0
        if parent == None:
            self.depth = 0
            
        else:
            self.depth = parent.depth + 1
            if self.depth > Node.MAX_DEPTH:
                Node.MAX_DEPTH = self.depth
                
        
        self.ishomog = 1
#        self.tlist = [] # contains the numbering of the nodes in the element
#        self.tpix = []       
        self.nsew = [0,0,0,0]
        self.hn = [Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1),Coordinate(-1,-1)]
        
        self.cube = cube 
        self.index = '-1'
        [p1,p2,p3,p4,p5,p6,p7,p8] = cube
#        
####        dx = abs(p1.x-p2.x)+1
####        dy = abs(p1.y-p4.y)+1
#####        ind_x = round(imageSize[0]/dx)
#####        ind_y = round(imageSize[1]/dy)
####        self.i = p2.x / dx
####        self.j = p4.y / dy
#
#
##        self.mat = 'Epoxy'
        self.enrichNodes = []
#        
        if self.parent == None:
            self.type = Node.ROOT
#        #elif abs(p1.x - p2.x) <= MIN_SIZE:
        elif ( abs(p1.x - p2.x) <= MIN_SIZE or 
            (self.children[0]==None and self.children[1]== None and self.children[2] == None and self.children[3] == None) ):
#            print self.cube
            self.type = Node.LEAF
        else:
            self.type = Node.BRANCH
        
        self.outImage = outImage
        self.inImage = inImage
        
        Lx = set_interval(imageSize[0],self.depth)
        Ly = set_interval(imageSize[1],self.depth)
        Lz = set_interval(imageSize[2],self.depth)

        self.i = list(Lx).index(p1.x)
        self.j = list(Ly).index(p1.y)
        self.k = list(Lz).index(p1.z)
        
    def subdivide(self): 
    # this method subdivides a node recursively if some
    # division criterion is satisfied
    
#        if self.type == Node.LEAF:
#            return
        
        p1,p2,p3,p4,p5,p6,p7,p8 = self.cube
        cMid12 = find_mid_point(p1,p2)
        cMid14 = find_mid_point(p1,p4)
        cMid23 = find_mid_point(p2,p3)
        cMid34 = find_mid_point(p3,p4)
        
        
        cMid15 = find_mid_point(p1,p5)
        cMid26 = find_mid_point(p2,p6)
        cMid37 = find_mid_point(p3,p7)
        cMid48 = find_mid_point(p4,p8)
        
        cMid56 = find_mid_point(p5,p6)
        cMid67 = find_mid_point(p6,p7)
        cMid87 = find_mid_point(p8,p7)
        cMid58 = find_mid_point(p5,p8)
        
        cMid1234 = find_mid_point(p2,p4)
        cMid5678 = find_mid_point(p6,p8)
        cMid3267 = find_mid_point(p2,p7)
        cMid4158 = find_mid_point(p1,p8)
        
        cMid1265 = find_mid_point(p2,p5)
        cMid4378 = find_mid_point(p3,p8)
        cMid = find_mid_point(p2,p8) # center of the cube
        
        cubes = []
        cubes.append((p1,cMid12,cMid1234,cMid14, cMid15, cMid1265, cMid, cMid4158)) #NW top child
        cubes.append((cMid12,p2,cMid23,cMid1234, cMid1265, cMid26, cMid3267, cMid)) #NE top child
        cubes.append((cMid14,cMid1234,cMid34,p4, cMid4158, cMid, cMid4378, cMid48)) #SW top child
        cubes.append((cMid1234,cMid23,p3,cMid34, cMid, cMid3267, cMid37, cMid4378)) #SE top child
        cubes.append((cMid15, cMid1265, cMid, cMid4158, p5, cMid56, cMid5678, cMid58)) # NW bottom child
        cubes.append(( cMid1265, cMid26, cMid3267, cMid, cMid56, p6, cMid67, cMid5678)) # NE bottom child
        cubes.append((cMid4158, cMid, cMid4378, cMid48, cMid58, cMid5678, cMid87, p8)) # SW bottom child
        cubes.append((cMid, cMid3267, cMid37, cMid4378, cMid5678, cMid67, p7, cMid87)) # SE bottom child
        
        for n in range(len(cubes)):
            span = self.division_criterion(cubes[n], self.inImage, self.outImage)

            if span == True:
#                print 'criterion is TRUE'
                self.children[n] = self.getinstance(cubes[n], self.inImage, self.outImage,imageSize)
                self.children[n].index = str(convert_to_base_8(tomorton(self.children[n].i, self.children[n].j)))
                diff_level = abs(len(self.children[n].index) - self.children[n].depth)
                if diff_level != 0:
                    self.children[n].index = '0'*diff_level + self.children[n].index
                
                p1r,p2r,p3r,p4r,p5r,p6r,p7r,p8r = cubes[n]
                L1 = linear_search(self.inImage,p1r,p2r)
                L2 = linear_search(self.inImage,p2r,p3r)
                L3 = linear_search(self.inImage,p4r,p3r)
                L4 = linear_search(self.inImage,p1r,p4r)
                
                L5 = linear_search(self.inImage,p5r,p6r)
                L6 = linear_search(self.inImage,p6r,p7r)
                L7 = linear_search(self.inImage,p8r,p7r)
                L8 = linear_search(self.inImage,p5r,p8r)
                
                L9 = linear_search(self.inImage,p1r,p5r);
                L10 = linear_search(self.inImage,p2r,p6r);
                L11 = linear_search(self.inImage,p3r,p7r);
                L12 = linear_search(self.inImage,p4r,p8r);
# 
#                if len(L1) == 1:
#                    L1 = L1[0]
#                    if in_child_k(cubes[n],L1) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L1]
#                if len(L2) == 1:
#                    L2 = L2[0]
#                    if in_child_k(cubes[n],L2) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L2]
#                if len(L3) == 1:
#                    L3 = L3[0]
#                    if in_child_k(cubes[n],L3) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L3]
#                if len(L4) == 1:
#                    L4 = L4[0]
#                    if in_child_k(cubes[n],L4) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L4]
#                        
#                if len(L5) == 1:
#                    L5 = L5[0]
#                    if in_child_k(cubes[n],L5) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L5]
#                if len(L6) == 1:
#                    L6 = L6[0]
#                    if in_child_k(cubes[n],L6) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L6]
#                if len(L7) == 1:
#                    L7 = L7[0]
#                    if in_child_k(cubes[n],L7) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L7]
#                if len(L8) == 1:
#                    L8 = L8[0]
#                    if in_child_k(cubes[n],L8) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L8]
#                        
#                if len(L9) == 1:
#                    L9 = L9[0]
#                    if in_child_k(cubes[n],L9) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L9]
#                if len(L10) == 1:
#                    L10 = L10[0]
#                    if in_child_k(cubes[n],L10) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L10]
#                if len(L11) == 1:
#                    L11 = L11[0]
#                    if in_child_k(cubes[n],L11) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L11]
#                if len(L12) == 1:
#                    L12 = L12[0]
#                    if in_child_k(cubes[n],L12) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L12]                                
            
                self.children[n].subdivide()
                
        if ( self.children[0] != None and
             self.children[1] != None and
             self.children[2] != None and
             self.children[3] != None and
             self.children[4] != None and
             self.children[5] != None and
             self.children[6] != None and
             self.children[7] != None             
             ):
            self.has_children = True
        else:
            self.has_children = False    
            
                        
    def divideOnce(self): 
    # this method divides once a given node
        
        p1,p2,p3,p4,p5,p6,p7,p8 = self.cube
        cMid12 = find_mid_point(p1,p2)
        cMid14 = find_mid_point(p1,p4)
        cMid23 = find_mid_point(p2,p3)
        cMid34 = find_mid_point(p3,p4)
        
        
        cMid15 = find_mid_point(p1,p5)
        cMid26 = find_mid_point(p2,p6)
        cMid37 = find_mid_point(p3,p7)
        cMid48 = find_mid_point(p4,p8)
        
        cMid56 = find_mid_point(p5,p6)
        cMid67 = find_mid_point(p6,p7)
        cMid87 = find_mid_point(p8,p7)
        cMid58 = find_mid_point(p5,p8)
        
        cMid1234 = find_mid_point(p2,p4)
        cMid5678 = find_mid_point(p6,p8)
        cMid3267 = find_mid_point(p2,p7)
        cMid4158 = find_mid_point(p1,p8)
        
        cMid1265 = find_mid_point(p2,p5)
        cMid4378 = find_mid_point(p3,p8)
        cMid = find_mid_point(p2,p8) # center of the cube
        
        cubes = []
        cubes.append((p1,cMid12,cMid1234,cMid14, cMid15, cMid1265, cMid, cMid4158)) #NW top child
        cubes.append((cMid12,p2,cMid23,cMid1234, cMid1265, cMid26, cMid3267, cMid)) #NE top child
        cubes.append((cMid14,cMid1234,cMid34,p4, cMid4158, cMid, cMid4378, cMid48)) #SW top child
        cubes.append((cMid1234,cMid23,p3,cMid34, cMid, cMid3267, cMid37, cMid4378)) #SE top child
        cubes.append((cMid15, cMid1265, cMid, cMid4158, p5, cMid56, cMid5678, cMid58)) # NW bottom child
        cubes.append(( cMid1265, cMid26, cMid3267, cMid, cMid56, p6, cMid67, cMid5678)) # NE bottom child
        cubes.append((cMid4158, cMid, cMid4378, cMid48, cMid58, cMid5678, cMid87, p8)) # SW bottom child
        cubes.append((cMid, cMid3267, cMid37, cMid4378, cMid5678, cMid67, p7, cMid87)) # SE bottom child
        
        for n in range(len(cubes)):
            span = self.division_criterion(cubes[n], self.inImage, self.outImage)

            if span == True:
#                print 'criterion is TRUE'
                self.children[n] = self.getinstance(cubes[n], self.inImage, self.outImage,imageSize)
                self.children[n].index = str(convert_to_base_4(tomorton(self.children[n].i, self.children[n].j)))
                diff_level = abs(len(self.children[n].index) - self.children[n].depth)
                if diff_level != 0:
                    self.children[n].index = '0'*diff_level + self.children[n].index
                
                p1r,p2r,p3r,p4r,p5r,p6r,p7r,p8r = cubes[n]
                L1 = linear_search(self.inImage,p1r,p2r)
                L2 = linear_search(self.inImage,p2r,p3r)
                L3 = linear_search(self.inImage,p4r,p3r)
                L4 = linear_search(self.inImage,p1r,p4r)
                
                L5 = linear_search(self.inImage,p5r,p6r)
                L6 = linear_search(self.inImage,p6r,p7r)
                L7 = linear_search(self.inImage,p8r,p7r)
                L8 = linear_search(self.inImage,p5r,p8r)
                
                L9 = linear_search(self.inImage,p1r,p5r);
                L10 = linear_search(self.inImage,p2r,p6r);
                L11 = linear_search(self.inImage,p3r,p7r);
                L12 = linear_search(self.inImage,p4r,p8r);
# 
#                if len(L1) == 1:
#                    L1 = L1[0]
#                    if in_child_k(cubes[n],L1) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L1]
#                if len(L2) == 1:
#                    L2 = L2[0]
#                    if in_child_k(cubes[n],L2) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L2]
#                if len(L3) == 1:
#                    L3 = L3[0]
#                    if in_child_k(cubes[n],L3) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L3]
#                if len(L4) == 1:
#                    L4 = L4[0]
#                    if in_child_k(cubes[n],L4) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L4]
#                        
#                if len(L5) == 1:
#                    L5 = L5[0]
#                    if in_child_k(cubes[n],L5) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L5]
#                if len(L6) == 1:
#                    L6 = L6[0]
#                    if in_child_k(cubes[n],L6) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L6]
#                if len(L7) == 1:
#                    L7 = L7[0]
#                    if in_child_k(cubes[n],L7) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L7]
#                if len(L8) == 1:
#                    L8 = L8[0]
#                    if in_child_k(cubes[n],L8) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L8]
#                        
#                if len(L9) == 1:
#                    L9 = L9[0]
#                    if in_child_k(cubes[n],L9) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L9]
#                if len(L10) == 1:
#                    L10 = L10[0]
#                    if in_child_k(cubes[n],L10) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L10]
#                if len(L11) == 1:
#                    L11 = L11[0]
#                    if in_child_k(cubes[n],L11) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L11]
#                if len(L12) == 1:
#                    L12 = L12[0]
#                    if in_child_k(cubes[n],L12) == True:
#                        self.children[n].enrichNodes = self.children[n].enrichNodes + [L12]                                
            
            
        if ( self.children[0] != None and
             self.children[1] != None and
             self.children[2] != None and
             self.children[3] != None and
             self.children[4] != None and
             self.children[5] != None and
             self.children[6] != None and
             self.children[7] != None             
             ):
            self.has_children = True
        else:
            self.has_children = False 
            
            
            
    def printcube(self):
        print [self.cube[0].x, self.cube[1].x, self.cube[0].y, self.cube[3].y, self.cube[0].z, self.cube[4].z]

    def getinstance(self, cube, inImage, outImage):
        return Node(self, cube, inImage, outImage)
    
    def division_criterion(self, cube, inImage, outImage):
        return False        
                
    def division_criterionOnce(self, cube, inImage, outImage):
        return False
    
class CNode(Node):
    
    def getinstance(self,cube,inImage,outImage,imageSize):
        return CNode(self,cube,inImage,outImage,imageSize)


    def division_criterion(self, cube, inImage, outImage):
        
        p1,p2,p3,p4,p5,p6,p7,p8 = self.cube
        cMid12 = find_mid_point(p1,p2)
        cMid14 = find_mid_point(p1,p4)
        cMid23 = find_mid_point(p2,p3)
        cMid34 = find_mid_point(p3,p4)
        
        
        cMid15 = find_mid_point(p1,p5)
        cMid26 = find_mid_point(p2,p6)
        cMid37 = find_mid_point(p3,p7)
        cMid48 = find_mid_point(p4,p8)
        
        cMid56 = find_mid_point(p5,p6)
        cMid67 = find_mid_point(p6,p7)
        cMid87 = find_mid_point(p8,p7)
        cMid58 = find_mid_point(p5,p8)
        
        cMid1234 = find_mid_point(p2,p4)
        cMid5678 = find_mid_point(p6,p8)
        cMid3267 = find_mid_point(p2,p7)
        cMid4158 = find_mid_point(p1,p8)
        
        cMid1265 = find_mid_point(p2,p5)
        cMid4378 = find_mid_point(p3,p8)
        cMid = find_mid_point(p2,p8) # center of the cube
        
        
#        print abs(p1.x-p2.x), (p1.y-p4.y), (p1.z-p5.z)
        
        if abs(p1.x - p2.x) >= MAX_SIZE_X and abs(p1.y - p4.y) >= MAX_SIZE_Y or abs(p1.z - p5.z) >= MAX_SIZE_Z:          
            
            draw_line(self.outImage, cMid12 , cMid1234 )
            draw_line(self.outImage, cMid1234 , cMid34 )
            draw_line(self.outImage, cMid14 , cMid1234 )
            draw_line(self.outImage, cMid1234 , cMid23 )
            
            draw_line(self.outImage, cMid23 , cMid3267 )
            draw_line(self.outImage, cMid3267 , cMid67 )
            draw_line(self.outImage, cMid37 , cMid3267 )
            draw_line(self.outImage, cMid3267 , cMid26 )
            
            draw_line(self.outImage, cMid12 , cMid1265 )
            draw_line(self.outImage, cMid1265 , cMid56 )
            draw_line(self.outImage, cMid15 , cMid1265 )
            draw_line(self.outImage, cMid1265 , cMid26 )
            
            draw_line(self.outImage, cMid14 , cMid4158 )
            draw_line(self.outImage, cMid4158 , cMid58 )
            draw_line(self.outImage, cMid48 , cMid4158 )
            draw_line(self.outImage, cMid4158 , cMid15 )

            draw_line(self.outImage, cMid34 , cMid4378 )
            draw_line(self.outImage, cMid4378 , cMid87 )
            draw_line(self.outImage, cMid48 , cMid4378 )
            draw_line(self.outImage, cMid4378 , cMid37 )
            
            draw_line(self.outImage, cMid58 , cMid5678 )
            draw_line(self.outImage, cMid5678 , cMid67 )
            draw_line(self.outImage, cMid87 , cMid5678 )
            draw_line(self.outImage, cMid5678 , cMid56 )
                        
            
            draw_line(self.outImage, cMid1234 , cMid )
            draw_line(self.outImage, cMid , cMid5678 )
            draw_line(self.outImage, cMid4158 , cMid )
            draw_line(self.outImage, cMid , cMid3267 )
            draw_line(self.outImage, cMid4378 , cMid )
            draw_line(self.outImage, cMid , cMid1265 )
            
            return True
        else:
            pxVal1 = self.inImage.GetPixel(p1.x,p1.y,p1.z)
            pxVal2 = self.inImage.GetPixel(p2.x,p2.y,p2.z)
            pxVal3 = self.inImage.GetPixel(p3.x,p3.y,p3.z)
            pxVal4 = self.inImage.GetPixel(p4.x,p4.y,p4.z)
            pxVal5 = self.inImage.GetPixel(p5.x,p5.y,p5.z)
            pxVal6 = self.inImage.GetPixel(p6.x,p6.y,p6.z)
            pxVal7 = self.inImage.GetPixel(p7.x,p7.y,p7.z)
            pxVal8 = self.inImage.GetPixel(p8.x,p8.y,p8.z)
            
            # are the 8 corners of the element in the same bin? i.e. homogeneous?
            isHomogeneous = eight_corners_test(pxVal1, pxVal2, pxVal3, pxVal4, pxVal5, pxVal6, pxVal7, pxVal8)

            # if the eight corners test fails
            if ( ( isHomogeneous == 0) and has_inclusions(self.inImage,p1,p2,p3,p4,p5) 
#                 and (abs(p1.x-p2.x) >= MIN_SIZE_X and abs(p1.y-p4.y) >= MIN_SIZE_Y and abs(p1.z - p5.z) >= MIN_SIZE_Z)
                 ):
                
                l1 = ends_in_same_bin(self.inImage,p1,p2)
                l2 = ends_in_same_bin(self.inImage,p2,p3)
                l3 = ends_in_same_bin(self.inImage,p4,p3)
                l4 = ends_in_same_bin(self.inImage,p1,p4)
                l5 = ends_in_same_bin(self.inImage,p5,p6)
                l6 = ends_in_same_bin(self.inImage,p6,p7)
                l7 = ends_in_same_bin(self.inImage,p8,p7)
                l8 = ends_in_same_bin(self.inImage,p5,p8)
                l9 = ends_in_same_bin(self.inImage,p1,p5)
                l10 = ends_in_same_bin(self.inImage,p2,p6)
                l11 = ends_in_same_bin(self.inImage,p3,p7)
                l12 = ends_in_same_bin(self.inImage,p4,p8)
            
                L1 = linear_search(self.inImage,p1,p2)
                L2 = linear_search(self.inImage,p2,p3)
                L3 = linear_search(self.inImage,p4,p3)
                L4 = linear_search(self.inImage,p1,p4)
                
                L5 = linear_search(self.inImage,p5,p6)
                L6 = linear_search(self.inImage,p6,p7)
                L7 = linear_search(self.inImage,p8,p7)
                L8 = linear_search(self.inImage,p5,p8)
                
                L9 = linear_search(self.inImage,p1,p5)
                L10 = linear_search(self.inImage,p2,p6)
                L11 = linear_search(self.inImage,p3,p7)
                L12 = linear_search(self.inImage,p4,p8)

                # list of coordinates:
                x_list_c = []
                y_list_c = []
                z_list_c = []
                

                if ( 
                    len(L2) > 1 or len(L4) > 1 or len(L1) > 1 or len(L3) > 1
                    or len(L5) > 1 or len(L6) > 1 or len(L7) > 1 or len(L8) > 1
                    or len(L9) > 1 or len(L10) > 1 or len(L11) > 1 or len(L12) > 1 
#                     or len(L2) < 1 or len(L4) < 1 or len(L1) < 1 or len(L3) < 1
#                     or len(L5) < 1 or len(L6) < 1 or len(L7) < 1 or len(L8) < 1
#                     or len(L9) < 1 or len(L10) < 1 or len(L11) < 1 or len(L12) < 1
                     ):
                    
                    # interface croses one edge multiple times
                    draw_line(self.outImage, cMid12 , cMid1234 )
                    draw_line(self.outImage, cMid1234 , cMid34 )
                    draw_line(self.outImage, cMid14 , cMid1234 )
                    draw_line(self.outImage, cMid1234 , cMid23 )
                    
                    draw_line(self.outImage, cMid23 , cMid3267 )
                    draw_line(self.outImage, cMid3267 , cMid67 )
                    draw_line(self.outImage, cMid37 , cMid3267 )
                    draw_line(self.outImage, cMid3267 , cMid26 )
                    
                    draw_line(self.outImage, cMid12 , cMid1265 )
                    draw_line(self.outImage, cMid1265 , cMid56 )
                    draw_line(self.outImage, cMid15 , cMid1265 )
                    draw_line(self.outImage, cMid1265 , cMid26 )
                    
                    draw_line(self.outImage, cMid14 , cMid4158 )
                    draw_line(self.outImage, cMid4158 , cMid58 )
                    draw_line(self.outImage, cMid48 , cMid4158 )
                    draw_line(self.outImage, cMid4158 , cMid15 )
        
                    draw_line(self.outImage, cMid34 , cMid4378 )
                    draw_line(self.outImage, cMid4378 , cMid87 )
                    draw_line(self.outImage, cMid48 , cMid4378 )
                    draw_line(self.outImage, cMid4378 , cMid37 )
                    
                    draw_line(self.outImage, cMid58 , cMid5678 )
                    draw_line(self.outImage, cMid5678 , cMid67 )
                    draw_line(self.outImage, cMid87 , cMid5678 )
                    draw_line(self.outImage, cMid5678 , cMid56 )
                                
                    
                    draw_line(self.outImage, cMid1234 , cMid )
                    draw_line(self.outImage, cMid , cMid5678 )
                    draw_line(self.outImage, cMid4158 , cMid )
                    draw_line(self.outImage, cMid , cMid3267 )
                    draw_line(self.outImage, cMid4378 , cMid )
                    draw_line(self.outImage, cMid , cMid1265 )
                
                    return True
                
                else:
                    if len(L1) == 1:
                        L1 = L1[0]
                        x_list_c.append(L1.x)
                        y_list_c.append(L1.y)
                        z_list_c.append(L1.z)
                    
                    if len(L2) == 1:
                        L2 = L2[0]
                        x_list_c.append(L2.x)
                        y_list_c.append(L2.y)
                        z_list_c.append(L2.z)
                    
                    if len(L3) == 1:
                        L3 = L3[0]
                        x_list_c.append(L3.x)
                        y_list_c.append(L3.y)
                        z_list_c.append(L3.z)
                    if len(L4) == 1:
                        L4 = L4[0]
                        x_list_c.append(L4.x)
                        y_list_c.append(L4.y)
                        z_list_c.append(L4.z)
                    if len(L5) == 1:
                        L5 = L5[0]
                        x_list_c.append(L5.x)
                        y_list_c.append(L5.y)
                        z_list_c.append(L5.z)
                    if len(L6) == 1:
                        L6 = L6[0]           
                        x_list_c.append(L6.x)
                        y_list_c.append(L6.y)
                        z_list_c.append(L6.z)                                     
                    if len(L7) == 1:
                        L7 = L7[0]
                        x_list_c.append(L7.x)
                        y_list_c.append(L7.y)
                        z_list_c.append(L7.z)
                    if len(L8) == 1:
                        L8 = L8[0]
                        x_list_c.append(L8.x)
                        y_list_c.append(L8.y)
                        z_list_c.append(L8.z)
                    if len(L9) == 1:
                        L9 = L9[0]
                        x_list_c.append(L9.x)
                        y_list_c.append(L9.y)
                        z_list_c.append(L9.z)
                    if len(L10) == 1:
                        L10 = L10[0]
                        x_list_c.append(L10.x)
                        y_list_c.append(L10.y)
                        z_list_c.append(L10.z)
                    if len(L11) == 1:
                        L11 = L11[0]
                        x_list_c.append(L11.x)
                        y_list_c.append(L11.y)
                        z_list_c.append(L11.z)
                    if len(L12) == 1:
                        L12 = L12[0]
                        x_list_c.append(L12.x)
                        y_list_c.append(L12.y)
                        z_list_c.append(L12.z)


                if len(x_list_c) > 0 and  len(y_list_c) > 0 and len(z_list_c) >0:
                        res = calc_plane_res(x_list_c, y_list_c, z_list_c)
                        if res <= 0.001:
                            draw_plane_connections(self.outImage, l1,l2,l3,l4, L1,L2,L3,L4) # 1234
                            draw_plane_connections(self.outImage, l1,l2,l6,l5, L1,L2,L6,L5) # 1265
                            draw_plane_connections(self.outImage, l3,l2,l6,l7, L3,L2,L6,L7) # 3267
                            draw_plane_connections(self.outImage, l4,l3,l7,l8, L4,L3,L7,L8) # 4378
                            draw_plane_connections(self.outImage, l4,l1,l5,l8, L4,L1,L5,L8) # 4158
                            draw_plane_connections(self.outImage, l5,l6,l7,l8, L5,L6,L7,L8) # 5678
                                
                            
                        else: # residual criterion approximation of the interface is not met
                            if (abs(p1.x-p2.x) >= 2*MIN_SIZE_X and abs(p1.y-p4.y) >= 2*MIN_SIZE_Y and abs(p1.z - p5.z)>= 2*MIN_SIZE_Z):
                                draw_line(self.outImage, cMid12 , cMid1234 )
                                draw_line(self.outImage, cMid1234 , cMid34 )
                                draw_line(self.outImage, cMid14 , cMid1234 )
                                draw_line(self.outImage, cMid1234 , cMid23 )
                                
                                draw_line(self.outImage, cMid23 , cMid3267 )
                                draw_line(self.outImage, cMid3267 , cMid67 )
                                draw_line(self.outImage, cMid37 , cMid3267 )
                                draw_line(self.outImage, cMid3267 , cMid26 )
                                
                                draw_line(self.outImage, cMid12 , cMid1265 )
                                draw_line(self.outImage, cMid1265 , cMid56 )
                                draw_line(self.outImage, cMid15 , cMid1265 )
                                draw_line(self.outImage, cMid1265 , cMid26 )
                                
                                draw_line(self.outImage, cMid14 , cMid4158 )
                                draw_line(self.outImage, cMid4158 , cMid58 )
                                draw_line(self.outImage, cMid48 , cMid4158 )
                                draw_line(self.outImage, cMid4158 , cMid15 )
                    
                                draw_line(self.outImage, cMid34 , cMid4378 )
                                draw_line(self.outImage, cMid4378 , cMid87 )
                                draw_line(self.outImage, cMid48 , cMid4378 )
                                draw_line(self.outImage, cMid4378 , cMid37 )
                                
                                draw_line(self.outImage, cMid58 , cMid5678 )
                                draw_line(self.outImage, cMid5678 , cMid67 )
                                draw_line(self.outImage, cMid87 , cMid5678 )
                                draw_line(self.outImage, cMid5678 , cMid56 )
                                            
                                
                                draw_line(self.outImage, cMid1234 , cMid )
                                draw_line(self.outImage, cMid , cMid5678 )
                                draw_line(self.outImage, cMid4158 , cMid )
                                draw_line(self.outImage, cMid , cMid3267 )
                                draw_line(self.outImage, cMid4378 , cMid )
                                draw_line(self.outImage, cMid , cMid1265 )
                                return True
                    
#                # NW
#                if (l1==0 and l2==1 and l3==1 and l4==0) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
#                    draw_line(self.outImage, L1, L4)
#                    self.ishomog = 0
#                    
#                # NE
#                if (l1==0 and l2==0 and l3==1 and l4==1) and (abs(p1.x-p2.x) < 2*MIN_SIZE):
#                    
#                    draw_line(self.outImage, L1, L2)
#                    self.ishomog = 0
#                # SE
#                if(l1==1 and l2==0 and l3==0 and l4==1) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
#                    draw_line(self.outImage, L2, L3)
#                    self.ishomog = 0
#                # SW
#                if (l1==1 and l2==1 and l3==0 and l4==0) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
#                    draw_line(self.outImage, L3, L4)
#                    self.ishomog = 0
#                # vertical
#                if (l1==0 and l2==1 and l3==0 and l4==1) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
#                    draw_line(self.outImage, L1, L3)
#                    self.ishomog = 0
#                # horizontal
#                if (l1==1 and l2==0 and l3==1 and l4==0) and (abs(p1.x-p2.x) < 2*MIN_SIZE) :
#                    draw_line(self.outImage, L4, L2)
#                    self.ishomog = 0


                
        return False
    
    def division_criterionOnce(self, cube, inImage, outImage):
        
        p1,p2,p3,p4,p5,p6,p7,p8 = self.cube
        cMid12 = find_mid_point(p1,p2)
        cMid14 = find_mid_point(p1,p4)
        cMid23 = find_mid_point(p2,p3)
        cMid34 = find_mid_point(p3,p4)
        
        
        cMid15 = find_mid_point(p1,p5)
        cMid26 = find_mid_point(p2,p6)
        cMid37 = find_mid_point(p3,p7)
        cMid48 = find_mid_point(p4,p8)
        
        cMid56 = find_mid_point(p5,p6)
        cMid67 = find_mid_point(p6,p7)
        cMid87 = find_mid_point(p8,p7)
        cMid58 = find_mid_point(p5,p8)
        
        cMid1234 = find_mid_point(p2,p4)
        cMid5678 = find_mid_point(p6,p8)
        cMid3267 = find_mid_point(p2,p7)
        cMid4158 = find_mid_point(p1,p8)
        
        cMid1265 = find_mid_point(p2,p5)
        cMid4378 = find_mid_point(p3,p8)
        cMid = find_mid_point(p2,p8) # center of the cube
            
        draw_line(self.outImage, cMid12 , cMid1234 )
        draw_line(self.outImage, cMid1234 , cMid34 )
        draw_line(self.outImage, cMid14 , cMid1234 )
        draw_line(self.outImage, cMid1234 , cMid23 )
        
        draw_line(self.outImage, cMid23 , cMid3267 )
        draw_line(self.outImage, cMid3267 , cMid67 )
        draw_line(self.outImage, cMid37 , cMid3267 )
        draw_line(self.outImage, cMid3267 , cMid26 )
        
        draw_line(self.outImage, cMid12 , cMid1265 )
        draw_line(self.outImage, cMid1265 , cMid56 )
        draw_line(self.outImage, cMid15 , cMid1265 )
        draw_line(self.outImage, cMid1265 , cMid26 )
        
        draw_line(self.outImage, cMid14 , cMid4158 )
        draw_line(self.outImage, cMid4158 , cMid58 )
        draw_line(self.outImage, cMid48 , cMid4158 )
        draw_line(self.outImage, cMid4158 , cMid15 )

        draw_line(self.outImage, cMid34 , cMid4378 )
        draw_line(self.outImage, cMid4378 , cMid87 )
        draw_line(self.outImage, cMid48 , cMid4378 )
        draw_line(self.outImage, cMid4378 , cMid37 )
        
        draw_line(self.outImage, cMid58 , cMid5678 )
        draw_line(self.outImage, cMid5678 , cMid67 )
        draw_line(self.outImage, cMid87 , cMid5678 )
        draw_line(self.outImage, cMid5678 , cMid56 )
                    
        
        draw_line(self.outImage, cMid1234 , cMid )
        draw_line(self.outImage, cMid , cMid5678 )
        draw_line(self.outImage, cMid4158 , cMid )
        draw_line(self.outImage, cMid , cMid3267 )
        draw_line(self.outImage, cMid4378 , cMid )
        draw_line(self.outImage, cMid , cMid1265 )

        return True

class OctoTree(Node):
    maxdepth = 1
    leaves = []
    allnodes = []
    
    def __init__(self,rootnode):
#        self = self
        rootnode.subdivide() # constructs the network or nodes
#        rootnode.division_criterionOnce(self, cube, inImage, outImage)visit
         
class COctoTree(OctoTree):
    def __init__(self,rootnode):
        OctoTree.__init__(self, rootnode)
    
            
if __name__ == "__main__":
    # two_channels.dcm contains the original data
    # img.vtk contains an empty dataset of the same dimensions with the original
    # orig_mesh.vtk is two_channels.dcm converted to VTK format
    # empty_mesh.vtk contains a mesh on an empty set
    print "Reading image in..."
    inputImage = sitk.ReadImage("dataset/channels_512.dcm")
    outputImage = sitk.ReadImage("dataset/channels_512.dcm")
#    outputImage = sitk.ReadImage("dataset/img.vtk")
#    inputImage = sitk.ReadImage((sys.argv[1]));
#    outputImage = sitk.ReadImage((sys.argv[1]));

#    sitk.WriteImage(inputImage,"dataset/orig_mesh.vtk");
    
    nameOutputImage = "dataset/out_mesh.vtk" 
    
    imageSize = inputImage.GetSize()
    print "Image size:", imageSize
 
#    dim = (imageSize[1], imageSize[2], imageSize[0])
#    origin = (0, 0, 0)
#    spacing = (1.0, 1.0, 1.0)
#    img = scipy.zeros(dim)*255
#    
#    
#    img_sitk = sitk.GetImageFromArray(img)
#    img_sitk.SetOrigin(origin)
#    img_sitk.SetSpacing(spacing)
#
##    sitk.Show(img_sitk)
#    sitk.WriteImage(img_sitk,"dataset/img.vtk")
#    print img_sitk.GetSize()
#    print inputImage.GetSize()
     
    # setting the 4 corners coordinates
    p1 = Coordinate(0,0,0)
    p2 = Coordinate(imageSize[0]-1,0,0)
    p3 = Coordinate(imageSize[0]-1,imageSize[1]-1,0)
    p4 = Coordinate(0,imageSize[1]-1,0)
    p5 = Coordinate(0,0,imageSize[2]-1)
    p6 = Coordinate(imageSize[0]-1,0,imageSize[2]-1)
    p7 = Coordinate(imageSize[0]-1,imageSize[1]-1,imageSize[2]-1)
    p8 = Coordinate(0,imageSize[1]-1,imageSize[2]-1)
    
#    print p1.x, p1.y, p1.z
#    print p2.x, p2.y, p2.z
#    print p3.x, p3.y, p3.z
#    print p4.x, p4.y, p4.z
#    print p5.x, p5.y, p5.z
#    print p6.x, p6.y, p6.z
#    print p7.x, p7.y, p7.z
#    print p8.x, p8.y, p8.z

#    draw_line(outputImage,p1,p3)
#    draw_line(outputImage,p4,p7)
#    draw_line(outputImage,p3,p6)
#    draw_line(outputImage,p1,p8)
#     
#    draw_line(outputImage,p1,p6)
#    cMid12  = find_mid_point(p1,p2)
#    print cMid12.x, cMid12.y, cMid12.z
#    
    cube = [p1,p2,p3,p4,p5,p6,p7,p8]
    rootNode = CNode(None,cube,inputImage,outputImage,imageSize)
    tree = COctoTree(rootNode)
    
    masterNode = CNode(None,cube,inputImage,outputImage,imageSize)
    

#    writer = sitk.ImageFileWriter()
#    writer.SetFileName ( nameOutputImage )
#    writer.Execute ( outputImage )

    print 'writing the image out'
    sitk.WriteImage(outputImage,nameOutputImage);
   