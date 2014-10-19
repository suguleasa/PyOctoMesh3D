# Developer: Elena Caraba, copyright 2014

# File setting up the experiment


# setting the boundaries for the pixel colors bins
binBnd = [0,150,256] # for SiBat_filtered.png


VOL_INCLUSION = 1048.0 # volume of significant inclusion
PROB = 0.3 # probability an inclusion of vol_inclusion size was missed


# steps size for plotting / drawing NURBS in element
T_STEP = 0.001

VAL = 0 # pixel color for drawing

# set element max and min size
MAX_SIZE_X = 127
MIN_SIZE_X = 127

MAX_SIZE_Y = 127
MIN_SIZE_Y = 127

MAX_SIZE_Z = 127
MIN_SIZE_Z = 127

ALT_MIN_SIZE = 40 # alternative minimum size of elements


# K-neighbor rule
k1_CONST = 5
k2_CONST = 2

# minimum number of elements between two interfaces
N_ELEMS = 4
STRESS_MIN = 4

# allow NURBS approximation of interface
NURBS_ON = 1
# pixel tolerance for NURBS approximation
TOL_NURBS = 3.0

# number of interior points: GRID_PTS x GRID_PTS for Coons NURBS surfaces
GRID_PTS = 2 # grid of 2 x 2