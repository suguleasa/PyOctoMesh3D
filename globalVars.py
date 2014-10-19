#binBnd = [0,200,256] # for real slice: microvascular
#binBnd = [0,110,256] # for fake images
#binBnd = [0,140,256] # for fake doughnut images
#binBnd = [0,90,256] # for Green Circles image: green1.png
#binBnd = [0,100,256] # for SiBat.png
binBnd = [0,150,256] # for SiBat_filtered.png
#binBnd = [0,170,256] # for SnBat.png
#binBnd = [0,200,256] # for SnBat_0054.png
#binBnd = [0,50,256] # for doughnut image
#binBnd = [0,110,256]
#binBnd = [0,200,256] # for SnBat
AREA_INCLUSION = 1048.0 # or a 4x4 small inclusion of pixels
PROB = 0.3
#MAX_SIZE = 80
#MIN_SIZE = 20#6
ALT_MIN_SIZE = 40#3
TOL_LINEARS = 2.0
TOL_QUADRATICS = -2.0 # cancel the cubic approximation by making a negative number
TOL_CUBICS = -2.0 # cancel the cubic approximation by making a negative number

TOL_error = 3
T_STEP = 0.001

VAL = 0

MAX_SIZE_X = 127
MIN_SIZE_X = 50

MAX_SIZE_Y = 127
MIN_SIZE_Y = 50

MAX_SIZE_Z = 127
MIN_SIZE_Z = 50

# K-neighbor rule
k1_CONST = 5
k2_CONST = 2

# minimum number of elements between two interfaces
N_ELEMS = 4
STRESS_MIN = 4

NURBS_ON = 1

TOL_NURBS = 3.0

GRID_PTS = 2 # grid of 2 x 2