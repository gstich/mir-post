##  Britton's Visit python script for 2d-nozzle case
##  Used for uniform plot generation and perhaps for movies

# This is needed for running in non-window mode,
# eg. "visit -cli -nowin -s script.py"
import sys

import math



# Import this module for launching a parallel engine
# Script must be issued on the machine which will be used
Source("visit_engine.py")
Source("utils.py")

# Launch/request the procs
e = Engine()
e.open(nprocs=128,part="pbatch",bank="views",rtime="120:00")


# Load the file
base = "/p/lscratchd/olson45/nozzle/nozzle";


## Coarse presets
toff = 600;       # Offset for the 2d data set
t1 = 600;         # Start index for viz
tf = 2258;        # End index for viz
t_viz = 2258;      # Start index for this render run...
resolution = "coarse3d";
odir = "/p/lscratchd/olson45/NOZ_VIZ/2d_cor/"
j_cut = 128-11+1;
zmid = 1.344;
i1 = 10;
iN = 450;

## Medium presets
#toff = 200;       # Offset for the 2d data set
#t1 = 200;         # Start index for viz
#tf = 1128;        # End index for viz
#t_viz = 1121;      # Start index for this render run...
#resolution = "medium3d";
#odir = "/p/lscratchd/olson45/NOZ_VIZ/2d_med/"
jcut = 12;  # y+=19
i1 = 10;
iN = 650;

## Fine presets
toff = 0;       # Offset for the 2d data set
t1 = 0;         # Start index for viz
tf = 230;        # End index for viz
t_viz = 230;      # Start index for this render run...
resolution = "fine3d";
odir = "/p/lscratchd/olson45/NOZ_VIZ/2d_fin/"
jcut = 372;#386 - 13 + 1;
zmid = 2.0;
i1 = 10;
iN = 920;
side = 1;

var = "Grho";
#var = "density";
#var = "vort";
#var = "w";

Nx = 1920;
Ny = 500;  # 1080/3
AR = Nx/Ny;
Ht = 1.78;

DefineScalarExpression("Grho","magnitude(gradient(density))");
DefineScalarExpression("vort","magnitude(curl(velocity))");
DefineScalarExpression("w","velocity[2]");

file = base + resolution + "/plot.mir"
OpenDatabase(file);
SetTimeSliderState(t_viz-toff);


# Add Pseudocolor plot
AddPlot("Pseudocolor",var);
psdo_atts = PseudocolorAttributes();
psdo_atts.centering = psdo_atts.Nodal;
psdo_atts.colorTableName = "orangehot";
psdo_atts.maxFlag = 1;
psdo_atts.minFlag = 1;
#psdo_atts.max = 5e6;#.01;  #.075 (sch)
psdo_atts.max = .05;
psdo_atts.min = 1e-6;
psdo_atts.scaling = psdo_atts.Log;
SetPlotOptions(psdo_atts);
#DrawPlots();

# Add slice using ijk-coordinates
islc_atts = IndexSelectAttributes();
islc_atts.dim = islc_atts.ThreeD;
islc_atts.xMin = i1;
islc_atts.xMax = iN;
islc_atts.yMax = j_cut;     # medium -  1-jcut
islc_atts.yMin = j_cut;     # fine-     jcut-Ny
AddOperator("IndexSelect");
SetOperatorOptions(islc_atts);
DrawPlots();

# Zoomed in of entire nozzle section
wide = 2*Ht;
x1 = -.1*Ht;
xn = 2*Ht*AR + x1;
v1 = GetView3D();
#v1.windowCoords = (x1, xn, -wide/2, wide/2)
#v1.viewportCoords = (0,1,0,1);
v1.viewNormal = ( 0 , side*1 , 0);
v1.focus = (6.75-2, side*1, zmid);
v1.viewUp = (0,0,-1*side);
v1.imageZoom = 35;
SetView3D(v1);


# Clean up the annotations
annot = AnnotationAttributes();
annot.legendInfoFlag = 0;
annot.userInfoFlag = 0;
annot.databaseInfoFlag = 0;
annot.triadFlag = 0;
annot.axes3D.visible = 0;
annot.axes3D.bboxFlag = 0;
SetAnnotationAttributes(annot)
#InvertBackgroundColor();

# Set the File Output type
win = SaveWindowAttributes();
win.width = Nx;
win.height = Ny;
win.family = 0;
win.outputToCurrentDirectory = 0;
win.outputDirectory = odir;
win.fileName = resolution + str(t_viz) + var;
win.format = win.TIFF;
win.resConstraint = win.NoConstraint;
SetSaveWindowAttributes(win);

Nt = tf - t_viz + 1
for tt in range(Nt):

    td = tt + t_viz - toff;  # Data/time slider index
    ti = tt + t_viz - t1
    SetTimeSliderState(td);
    win.fileName = resolution + var + '_' + str(ti).zfill(4);
    SetSaveWindowAttributes(win);
    DrawPlots();
    SaveWindow();



# Close the compute engine
e.close

# Again, needed when issueing -nowin option
sys.exit()





