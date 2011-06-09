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
e.open(nprocs=32,part="pbatch",bank="views",rtime="120:00")


# Load the file
base = "/p/lscratchd/olson45/nozzle/post_proc";


## Coarse presets
toff = 600;       # Offset for the 2d data set
t1 = 600;         # Start index for viz
tf = 2258;        # End index for viz
t_viz = 600;      # Start index for this render run...
resolution = "coarse3d";
odir = "/p/lscratchd/olson45/NOZ_VIZ/2d_cor/"

## Medium presets
#toff = 200;       # Offset for the 2d data set
#t1 = 200;         # Start index for viz
#tf = 1128;        # End index for viz
#t_viz = 1121;      # Start index for this render run...
#resolution = "medium3d";
#odir = "/p/lscratchd/olson45/NOZ_VIZ/2d_med/"

var = "gradRHO";


Nx = 1920;
Ny = 360;  # 1080/3
AR = Nx/Ny;
Ht = 1.78;

# Zoomed in of entire nozzle section
wide = 2*Ht;
x1 = -.1*Ht;
xn = 2*Ht*AR + x1;
v1 = GetView2D();
v1.windowCoords = (x1, xn, -wide/2, wide/2)
v1.viewportCoords = (0,1,0,1);

# Bit more zoomed in
#v2 = GetView2D();
#v2.windowCoords = (0, 15.6229, -3.92, 3.92)
#v2.viewportCoords = (0,1,0,1);

SetView2D(v1);

file = base + resolution + "/post.mir"
OpenDatabase(file);
SetTimeSliderState(t_viz-toff);


DefineScalarExpression("sch","0.8*exp( -15*gradRHO/.015)");


# Add Pseudocolor plot
AddPlot("Pseudocolor",var);
psdo_atts = PseudocolorAttributes();
psdo_atts.centering = psdo_atts.Nodal;
psdo_atts.colorTableName = "xray";
psdo_atts.maxFlag = 1;
psdo_atts.max = .005;  #.075 (sch)x
SetPlotOptions(psdo_atts);
DrawPlots();


# Clean up the annotations
annot = AnnotationAttributes();
annot.legendInfoFlag = 0;
annot.userInfoFlag = 0;
annot.databaseInfoFlag = 0;
annot.axes2D.visible = 0;
SetAnnotationAttributes(annot)
InvertBackgroundColor();

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





