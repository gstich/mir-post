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
resolution = "medium3d";
t_viz = 350;
var = "gradRHO";
#var = "sch";

# Zoomed in of entire nozzle section
v1 = GetView2D();
v1.windowCoords = (-2.5, 20.6229, -5.72, 5.72)
v1.viewportCoords = (0,1,0,1);

# Bit more zoomed in
v2 = GetView2D();
v2.windowCoords = (0, 15.6229, -3.92, 3.92)
v2.viewportCoords = (0,1,0,1);

SetView2D(v2);

file = base + resolution + "/post.mir"
OpenDatabase(file);
SetTimeSliderState(t_viz);


DefineScalarExpression("sch","0.8*exp( -15*gradRHO/.015)");


# Add Pseudocolor plot
AddPlot("Pseudocolor",var);
psdo_atts = PseudocolorAttributes();
psdo_atts.centering = psdo_atts.Nodal;
psdo_atts.colorTableName = "xray";
psdo_atts.maxFlag = 1;
psdo_atts.max = .075;
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
win.width = 2048;
win.height = 1024;
win.family = 0;
win.fileName = resolution + str(t_viz) + var;
win.format = win.TIFF;
win.resConstraint = win.NoConstraint;
SetSaveWindowAttributes(win);


DrawPlots();
SaveWindow();






# Close the compute engine
e.close

# Again, needed when issueing -nowin option
sys.exit()





