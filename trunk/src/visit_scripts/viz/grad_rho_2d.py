##  Britton's Visit python script for 2d-nozzle case
##  Used for uniform plot generation and perhaps for movies

# This is needed for running in non-window mode,
# eg. "visit -cli -nowin -s script.py"
import sys

import math


# Sets the view based on interpolating a spline
def FlyView(i, N):
  ## Pan-back- Starting view of whole domain
  v1 = GetView2D();
  v1.windowCoords = (-60.5, 55.6229, -43.72, 43.72)
  v1.viewportCoords = (0,1,0,1);
  
  # Zoomed in of entire nozzle section
  v2 = GetView2D();
  v2.windowCoords = (-2.5, 20.6229, -5.72, 5.72)
  v2.viewportCoords = (0,1,0,1);

  vpts = (v1, v2)
  x = [0, 1.0]
  ramp = 32;
  pos=SineParameterize(N, i, ramp)
  v = EvalCubicSpline(pos, x, vpts)
  v.fullFrameActivationMode = v.Off;
  SetView2D(v)

def FlyView2(i, N):

  # Zoomed in of entire nozzle section
  v1 = GetView2D();
  v1.windowCoords = (-2.5, 20.6229, -5.72, 5.72)
  v1.viewportCoords = (0,1,0,1);

  # Bit more zoomed in
  v2 = GetView2D();
  v2.windowCoords = (0, 15.6229, -3.92, 3.92)
  v2.viewportCoords = (0,1,0,1);
  
  vpts = (v1, v2)
  x = [0, 1.0]
  ramp = 32;
  pos=SineParameterize(N, i, ramp)
  v = EvalCubicSpline(pos, x, vpts)
  v.fullFrameActivationMode = v.Off;
  SetView2D(v)


# Import this module for launching a parallel engine
# Script must be issued on the machine which will be used
Source("visit_engine.py")
Source("utils.py")

# Launch/request the procs
#e = Engine()
#e.open(nprocs=32,part="pdebug",bank="views",rtime="120:00")


# Load the file
base = "/p/lscratchd/olson45/nozzle/post_proc";
resolution = "medium3d";
t_start = 0;
t_end = 100;
#resolution = "coarse3d";
#t_start = 500;
#t_end = 800;

file = base + resolution + "/mean.mir"
OpenDatabase(file);
SetTimeSliderState(t_start);


# Add Pseudocolor plot
AddPlot("Pseudocolor","gradRHO");
psdo_atts = PseudocolorAttributes();
psdo_atts.centering = psdo_atts.Nodal;
psdo_atts.colorTableName = "xray";
psdo_atts.maxFlag = 1;
psdo_atts.max = .01;
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
win.fileName = "grad";
win.format = win.TIFF;
win.resConstraint = win.NoConstraint;
SetSaveWindowAttributes(win);


## Stage-1: Zoom in from full mesh to nozzle local -- 5sec -> 5*20 = 100 frames
nF = 100
for i in range(nF):
    FlyView(i,nF);
    DrawPlots();
    SaveWindow();


## Stage-2: Loop through the time series
SS = 5;
nS = math.floor((t_end-t_start+1)/SS);
for t in range( nS ):
    tt = t*SS
    SetTimeSliderState(tt);
    DrawPlots();
    SaveWindow();


## Stage-3 Turn on U and Contour (sep. pt
## and turn off grad_rho using opacities

# Hide the first plot
SetActivePlots(0);
HideActivePlots();

# Add Pseudocolor plot
AddPlot("Pseudocolor","uu");
psdo_atts1 = PseudocolorAttributes();
psdo_atts1.centering = psdo_atts.Nodal;
psdo_atts1.colorTableName = "hot_desaturated";
SetActivePlots(1)
SetPlotOptions(psdo_atts1);
DrawPlots();
#nF = 100  ## Not working to dynamically chnage opacity...?
#for i in range(nF):
#    opaque = (i-1.0)/(nF-1.0);
#    SetActivePlots(0);
#    psdo_atts.opacity = opaque;
#    SetPlotOptions(psdo_atts);
#    DrawPlots();
#    SetActivePlots(1);
#    psdo_atts1.opacity = 1-opaque;
#    SetPlotOptions(psdo_atts1);
#    DrawPlots();
#    SaveWindow();


## Stage-3: Loop again at (SS)x speed
## Maybe plot Mach + seperation contour
SS = 1;
nS =  math.floor((t_end-t_start+1)/SS);
for t in range( nS ):
    tt = t*SS;
    FlyView2(t,nS);
    SetTimeSliderState(tt);
    DrawPlots();
    SaveWindow();

# Close the compute engine
#e.close

# Again, needed when issueing -nowin option
#sys.exit()





