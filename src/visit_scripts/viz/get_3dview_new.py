##  Britton's Visit python script for 3d-nozzle case
##  Used for uniform plot generation and perhaps for movies

# This is needed for running in non-window mode,
# eg. "visit -cli -nowin -s script.py"
import sys
import math


# Sets the view based on interpolating a spline
def FlyView(i, N, side):
  ## Pan-back- Starting view of whole domain
  v1 = GetView3D();
  v1.viewNormal = (0, 0, 1*side)
  #v1.viewNormal = (0, 0, -1)
  v1.focus = (6.0, 0.0, 1.7)
  v1.viewUp = (0, 1*side, 0)
  #v1.viewUp = (0, -1, 0)
  v1.viewAngle = 30
  v1.parallelScale = 51.4056
  v1.nearPlane = -102.811
  v1.farPlane = 102.811
  v1.imagePan = (0,0)
  v1.imageZoom = 10
  v1.perspective = 1
  v1.eyeAngle = 2
  v1.axis3DScaleFlag = 0
  v1.axis3DScales = (1, 1, 1)
  v1.shear = (0, 0, 1)

  # Zoom and turn
  v2 = GetView3D();
  v2.viewNormal = (-1.00848, 0.256162*side, 0.705451*side)
  #v2.viewNormal = (-0.660848, -0.256162, -0.705451)
  v2.focus = (2.2244, -1.0*side, 1.792)
  #v2.focus = (2.2244, 1.0, 1.792)
  v2.viewUp = (0, 1*side, 0)
  #v2.viewUp = (0, -1, 0)
  v2.viewAngle = 30
  v2.parallelScale = 51.4056
  v2.nearPlane = -102.811
  v2.farPlane = 102.811
  v2.imagePan = ( 0 , 0)
  v2.imageZoom = 55.5477
  v2.perspective = 1
  v2.eyeAngle = 2
  v2.axis3DScaleFlag = 0
  v2.axis3DScales = (1, 1, 1)
  v2.shear = (0, 0, 1)

  # Turn around lower BL ... zoom in direct in line
  v3 = GetView3D();
  v3.viewNormal = (-0.95412, 0.206522, 0)
  v3.focus = (3.2244, -1.0, 1.792)
  v3.viewUp = (0, 1 , 0)
  v3.viewAngle = 30
  v3.parallelScale = 51.4056
  v3.nearPlane = -102.811
  v3.farPlane = 102.811
  v3.imagePan = (0, 0)
  v3.imageZoom = 15.2593
  v3.perspective = 1
  v3.eyeAngle = 2
  v3.axis3DScaleFlag = 0
  v3.axis3DScales = (1, 1, 1)
  v3.shear = (0, 0, 1)

  ## Run along the length of the nozzle
  v4 = GetView3D();
  v4.viewNormal = (-0.60848, 0.56162*side, 0.705451*side)
  v4.focus = (7.2244, 0.0, 1.792)
  v4.viewUp = (0, 1*side, 0)
  v4.viewAngle = 30
  v4.parallelScale = 51.4056
  v4.nearPlane = -102.811
  v4.farPlane = 102.811
  v4.imagePan = ( 0 , 0)
  v4.imageZoom = 40.5477
  v4.perspective = 1
  v4.eyeAngle = 2
  v4.axis3DScaleFlag = 0
  v4.axis3DScales = (1, 1, 1)
  v4.shear = (0, 0, 1)

  ## View 5
  v5 = GetView3D();
  v5.viewNormal = (-0.60848, 0.26162*side, -0.705451*side)
  v5.focus = (8.2244, 0.0, 1.792)
  v5.viewUp = (0, 1*side, 0)
  v5.viewAngle = 30
  v5.parallelScale = 51.4056
  v5.nearPlane = -102.811
  v5.farPlane = 102.811
  v5.imagePan = ( 0 , 0)
  v5.imageZoom = 40.5477
  v5.perspective = 1
  v5.eyeAngle = 2
  v5.axis3DScaleFlag = 0
  v5.axis3DScales = (1, 1, 1)
  v5.shear = (0, 0, 1)



  v6 = GetView3D();
  v6.viewNormal = (-0.660848, 0.256162*side, -0.705451*side)
  v6.focus = (2.2244, -1.0*side, 1.792)
  v6.viewUp = (0, 1*side, 0)
  v6.viewAngle = 30
  v6.parallelScale = 51.4056
  v6.nearPlane = -102.811
  v6.farPlane = 102.811
  v6.imagePan = ( 0 , 0)
  v6.imageZoom = 45.5477
  v6.perspective = 1
  v6.eyeAngle = 2
  v6.axis3DScaleFlag = 0
  v6.axis3DScales = (1, 1, 1)
  v6.shear = (0, 0, 1)



  v7 = GetView3D();
  v7.viewNormal = (-0.40848, 0.256162*side, 0.805451*side)
  v7.focus = (5.2, -1.0*side, 2.0)
  v7.viewUp = (0, 1*side, 0)
  v7.viewAngle = 30
  v7.parallelScale = 51.4056
  v7.nearPlane = -102.811
  v7.farPlane = 102.811
  v7.imagePan = ( 0 , 0)
  v7.imageZoom = 12.5477
  v7.perspective = 1
  v7.eyeAngle = 2
  v7.axis3DScaleFlag = 0
  v7.axis3DScales = (1, 1, 1)
  v7.shear = (0, 0, 1)  



  
  vpts = (v1, v2, v4, v5, v6, v7)
  x = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
  ramp = 32;
  pos=SineParameterize(N, i, ramp)
  v = EvalCubicSpline(pos, x, vpts)
  SetView3D(v)
  
def FlyView2(i, N):
  ## Pan-back- Starting view of whole domain
  v1 = GetView3D();
  
  # Zoomed in of entire nozzle section
  v2 = GetView3D();


  vpts = (v1, v2)
  x = [0, 1.0]
  ramp = 1;
  pos=SineParameterize(N, i, ramp)
  v = EvalCubicSpline(pos, x, vpts)
  SetView2D(v)



# Import this module for launching a parallel engine
# Script must be issued on the machine which will be used
Source("visit_engine.py")
Source("utils.py")

# Launch/request the procs
e = Engine()
e.open(nprocs=100,part="pbatch",bank="views",rtime="480:00")
#e.open(nprocs=128,part="pdebug",rtime="30:00")


# Load the file
base = "/p/lscratchd/olson45/nozzle/nozzle";

del_bl = .2;
Ht = 1.78;


## Fine case fly
resolution = "fine3d";
j_cut = 159;     # fine
t_start = 233;
t_end = 233;
side = -1;     # Side= top(-1), bot(1)
zoff = 4.0215

## Medium case fly
#resolution = "medium3d";
#j_cut = 150;    # medium
#t_start = 233;
#t_end = 233;
#side = 1;     # Side= top(-1), bot(1)
#zoff = 0;

## Coarse case fly
#resolution = "coarse3d";
#j_cut = 60;
#t_start = 1800;
#t_end = 1800;
#side = -1;
#zoff = 2.667;

file = base + resolution + "/plot.mir"
OpenDatabase(file);
SetTimeSliderState(t_start);


# Define Q-citerion
DefineVectorExpression("gU","gradient(velocity[0])");
DefineVectorExpression("gV","gradient(velocity[1])");
DefineVectorExpression("gW","gradient(velocity[2])");
DefineTensorExpression("gradU","{gU,gV,gW}");
DefineTensorExpression("S","1/2*(gradU + transpose(gradU)) ");
DefineTensorExpression("Omega","1/2*(gradU - transpose(gradU)) ");
DefineScalarExpression("Q","1/2*( trace(Omega*transpose(Omega))  - trace(S*transpose(S)) )")

# Define Div / Mach
DefineScalarExpression("div","divergence(velocity)");
DefineScalarExpression("Mach","velocity_magnitude/sound_speed");

# Mag. of Grad(rho)
DefineScalarExpression("Grho","magnitude(gradient(density))");

# Add Pseudocolor plot - Mach
AddPlot("Pseudocolor","Mach");
psdo_atts = PseudocolorAttributes();
psdo_atts.centering = psdo_atts.Nodal;
psdo_atts.colorTableName = "hot_desaturated";
psdo_atts.maxFlag = 1;
#psdo_atts.max = 1.5;         ## Medium Case
psdo_atts.max = 1.75;        ## Fine Case
SetPlotOptions(psdo_atts);
#DrawPlots();

# Add the Isovolume of Q to Mach Contours
SetActivePlots(0);
isoV_atts = IsovolumeAttributes();
isoV_atts.variable = "Q"
isoV_atts.lbound = 5e9;
isoV_atts.ubound = 4e11;    ## Coarse case
#isoV_atts.ubound = 4e11;    ## Medium case
AddOperator("Isovolume");
SetOperatorOptions(isoV_atts);
DrawPlots();

# Add cut to remove one side
## Make this just a cut at
#SetActivePlots(0);
#islc_atts = IndexSelectAttributes();
#islc_atts.dim = islc_atts.ThreeD;
#islc_atts.yMax = j_cut;     # medium -  1-jcut
#islc_atts.yMin = j_cut;     # fine-     jcut-Ny
#AddOperator("IndexSelect");
#SetOperatorOptions(islc_atts);
#DrawPlots();





## VISUALIZE THE SHOCK

# Add Pseudocolor plot - Divergence
AddPlot("Pseudocolor","div");
psdo1_atts = PseudocolorAttributes();
psdo1_atts.centering = psdo1_atts.Nodal;
psdo1_atts.colorTableName = "bluehot";
psdo1_atts.opacity = 0.5;
psdo1_atts.maxFlag = 1;
psdo1_atts.max = 0;
SetPlotOptions(psdo1_atts);
#DrawPlots()




# Add the Isovolume Operator
SetActivePlots(1);
isoV1_atts = IsovolumeAttributes();
isoV1_atts.ubound = -50000;   #-8000 med
AddOperator("Isovolume");
SetOperatorOptions(isoV1_atts);
#DrawPlots();


# Clip the first bit to remove the noise
SetActivePlots(1);
clp_atts = ClipAttributes();
clp_atts.plane1Origin = (-2.0,0,0);
clp_atts.plane1Normal = (-1,0,0);
AddOperator("Clip");
SetOperatorOptions(clp_atts);
#DrawPlots();

# Clip the to to clear view
#SetActivePlots(1);
#clp1_atts = ClipAttributes();
#clp1_atts.plane1Origin = (0,0.7,0);
#clp1_atts.plane1Normal = (0,1,0);
#AddOperator("Clip");
#SetOperatorOptions(clp1_atts);
#DrawPlots();

SetActivePlots((0,1));
clp2_atts = ClipAttributes();
clp2_atts.plane1Origin = (0, side*(Ht/2.0-del_bl*2.0) , 0);
clp2_atts.plane1Normal = (0,side,0);
AddOperator("Clip");
SetOperatorOptions(clp2_atts);
DrawPlots();

# Add back plane of bw |grad(rho)|
AddPlot("Pseudocolor","Grho");
psdo3_atts = PseudocolorAttributes();
psdo3_atts.centering = psdo3_atts.Nodal;
psdo3_atts.colorTableName = "xray";
psdo3_atts.maxFlag = 1;
psdo3_atts.max = 0.005;
SetPlotOptions(psdo3_atts);

SetActivePlots(2)
slc_atts = SliceAttributes();
slc_atts.project2d = 0;
slc_atts.originType = slc_atts.Intercept;
slc_atts.originIntercept = zoff;
slc_atts.axisType = slc_atts.ZAxis;
AddOperator("Slice");
SetOperatorOptions(slc_atts);


# Clean up the annotations
annot = AnnotationAttributes();
annot.legendInfoFlag = 0;
annot.userInfoFlag = 0;
annot.databaseInfoFlag = 0;
annot.axes3D.visible = 0;
annot.axes3D.bboxFlag = 0;
#annot.backgroundColor = (0, 0, 0, 255);
#annot.foregroundColor = (255, 255, 255, 255);
#annot.backgroundMode = annot.Gradient;
#annot.gradientColor1 = (0, 0, 127, 255);
#annot.gradientColor2 = (0, 0, 0, 255);
SetAnnotationAttributes(annot)


# Set the File Output type
win = SaveWindowAttributes();
win.width = 1920*5;        ## High-Def
win.height = 1080*5;       ## High-Def
win.fileName = "HIGH_" + resolution;
win.outputToCurrentDirectory = 0;
win.outputDirectory = "/p/lscratchd/olson45/NOZ_VIZ/"
win.format = win.TIFF;
win.quality = 100;
win.resConstraint = win.NoConstraint;
SetSaveWindowAttributes(win);


## Stage-1: Zoom in from full mesh to nozzle local -- 5sec -> 5*20 = 100 frames
nF = 2000
off = 0

FlyView(1999,2000,-1)
DrawPlots();
SaveWindow();

#for i in range(nF-off):
#  ii = i + off
#  FlyView(ii,nF,side);
#  DrawPlots();
#  SaveWindow();
  
## Run through the data
#for i in range(t_end-t_start+1):
#  SetTimeSliderState(t_start+i);
#  DrawPlots();
#  SaveWindow();
  



# Close the compute engine
e.close

# Again, needed when issueing -nowin option
sys.exit()





