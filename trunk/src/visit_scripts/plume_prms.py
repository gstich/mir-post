# This is needed for running in non-window mode,
# eg. "visit -cli -nowin -s script.py"
import sys


import math, random
from math import sqrt;

# Import this module for launching a parallel engine
# Script must be issued on the machine which will be used
Source("visit_engine.py")

# Launch/request the procs
e = Engine()
e.open(nprocs=32,part="pbatch",bank="views",rtime="30:00")
 


def get_data(VAR,time):
    
    ofile = "%04d.dat" % time
    ofile = VAR + "_" + ofile
    
    # Get the data
    SetActiveWindow(2)
    SetActivePlots(0)
    vals = GetPlotInformation()["Curve"]
    SetActivePlots(1)
    xc = GetPlotInformation()["Curve"]
    SetActivePlots(2)
    yc = GetPlotInformation()["Curve"]
    SetActivePlots(3)
    zc = GetPlotInformation()["Curve"]
    
    # Write it as "x y z val"
    f = open(ofile, "wt")
    #f.write("# samples\n")
    for i in range(len(vals) / 2):
        idx = i*2+1
        f.write("%g %g %g %g\n" % (xc[idx], yc[idx], zc[idx], vals[idx]))
        
    f.close()


############################################################################
############################################################################
############################################################################
## Actual script.....

DB = "/p/lscratchd/olson45/nozzle/nozzlemedium3d/plot.mir";
t0 = 350;
tf = 350;
Var = "pressure";
Mesh = "mesh";

Ht = 2.23; #cm
x0 = 12.23; #cm
y0 = 0.0;  #cm
z0 = 0.0;  #cm

xp = 0.5 * Ht;
yp = 2.0 * Ht;

p0 = (x0 + xp, y0 - yp, z0);
p1 = (x0 + xp, y0 + yp, z0);

npts = 100;

LOA = GetGlobalLineoutAttributes();
LOA.numSamples = npts;
SetGlobalLineoutAttributes(LOA);


OpenDatabase(DB)
SetTimeSliderState(t0)
AddPlot("Pseudocolor", Var)
DrawPlots()
DefineScalarExpression("xc", "coord(%s)[0]" % Mesh)
DefineScalarExpression("yc", "coord(%s)[1]" % Mesh)
DefineScalarExpression("zc", "coord(%s)[2]" % Mesh)


# Do a lineout on all 4 variables to produce 4 curves.
Lineout(p0, p1, ("default", "xc", "yc", "zc"))

nSteps = tf - t0 + 1
for step in range(nSteps):

   time = step + t0
   SetTimeSliderState(time)
   get_data(Var,time)




# Again, needed when issueing -nowin option
sys.exit()
    
