## Script to take the experimental images, rotate them, make them B&W
## and put them in an arbitrary order
import os

Nx = 1920
Ny = 360
base='medium3dgradRHO_';
fmt='.tiff'
Iname = []
#Iname.append('1100')
#Iname.append('1125')
Iname.append('0600')
Iname.append('0625')
Iname.append('0650')
#Iname.append('1225')

for i in range(len(Iname)):
    Iname[i] = base + Iname[i] + fmt;
    

## Make tmp workspace
cmd = 'mkdir tmp'
os.system(cmd)

## Copy all the images
path = '/p/lscratchd/olson45/NOZ_VIZ/2d_med/'
for i in range(len(Iname)):
    cmd = "cp " + path + Iname[i] + ' tmp/' + Iname[i]
    os.system(cmd)

## Crop to remove black
Nxx = int(Nx/1.6) - 124;
xoff = 124
Nyy = int(Ny/1.25);
yoff = int((Ny-Nyy)/2);
crop_args = " -crop " + str(Nxx)+"x" + str(Nyy)+ "+"+str(xoff)+"+" + str(yoff) + " ";
for i in range(len(Iname)):
    cmd = 'cd tmp; convert ' + Iname[i] + crop_args + Iname[i]
    os.system(cmd)


## Flip all images to get large lambda on top
for i in range(len(Iname)):
    cmd = "cd tmp; mogrify -flip " + Iname[i]
    os.system(cmd);

## Make a montage of the images
ii = [0]*len(Iname)
ii[0] = 0;
ii[1] = 1;
ii[2] = 2;
#ii[3] = 3;
#ii[4] = 4;
#ii[5] = 5;

mN = 3;                    # Number of images in montage
bdr = 5;                   # Border size
mNx = Nxx                  # Nx
mNy = Nyy*mN + bdr*(mN+1)  # Ny
ofile = 'test.jpg'

files = ' ';
for i in range(len(Iname)):
    j = ii[i]
    files = files + ' ' + Iname[j] + " "

cmd = "cd tmp; montage -border " + str(bdr) + " -bordercolor white " +  " -geometry "+ str(mNx)+"x"+ " -tile 1x"+str(mN)+files+ofile
os.system(cmd)


## Make a grayscale image
cmd = "cd tmp; mogrify -colorspace Gray " + ofile
os.system(cmd)

## Clean up stuff
cmd = "cp tmp/"+ofile+" sim_montage.jpg"
os.system(cmd)

cmd = 'rm -rf tmp'
os.system(cmd)

#cmd = "mogrify -resize 1000x sim_montage.jpg"
#os.system(cmd)

print "Done with montage:"
print str(mN) + " images written"
print str(mNx)+"x"+str(Nyy)+" (*"+str(mN)+")"
