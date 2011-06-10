## Script to take the experimental images, rotate them, make them B&W
## and put them in an arbitrary order
import os

Nx = 2048
Ny = 1536
base='NWM059';
fmt='.JPG'
Iname = []
Iname.append('15')
Iname.append('16')
Iname.append('18')
Iname.append('19')
Iname.append('20')
Iname.append('21')

for i in range(len(Iname)):
    Iname[i] = base + Iname[i] + fmt;
    

## Make tmp workspace
cmd = 'mkdir tmp'
os.system(cmd)

## Copy all the images
for i in range(len(Iname)):
    cmd = 'cp WALL/' + Iname[i] + ' tmp/' + Iname[i]
    os.system(cmd)

## Rotate all the images
#for i in range(len(Iname)):
#    cmd = 'cd tmp; mogrify -rotate 180 ' + Iname[i]
#    os.system(cmd)


## Crop to remove black
Nyy = int(Ny/2.8);
yoff = int(Ny/7.25);
crop_args = " -crop " + str(Nx)+"x" + str(Nyy)+ "+0+" + str(yoff) + " ";
for i in range(len(Iname)):
    cmd = 'cd tmp; convert ' + Iname[i] + crop_args + Iname[i]
    os.system(cmd)

## Resize for paper (to fit sim_data
Nxx = 1076;
Nyy = 288;
crop_args = " -resize " + str(Nxx)+"x" + str(Nyy)+" ";
for i in range(len(Iname)):
    cmd = 'cd tmp; mogrify ' + crop_args + Iname[i]
    os.system(cmd)

## Flip some images to get large lambda on top
cmd = "cd tmp; mogrify -flip " + Iname[1]
os.system(cmd);

## Make a montage of the images
ii = [0]*len(Iname)
j=0;
#ii[j] = 2;j=j+1
ii[j] = 5;j=j+1
ii[j] = 4;j=j+1
ii[j] = 3;j=j+1
#ii[j] = 0;j=j+1
#ii[j] = 1;

mN = 3;                    # Number of images in montage
bdr = 5;                   # Border size
mNx = Nxx                  # Nx
mNy = Nyy*mN + bdr*(mN+1)  # Ny
ofile = 'test.jpg'

files = ' ';
for i in range(mN):
    j = ii[i]
    files = files + ' ' + Iname[j] + " "

cmd = "cd tmp; montage -border " + str(bdr) + " -bordercolor white " +  " -geometry "+ str(mNx)+"x"+ " -tile 1x"+str(mN)+files+ofile
os.system(cmd)


## Make a grayscale image
cmd = "cd tmp; mogrify -colorspace Gray " + ofile
os.system(cmd)

## Vary the brightness
cmd = "cd tmp; mogrify -modulate 150 " + ofile
os.system(cmd)

## Clean up stuff
cmd = "cp tmp/"+ofile+" montage.jpg"
os.system(cmd)

cmd = 'rm -rf tmp'
os.system(cmd)

#cmd = "mogrify -resize 1000x montage.jpg"
#os.system(cmd)

print "Done with montage:"
print str(mN) + " images written"
print str(mNx)+"x"+str(Nyy)+" (*"+str(mN)+")"
