## Britton's script to make Enve31 stitch together image files into
## one single .mp4 file

import dircache

#############################
## Get the files in the    ##
## current directory       ## 
#############################
def get_frames( prefix , dir ):
    files = dircache.listdir(dir);
    npre = len(prefix);
    count = 0;
    images = list();
    for file in files:
        test = file[0:npre];
        if test == prefix:
            images.append(dir + '/' + file);
            
    return  images
#############################
## END get_frames function ##
#############################

###############################
###### Movie Parameters #######
###############################

FPS = 20;            ## Frame rate
mov_name = "3d_fly_med";  ## Name of movie- appended with file type
dir = '.';          ## Directory location of files
prefix = 'grad0'; ## Prefix of the image frames
nx = 2048;           ## X dimensions of the movie
ny = 1024;           ## Y dimensions of the movie

        
# Get the image file names and count them
images = get_frames(prefix,dir);
nframes = len(images);

# Create the movie file and set the selected
# attributes
x = enve.movie(enve.MOVIE_WRITE);
x.fps = FPS;
x.filename = mov_name;
x.format = "MPEG";
x.options = "COMPRESSION MPEG4";
x.dims = nx,ny;
err = x.addcount(nframes);

# Open the movie file and write each frame
# output the progress and close the movie
err1 = x.open();
img = enve.image();
count = 0;
for image in images:
    file = image;
    err2 = img.load(file);
    err3 = x.append(img)
    count = count + 1;
    print str(100.0*float(count)/float(nframes))+'%' + ' file ' + file + str(err2)
 
err = x.close();
