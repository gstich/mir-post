f90=	ifort
#f90=	gfortran

flags = -O3

SOURCE=	globals.o interfaces.o lines.o driver.o block_ops.o post_miranda.o
EXE=	mir-post

$(EXE) :$(SOURCE) 
	$(f90)  $(SOURCE) $(flags) -o $(EXE) 

%.o: %.f90
	$(f90) $(flags) -c $<

clean : 
	rm -f *.log $(EXE) *.o *.mod