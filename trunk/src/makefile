CC=	ifort
#CC=	gfortran
 
SOURCE=	post_miranda.f90
EXE=	mir-post

$(EXE) :$(SOURCE) 
	$(CC)  $(SOURCE) $(CFLAGS) -o $(EXE) 

main.o: post_miranda.f90
	$(CC)  $(CFLAGS) -c post_miranda.f90

clean : 
	rm -f *.log $(EXE) *~ *.o *.mod