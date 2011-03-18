f90=	ifort
#f90=	gfortran

flags = 

SOURCE=	globals.o interfaces.o lines.o driver.o block_ops.o #post_miranda.o
EXE=	mir-post

### $(EXE) :$(SOURCE) 
###	$(f90)  $(SOURCE) $(flags) -o $(EXE) 

warn :
	@echo "Please select a case to compile:"

%.o: %.f90
	$(f90) $(flags) -c $<

clean : 
	rm -f *.log $(EXE) *.o *.mod

post :$(SOURCE) post_miranda.o
	$(f90) $(SOURCE) post_miranda.o $(flags) -o post

yplus :$(SOURCE) yplus.o
	$(f90) $(SOURCE) yplus.o $(flags) -o yplus

length :$(SOURCE) length.o
	$(f90) $(SOURCE) length.o $(flags) -o length

span :$(SOURCE) span.o
	$(f90) $(SOURCE) span.o $(flags) -o span
