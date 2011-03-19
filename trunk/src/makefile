f90=	ifort
#f90=	gfortran

flags = 

SOURCE=	globals.o interfaces.o lines.o driver.o block_ops.o 
PROBS= yplus length span

warn :
	@echo "Please select a case to compile:"
	@echo "  - $(PROBS)"

all : $(PROBS)

%.o: %.f90
	$(f90) $(flags) -c $<

clean : 
	rm -f *.log *.o *.mod

cleanall : 
	rm -f *.log *.o *.mod *_ex

post :$(SOURCE) post_miranda.o
	$(f90) $(SOURCE) post_miranda.o $(flags) -o post_ex

yplus :$(SOURCE) yplus.o
	$(f90) $(SOURCE) yplus.o $(flags) -o yplus_ex

length :$(SOURCE) length.o
	$(f90) $(SOURCE) length.o $(flags) -o length_ex

span :$(SOURCE) span.o
	$(f90) $(SOURCE) span.o $(flags) -o span_ex
