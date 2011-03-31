f90=	ifort
#f90=	gfortran

flags = 

SOURCE=	globals.o interfaces.o lines.o planes.o driver.o block_ops.o 
PROBS= yplus length span span_flc writeplanes

warn :
	@echo "Please select a case to compile:"
	@echo "  "
	@echo "Case are: $(PROBS)"
	@echo "  "
	@echo "Also : all, clean, cleanall"

all : $(PROBS)

%.o: %.f90
	$(f90) $(flags) -c $<

clean : 
	rm -f *.log *.o *.mod

cleanall : 
	rm -f *.log *.o *.mod *_ex

post :$(SOURCE) post_miranda.o
	$(f90) $(SOURCE) post_miranda.o $(flags) -o post_ex

yplus :$(SOURCE) cases/yplus.o
	$(f90) $(SOURCE) yplus.o $(flags) -o yplus_ex

length :$(SOURCE) cases/length.o
	$(f90) $(SOURCE) length.o $(flags) -o length_ex

span :$(SOURCE) cases/span.o
	$(f90) $(SOURCE) span.o $(flags) -o span_ex

span_flc :$(SOURCE) cases/span_flc.o
	$(f90) $(SOURCE) span_flc.o $(flags) -o span_flc_ex

writeplanes :$(SOURCE) cases/writeplanes.o
	$(f90) $(SOURCE) writeplanes.o $(flags) -o writeplanes_ex