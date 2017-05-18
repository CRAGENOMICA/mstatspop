# Command lines must start with a '\t' otherwise make shows: "make: *** falta un separador: Alto"
# If Eclipse debug shows: "No source available for main()", then set this to the gcc command: -g -O0
# If make shows "make: Error 1" then be sure that you set the gcc arguments correctly

all: clean mstatspop

mstatspop: 
	$(info ************  mstatspop ************)
	gcc *.c -lm -g -O0 -o ../bin/mstatspop -I/usr/local/include/gsl /usr/local/lib/libgsl.a -DinGSL=1 -lgsl -lgslcblas -lz

clean:
	$(info ************  clean ************)
	@rm -f *.o ../bin/mstatspop
