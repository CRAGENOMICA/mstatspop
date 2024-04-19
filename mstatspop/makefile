
all: mstatspop

mstatspop: 
    gcc *.c -lm -o ../bin/mstatspop -Wall -O3 -I/usr/local/include/gsl /usr/local/lib/libgsl.a -DinGSL=1 -lgsl -lgslcblas -lz

clean:
	@rm -f *.o ../bin/mstatspop

