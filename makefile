cc=gcc -std=gnu99 -mmmx -msse -msse2 -msse3 -g -O3 -Wall -Ibitbucket/m4ri 
#mmio=../../decoder/mmio
mmio=mmio

## m4ri/ple.lo -MD -MP -MF m4ri/.deps/ple.Tpo -c m4ri/ple.c  -DDLL_EXPORT -DPIC -o

all: qdist_m4ri qdec_m4ri

test: test_elimination.c
	${cc} -o test_elimination test_elimination.c -lm4ri -lm 

qdec_m4ri: qdec_m4ri.c util_io.o util_m4ri.o ${mmio}.o  makefile 
	${cc} -DDEBUG -o qdec_m4ri $< ${mmio}.o util_m4ri.o util_io.o -lm4ri -lm

qdist_rec: qdist_m4ri_rec.c util_m4ri.h makefile ${mmio}.o util_m4ri.o 
	${cc} -DDEBUG -o qdist_rec $< ${mmio}.o util_m4ri.o -lm4ri -lm


qdist_m4ri: qdist_m4ri.c util_io.o util_m4ri.o ${mmio}.o  makefile 
	${cc} -DDEBUG -o qdist_m4ri $< ${mmio}.o util_m4ri.o util_io.o -l:libm4ri.a -lm

qdist_m4ri_old: qdist_m4ri_old.c util_io.o util_m4ri.o ${mmio}.o  makefile 
	${cc} -DDEBUG -o qdist_m4ri_old $< ${mmio}.o util_m4ri.o util_io.o -lm4ri -lm

#${mmio}.o: ${mmio}.c
#	cd ../../decoder 
#	make mmio.o

util_m4ri.o: util_m4ri.c util_m4ri.h makefile 
	${cc} -c -o util_m4ri.o $<  

util_io.o: util_io.c util_io.h util_m4ri.h makefile 
	${cc} -c -o util_io.o $<  

mmio.o: mmio.c mmio.h makefile
	${CC} -c $< 


world.o:hello.cpp
	g++ hello.cpp -o world.o

compiler=g++
hello.o:hello.cpp
	$(compiler) $< -o $@