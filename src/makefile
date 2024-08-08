WARN = -Wall -Wsign-compare -Wextra
# use `valgrind` for memory / array bounds check
# use `-ggdb3` for use with `valgrind`
# use `-pg` switch for use with `gprof` (profiling/performance analyzer)
## EXTRA may be set on the command line, e.g., "EXTRA=''" or "EXTRA='-pg -DNDEBUG'"
# EXTRA = -DNDEBUG -Wno-unused-variable -Wno-unused-parameter  
# OPT  = -g -march=native -mtune=native -O3 ${PRODUCTION}
EXTRA=-DNDEBUG
# use -DNDEBUG to disable assert, -Wfatal-errors to stop at 1st error
# to compile on Mac you may need to disable '-march=native` 
OPT = -g -march=native -mtune=native -O3
#OPT = -g -mtune=native -O3 ${EXTRA}
CC = gcc
CFLAGS = ${OPT} ${VER} ${WARN} ${EXTRA}
CPP = g++

default: give_help
all: dist_m4ri 

dist_m4ri: dist_m4ri.c dist_m4ri.h util_io.o util_m4ri.o mmio.o dist_cc.o makefile 
	${CC} ${CFLAGS} -DSTANDALONE -o dist_m4ri $< dist_cc.o mmio.o util_m4ri.o util_io.o -lm4ri -lm

util_m4ri.o: util_m4ri.c util_m4ri.h mmio.h makefile 
	${CC} ${CFLAGS} -c -o util_m4ri.o $<  

util_io.o: util_io.c util_io.h util_m4ri.h mmio.h makefile 
	${CC} ${CFLAGS} -c -o util_io.o $<  

dist_cc.o: dist_cc.c util_io.h util_m4ri.h mmio.h makefile 
	${CC} ${CFLAGS} -c -o dist_cc.o $<  

mmio.o: mmio.c mmio.h makefile
	${CC} ${CFLAGS} -c $< 


clean:
	rm -f *~ *.tmp *.out *.o ../input/*~ ../examples/*~ TAGS 

#dox: 
#	doxygen ../Doxyfile mainpage.h

tags:
	etags *.c *.h 


veryclean: clean
	rm -f dist_m4ri *.exe 

zip: ${src}
	zip -9 -o dist_m4ri.zip ${src} ../Doxyfile ../README.md

give_help:
	@echo  "\e[36m make -j dist_m4ri \e[0m" "\t\t compile \`dist_m4ri\` \e[32m\e[0m"
	@echo  "\e[36m make -j dist_m4ri EXTRA=\"\"\e[0m" "\t compile with add'l checks / assert statements\e[32m\e[0m"
	@echo  "\e[36m make -j dist_m4ri OPT=\"-pg -O2\"\e[0m" "\t compile for use with \`gprof\`\e[32m\e[0m"
	@echo  "\e[36m make clean \e[0m   " "\t to\e[31m remove\e[0m object files and such"
	@echo  "\e[36m make veryclean \e[0m" "\t also\e[31m remove\e[0m executables and zip file"
	@echo  "\e[36m make zip \e[0m    " "\t to make \e[35m \`dist_m4ri.zip\`\e[0m file with all sources"
	@echo  "\e[36m make tags \e[0m   " "\t to create \e[35m TAGS\e[0m file"
	@echo
	@echo  "\e[31m The program does not compile?\e[0m"
	@echo  "  - Make sure you have '\e[32mlibm4ri_dev\e[0m' on your system."
	@echo  "    Please see '../README.md' for installation instructions."
	@echo  "  - On a Mac, you may want to disable '\e[32m-march=native\e[0m' compiler switch by running" 
	@echo  "\t\e[36m make dist_m4ri OPT='-g -mtune=native -O3'\e[0m" 

test:
	echo "start test"
	./dist_m4ri method=1
	./dist_m4ri method=2
