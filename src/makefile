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

.PHONY: default all test clean veryclean dox tags zip help give_help

default: all
all: dist_m4ri dist_m4ri_old

test: all
	./tests/test_run.sh

dist_m4ri: dist_m4ri.c dist_m4ri.h dist_cc.h rw_triplet.h util_io.o util_m4ri.o mmio.o dist_cc.o rw_triplet.o makefile
	${CC} ${CFLAGS} -pthread -o dist_m4ri $< dist_cc.o mmio.o util_m4ri.o util_io.o rw_triplet.o -lm4ri -lm -lpthread

dist_m4ri_old: dist_rw.c dist_m4ri.h dist_cc.h util_io.o util_m4ri.o mmio.o dist_cc.o makefile 
	${CC} ${CFLAGS} -DSTANDALONE -o dist_m4ri_old $< dist_cc.o mmio.o util_m4ri.o util_io.o -lm4ri -lm

distfork: dist_m4ri
	ln -sf dist_m4ri distfork

util_m4ri.o: util_m4ri.c util_m4ri.h mmio.h makefile 
	${CC} ${CFLAGS} -c -o util_m4ri.o $<  

util_io.o: util_io.c util_io.h util_m4ri.h mmio.h util_hash.h makefile 
	${CC} ${CFLAGS} -c -o util_io.o $<  

rw_triplet.o: rw_triplet.c rw_triplet.h util_io.h util_m4ri.h mmio.h util_hash.h makefile
	${CC} ${CFLAGS} -c -o rw_triplet.o $<

dist_cc.o: dist_cc.c dist_cc.h util_io.h util_m4ri.h mmio.h util_hash.h makefile 
	${CC} ${CFLAGS} -c -o dist_cc.o $<  

mmio.o: mmio.c mmio.h makefile
	${CC} ${CFLAGS} -c $< 


clean:
	rm -f *~ *.tmp *.out *.o ../input/*~ ../examples/*~ TAGS 

dox: 
	cd .. && doxygen Doxyfile 

tags:
	etags *.c *.h 


veryclean: clean
	rm -f dist_m4ri dist_m4ri_old distfork *.exe dist_m4ri.zip

src = dist_m4ri.c dist_rw.c dist_cc.c util_io.c util_m4ri.c mmio.c rw_triplet.c dist_m4ri.h dist_cc.h util_io.h \
      util_m4ri.h rw_triplet.h util_hash.h mmio.h sort.h uthash.h makefile
zip: ${src}
	zip -9 -o dist_m4ri.zip ${src} ../Doxyfile ../README.md

give_help: help

help:
	@printf '%b\n' '\033[36m make all \033[0m\t\t compile `dist_m4ri` and `dist_m4ri_old`'
	@printf '%b\n' '\033[36m make dist_m4ri \033[0m\t compile multithreaded `dist_m4ri`'
	@printf '%b\n' '\033[36m make dist_m4ri_old \033[0m\t compile single-threaded `dist_m4ri_old`'
	@printf '%b\n' '\033[36m make distfork \033[0m\t\t create \033[35m`distfork`\033[0m symlink to `dist_m4ri`'
	@printf '%b\n' "\033[36m make -j dist_m4ri EXTRA=\"\"\033[0m\t compile with add'l checks / assert statements"
	@printf '%b\n' '\033[36m make -j dist_m4ri OPT="-pg -O2"\033[0m compile for use with `gprof`'
	@printf '%b\n' '\033[36m make test \033[0m\t\t build all, then run `tests/test_run.sh` (about 2 minutes)'
	@printf '%b\n' '\033[36m make dox \033[0m\t\t build \033[35mDoxygen\033[0m documentation (`../Doxyfile`)'
	@printf '%b\n' '\033[36m make clean \033[0m\t\t to \033[31mremove\033[0m object files and such'
	@printf '%b\n' '\033[36m make veryclean \033[0m\t also \033[31mremove\033[0m executables, symlink, and zip file'
	@printf '%b\n' '\033[36m make zip \033[0m\t\t to make \033[35m`dist_m4ri.zip`\033[0m file with all sources'
	@printf '%b\n' '\033[36m make tags \033[0m\t\t to create \033[35mTAGS\033[0m file'
	@printf '%b\n' ''
	@printf '%b\n' '\033[32m Python tests (from the repository root):\033[0m'
	@printf '%b\n' '  \033[36mpython3 -m pytest -q tests/test_dist_m4ri.py\033[0m'
	@printf '%b\n' '  - `test_auto_stage_b_real_qec_cert` (needs `\033[32mqec_cert\033[0m`) is skipped unless'
	@printf '%b\n' '    `qec_cert` and `\033[32mlibqeccert.so\033[0m` load: build a `qec-cert` checkout (`make`'
	@printf '%b\n' '    there), then run'
	@printf '%b\n' '  \033[36mPYTHONPATH=<path to qec-cert> python3 -m pytest -q -rs tests/test_dist_m4ri.py\033[0m'
	@printf '%b\n' '    (alternatives: `\033[36mpip install -e <path to qec-cert>\033[0m`, or'
	@printf '%b\n' '    `\033[36mQEC_CERT_LIB=<path to libqeccert.so>\033[0m` with an installed `qec_cert`).'
	@printf '%b\n' ''
	@printf '%b\n' '\033[31m The program does not compile?\033[0m'
	@printf '%b\n' "  - Make sure you have '\033[32mlibm4ri-dev\033[0m' on your system."
	@printf '%b\n' "    Please see '../README.md' for installation instructions."
	@printf '%b\n' "  - On a Mac, you may want to disable '\033[32m-march=native\033[0m' compiler switch by running"
	@printf '%b\n' "\t\033[36m make dist_m4ri OPT='-g -mtune=native -O3'\033[0m"


