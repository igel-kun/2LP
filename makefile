include makefile_common

TARGET=main
PROG_NAME=cr
SUBDIRS=util reduction solv
LIB_CPPS=$(shell ls -f $(addsuffix /*.cpp,$(SUBDIRS)))
LIB_OS=$(LIB_CPPS:.cpp=.o)

all: $(SUBDIRS)
	g++ $(CFLAGS) -std=c++0x -Wall -static  ${LIB_OS} ${TARGET}.cpp -o ${PROG_NAME}  2>&1 | tee error.log || sleep 1


$(SUBDIRS):
	make -s -C $@

tests:
	cd tests && { ./testing ; cd .. ; }

clean:
	rm -f $(shell find -name "*.o") ${PROG_NAME}

.PHONY: $(SUBDIRS) tests clean
