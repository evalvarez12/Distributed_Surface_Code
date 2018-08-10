DIRS := . MinCost GEOM

SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp)) 
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

#CFLAGS := -O3 -D_NDEBUG	# flags for .o compilation
CFLAGS:= -c -fPIC
CXX ?= icc		# cpp compiler
LIBS := -lrt
INCLUDES :=
LIBDIR :=

all: PMlib.so

PMlib.so: ${OBJS}
	$(CXX) -shared ${OBJS} -o $@  ${LIBS}

#blossom5: ${OBJS}
#	$(CXX) $(CFLAGS) ${LIBDIR} -o $@ ${OBJS} ${LIBS}

.cpp.o:
	$(CXX) $(CFLAGS) ${INCLUDES} $< -o $@

clean:
	rm -f ${OBJS} blossom5
