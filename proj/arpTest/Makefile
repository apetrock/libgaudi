include ../include.mk

# calls:
CC         = g++ -DDEBUG -std=c++0x -g -ggdb -rdynamic
#CC         = g++
#CC         = clang++ -DDEBUG -g -ggdb -fsanitize=address -fno-omit-frame-pointer
CFLAGS     = ${CFLAGS_COMMON} ${ARPACK_INC} -O0
LDFLAGS    = ${LDFLAGS_COMMON} 
EXECUTABLE = ../../bin/arpTest

SOURCES    = arpTest.cpp \
	 	../../src/TIMER.cpp \
		../../src/manifold/manifold_singleton.cpp \
		../../src/manifold/debugger.cpp \

OBJECTS    = $(NANO_SRC:.cpp=.o) $(SOURCES:.cpp=.o)


#all: $(NANO_SRC) $(SOURCES) $(EXECUTABLE)
all:  $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
