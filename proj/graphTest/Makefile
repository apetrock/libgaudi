include ../include.mk

# calls:
CC         = g++ -DDEBUG -std=c++0x -g -ggdb -rdynamic
#CC         = g++
#CC         = clang++ -DDEBUG -g -ggdb -fsanitize=address -fno-omit-frame-pointer
CFLAGS     = ${CFLAGS_COMMON} -O0
LDFLAGS    = ${LDFLAGS_COMMON}
EXECUTABLE = ../../bin/graphTest

SOURCES    = main.cpp \

OBJECTS    = $(NANO_SRC:.cpp=.o) $(SOURCES:.cpp=.o)


#all: $(NANO_SRC) $(SOURCES) $(EXECUTABLE)
all:  $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
