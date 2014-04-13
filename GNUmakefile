VPATH= . ./src ./lib/timer

CPP= g++
CPPFLAGS:=-g -O2 -Wall -Wno-sign-compare -std=c++11 -DDEBUG_ON
LIBFLAG:= -lconfig++ -lboost_filesystem -lboost_system
INC:= $(patsubst %, -I%, $(VPATH))
CPPFLAGS += $(INC)

# Only builds objects into SRC for files with headers that can be included (aka not for main() code)
HEADER:= $(wildcard ./src/*.h)
SRC:= $(patsubst %.h, %.cpp, $(HEADER)) $(wildcard ./lib/timer/*.cpp)
OBJ:= $(patsubst %.cpp, %.o, $(SRC))

.PHONY: all clean

all: src/traj2dEcsv src/dEcsv2corr

src/traj2dEcsv: FMOparse.cpp $(OBJ) GNUmakefile
	$(CPP) $(CPPFLAGS) -L/usr/lib -L/usr/lib/x86_64-linux-gnu/ -o $@ $< $(OBJ) $(LIBFLAG)

src/dEcsv2corr: MergeDE.cpp $(OBJ) GNUmakefile
	$(CPP) $(CPPFLAGS) -L/usr/lib -L/usr/lib/x86_64-linux-gnu/ -o $@ $< $(OBJ) $(LIBFLAG)

src/unittest: _test_grom2atomFMO.cpp $(OBJ) GNUmakefile
	$(CPP) $(CPPFLAGS) -o $@ $< $(OBJ) $(LIBFLAG)


clean: .
	rm *.o src/*.o
