CC = g++
CPFLAGS = -O3 -std=c++11
CFLAGS = -O3 

ifneq "$(NETCDF)" ""
        INCLUDES += -I$(NETCDF)/include
        LIBS += -L$(NETCDF)/lib
        NCLIB = -lnetcdf
        NCLIBF = -lnetcdf_c++4
        ifneq ($(wildcard $(NETCDF)/lib/libnetcdf_c++.*), ) # CHECK FOR NETCDF4
                LIBS += $(NCLIBF)
        endif # CHECK FOR NETCDF4
        LIBS += $(NCLIB)
endif

all: meshmaker

meshmaker: triangle.o triangulator.cpp
	$(CC) $(CPFLAGS) triangle.o triangulator.cpp -o meshmaker

triangle.o: triangle.h
	$(CC) $(CFLAGS) -c triangle.c

#triangulator.o: triangulator.cpp triangle.o
#	$(CC) $(CPFLAGS) -c triangle.o triangulator.cpp

density: addMeshDensity.cpp
	$(CC) $(CPFLAGS) addMeshDensity.cpp -o addDensity $(INCLUDES) $(LIBS) 

clean: 
	rm *.o meshmaker

icpc:
	( $(MAKE) all \
    "CC = icpc" )
