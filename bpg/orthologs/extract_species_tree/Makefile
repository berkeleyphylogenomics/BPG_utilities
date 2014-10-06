CXX = g++
CXXLD = g++
CXXFLAGS = -O2
CXXCOMPILE = $(CXX) $(CXXFLAGS)
CXXLINK = $(CXXLD) $(CXXFLAGS) -o $@
LIBS = -lmysqlpp
INCLUDE_PATH = -I/usr/include/mysql -I/usr/include/mysql++


find_species_SOURCES = find_species.cpp

find_species_OBJECTS = find_species.o

all: find_species

find_species.o: find_species.cpp

find_species: $(find_species_OBJECTS)
	@rm -f find_species
	$(CXXLINK) $(find_species_OBJECTS) $(INCLUDE_PATH) $(LIBS)
.cpp.o:
	$(CXXCOMPILE) $(INCLUDE_PATH) -c -o $@ $<
