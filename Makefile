IDIR=include/
ODIR=obj/
SDIR=src/
CXX=g++
CFLAGS=-O2 -I $(IDIR) -fPIC

LIBS=-fopenmp

HEADS=aggregate.h column.h dendrite.h geom_lib.h math_lib.h hex_prism.h plate.h population.h pristine.h rosette.h
HEADERS=$(addprefix $(IDIR),$(HEADS)) # append IDIR to every HEAD

OBJS=aggregate.o column.o dendrite.o geom_lib.o math_lib.o hex_prism.o plate.o population.o pristine.o rosette.o
OBJECTS=$(addprefix $(ODIR),$(OBJS)) # append ODIR to every OBJ

# By default make calls the first rule
all: aggregation collection melting

$(ODIR):
	mkdir -p $(ODIR)
	mkdir -p dat
	mkdir -p vtk

$(ODIR)%.o:  $(SDIR)%.cpp | $(ODIR)
	$(CXX) $(CFLAGS) $(LIBS) -c $< -o $@

aggregation: $(OBJECTS)
	@echo "Recipe for aggregation"
	$(CXX) aggregation.cpp $(OBJECTS) $(CFLAGS) $(LIBS) -o $@

collection: $(OBJECTS)
	$(CXX) collection.cpp $(OBJECTS) $(CFLAGS) $(LIBS) -o $@

melting: $(meltingOBJ)
	@echo "No recipe for melting"

.PHONY: all clean cleanoutput cleanall aggregation collection melting

cleanall: clean cleanoutput

clean:
	rm -f $(ODIR)*.o

cleanoutput:
	rm -f *.dat *.vtk dat/*.dat vtk/*.vtk