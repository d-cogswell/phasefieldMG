TARGET = phasefieldMG

CXX = g++ 
CXXFLAGS = -fopenmp -O3 -march=native
OBJECTS = main.o phasefield3DMG.o systm.o grid3D.o
LIBS = `Magick++-config --cxxflags --libs` -lnetcdf

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -o $(TARGET)

main.o: main.cpp multigrid3D.h
	$(CXX) $(CXXFLAGS) $(LIBS) -c $<

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $<

lib%.a: %.o
	ar rc $@ $?
	ranlib $@

clean:
	rm -f $(TARGET) *.o *.a
