TARGET = phasefieldMG

CXX = g++ 
CXXFLAGS = -DMAGICKCORE_HDRI_ENABLE=0 -DMAGICKCORE_QUANTUM_DEPTH=16 -O3 -ffast-math
OBJECTS = main.o multigrid3D.o phasefield3DMG.o libgrid3D.a
LIBS = `Magick++-config --libs` -L. -lgrid3D

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $?

lib%.a: %.o
	ar rc $@ $?
	ranlib $@

clean:
	rm -f $(TARGET) *.o *.a
