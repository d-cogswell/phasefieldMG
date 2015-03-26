TARGET = phasefieldMG

CXX = g++ 
CXXFLAGS = -DMAGICKCORE_HDRI_ENABLE=0 -DMAGICKCORE_QUANTUM_DEPTH=16 -O3 -ffast-math -march=native
OBJECTS = main.cpp multigrid3D.cpp phasefield3DMG.cpp grid3D.cpp
LIBS = `Magick++-config --libs`

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $?

lib%.a: %.o
	ar rc $@ $?
	ranlib $@

clean:
	rm -f $(TARGET) *.o *.a
