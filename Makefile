TARGET = phasefieldMG

CXX = g++ 
CXXFLAGS = -fopenmp -O3 -march=native
OBJECTS = main.cpp phasefield3DMG.cpp systm.cpp grid3D.cpp
LIBS = `Magick++-config --cxxflags --libs`

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $?

lib%.a: %.o
	ar rc $@ $?
	ranlib $@

clean:
	rm -f $(TARGET) *.o *.a
