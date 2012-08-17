TARGET = mse3D

CXX = g++ 
CXXFLAGS = -O3
OBJECTS = main.o multigrid3D.o electrochem.o electrodep3DMG.o libgrid3D.a
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
