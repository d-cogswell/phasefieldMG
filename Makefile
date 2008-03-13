TARGET = mse3D

CXX = icpc
CXXFLAGS = -O3 -xT -parallel -pthread
OBJECTS = main.o phasefield3D.o libgrid3D.a
LIBS = `Magick++-config --libs` -L. -lgrid3D

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $?

lib%.a: %.o
	ar rc $@ $?
	ranlib $@

clean:
	rm $(TARGET) *.o *.a output/*
