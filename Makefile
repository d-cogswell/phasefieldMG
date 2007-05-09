TARGET = mse3D 
SOURCE = main.cpp phasefield3D.cpp
SOURCEMPI = main.cpp convexificationMPI.cpp
#CC = icc
CC = g++ 
CCMPI = mpic++
CLIBS = -L. -lgrid3D -lpthread -lX11
CLIBSMPI = -L. -lgrid3DMPI
#CFLAGS = -static -O3 -ipo -xN
CFLAGS = -O3
CFLAGSMPI = -O3 #-Wno-long-double

$(TARGET): $(SOURCE)
	#make grid3D
	$(CC) $(CFLAGS) $(SOURCE) $(CLIBS) -o $(TARGET)

debug: $(SOURCE)
	make grid3D
	$(CC) $(CFLAGS) -g $(SOURCE) $(CLIBS) -o $(TARGET)

MPI: $(SOURCEMPI) 
	make grid3D
	make grid3DMPI
	$(CCMPI) $(CFLAGSMPI) $(SOURCEMPI) $(CLIBSMPI) -o $(TARGET)

grid3D:
	$(CC) -c $(CFLAGS) grid3D.cpp
	ar rc libgrid3D.a grid3D.o
	ranlib libgrid3D.a

grid3DMPI:
	$(CCMPI) -c $(CFLAGS) grid3DMPI.cpp
	ar rc libgrid3DMPI.a grid3DMPI.o
	ranlib libgrid3DMPI.a

merge:
	make grid3D
	$(CC) merge.cpp -o merge -L. -lgrid3D

dxconvert:
	make grid3D
	$(CC) dxconvert.cpp -o dxconvert -L. -lgrid3D

clean:
	rm $(TARGET)
	rm output/*
