all: sphereToRaw voxelMesh.so imageFileConvert   rawToFoam  clean



voxelMesh.so:
	g++ -c -fPIC imageReader.cpp -o imageReader.o
	g++ -c -fPIC voxelMesh.cpp -o voxelMesh.o
	g++ -shared voxelMesh.o imageReader.o -o $(I_DIR)/bin/libvoxelMesh.so


sphereToRaw.o:	sphereToRaw.cpp
	g++ -std=c++0x -c  sphereToRaw.cpp
    
sphereToRaw:	sphereToRaw.o 
	g++	-o	$(I_DIR)/bin/sphereToRaw	sphereToRaw.o	


imageFileConvert.o:	imageFileConvert.cpp
	g++ -c imageFileConvert.cpp
imageFileConvert:	imageFileConvert.o
	g++	  -fPIC  -L$(I_DIR)/bin/  imageFileConvert.o -lvoxelMesh -o	$(I_DIR)/bin/imageFileConvert 



rawToFoam.o:	rawToFoam.cpp
	g++   -fPIC -c rawToFoam.cpp
rawToFoam:	rawToFoam.o
	g++	  -fPIC  -o	$(I_DIR)/bin/rawToFoam -L$(I_DIR)/bin rawToFoam.o -lvoxelMesh



clean:
	rm  -f *.o
