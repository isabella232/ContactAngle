all: writeRaCLInds writeRaPoreInds writepAreaPoreInds writeCAThroatInds writeRcThroatInds writeCAPoreInds writeRcPoreInds writeCACornerInds writeRcCornerInds    
# warning type clean, dependencies are not handelled yet

myBinDir ?= ../../../bin


rawToFoam.o:	rawToFoam.cpp
	g++   -O3 -Wall   -std=c++0x  -pedantic -c rawToFoam.cpp
rawToFoam:	rawToFoam.o
	g++	   -std=c++0x  -o	$(myBinDir)/rawToFoam rawToFoam.o 


rawToFoamPar.o:	rawToFoamPar.cpp
	g++    -O3 -Wall   -std=c++0x -pedantic -c rawToFoamPar.cpp
rawToFoamPar:	rawToFoamPar.o
	g++	 -O3 -Wall  -std=c++0x  -o	$(myBinDir)/rawToFoamPar rawToFoamPar.o 



voxelImageConvert.o:	voxelImageConvert.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c voxelImageConvert.cpp
voxelImageConvert:	voxelImageConvert.o
	g++ -O3 -Wall    voxelImageConvert.o  -o $(myBinDir)/voxelImageConvert 
	ln -f -s $(myBinDir)/voxelImageConvert $(myBinDir)/imageFileConvert


Ufraw2Ucdat.o:	Ufraw2Ucdat.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c Ufraw2Ucdat.cpp
Ufraw2Ucdat:	Ufraw2Ucdat.o
	g++ -O3 -Wall    Ufraw2Ucdat.o  -o $(myBinDir)/Ufraw2Ucdat 



clean:
	rm  -f *.o


writepAreaPoreInds.o:	writepAreaPoreInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writepAreaPoreInds.cpp
writepAreaPoreInds:	writepAreaPoreInds.o
	g++ -O3 -Wall    writepAreaPoreInds.o  -o $(myBinDir)/writepAreaPoreInds 
	ln -f -s $(myBinDir)/writepAreaPoreInds $(myBinDir)/imageFileConvert
	
writeCAPoreInds.o:	writeCAPoreInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeCAPoreInds.cpp
writeCAPoreInds:	writeCAPoreInds.o
	g++ -O3 -Wall    writeCAPoreInds.o  -o $(myBinDir)/writeCAPoreInds 
	ln -f -s $(myBinDir)/writeCAPoreInds $(myBinDir)/imageFileConvert
		
writeRcPoreInds.o:	writeRcPoreInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRcPoreInds.cpp
writeRcPoreInds:	writeRcPoreInds.o
	g++ -O3 -Wall    writeRcPoreInds.o  -o $(myBinDir)/writeRcPoreInds 
	ln -f -s $(myBinDir)/writeRcPoreInds $(myBinDir)/imageFileConvert

writeRaPoreInds.o:	writeRaPoreInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRaPoreInds.cpp
writeRaPoreInds:	writeRaPoreInds.o
	g++ -O3 -Wall    writeRaPoreInds.o  -o $(myBinDir)/writeRaPoreInds 
	ln -f -s $(myBinDir)/writeRaPoreInds $(myBinDir)/imageFileConvert


writeRaCLInds.o:	writeRaCLInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRaCLInds.cpp
writeRaCLInds:	writeRaCLInds.o
	g++ -O3 -Wall    writeRaCLInds.o  -o $(myBinDir)/writeRaCLInds 
	ln -f -s $(myBinDir)/writeRaCLInds $(myBinDir)/imageFileConvert

	
writeCACornerInds.o:	writeCACornerInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeCACornerInds.cpp
writeCACornerInds:	writeCACornerInds.o
	g++ -O3 -Wall    writeCACornerInds.o  -o $(myBinDir)/writeCACornerInds 
	ln -f -s $(myBinDir)/writeCACornerInds $(myBinDir)/imageFileConvert

writeRcCornerInds.o:	writeRcCornerInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRcCornerInds.cpp
writeRcCornerInds:	writeRcCornerInds.o
	g++ -O3 -Wall    writeRcCornerInds.o  -o $(myBinDir)/writeRcCornerInds 
	ln -f -s $(myBinDir)/writeRcCornerInds $(myBinDir)/imageFileConvert



writeCAThroatInds.o:	writeCAThroatInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeCAThroatInds.cpp
writeCAThroatInds:	writeCAThroatInds.o
	g++ -O3 -Wall    writeCAThroatInds.o  -o $(myBinDir)/writeCAThroatInds 
	ln -f -s $(myBinDir)/writeCAThroatInds $(myBinDir)/imageFileConvert
	
writeRcThroatInds.o:	writeRcThroatInds.cpp
	g++ -O3 -Wall -std=c++0x -pedantic -c writeRcThroatInds.cpp
writeRcThroatInds:	writeRcThroatInds.o
	g++ -O3 -Wall    writeRcThroatInds.o  -o $(myBinDir)/writeRcThroatInds 
	ln -f -s $(myBinDir)/writeRcThroatInds $(myBinDir)/imageFileConvert


