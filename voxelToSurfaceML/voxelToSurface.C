/*-------------------------------------------------------------------------*\
This code is part of poreFOAM, a suite of codes written using OpenFOAM
for direct simulation of flow at the pore scale. 	
You can redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.



The code has been developed by Ali Qaseminejad Raeini as a part his PhD 
at Imperial College London, under the supervision of Branko Bijeljic 
and Martin Blunt. 
Please see our website for relavant literature:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:    a.qaseminejad-raeini09@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
 
 Description:
	creates a surface between the pore and the solid from a 3D rock image
\*-------------------------------------------------------------------------*/

    #include <fstream>
    #include <iostream>
    #include <vector>

    #include <assert.h>

#include "fvCFD.H"

#include "argList.H"
#include "timeSelector.H"
#include "graph.H"
#include "mathematicalConstants.H"
#include "labelVector.H"

#include "OFstream.H"
#include "triFaceList.H"
#include "triSurface.H"

#include "DynamicField.H"
#include "MeshedSurfaces.H"

#include "voxelMesh.h"
#include "createSurface.h"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:



int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

    IOdictionary meshingDict
    (
        IOobject
        (
            "meshingDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    word headerName(meshingDict.lookup("headerName"));
    word inputName(meshingDict.lookup("inputName"));
    //word outputName(meshingDict.lookup("outputName"));
    word outputSurface(meshingDict.lookup("outputSurface"));





	voxelMesh vximage;

	readFromHeader(vximage, headerName, inputName);

	vximage.growBox(0);


	
    unsigned int n[3];
	vximage.getSize(n[0],n[1],n[2]);
    //vximage.crop(0,n[0]-1,0,n[1]-1,0,n[2]-1, 2 ,0);//         XXXXXXXXXXXXXXXXXXXXXXXXXXXX
    vximage.crop(0,n[0]-1,0,n[1]-1,0,n[2]-1, 2 ,1);//         XXXXXXXXXXXXXXXXXXXXXXXXXXXX

	 vximage.median(6);		
	 vximage.median(4);	
	 vximage.median(3);
	 vximage.median(2);	

	 vximage.median(4);
	 vximage.median(6);	
	 vximage.median(3);
	 vximage.median(2);
	 vximage.median(3);
	 vximage.median(4);
	 vximage.median(2);
	 vximage.median(3);
	 vximage.median(4);
	 vximage.median(2);
	 vximage.median(3);
	 vximage.median(4);
	 vximage.median(3);
	 vximage.median(2);
	 vximage.median(4);
	 vximage.median(3);
	 vximage.median(2);
	 vximage.median(4);
	 vximage.median(6);
	
			
    vximage.printInfo();
 
    writeSTLBINARY(vximage, outputSurface);//         XXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    Info<<"finished  writeSTLBINARY, "<<endl;

    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
