/*-------------------------------------------------------------------------*\	
You can redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.


The code has been developed by Ali Qaseminejad Raeini as a part his PhD 
at Imperial College London, under the supervision of Branko Bijeljic 
and Martin Blunt. 
* 
Please see our website for relavant literature:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:    a.qaseminejad-raeini09@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
\*-------------------------------------------------------------------------*/

    #include <fstream>
    #include <iostream>
    #include <vector>

    #include <assert.h>

#include "voxelMesh.h"




int main(int argc, char *argv[])
{

    if(argc<=3)
    {
        std::cout<<"   imageFileConvert utility to convert files from binary (if suffix is .raw) to ascii format\n"
                 <<"   and vice versa with optional trimming and re-sampling \n"
                 <<"usage: \n"
                 <<"\n    imageFileConvert headerFile inputFile outPutFile \n\n"
                 <<"header file contents should be like: \n"
                 <<"  Nxyz  \n"
                 <<"  dxXo \n"
                 <<"  500  300  400 \n"
                 <<"  5.2  5.2  5.2 \n"
                 <<"  0.0  0.0  0.0 \n"
                 <<"  crop  0 300  0 300  0 300 \n"
                 <<"  resample    2 \n"
                 <<"  pore     0 0  \n"
                 <<"only the numbers can be changed in this file. crop, pore and resample keywords are optional\n"               
                 << std::endl;
        return 1;
    }







	std::string headerName(argv[1]);
	std::string inputName(argv[2]);
	std::string outputName(argv[3]);
	//~ std::string trimName;


    //~ unsigned int n[3];
    //~ double  xmin[3];
    //~ double dx[3];    
    voxelMesh voxelImage;

	readFromHeader(voxelImage, headerName, inputName);
    unsigned int n[3];	voxelImage.getSize(n[0],n[1],n[2]);
    Double3D xmin=voxelImage.X0();
    Double3D dx=voxelImage.dx();
    //~ voxelImage.crop(0,n[0]-1,0,n[1]-1,0,n[2]-1,2,1);//         XXXXXXXXXXXXXXXXXXXXXXXXXXXX



	voxelImage.write(outputName);



    //~ voxelImage.crop(0,0.5*(n[0]-1), 0,0.5*(n[1]-1), 0,0.5*(n[2]-1), 0,1);
    //~ voxelImage.writeAConnectedPoreVoxel(outputName+"_aPoreVoxel");
    
    std::cout<< "end" << std::endl;
    
	voxelImage.printInfo();


    return 0;
}


