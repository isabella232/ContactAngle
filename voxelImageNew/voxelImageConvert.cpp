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

#include "voxelImage.h"
//using namespace std;

int usage()
{
	std::cout<<"   voxelImageConvert utility to convert files from binary (if suffix is .raw) to ascii format\n"
		<<"   and vice versa with optional trimming and re-sampling \n"
		<<"usage: \n"
		<<"\n    voxelImageConvert headerFile inputFile outPutFile \n"
		<<"\n    or (works only with simple mhd files):"
		<<"\n    voxelImageConvert inputMetaImage.mhd outPutFile \n\n"
		<<"header file contents should be like: \n"
		<<"  NxyzdXXo \n"
		<<"  500  300  400 \n" 				<<"  5.2  5.2  5.2 \n"		<<"  0.0  0.0  0.0 \n\n"         
		<<"  crop  0 299  0 99  0 199 \n"		<<"  resample    2 \n"		<<"  pore     0 0 \n"		<<"  direction  z \n"
		<<"only the numbers can be changed in this file. \n"				<<"crop, pore, resample and direction keywords are optional\n"               
		<< std::endl;
	return -1;
}


int main(int argc, char *argv[])
{


	if(argc!=3) usage();


	std::string headerName(argv[1]);
	if(headerName.size()<4) return usage();
	std::string outputName(argv[2]);

	std::cout<<"voxelImageConvert "<<headerName<<"  "<<outputName<<std::endl;

	std::unique_ptr<voxelImageTBase> vxlImage = readImage(headerName);


	vxlImage->write(outputName);


	std::cout<< "end" << std::endl;

	vxlImage->printInfo();


	return 0;
}


