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
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

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
int usage()
{
	std::cout<<"Ufraw2Ucdat "<<std::endl;
		std::cout
		<<" converts vxlImage.dat to vxlImage.raw, p.raw to p.dat, \n"
		<<" Ufx.raw to Uxcc.dat, Ufy.raw to Uycc.dat and, Ufz.raw to Uzcc.dat\n\n"
		<<"To continue, type:  \n"
		<<" Ufraw2Ucdat  vxlImage.mhd"<< std::endl;
	return 1;
}

int main(int argc, char** argv)
{

	if(argc!=2)		return usage();
	std::string headerName(argv[1]);
	if(headerName.size()<4 || headerName.compare(0,headerName.size(),"vxlImage.mhd") != 0) return usage();

	voxelImage vimage("vxlImage.mhd");
	if(!vimage.size()) {std::cout<<"Error: vxlImage.mhd not read"<<std::endl; return 1;}
	int3 n=vimage.size3();
	vimage.write("vxlImage.dat");
	vimage.resize(0);

	{
		voxelImageT<float> fField(n[0]+1,n[1],n[2],0.0);
		fField.readBin("Ufx.raw");

		for (int k = 0; k<int(fField.size()) ; k++ )
		 for ( int j = 0; j<int(fField[k].size()) ; j++ )
		  for ( int i = 0; i<int(fField[0][0].size())-1 ; i++ )
			fField[k][j][i]=0.5*(fField[k][j][i]+fField[k][j][i+1]);
		
		fField.writeAscii("Uccx.dat", 0,n[0],0,n[1],0,n[2]);
	}
	{
		voxelImageT<float> fField(n[0],n[1]+1,n[2],0.0);
		fField.readBin("Ufy.raw");

		for (int k = 0; k<int(fField.size()) ; k++ )
		 for ( int j = 0; j<int(fField[k].size())-1 ; j++ )
		  for ( int i = 0; i<int(fField[0][0].size()) ; i++ )
			fField[k][j][i]=0.5*(fField[k][j][i]+fField[k][j+1][i]);
		fField.writeAscii("Uccy.dat", 0,n[0],0,n[1],0,n[2]);
	}
	{
		voxelImageT<float> fField(n[0],n[1],n[2]+1,0.0);
		fField.readBin("Ufz.raw");
		for (int k = 0; k<int(fField.size())-1 ; k++ )
		 for ( int j = 0; j<int(fField[k].size()) ; j++ )
		  for ( int i = 0; i<int(fField[0][0].size()) ; i++ )
			fField[k][j][i]=0.5*(fField[k][j][i]+fField[k+1][j][i]);
		fField.writeAscii("Uccz.dat", 0,n[0],0,n[1],0,n[2]);
	}
	{
		voxelImageT<float> pField(n[0],n[1],n[2],0.0);
		pField.readBin("Ufz.raw");
		pField.writeAscii("p.dat");
	}

	std::cout<< "end" << std::endl;


	return 0;
}


