/*-------------------------------------------------------------------------*\	
You pKcn redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.


The code has been developed by Ahmed AlRatrout as a part his PhD 
at Imperial College London, under the supervision of Branko Bijeljic 
and Martin Blunt. 


For further information please contact us by email:
Ahmed AlRarout:  a.alratrout14@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
\*-------------------------------------------------------------------------*/

    #include <fstream>
    #include <iostream>
    #include <vector>
    #include <algorithm>

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
		<<"only the numbers pKcn be changed in this file. \n"				<<"crop, pore, resample and direction keywords are optional\n"               
		<< std::endl;
	return -1;
}


int main(int argc, char *argv[])
{


	if(argc!=3) usage();


	std::string headerName(argv[1]);	///(fileName_VSubElems0.mhd) 
	if(headerName.size()<4) return usage();
	//~ std::string outputName(argv[2]);

	//~ std::cout<<"voxelImageConvert "<<headerName<<"  "<<outputName<<std::endl;

	voxelImageT<int> vimage(headerName);
	
	int3 n=vimage.size3();
	vec3 X0=vimage.X0();
	vec3 dx=vimage.dx();

	//~ dx*=1.0e-6;
	//~ X0*=1.0e-6;
	
	//char* fileGNM = argv[2];  // file name
	//std::ifstream inGN;          ///(fileName.gnm) 
	//std::cout<<"\n\n Openning VElems file: "<<fileGNM<<std::endl;
	//inGN.open(fileGNM);
	//assert(inGN);

	//std::string lline;
	
	
	char* fileRc = argv[2];  // file name
	std::ifstream inRc;     ///Rc_x.txt
	std::ofstream outRcT("outRcCorner.txt");     
	std::cout<<"\nOpenning Rc file: "<<fileRc<<std::endl;
	inRc.open(fileRc);
	assert(inRc);

	std::string tmpc;
	std::string ttmpc;
	std::string line;
	int countRc = 0;
	double xx,yy,zz, pKc;
	int  nThroats;
	
	
	
	for ( unsigned int i=0; i<4;i++)   inRc>>tmpc, std::cout<<" "<<tmpc;
	//for ( unsigned int i=0; i<5;i++)   inGN>>ttmpc, std::cout<<" "<<ttmpc; 
	//inGN>>nPors, std::cout<<"\n  nPores = "<<nPors<<std::endl;
	//inGN>>nThroats, std::cout<<"\n  nThroats = "<<nThroats<<std::endl;

	//std::vector<double> PRs(nPors,0.0);
	//std::vector<int> i(nPors,0.0);
	//std::vector<int> j(nPors,0.0);
	//std::vector<int> k(nPors,0.0);
	
	
	//~ int i,j,k;
	int ii,jj,kk;
	//int countPI = 0;
	//for (;countPI<nPors; ++countPI)
	//{
		//getline(inGN, lline);
		//char g;
		//inGN >> i[countPI]>>j[countPI]>>k[countPI]>>g>>PRs[countPI]>>tmpc;
	//}
	
	//std::vector<int> PI1(nThroats,0.0);
	//std::vector<int> PI2(nThroats,0.0);
	//std::vector<int> TCs(nThroats,0.0);
	//std::vector<double> TRs(nThroats,0.0);
	//std::vector<double> I(nThroats,0.0);
	//std::vector<double> J(nThroats,0.0);
	//std::vector<double> K(nThroats,0.0);	
	//int countTI = 0;
	//for (;countTI<nThroats; ++countTI)
	//{
		
		//getline(inGN, lline);
		//char g;
		//inGN >> PI1[countTI]>> PI2[countTI]>>g>>I[countTI]>>J[countTI]>>K[countTI]>>TRs[countTI]>>TCs[countTI];
		//std::cout<< "\n"; 
		//std::cout<< " " <<countTI<< " " <<PI1[countTI]<< " " <<PI2[countTI]<< " " <<I[countTI]<< " " <<J[countTI]<< " " <<K[countTI]<<" "<<g<<" "<<TRs[countTI]<<" "<<TCs[countTI]<<"  \n";
	//}

	outRcT<<"Curvature in Coners: "<<nThroats<<std::endl;	
	for ( ; getline( inRc, line );)
	{
		++countRc;
		inRc >> xx>>yy>>zz>>pKc; 
		//~ std::cout<< "\n " <<i<< " " <<j<< " " <<k<<" "<<cosCA<<" "<<CA<<"  \n";
		ii = std::min(int((xx/dx[0])+0.5-X0[0]),n[0]-1);
		jj = std::min(int((yy/dx[1])+0.5-X0[0]),n[1]-1);
		kk = std::min(int((zz/dx[2])+0.5-X0[0]),n[2]-1);
		//~ ii = xx/6.40002;
		//~ jj = yy/6.40002;
		//~ kk = zz/6.40002;
		std::cout<< " " <<kk<< " " <<jj<< " " <<ii<< " : " <<std::endl;
		int TI = vimage[kk][jj][ii];

		if(TI<=0)	{TI = vimage[std::min(kk+1,n[2]-1)][jj][ii]; 										///1,0,0
		if(TI<=0)	TI = vimage[std::min(kk+1,n[2]-1)][std::min(jj+1,n[1]-1)][ii];						///1,1,0
		if(TI<=0)	TI = vimage[std::min(kk+1,n[2]-1)][std::min(jj+1,n[1]-1)][std::min(ii+1,n[0]-1)];	///1,1,1
		if(TI<=0)	TI = vimage[kk][std::min(jj+1,n[1]-1)][ii];											///0,1,0
		if(TI<=0)	TI = vimage[kk][std::min(jj+1,n[1]-1)][std::min(ii+1,n[0]-1)];						///0,1,1
		if(TI<=0)	TI = vimage[kk][jj][std::min(ii+1,n[0]-1)];											///0,0,1
		if(TI<=0)	TI = vimage[std::min(kk+1,n[2]-1)][jj][std::min(ii+1,n[0]-1)];						///1,0,1
		if(TI<=0)	TI = vimage[std::max(kk-1,0)][jj][ii];												///-1,0,0
		if(TI<=0)	TI = vimage[std::max(kk-1,0)][std::max(jj-1,0)][ii];								///-1,-1,0
		if(TI<=0)	TI = vimage[std::max(kk-1,0)][std::max(jj-1,0)][std::max(ii-1,0)];					///-1,-1,-1
		if(TI<=0)	TI = vimage[kk][std::max(jj-1,0)][ii];												///0,-1,0
		if(TI<=0)	TI = vimage[kk][std::max(jj-1,0)][std::max(ii-1,0)];								///0,-1,-1
		if(TI<=0)	TI = vimage[kk][jj][std::max(ii-1,0)];												///0,0,-1
		if(TI<=0)	TI = vimage[std::max(kk-1,0)][jj][std::max(ii-1,0)];								///-1,0,-1
		
		if(TI<=0)	TI = vimage[std::min(kk+1,n[2]-1)][std::max(jj-1,0)][ii];							///1,-1,0
		if(TI<=0)	TI = vimage[std::max(kk-1,0)][std::min(jj+1,n[1]-1)][ii];							///-1,1,-1
		if(TI<=0)	TI = vimage[std::min(kk+1,n[2]-1)][std::min(jj+1,n[1]-1)][std::max(ii-1,0)];		///1,1,-1
		if(TI<=0)	TI = vimage[std::min(kk+1,n[2]-1)][std::max(jj-1,0)][std::min(ii+1,n[0]-1)];		///1,-1,1
		if(TI<=0)	TI = vimage[std::max(kk-1,0)][std::max(jj-1,0)][std::min(ii+1,n[0]-1)];				///-1,-1,1
		if(TI<=0)	TI = vimage[std::min(kk+1,n[2]-1)][std::max(jj-1,0)][std::max(ii-1,0)];				///1,-1,-1
		if(TI<=0)	TI = vimage[std::max(kk-1,0)][std::min(jj+1,n[1]-1)][std::min(ii+1,n[0]-1)];		///-1,1,1
		if(TI<=0)	TI = vimage[std::max(kk-1,0)][std::min(jj+1,n[1]-1)][std::max(ii-1,0)];				///-1,1,-1
		if(TI<=0)	TI = vimage[std::max(kk-1,0)][jj][std::min(ii+1,n[0]-1)];							///-1,0,1
		if(TI<=0)	TI = vimage[std::min(kk+1,n[2]-1)][jj][std::max(ii-1,0)];							///1,0,-1
		if(TI<=0)	TI = vimage[kk][std::max(jj-1,0)][std::min(ii+1,n[0]-1)];							///0,-1,1
		if(TI<=0)	TI = vimage[kk][std::min(jj+1,n[1]-1)][std::max(ii-1,0)];}							///0,1,-1
		
		//if(TI<=0)	{TI = vimage[std::min(kk+2,n[2]-1)][jj][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+2,n[2]-1)][std::min(jj+2,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+2,n[2]-1)][std::min(jj+2,n[1]-1)][std::min(ii+2,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+2,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+2,n[1]-1)][std::min(ii+2,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][jj][std::min(ii+2,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+2,n[2]-1)][jj][std::min(ii+2,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-2,0)][jj][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-2,0)][std::max(jj-2,0)][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-2,0)][std::max(jj-2,0)][std::max(ii-2,0)];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-2,0)][ii];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-2,0)][std::max(ii-2,0)];
		//if(TI<=0)	TI = vimage[kk][jj][std::max(ii-2,0)];
		//if(TI<=0)	TI = vimage[std::max(kk-2,0)][jj][std::max(ii-2,0)];
		
		//if(TI<=0)	TI = vimage[std::min(kk+2,n[2]-1)][std::max(jj-2,0)][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-2,0)][std::min(jj+2,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+2,n[2]-1)][std::min(jj+2,n[1]-1)][std::max(ii-2,0)];
		//if(TI<=0)	TI = vimage[std::min(kk+2,n[2]-1)][std::max(jj-2,0)][std::min(ii+2,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-2,0)][std::max(jj-2,0)][std::min(ii+2,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+2,n[2]-1)][std::max(jj-2,0)][std::max(ii-2,0)];
		//if(TI<=0)	TI = vimage[std::max(kk-2,0)][std::min(jj+2,n[1]-1)][std::min(ii+2,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-2,0)][std::min(jj+2,n[1]-1)][std::max(ii-2,0)];	
		//if(TI<=0)	TI = vimage[std::max(kk-2,0)][jj][std::min(ii+2,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+2,n[2]-1)][jj][std::max(ii-2,0)];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-2,0)][std::min(ii+2,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+2,n[1]-1)][std::max(ii-2,0)];}}


		//if(TI<=0)	{TI = vimage[std::min(kk+3,n[2]-1)][jj][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+3,n[2]-1)][std::min(jj+3,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+3,n[2]-1)][std::min(jj+3,n[1]-1)][std::min(ii+3,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+3,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+3,n[1]-1)][std::min(ii+3,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][jj][std::min(ii+3,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+3,n[2]-1)][jj][std::min(ii+3,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-3,0)][jj][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-3,0)][std::max(jj-3,0)][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-3,0)][std::max(jj-3,0)][std::max(ii-3,0)];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-3,0)][ii];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-3,0)][std::max(ii-3,0)];
		//if(TI<=0)	TI = vimage[kk][jj][std::max(ii-3,0)];
		//if(TI<=0)	TI = vimage[std::max(kk-3,0)][jj][std::max(ii-3,0)];
		
		//if(TI<=0)	TI = vimage[std::min(kk+3,n[2]-1)][std::max(jj-3,0)][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-3,0)][std::min(jj+3,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+3,n[2]-1)][std::min(jj+3,n[1]-1)][std::max(ii-3,0)];
		//if(TI<=0)	TI = vimage[std::min(kk+3,n[2]-1)][std::max(jj-3,0)][std::min(ii+3,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-3,0)][std::max(jj-3,0)][std::min(ii+3,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+3,n[2]-1)][std::max(jj-3,0)][std::max(ii-3,0)];
		//if(TI<=0)	TI = vimage[std::max(kk-3,0)][std::min(jj+3,n[1]-1)][std::min(ii+3,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-3,0)][std::min(jj+3,n[1]-1)][std::max(ii-3,0)];	
		//if(TI<=0)	TI = vimage[std::max(kk-3,0)][jj][std::min(ii+3,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+3,n[2]-1)][jj][std::max(ii-3,0)];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-3,0)][std::min(ii+3,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+3,n[1]-1)][std::max(ii-3,0)];
		
		
		//if(TI<=0)	{TI = vimage[std::min(kk+4,n[2]-1)][jj][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+4,n[2]-1)][std::min(jj+4,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+4,n[2]-1)][std::min(jj+4,n[1]-1)][std::min(ii+4,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+4,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+4,n[1]-1)][std::min(ii+4,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][jj][std::min(ii+4,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+4,n[2]-1)][jj][std::min(ii+4,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-4,0)][jj][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-4,0)][std::max(jj-4,0)][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-4,0)][std::max(jj-4,0)][std::max(ii-4,0)];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-4,0)][ii];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-4,0)][std::max(ii-4,0)];
		//if(TI<=0)	TI = vimage[kk][jj][std::max(ii-4,0)];
		//if(TI<=0)	TI = vimage[std::max(kk-4,0)][jj][std::max(ii-4,0)];
		
		//if(TI<=0)	TI = vimage[std::min(kk+4,n[2]-1)][std::max(jj-4,0)][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-4,0)][std::min(jj+4,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+4,n[2]-1)][std::min(jj+4,n[1]-1)][std::max(ii-4,0)];
		//if(TI<=0)	TI = vimage[std::min(kk+4,n[2]-1)][std::max(jj-4,0)][std::min(ii+4,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-4,0)][std::max(jj-4,0)][std::min(ii+4,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+4,n[2]-1)][std::max(jj-4,0)][std::max(ii-4,0)];
		//if(TI<=0)	TI = vimage[std::max(kk-4,0)][std::min(jj+4,n[1]-1)][std::min(ii+4,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-4,0)][std::min(jj+4,n[1]-1)][std::max(ii-4,0)];	
		//if(TI<=0)	TI = vimage[std::max(kk-4,0)][jj][std::min(ii+4,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+4,n[2]-1)][jj][std::max(ii-4,0)];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-4,0)][std::min(ii+4,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+4,n[1]-1)][std::max(ii-4,0)];
		
		
		//if(TI<=0)	{TI = vimage[std::min(kk+5,n[2]-1)][jj][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+5,n[2]-1)][std::min(jj+5,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+5,n[2]-1)][std::min(jj+5,n[1]-1)][std::min(ii+5,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+5,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+5,n[1]-1)][std::min(ii+5,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][jj][std::min(ii+5,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+5,n[2]-1)][jj][std::min(ii+5,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-5,0)][jj][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-5,0)][std::max(jj-5,0)][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-5,0)][std::max(jj-5,0)][std::max(ii-5,0)];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-5,0)][ii];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-5,0)][std::max(ii-5,0)];
		//if(TI<=0)	TI = vimage[kk][jj][std::max(ii-5,0)];
		//if(TI<=0)	TI = vimage[std::max(kk-5,0)][jj][std::max(ii-5,0)];
		
		//if(TI<=0)	TI = vimage[std::min(kk+5,n[2]-1)][std::max(jj-5,0)][ii];
		//if(TI<=0)	TI = vimage[std::max(kk-5,0)][std::min(jj+5,n[1]-1)][ii];
		//if(TI<=0)	TI = vimage[std::min(kk+5,n[2]-1)][std::min(jj+5,n[1]-1)][std::max(ii-5,0)];
		//if(TI<=0)	TI = vimage[std::min(kk+5,n[2]-1)][std::max(jj-5,0)][std::min(ii+5,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-5,0)][std::max(jj-5,0)][std::min(ii+5,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+5,n[2]-1)][std::max(jj-5,0)][std::max(ii-5,0)];
		//if(TI<=0)	TI = vimage[std::max(kk-5,0)][std::min(jj+5,n[1]-1)][std::min(ii+5,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::max(kk-5,0)][std::min(jj+5,n[1]-1)][std::max(ii-5,0)];	
		//if(TI<=0)	TI = vimage[std::max(kk-5,0)][jj][std::min(ii+5,n[0]-1)];
		//if(TI<=0)	TI = vimage[std::min(kk+5,n[2]-1)][jj][std::max(ii-5,0)];
		//if(TI<=0)	TI = vimage[kk][std::max(jj-5,0)][std::min(ii+5,n[0]-1)];
		//if(TI<=0)	TI = vimage[kk][std::min(jj+5,n[1]-1)][std::max(ii-5,0)];}}}}}
		
		if(TI>0)
		{
			outRcT << xx<<" "<<yy<<" "<<zz<<" "<<pKc<<" "<<TI<<"\n"; 
		}
		
	}
	std::cout<< "\nTotal Curvatures in corners = " << countRc;

	
	
	
	//~ double i_ID_max=-1,j_ID_max=-1,k_ID_max=-1;
	
	
	//~ std::cout<<"i_ID_max = "<<i_ID_max<<", "<<"j_ID_max = "<<j_ID_max<<", "<<"k_ID_max = "<<k_ID_max<<", \n";
	std::cout<< "\nTotal pKc = " << countRc;
	
	
	inRc.close();

	//~ if (!ifContAngle) {
    //~ std::cerr << "Unable to open file datafile.txt";
    //~ exit(1);   // pKcll system to stop
	//~ } else{
		//~ for ( unsigned int i=0; i<2;i++)   std::cin>>tmpc, std::cout<<" "<<tmpc;
		
	




//~ writeRcPoreInds (ImageName)_VElems.mhd contactAngles_x_PI.txt













		//~ int PI = vimage[k][j][i];

	std::cout<< "\n ***** END ***** \n\n\n" << std::endl;

	//~ vimage.printInfo();


	return 0;
}


