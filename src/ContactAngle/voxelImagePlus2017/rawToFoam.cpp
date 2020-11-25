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
Ali Q Raeini:	a.qaseminejad-raeini09@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
\*-------------------------------------------------------------------------*/

//~ #define _2D_

#include <sys/stat.h>

#include <fstream>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <valarray>
#include <iostream>

#include "voxelImage.h"

using namespace std;


int usage()
{
	cout<<"convert micro-Ct images to OpenFOAM mesh"<<endl;
	cout<<"usage: example:"<<endl;
	cout<<"	   rawToFoam imageName.mhd"<<endl;
	return 1;
}

void fixImage(voxelImage& vimage);
int main(int argc, char** argv)
{

	if(argc!=2)		return usage();
	std::string headerName(argv[1]);
	if(headerName.size()<4 || headerName.compare(headerName.size()-4,4,".mhd") != 0) return usage();

	voxelImage vimage(headerName);

	
	int3 n=vimage.size3();
	vec3 X0=vimage.X0();
	vec3 dx=vimage.dx();

	dx*=1.0e-6;
	X0*=1.0e-6;


	vimage.printInfo();

	vimage.writeHeader("vxlImage.mhd");

	cout <<"converting to OpenFOAM format "<<endl;




	vimage.threshold101(0,0);
	
	vimage.growBox(1);
	vimage.FaceMedian06(1,5);
	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);

	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);
	vimage.FaceMedian06(2,4);

	vimage.FaceMedian06(1,5);

	vimage.crop(1,n[0],1,n[1],1,n[2],1,1);	//		 XXXXXXXXXXXXXXXXXXXXXXXXXXXX


	fixImage(vimage);


	vimage.printInfo();

	::mkdir("constant",0777);
	::mkdir("constant/polyMesh",0777);

	if (!::mkdir("constant",0777))
	{
		cout  << "Couldn't create 'constant' directory ";
		exit(0);
	}
	if (!::mkdir("constant/polyMesh",0777))
	{
		cout  << "Couldn't create 'constant/polyMesh' directory ";
		exit(0);
	}
	ofstream pointsf("constant/polyMesh/points");
	assert(pointsf);


	pointsf<<
	"FoamFile\n"
	"{\n"
	"	version	 2.0;\n"
	"	format	  ascii;\n"
	"	class	   vectorField;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  points;\n"
	"}\n\n";


//=======================================================================

	voxelField<int> point_mapper(n[0]+1,n[1]+1,n[2]+1,-1);



	int nPoints=0;

{
	std::vector<vector<vector<unsigned char> > >::iterator ddd=vimage.begin();//,d1p1=vimage.begin()+1;//		int iz=ppp-point_mapper.begin();
	register int iz=0;
	for (;ddd<vimage.end()-1;ddd++,iz++)
	{
		int iy=0;
		std::vector<vector<unsigned char> >::iterator dd1=ddd->begin(), dd2=(ddd+1)->begin();
		for (;dd1<ddd->end()-1;dd1++,dd2++,iy++)
		{
			int ix=0;
			std::vector<unsigned char>::iterator d11=dd1->begin(), d12=(dd1+1)->begin();
			std::vector<unsigned char>::iterator d21=dd2->begin(), d22=(dd2+1)->begin();
			for (;d11<dd1->end()-1;d11++,d12++,d21++,d22++,ix++)
			{
				if (*d11==0 || *(d11+1)==0 ||
					*d12==0 || *(d12+1)==0 ||
					*d21==0 || *(d21+1)==0 ||
					*d22==0 || *(d22+1)==0 )
				{
					++nPoints;
				}
			}
		}
	}
}
	cout<<"nPoints: "<<nPoints<<endl;
	cout <<"writing points";cout.flush();


	pointsf<<nPoints/*(n[0]*n[1]*n[2])*/<<endl<<"("<<endl;
	register unsigned int iPoints=0;
	std::vector<vector<vector<unsigned char> > >::iterator ddd=vimage.begin();//,d1p1=vimage.begin()+1;//		int iz=ppp-point_mapper.begin();
	register unsigned int iz=0;
	for (;ddd<vimage.end()-1;ddd++,iz++)
	{
		register unsigned int iy=0;
		double z=iz*dx[2]+X0[2];
		std::vector<vector<unsigned char> >::iterator dd1=ddd->begin(), dd2=(ddd+1)->begin();
		for (;dd1<ddd->end()-1;dd1++,dd2++,iy++)
		{
			register unsigned int ix=0;
			double y=iy*dx[1]+X0[1];
			std::vector<unsigned char>::iterator d11=dd1->begin(), d12=(dd1+1)->begin();
			std::vector<unsigned char>::iterator d21=dd2->begin(), d22=(dd2+1)->begin();
			for (;d11<dd1->end()-1;d11++,d12++,d21++,d22++,ix++)
			{
				if (*d11==0 || *(d11+1)==0 ||
					*d12==0 || *(d12+1)==0 ||
					*d21==0 || *(d21+1)==0 ||
					*d22==0 || *(d22+1)==0 )
				{
					double x=ix*dx[0]+X0[0];
					pointsf<< "("<<x<< ' '<<y<<' '<<z<<")\n";
					point_mapper[iz][iy][ix]=iPoints;iPoints++;
				}
			}
		}
	}
	pointsf<<endl<<")"<<endl;

	pointsf.close();

	cout <<"   done"<<endl;	cout.flush();
//=======================================================================



//______________________________________________________________________________________________________________________

	ofstream faces("constant/polyMesh/faces");
	assert(faces);
	faces<<
	"FoamFile\n"
	"{\n"
	"	version	 2.0;\n"
	"	format	  ascii;\n"
	"	class	   faceList;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  faces;\n"
	"}\n\n";


	ofstream owner("constant/polyMesh/owner");
	assert(owner);
	owner<<
	"FoamFile\n"
	"{\n"
	"	version	 2.0;\n"
	"	format	  ascii;\n"
	"	class	   labelList;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  owner;\n"
	"}\n\n";
 
	ofstream neighbour("constant/polyMesh/neighbour");
	assert(neighbour);
	neighbour<<
	"FoamFile\n"
	"{\n"
	"	version	 2.0;\n"
	"	format	  ascii;\n"
	"	class	   labelList;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  neighbour;\n"
	"}\n\n";


	ofstream boundary("constant/polyMesh/boundary");
	assert(boundary);
	boundary<<
	"FoamFile\n"
	"{\n"
	"	version	 2.0;\n"
	"	format	  ascii;\n"
	"	class	   polyBoundaryMesh;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  boundary;\n"
	"}\n\n";





	size_t nCells=0;
	size_t nInternalFaces=0;
	size_t nGrainwallsFaces=0;
	size_t nBackFaces=0;
	size_t nFrontFaces=0;
	size_t nBottomFaces=0;
	size_t nTopFaces=0;
	size_t nLeftFaces=0;
	size_t nRightFaces=0;

#define  mcountFace(type) n##type##Faces++;


	for (register int iz=1;iz<=n[2];iz++)
	{
		for (register int iy=1;iy<=n[1];iy++)

		{
			for (register int ix=1;ix<=n[0];ix++)
			{
					if (!vimage[iz][iy][ix])
					{
						nCells++;
						if (iz!=1)
						{
						  if (!vimage[iz-1][iy][ix])  {}
						  else					mcountFace( Grainwalls)
						}else mcountFace( Back);
						if (iz!=n[2])
						{
						  if (!vimage[iz+1][iy][ix])	mcountFace( Internal)
						  else					mcountFace( Grainwalls)
						}else mcountFace( Front);

						if (iy!=1)
						{
						  if (!vimage[iz][iy-1][ix])  {}
						  else					mcountFace( Grainwalls)
						}else mcountFace( Bottom);
						if (iy!=n[1])
						{
						  if (!vimage[iz][iy+1][ix])	mcountFace( Internal)
						  else					mcountFace( Grainwalls)
						}else mcountFace( Top);

						if (ix!=1)
						{
						  if (!vimage[iz][iy][ix-1])  {}
						  else					mcountFace( Grainwalls)
						}else mcountFace( Left);
						if (ix!=n[0])
						{
						  if (!vimage[iz][iy][ix+1])	mcountFace( Internal)
						  else					mcountFace( Grainwalls)
						}else mcountFace( Right);

			   }
			}
		}
	}

	cout<<"nCells: "<<nCells<<endl;
	cout<<"nInternalFaces: "<<nInternalFaces<<endl;
	cout<<"nGrainwallsFaces: "<<nGrainwallsFaces<<endl;
	cout<<"nBackFaces: "<<nBackFaces<<endl;
	cout<<"nFrontFaces: "<<nFrontFaces<<endl;
	cout<<"nBottomFaces: "<<nBottomFaces<<endl;
	cout<<"nTopFaces: "<<nTopFaces<<endl;
	cout<<"nLeftFaces: "<<nLeftFaces<<endl;
	cout<<"nRightFaces: "<<nRightFaces<<endl;



	int iInternalStartFace=0;
	int iGrainwallsStartFace=nInternalFaces;

	int iLeftStartFace=iGrainwallsStartFace+nGrainwallsFaces;
	int iRightStartFace=iLeftStartFace+nLeftFaces;

	int iBottomStartFace=iRightStartFace+nRightFaces;
	int iTopStartFace=iBottomStartFace+nBottomFaces;

	int iBackStartFace=iTopStartFace+nTopFaces;
	int iFrontStartFace=iBackStartFace+nBackFaces;

	int nFaces=iFrontStartFace+nFrontFaces;



	{
		#define write_boundary(name,type)                              \
		if (n ##name ##Faces >0)                                       \
		boundary<<		                                                \
		"	"<<  #name										 <<endl<<			   \
		"	{"												 <<endl<<		         \
		"		type			"<<#type<<";"				  <<endl<<			   \
		"		nFaces		  "<<n ##name ##Faces<<';'	   <<endl<<	      \
		"		startFace	   "<<i ##name ##StartFace<<';'   <<endl<<   \
		"	}"												 <<endl;

			boundary<<int(nGrainwallsFaces>0)+int(nLeftFaces>0)+int(nRightFaces>0)+
				int(nBottomFaces>0)+int(nTopFaces>0)+int(nBackFaces>0)+int(nFrontFaces>0)  <<endl
					<<'('<<endl;
			write_boundary(Grainwalls,patch);
			write_boundary(Left,patch);
			write_boundary(Right,patch);
			write_boundary(Bottom,patch);
			write_boundary(Top,patch);
		#ifdef _2D_
			write_boundary(Back,empty);
			write_boundary(Front,empty);
		#else
			write_boundary(Back,patch);
			write_boundary(Front,patch);
		#endif

			boundary<<")"   <<endl;
	}







	cout<<"creating faces"<<endl;


	std::vector<vector<int> > faces_Internal(nInternalFaces, vector<int>(6,-1.0));
	std::vector<vector<int> > faces_Grainwalls(nGrainwallsFaces, vector<int>(5,-1.0));
	std::vector<vector<int> > faces_Back(nBackFaces, vector<int>(5,-1.0));
	std::vector<vector<int> > faces_Front(nFrontFaces , vector<int>(5,-1.0));
	std::vector<vector<int> > faces_Bottom(nBottomFaces,  vector<int>(5,-1.0));
	std::vector<vector<int> > faces_Top(nTopFaces,   vector<int>(5,-1.0));
	std::vector<vector<int> > faces_Left(nLeftFaces,   vector<int>(5,-1.0));
	std::vector<vector<int> > faces_Right(nRightFaces,  vector<int>(5,-1.0));


	voxelField<vector<int> > ownerMapper(n[0]+1,n[1]+1,n[2]+1,vector<int>(3,-1));


	cout<<"collecting faces"<<endl;


#define recordF_m( l10,l11,l20,l21,l30,l31,dir,ii,jj,kk,type )				   \
  {																	 \
	if (ownerMapper[ii][jj][kk][dir]<0)										 \
	{		   ++i##type ##Faces ;																		 \
			ownerMapper[ii][jj][kk][dir]=i##type ##Faces + i##type ##StartFace;							 \
			faces_##type[i##type ##Faces][4]=iCells;	/* cell number (for owners) */					 \
			faces_##type[i##type ##Faces][0]=point_mapper[ii][jj][kk];								  \
			faces_##type[i##type ##Faces][1]=point_mapper[ii+l11][jj+l21][kk+l31];					  \
			faces_##type[i##type ##Faces][2]=point_mapper[ii+l10+l11][jj+l20+l21][kk+l30+l31];		  \
			faces_##type[i##type ##Faces][3]=point_mapper[ii+l10][jj+l20][kk+l30];					  \
	}																			   \
	else																			\
	{																			   \
			faces_##type[ownerMapper[ii][jj][kk][dir]][5]=iCells;		/* cell number (for neighbours) */		\
	}																			   \
  }


#define  iclockwiserecordF(type) recordF_m( 0,1,1,0,0,0, 2, iz-1,iy-1,ix-1, type)
#define iuclockwiserecordF(type) recordF_m( 1,0,0,1,0,0, 2, iz-1,iy-1,ix  , type)
#define  jclockwiserecordF(type) recordF_m( 1,0,0,0,0,1, 1, iz-1,iy-1,ix-1, type)
#define juclockwiserecordF(type) recordF_m( 0,1,0,0,1,0, 1, iz-1,iy  ,ix-1, type)
#define  kclockwiserecordF(type) recordF_m( 0,0,0,1,1,0, 0, iz-1,iy-1,ix-1, type)
#define kuclockwiserecordF(type) recordF_m( 0,0,1,0,0,1, 0, iz  ,iy-1,ix-1, type)




	int iCells=-1;

	int iInternalFaces=-1;
	int iGrainwallsFaces=-1;
	int iBackFaces=-1;
	int iFrontFaces=-1;
	int iBottomFaces=-1;
	int iTopFaces=-1;
	int iLeftFaces=-1;
	int iRightFaces=-1;

	for (register int iz=1;iz<=n[2];iz++)
	{cout<<'.';cout.flush();
		for (register int iy=1;iy<=n[1];iy++)
		{
			for (register int ix=1;ix<=n[0];ix++)
			{
				if (!vimage[iz][iy][ix])
				{

					iCells++;

					//int iface=0;
					if (ix!=1)
					{
					  if (!vimage[iz][iy][ix-1])  {iclockwiserecordF( Internal);}
					  else					{iclockwiserecordF( Grainwalls)}
					}else {iclockwiserecordF( Left);}
					if (ix!=n[0])
					{
					  if (!vimage[iz][iy][ix+1])  {iuclockwiserecordF( Internal);}
					  else					{iuclockwiserecordF( Grainwalls)}
					}else {iuclockwiserecordF( Right); }

					if (iy!=1)
					{
					  if (!vimage[iz][iy-1][ix])  {jclockwiserecordF( Internal);}
					  else					{jclockwiserecordF( Grainwalls)}
					}else {jclockwiserecordF( Bottom)}
					if (iy!=n[1])
					{
					  if (!vimage[iz][iy+1][ix])  {juclockwiserecordF( Internal);}
					  else					{juclockwiserecordF( Grainwalls)}
					}else {juclockwiserecordF( Top)}

					if (iz!=1)
					{
					  if (!vimage[iz-1][iy][ix])  {kclockwiserecordF( Internal);}
					  else					{kclockwiserecordF( Grainwalls)}
					}else {kclockwiserecordF( Back)}
					if (iz!=n[2])
					{
					  if (!vimage[iz+1][iy][ix])  {kuclockwiserecordF( Internal);}
					  else					{kuclockwiserecordF( Grainwalls)}
					}else {kuclockwiserecordF( Front)}


				}
			}
		}
	}




point_mapper.resize(0);



//______________________________________________________________________________________________________________________
#define write_faces_owners(type) \
  for (std::vector<vector<int> >::iterator ff=faces_##type .begin();ff<faces_##type .end();ff++)	\
	{ \
		faces<<'('<<((*ff))[0]<<' '<<(*ff)[1]<<' '<<(*ff)[2]<<' '<<(*ff)[3]<<")\n"; \
		owner<<(*ff)[4]<<"\n";  \
	}


	owner<<nFaces<<endl;
	owner<<"("<<endl;
	faces<<nFaces<<endl
		 <<"("<<endl;

		write_faces_owners(Internal)
		write_faces_owners(Grainwalls)
		write_faces_owners(Left)
		write_faces_owners(Right)
		write_faces_owners(Bottom)
		write_faces_owners(Top)
		write_faces_owners(Back)
		write_faces_owners(Front)

	faces<<")"<<endl;
	faces.close();
	owner<<")"<<endl;
	owner.close();


	neighbour<<nInternalFaces<<endl;
	neighbour<<"("<<endl;
	  for (std::vector<vector<int> >::iterator ff=faces_Internal.begin();ff<faces_Internal.end();ff++)
		{
			neighbour<<(*ff)[5]<<"\n";
		}

	neighbour<<")"<<endl;
	neighbour.close();




   return 0;
}


template<typename T>
void replaceValue(voxelImageT<T>& vImage, T v, T newv)
{
	{
	 (cout<<int(v)<<"->"<<int(newv)<<"    ").flush();
    for ( unsigned int k=0; k<vImage.size() ; k++ )
        for ( unsigned int j=0; j<vImage[k].size() ; j++ )
            for ( unsigned int i=0; i<vImage[k][j].size() ; i++ )
                if (vImage[k][j][i]==v)  vImage[k][j][i]=newv;
	}
}


void fixImage(voxelImage& voxels)
{
	//voxels.write("dump1.mhd");
	int3 n = voxels.size3();
	const unsigned int bigN=255*255*255*127; 


	cout<<"removing disconected parts of the image "<<endl;
	int nxmid=n[0]/2;
	unsigned int vmax=n[1]*n[2]+1;
	voxelImageT<unsigned int> vxlsMids(1,n[1], n[2],bigN);
	voxelImageT<unsigned int> vxlsMidMap(1,n[1], n[2],vmax);
	vector<unsigned int> vxlsMidCompresdReg(n[1]*n[2],bigN);
	//std::valarray<unsigned int> vxlsMidMap(0,n[1]*n[2]);
	//std::valarray<unsigned int> vxlsMidMapCount(0,n[1]*n[2]);
	//for ( unsigned int i=0; i<vxlsMidMap->size() ; i++ )	vxlsMidMap[i]=i;

	for ( unsigned int k=0; k<vxlsMidMap.size() ; k++ )
	 for ( unsigned int j=0; j<vxlsMidMap[k].size() ; ++j )
	  //if(voxels[k][j][nxmid]==0)
		vxlsMidMap[k][j][0]=k*n[1]+j;
	for ( unsigned int k=0; k<vxlsMids.size() ; k++ )
	 for ( unsigned int j=0; j<vxlsMids[k].size() ; ++j )
		vxlsMids[k][j][0]=voxels[k][j][nxmid];
	long long nchanges=1;
	while(nchanges)
	{	nchanges = 0;
		for ( unsigned int k=1; k<vxlsMids.size() ; k++ )
		 for ( unsigned int j=1; j<vxlsMids[k].size() ; ++j )
			if (vxlsMids[k][j][0]==0)
			{
				if (vxlsMids[k-1][j][0]==0 && vxlsMidMap[k][j][0]>vxlsMidMap[k-1][j][0])
				{	vxlsMidMap[k][j][0]=vxlsMidMap[k-1][j][0]; ++nchanges;	}
				if (vxlsMids[k][j-1][0]==0 && vxlsMidMap[k][j][0]>vxlsMidMap[k][j-1][0])
				{	vxlsMidMap[k][j][0]=vxlsMidMap[k][j-1][0]; ++nchanges;	}
			}
		for ( unsigned int k=0; k<vxlsMids.size()-1 ; k++ )
		 for ( unsigned int j=0; j<vxlsMids[k].size()-1 ; ++j )
			if (vxlsMids[k][j][0]==0)
			{
				if (vxlsMids[k+1][j][0]==0 && vxlsMidMap[k][j][0]>vxlsMidMap[k+1][j][0])
				{	vxlsMidMap[k][j][0]=vxlsMidMap[k+1][j][0]; ++nchanges;	}
				if (vxlsMids[k][j+1][0]==0 && vxlsMidMap[k][j][0]>vxlsMidMap[k][j+1][0])
				{	vxlsMidMap[k][j][0]=vxlsMidMap[k][j+1][0]; ++nchanges;	}
			}
		//cout<<nchanges<<endl;
	}
	//cout<<endl;
	
	unsigned int nRegs=1;
	for ( unsigned int k=0; k<vxlsMids.size() ; k++ )
	 for ( unsigned int j=0; j<vxlsMids[k].size() ; ++j )
		if (vxlsMids[k][j][0]==0)
		{
			int jj= vxlsMidMap[k][j][0] % n[1];
			int kk= vxlsMidMap[k][j][0] / n[1];
			if(vxlsMids[kk][jj][0]!=0) cout<<"!"<<"  "<<k<<":"<<kk<<"    "<<j<<":"<<jj<<"    "<<vxlsMidMap[k][j][0]<<"  "<<endl;;
			if (vxlsMidCompresdReg[vxlsMidMap[kk][jj][0]]==bigN) { vxlsMidCompresdReg[vxlsMidMap[kk][jj][0]]=nRegs; ++nRegs; };
			vxlsMidCompresdReg[vxlsMidMap[k][j][0]]=vxlsMidCompresdReg[vxlsMidMap[kk][jj][0]];
		}
		//else vxlsMidMap[k][j][0]=bigN+1;
	if (nRegs>bigN) {cout<<"Error: nRegs >bigN"<<endl; exit(-1);}
	//for ( unsigned int k=0; k<vxlsMids.size() ; k++ )
	 //for ( unsigned int j=0; j<vxlsMids[k].size() ; ++j )
		//if (vxlsMids[k][j][0]==0)
		//{
			//int jj= vxlsMidMap[k][j][0] % n[1];
			//int kk= vxlsMidMap[k][j][0] / n[1];
			//if(vxlsMids[kk][jj][0]!=0) cout<<"!"<<"  "<<k<<":"<<kk<<"    "<<j<<":"<<jj<<"    "<<vxlsMidMap[k][j][0]<<"  "<<endl;;
			//vxlsMidMap[k][j][0]=vxlsMidMap[kk][jj][0];
		//}



	voxelImageT<unsigned int> vxlImg(n[0],n[1],n[2],bigN);
	//for ( unsigned int k=0; k<vxlImg.size() ; k++ )
	 //for ( unsigned int j=0; j<vxlImg[k].size() ; ++j )
		//for ( unsigned int i=0; i<vxlImg[k][j].size() ; ++i )
			//vxlImg[k][j][i]=voxels[k][j][i];


	//for ( unsigned int k=0; k<vxlImg.size() ; k++ )
	 //for ( unsigned int j=0; j<vxlImg[k].size() ; ++j )
		//for ( unsigned int i=0; i<vxlImg[k][j].size() ; ++i )
			//vxlImg[k][j][i]=bigN;

	for ( unsigned int k=0; k<vxlImg.size() ; k++ )
	 for ( unsigned int j=0; j<vxlImg[k].size() ; ++j )
			vxlImg[k][j][nxmid]=vxlsMidCompresdReg[vxlsMidMap[k][j][0]];


	  for ( unsigned int k=2; k<vxlImg.size()-2 ; k++ )
	   for ( unsigned int j=2; j<vxlImg[k].size()-2 ; ++j )
		{
		 for (int iter=0; iter<2; ++iter)
		 {
		  for ( unsigned int i=nxmid-1; i<vxlImg[k][j].size()-1 ; ++i )
		   if ( voxels[k][j][i]==0 )
		   {
			  if (voxels[k][j][i+1]==0 &&  vxlImg[k][j][i]<vxlImg[k][j][i+1]) vxlImg[k][j][i+1]=vxlImg[k][j][i];
		   }
		  for ( unsigned int i=nxmid+1; i>0 ; --i )
		   if ( voxels[k][j][i]==0 )
		   {
			  if (voxels[k][j][i-1]==0 && vxlImg[k][j][i]<vxlImg[k][j][i-1]) vxlImg[k][j][i-1]=vxlImg[k][j][i];
		   }
		 }
		}
	nchanges=1;
	while(nchanges)
	{ nchanges = 0;
	  for ( unsigned int k=1; k<vxlImg.size()-1 ; ++k )
	    for ( unsigned int j=1; j<vxlImg[k].size()-1 ; ++j )
			for ( unsigned int i=1; i<vxlImg[k][j].size()-1 ; ++i )
				 if (voxels[k][j][i]==0)
				 {
					 if(vxlImg[k][j][i]==bigN)
					 {
						unsigned int minv = vxlImg[k][j][i];
						minv=min(minv,vxlImg[k][j][i-1]);
						minv=min(minv,vxlImg[k][j][i+1]);
						minv=min(minv,vxlImg[k][j-1][i]);
						minv=min(minv,vxlImg[k][j+1][i]);
						minv=min(minv,vxlImg[k-1][j][i]);
						minv=min(minv,vxlImg[k+1][j][i]);
						if(minv<vxlImg[k][j][i]) {vxlImg[k][j][i]=minv; ++nchanges;}
					 }
					 else
					 {
						unsigned int minv = vxlImg[k][j][i];
						minv=min(minv,vxlImg[k][j][i-1]);
						minv=min(minv,vxlImg[k][j][i+1]);
						minv=min(minv,vxlImg[k][j-1][i]);
						minv=min(minv,vxlImg[k][j+1][i]);
						minv=min(minv,vxlImg[k-1][j][i]);
						minv=min(minv,vxlImg[k+1][j][i]);
						
						if(minv < vxlImg[k][j][i]) 
							{ cout<<"s1 "; replaceValue( vxlImg, vxlImg[k][j][i], minv );  ++nchanges; }
					 }
				 }
	  for ( unsigned int k=vxlImg.size()-2; k>0 ; --k )
	    for ( unsigned int j=vxlImg[k].size()-2; j>0 ; --j )
			for ( unsigned int i=vxlImg[k][j].size()-2; i>0 ; --i )
				 if (voxels[k][j][i]==0)
				 {
					 if(vxlImg[k][j][i]==bigN)
					 {
						unsigned int minv = vxlImg[k][j][i];
						minv=min(minv,vxlImg[k][j][i-1]);
						minv=min(minv,vxlImg[k][j][i+1]);
						minv=min(minv,vxlImg[k][j-1][i]);
						minv=min(minv,vxlImg[k][j+1][i]);
						minv=min(minv,vxlImg[k-1][j][i]);
						minv=min(minv,vxlImg[k+1][j][i]);
						if(minv<vxlImg[k][j][i]) {vxlImg[k][j][i]=minv; ++nchanges;}
					 }
					 else
					 {
						unsigned int minv = vxlImg[k][j][i];
						minv=min(minv,vxlImg[k][j][i-1]);
						minv=min(minv,vxlImg[k][j][i+1]);
						minv=min(minv,vxlImg[k][j-1][i]);
						minv=min(minv,vxlImg[k][j+1][i]);
						minv=min(minv,vxlImg[k-1][j][i]);
						minv=min(minv,vxlImg[k+1][j][i]);

						if(minv < vxlImg[k][j][i]) 
							{ cout<<"s2 "; replaceValue( vxlImg, vxlImg[k][j][i], minv );  ++nchanges; }
					 }
				 }
                
	 cout<<": "<<nchanges<<"   "; cout.flush();

	}
	
	std::valarray<unsigned int> vxlsMidMapCount(0u,nRegs);
	for ( unsigned int k=0; k<vxlImg.size() ; k++ )
	 for ( unsigned int j=0; j<vxlImg[k].size() ; ++j )
		for ( unsigned int i=0; i<vxlImg[k][j].size() ; ++i )
		 if(vxlImg[k][j][i]<nRegs)
			++vxlsMidMapCount[vxlImg[k][j][i]];
	
	unsigned int maxReg=0; unsigned int maxRegCount=0;
	for(unsigned int i=0; i < vxlsMidMapCount.size();++i) if(vxlsMidMapCount[i] > maxRegCount) {maxRegCount=vxlsMidMapCount[i]; maxReg=i; };
	
	cout<<"maxReg  "<<maxReg<<"    maxRegCount:  "<<maxRegCount<<endl;
	//.threshold101(maxReg,maxReg);
	for ( unsigned int k=0; k<vxlImg.size() ; k++ )
	 for ( unsigned int j=0; j<vxlImg[k].size() ; ++j )
		for ( unsigned int i=0; i<vxlImg[k][j].size() ; ++i )
		if(vxlImg[k][j][i]==maxReg)
			voxels[k][j][i]=0;
		else
			voxels[k][j][i]=1;

	//voxels.write("dump2.mhd");

}

