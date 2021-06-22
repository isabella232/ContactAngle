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
#include <array>
#include <valarray>
#include <iostream>

#include "voxelImage.h"

using namespace std;


int usage()
{
	cout<<"convert micro-Ct images to OpenFOAM simple parallel mesh"<<endl;
	cout<<"usage: example:"<<endl;
	cout<<"	   rawToFoamPar imageName.mhd 3 2 2 resetX0 unit"<<endl;
	return 1;
}

void fixImage(voxelImage& vimage);
void toFoam(voxelImage& vxlImg, const voxelField<int>& procIsijk, int iProc, int jProc, int kProc);

int main(int argc, char** argv)
{
	if(argc<5)		return usage();
	std::string headerName(argv[1]);
	if(headerName.size()<4 || headerName.compare(headerName.size()-4,4,".mhd") != 0) return usage();
	int nPar[3]={1,1,1};
	nPar[0] = atoi(argv[2]);
	nPar[1] = atoi(argv[3]);
	nPar[2] = atoi(argv[4]);
	char resetX0 = argc>5 ? argv[5][0] : 'F';
	char unit = argc>6 ? argv[6][0] : 'u';
	if (nPar[2]*nPar[1]*nPar[0]<1) {cout<<"\nError: nPar[2]*nPar[1]*nPar[0]<1\n"; return usage();}
	voxelImage vimage(headerName);
	int3 n = vimage.size3(); 


	vimage.printInfo();

	vimage.writeHeader("vxlImage.mhd");

	cout <<"finding connected parts of the image "<<endl;
	vimage.threshold101(0,0);

	vimage.crop(0,n[0]-1,0,n[1]-1,0,n[2]-1,1,2);	//		 XXXXXXXXXXXXXXXXXXXXXXXXXXXX
	fixImage(vimage);

	vimage.crop(1,n[0],1,n[1],1,n[2],1,2);	//		 XXXXXXXXXXXXXXXXXXXXXXXXXXXX

	cout <<"converting to OpenFOAM format "<<endl;
	vimage.printInfo();

	if(resetX0=='T' || resetX0=='t')		vimage.X0Ch()=-vimage.dx();
	if(unit=='u')	{ vimage.X0Ch()*=1.0e-6;  vimage.dxCh()*=1.0e-6; }


	voxelField<int> procIsijk(nPar[2]+2,nPar[1]+2,nPar[0]+2,-1);
	int iProc=-1;
	if (nPar[2]*nPar[1]*nPar[0]==1)
		toFoam(vimage, procIsijk, 1, 1, 1);
	else if (nPar[2]*nPar[1]*nPar[0]>1)
	{
		voxelField<voxelImage>  vimages(nPar[0],nPar[1],nPar[2],voxelImage());

		vector<int> iBs(nPar[0]+1,n[0]);
		vector<int> jBs(nPar[1]+1,n[1]);
		vector<int> kBs(nPar[2]+1,n[2]);
		for (int iz=0;iz<nPar[2];iz++)	kBs[iz]=int(n[2]/nPar[2])*iz;
		for (int iy=0;iy<nPar[1];iy++)	jBs[iy]=int(n[1]/nPar[1])*iy;
		for (int ix=0;ix<nPar[0];ix++)	iBs[ix]=int(n[0]/nPar[0])*ix;


		cout<<"iBs: "<<*iBs.begin()<<" ... "<<*iBs.rbegin()<<endl;
		cout<<"jBs: "<<*jBs.begin()<<" ... "<<*jBs.rbegin()<<endl;
		cout<<"kBs: "<<*kBs.begin()<<" ... "<<*kBs.rbegin()<<endl;

		 for (int iz=0;iz<nPar[2];iz++)
		  for (int iy=0;iy<nPar[1];iy++)
			for (int ix=0;ix<nPar[0];ix++)
			{
				vimages[iz][iy][ix].reset(iBs[ix+1]-iBs[ix]+2, jBs[iy+1]-jBs[iy]+2, kBs[iz+1]-kBs[iz]+2,0);
				vimages[iz][iy][ix].setFrom(vimage, iBs[ix], jBs[iy], kBs[iz]);
				if(vimages[iz][iy][ix].volFraction(0,0)>1.0e-12)		procIsijk[ix+1][iy+1][iz+1]=++iProc;
				cout<<"poro_"<<ix<<iy<<iz<<":"<<vimages[iz][iy][ix].volFraction(0,0)<<"  iBx"<<iBs[ix+1]-iBs[ix]+2<<endl;
			};

		 vimage.resize(0);
		 for (int iz=0;iz<nPar[2];iz++)
		  for (int iy=0;iy<nPar[1];iy++)
			for (int ix=0;ix<nPar[0];ix++)
			 if(procIsijk[ix+1][iy+1][iz+1]>=0)
			 {
				cout<<"************* processor: "<<ix<<" "<<iy<<" "<<iz<<", Phi="<<100*vimages[iz][iy][ix].volFraction(0,0)<<" *************"<<endl;
				toFoam(vimages[iz][iy][ix], procIsijk, ix+1, iy+1, iz+1);
				cout<<endl;
			 };;;;;
	} else cout<<"!!! npx x npy x npz = "<<nPar[2]*nPar[1]*nPar[0]<<endl;
	cout<<":/"<<endl;
   return 0;
}

void toFoam(voxelImage& vxlImg, const voxelField<int>& procIsijk, int iProc, int jProc, int kProc)
{
	int myprocI=procIsijk[iProc][jProc][kProc];
	int3 n=vxlImg.size3();n[0]-=2;n[1]-=2;n[2]-=2;
	vec3 X0=vxlImg.X0(); 
	vec3 dx=vxlImg.dx();
	X0+=dx;



	string Folder = myprocI>=0 ?  "processor"+toStr(myprocI)  :  ".";
	cout<<"procIs, size: "<<procIsijk.size()<<" "<<procIsijk[0].size()<<" "<<procIsijk[0][0].size()<<",  ijk: "<<iProc<<"  "<<jProc<<"  "<<kProc<<endl;
	cout<<Folder<<"  N: "<<n[0] <<" "<<n[1] <<" "<< n[2]<<endl;
	::mkdir(Folder.c_str(),0777);
	::mkdir((Folder+"/constant").c_str(),0777);
	Folder=Folder+"/constant/polyMesh";
	::mkdir((Folder).c_str(),0777);
	if (!::mkdir((Folder).c_str(),0777))
	{  cout << "Couldn't create directory "+Folder;  exit(0);  }

	ofstream pointsf((Folder+"/points").c_str());
	assert(pointsf);


	pointsf<<
	"FoamFile\n"
	"{\n"
	"	version	 2.0;\n"
	"	format	  ascii;\n"
	"	class	   vectorField;\n"
	"	location	\""<<Folder+"\";\n"
	"	object	  points;\n"
	"}\n\n";


//=======================================================================

	voxelField<int> point_mapper(n[0]+1,n[1]+1,n[2]+1,-1);



	int nPoints=0;

{
	std::vector<vector<vector<unsigned char> > >::iterator ddd=vxlImg.begin();//,d1p1=vxlImg.begin()+1;//		int iz=ppp-point_mapper.begin();
	register int iz=0;
	for (;ddd<vxlImg.end()-1;ddd++,iz++)
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
	std::vector<vector<vector<unsigned char> > >::iterator ddd=vxlImg.begin();//,d1p1=vxlImg.begin()+1;//		int iz=ppp-point_mapper.begin();
	register unsigned int iz=0;
	for (;ddd<vxlImg.end()-1;ddd++,iz++)
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



	size_t nCells=0;
	size_t nInternalFaces=0;
	size_t nGrainwallsFaces=0;
	size_t nBackFaces=0;
	size_t nFrontFaces=0;
	size_t nBottomFaces=0;
	size_t nTopFaces=0;
	size_t nLeftFaces=0;
	size_t nRightFaces=0;

#define  mcountFace(type) n##type##Faces++


	for (register int iz=1;iz<=n[2];iz++)
	 for (register int iy=1;iy<=n[1];iy++)
	  for (register int ix=1;ix<=n[0];ix++)
	  {
			if (!vxlImg[iz][iy][ix])
			{
				nCells++;
				if (iz!=1)
				{
				  if (vxlImg[iz-1][iy][ix])    mcountFace( Grainwalls);
				  //else
				}
				else if(vxlImg[iz-1][iy][ix]==1)  mcountFace( Grainwalls);
				else                           mcountFace( Back);

				if (iz!=n[2])
				{
				  if (vxlImg[iz+1][iy][ix])	 mcountFace( Grainwalls);
				  else                         mcountFace( Internal);
				}else if(vxlImg[iz+1][iy][ix]==1) mcountFace( Grainwalls);
				else                           mcountFace( Front);


				if (iy!=1)
				{
				  if (vxlImg[iz][iy-1][ix])    mcountFace( Grainwalls);
				  //else
				}
				else if(vxlImg[iz][iy-1][ix]==1)  mcountFace( Grainwalls);
				else                           mcountFace( Bottom);

				if (iy!=n[1])
				{
				  if (vxlImg[iz][iy+1][ix])	 mcountFace( Grainwalls);
				  else                         mcountFace( Internal);
				}else if(vxlImg[iz][iy+1][ix]==1) mcountFace( Grainwalls);
				else                           mcountFace( Top);


				if (ix!=1)
				{
				  if (vxlImg[iz][iy][ix-1])    mcountFace( Grainwalls);
				  //else
				}else if(vxlImg[iz][iy][ix-1]==1)  mcountFace( Grainwalls);
				else                            mcountFace( Left);

				if (ix!=n[0])
				{
				  if (vxlImg[iz][iy][ix+1])	  mcountFace( Grainwalls);
				  else                          mcountFace( Internal);
				}else if(vxlImg[iz][iy][ix+1]==1)  mcountFace( Grainwalls);
				else                            mcountFace( Right);

		}
	  }

	cout<<"nCells: "<<nCells<<endl;
	cout<<"nInternalFaces: "<<nInternalFaces<<endl;
	cout<<"nGrainwallsFaces: "<<nGrainwallsFaces<<endl;
	cout<<"nLeftFaces: "<<nLeftFaces<<endl;
	cout<<"nRightFaces: "<<nRightFaces<<endl;
	cout<<"nBottomFaces: "<<nBottomFaces<<endl;
	cout<<"nTopFaces: "<<nTopFaces<<endl;
	cout<<"nBackFaces: "<<nBackFaces<<endl;
	cout<<"nFrontFaces: "<<nFrontFaces<<endl;


	int iInternalStartFace=0;
	int iGrainwallsStartFace=0;

	int iLeftStartFace=0;
	int iRightStartFace=0;
	int iBottomStartFace=0;
	int iTopStartFace=0;
	int iBackStartFace=0;
	int iFrontStartFace=0;


	//cout<<"iInternalStartFace: "<<iInternalStartFace<<endl;
	//cout<<"iGrainwallsStartFace: "<<iGrainwallsStartFace<<endl;
	//cout<<"iLeftStartFace: "<<iLeftStartFace<<endl;
	//cout<<"iRightStartFace: "<<iRightStartFace<<endl;
	//cout<<"iBottomStartFace: "<<iBottomStartFace<<endl;
	//cout<<"iTopStartFace: "<<iTopStartFace<<endl;
	//cout<<"iBackStartFace: "<<iBackStartFace<<endl;
	//cout<<"iFrontStartFace: "<<iFrontStartFace<<endl;








	{	ofstream boundary((Folder+"/boundary").c_str());
		assert(boundary);
		boundary<<
		"FoamFile\n"
		"{\n"
		"	version	 2.0;\n"
		"	format	  ascii;\n"
		"	class	   polyBoundaryMesh;\n"
		"	location	\""<<Folder+"\";\n"
		"	object	  boundary;\n"
		"}\n\n";


		#define write_boundary(nwrite, name,type)          {          \
			i ##name ##StartFace = iLastFace;										\
			boundary<<		                                                \
			"	"<<  #name										 <<endl<<			   \
			"	{"												 <<endl<<		         \
			"		type			"<<#type<<";"				  <<endl<<			   \
			"		nFaces		  "<<(nwrite)* n ##name ##Faces<<';'<<endl<<	\
			"		startFace	   "<<iLastFace<<';'   <<endl<<  			 \
			"	}"												 <<endl; 					\
			iLastFace += (nwrite)* n ##name ##Faces;			}

		#define write_procBoundary(name,myId,neiId)     {              \
			i ##name ##StartFace = iLastFace;										\
			if (n ##name ##Faces >0)                                       \
			boundary<<		                                                \
			"	"<<  "processor"+	toStr(myId)+"to"+toStr(neiId)	<<endl<<	   \
			"	{"												 <<endl<<		         \
			"		type			processor;"				  <<endl<<			      \
			"		inGroups        1(processor);"	  <<endl<<			      \
			"		nFaces		  "<<n ##name ##Faces<<';'	   <<endl<<	      \
			"		startFace	   "<<iLastFace<<';'   <<endl<<   \
			"		matchTolerance  0.0001;"	  <<endl<<			            \
			"		transform       unknown;"	  <<endl<<			            \
			"		myProcNo        "<<myId<<";"	  <<endl<<			         \
			"		neighbProcNo    "<<neiId<<";"	  <<endl<<			         \
			"	}"												 <<endl;						\
			iLastFace += n ##name ##Faces;					  }
        
        
		boundary<<7+ 
						int(procIsijk[iProc-1][jProc][kProc]>=0 && nLeftFaces)+ 
						int(procIsijk[iProc+1][jProc][kProc]>=0 && nRightFaces)+
						 int(procIsijk[iProc][jProc-1][kProc]>=0 && nBottomFaces)+ 
						 int(procIsijk[iProc][jProc+1][kProc]>=0 && nTopFaces)+
						 int(procIsijk[iProc][jProc][kProc-1]>=0 && nBackFaces)+ 
						 int(procIsijk[iProc][jProc][kProc+1]>=0 && nFrontFaces)  <<endl
				<<'('<<endl;






	int iLastFace = nInternalFaces;

		write_boundary(true, Grainwalls,patch);
		write_boundary(procIsijk[iProc-1][jProc][kProc]<0, Left,patch);
		write_boundary(procIsijk[iProc+1][jProc][kProc]<0, Right,patch);
		write_boundary(procIsijk[iProc][jProc-1][kProc]<0, Bottom,patch);
		write_boundary(procIsijk[iProc][jProc+1][kProc]<0, Top,patch);
	 #ifdef _2D_
		write_boundary(true, Back,empty);
		write_boundary(true, Front,empty);
	 #else
		write_boundary(procIsijk[iProc][jProc][kProc-1]<0, Back,patch);
		write_boundary(procIsijk[iProc][jProc][kProc+1]<0, Front,patch);
	 #endif

	//if (procIsijk[iProc-1][jProc][kProc]>=0)  { iLeftStartFace=iLastFace; iLastFace+=nLeftFaces; }
	//if (procIsijk[iProc+1][jProc][kProc]>=0)  { iRightStartFace=iLastFace; iLastFace+=nRightFaces; }
	//if (procIsijk[iProc][jProc-1][kProc]>=0)  { iBottomStartFace=iLastFace; iLastFace+=nBottomFaces; }
	//if (procIsijk[iProc][jProc+1][kProc]>=0)  { iTopStartFace=iLastFace; iLastFace+=nTopFaces; }
	//if (procIsijk[iProc][jProc][kProc-1]>=0)  { iBackStartFace=iLastFace; iLastFace+=nBackFaces; }
	//if (procIsijk[iProc][jProc][kProc+1]>=0)  { iFrontStartFace=iLastFace; iLastFace+=nFrontFaces; }

		if (procIsijk[iProc-1][jProc][kProc]>=0) write_procBoundary(Left,myprocI, procIsijk[iProc-1][jProc][kProc]);
		if (procIsijk[iProc+1][jProc][kProc]>=0) write_procBoundary(Right,myprocI, procIsijk[iProc+1][jProc][kProc]);
		if (procIsijk[iProc][jProc-1][kProc]>=0) write_procBoundary(Bottom,myprocI, procIsijk[iProc][jProc-1][kProc]);
		if (procIsijk[iProc][jProc+1][kProc]>=0) write_procBoundary(Top,myprocI, procIsijk[iProc][jProc+1][kProc]);
		if (procIsijk[iProc][jProc][kProc-1]>=0) write_procBoundary(Back,myprocI, procIsijk[iProc][jProc][kProc-1]);
		if (procIsijk[iProc][jProc][kProc+1]>=0) write_procBoundary(Front,myprocI, procIsijk[iProc][jProc][kProc+1]);

		boundary<<")"   <<endl;
		boundary.close();
	}







	cout<<"creating faces"<<endl;


	std::vector<array<int,6> > faces_Internal(nInternalFaces, array<int,6>{{-1,-1,-1,-1,-1,-1}});
	std::vector<array<int,6> > faces_Grainwalls(nGrainwallsFaces, array<int,6>{{-1,-1,-1,-1,-1,-1}});
	std::vector<array<int,6> > faces_Back(nBackFaces, array<int,6>{{-1,-1,-1,-1,-1,-1}});
	std::vector<array<int,6> > faces_Front(nFrontFaces , array<int,6>{{-1,-1,-1,-1,-1,-1}});
	std::vector<array<int,6> > faces_Bottom(nBottomFaces,  array<int,6>{{-1,-1,-1,-1,-1,-1}});
	std::vector<array<int,6> > faces_Top(nTopFaces,   array<int,6>{{-1,-1,-1,-1,-1,-1}});
	std::vector<array<int,6> > faces_Left(nLeftFaces,   array<int,6>{{-1,-1,-1,-1,-1,-1}});
	std::vector<array<int,6> > faces_Right(nRightFaces,  array<int,6>{{-1,-1,-1,-1,-1,-1}});


	voxelField<array<int,3> > ownerMapper(n[0]+1,n[1]+1,n[2]+1,array<int,3>{{-1,-1,-1}});


	cout<<"collecting faces"<<endl;


#define recordF_m( l10,l11,l20,l21,l30,l31,dir,ii,jj,kk,type )		            		   \
  {																	                                   \
	if (ownerMapper[ii][jj][kk][dir]<0)										                       \
	{		   ++i##type ##Faces ;																		     \
			ownerMapper[ii][jj][kk][dir]=i##type ##Faces + i##type ##StartFace;	    	      \
			faces_##type[i##type ##Faces][4]=iCells;	/* cell number (for owners) */			 \
			faces_##type[i##type ##Faces][0]=point_mapper[ii][jj][kk];								   \
			faces_##type[i##type ##Faces][1]=point_mapper[ii+l11][jj+l21][kk+l31];					\
			faces_##type[i##type ##Faces][2]=point_mapper[ii+l10+l11][jj+l20+l21][kk+l30+l31];  \
			faces_##type[i##type ##Faces][3]=point_mapper[ii+l10][jj+l20][kk+l30];					\
	}																			   \
	else																			\
	{																			   \
			faces_##type[ownerMapper[ii][jj][kk][dir]][5]=iCells; /* cell number (for neighbours) */ \
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
	{  cout<<(iz%50 ? '.' : '\n');cout.flush();
		for (register int iy=1;iy<=n[1];iy++)
		 for (register int ix=1;ix<=n[0];ix++)
		  if (!vxlImg[iz][iy][ix])
		  {

				iCells++;

				//int iface=0;
				if (iz!=1)
				{
				  if (vxlImg[iz-1][iy][ix])    {kclockwiserecordF( Grainwalls)}
				  else                         {kclockwiserecordF( Internal);}
				}else if(vxlImg[iz-1][iy][ix]==1) {kclockwiserecordF( Grainwalls)}
				else                           {kclockwiserecordF( Back)}

				if (iz!=n[2])
				{
				  if (vxlImg[iz+1][iy][ix])    {kuclockwiserecordF( Grainwalls)}
				  else                         {kuclockwiserecordF( Internal);}
				}else if(vxlImg[iz+1][iy][ix]==1) {kuclockwiserecordF( Grainwalls)}
				else                           {kuclockwiserecordF( Front)}

				if (iy!=1)
				{
				  if (vxlImg[iz][iy-1][ix])    {jclockwiserecordF( Grainwalls)}
				  else                         {jclockwiserecordF( Internal);}
				}else if(vxlImg[iz][iy-1][ix]==1) {jclockwiserecordF( Grainwalls)}
				else                           {jclockwiserecordF( Bottom)}

				if (iy!=n[1])
				{
				  if (vxlImg[iz][iy+1][ix])    {juclockwiserecordF( Grainwalls)} 
				  else                         {juclockwiserecordF( Internal);}
				}else if(vxlImg[iz][iy+1][ix]==1) {juclockwiserecordF( Grainwalls)} 
				else                           {juclockwiserecordF( Top)}
				
				if (ix!=1)
				{
				  if (vxlImg[iz][iy][ix-1])    {iclockwiserecordF( Grainwalls)}
				  else                         {iclockwiserecordF( Internal);}
				}else if(vxlImg[iz][iy][ix-1]==1) {iclockwiserecordF( Grainwalls)}
				else                           {iclockwiserecordF( Left);}
				if (ix!=n[0])
				{
				  if (vxlImg[iz][iy][ix+1])    {iuclockwiserecordF( Grainwalls)}
				  else                         {iuclockwiserecordF( Internal);}
				}else if(vxlImg[iz][iy][ix+1]==1) {iuclockwiserecordF( Grainwalls)}
				else                           {iuclockwiserecordF( Right); }

		  }
	}




point_mapper.resize(0);




	ofstream faces((Folder+"/faces").c_str());
	assert(faces);
	faces<<
	"FoamFile\n"
	"{\n"
	"	version	 2.0;\n"
	"	format	  ascii;\n"
	"	class	   faceList;\n"
	"	location	\""<<Folder+"\";\n"
	"	object	  faces;\n"
	"}\n\n";


	ofstream owner((Folder+"/owner").c_str());
	assert(owner);
	owner<<
	"FoamFile\n"
	"{\n"
	"	version	 2.0;\n"
	"	format	  ascii;\n"
	"	class	   labelList;\n"
	"	location	\""<<Folder+"\";\n"
	"	object	  owner;\n"
	"}\n\n";
 
	ofstream neighbour((Folder+"/neighbour").c_str());
	assert(neighbour);
	neighbour<<
	"FoamFile\n"
	"{\n"
	"	version	 2.0;\n"
	"	format	  ascii;\n"
	"	class	   labelList;\n"
	"	location	\""<<Folder+"\";\n"
	"	object	  neighbour;\n"
	"}\n\n";

//______________________________________________________________________________________________________________________
#define write_faces_owners(type) \
  for (std::vector<array<int,6> >::iterator ff=faces_##type .begin();ff<faces_##type .end();ff++)	\
	{ \
		faces<<'('<<((*ff))[0]<<' '<<(*ff)[1]<<' '<<(*ff)[2]<<' '<<(*ff)[3]<<")\n"; \
		owner<<(*ff)[4]<<"\n";  \
	}

	int nFaces=nInternalFaces+nGrainwallsFaces+nLeftFaces+nRightFaces+nBottomFaces+nTopFaces+nBackFaces+nFrontFaces;

	owner<<nFaces<<endl;
	owner<<"("<<endl;
	faces<<nFaces<<endl
		 <<"("<<endl;

		//write_boundary(Grainwalls,patch);
		//if (procIsijk[iProc][jProc][kProc]<0) write_boundary(Left,patch);
		//if (procIsijk[iProc][jProc][kProc]<0) write_boundary(Right,patch);
		//if (procIsijk[iProc][jProc][kProc]<0) write_boundary(Bottom,patch);
		//if (procIsijk[iProc][jProc][kProc]<0) write_boundary(Top,patch);
	 //#ifdef _2D_
		//write_boundary(Back,empty);
		//write_boundary(Front,empty);
	 //#else
		//if (procIsijk[iProc][jProc][kProc]<0) write_boundary(Back,patch);
		//if (procIsijk[iProc][jProc][kProc]<0) write_boundary(Front,patch);
	 //#endif
		//int myprocI=procIsijk[iProc][jProc][kProc];
//
		//if (procIsijk[iProc-1][jProc][kProc]>=0) write_procBoundary(Left,myprocI, procIsijk[iProc-1][jProc][kProc]);
		//if (procIsijk[iProc+1][jProc][kProc]>=0) write_procBoundary(Right,myprocI, procIsijk[iProc+1][jProc][kProc]);
		//if (procIsijk[iProc][jProc-1][kProc]>=0) write_procBoundary(Bottom,myprocI, procIsijk[iProc][jProc-1][kProc]);
		//write_procBoundary(Top,myprocI, procIsijk[iProc][jProc+1][kProc]);
		//write_procBoundary(Back,myprocI, procIsijk[iProc][jProc][kProc-1]);
		//if (procIsijk[iProc][jProc][kProc+1]>=0) write_procBoundary(Front,myprocI, procIsijk[iProc][jProc][kProc+1]);



		write_faces_owners(Internal)
		write_faces_owners(Grainwalls)
		if (procIsijk[iProc-1][jProc][kProc]<0) write_faces_owners(Left)
		if (procIsijk[iProc+1][jProc][kProc]<0) write_faces_owners(Right)
		if (procIsijk[iProc][jProc-1][kProc]<0) write_faces_owners(Bottom)
		if (procIsijk[iProc][jProc+1][kProc]<0) write_faces_owners(Top)
		if (procIsijk[iProc][jProc][kProc-1]<0) write_faces_owners(Back)
		if (procIsijk[iProc][jProc][kProc+1]<0) write_faces_owners(Front)

		if (procIsijk[iProc-1][jProc][kProc]>=0) write_faces_owners(Left)
		if (procIsijk[iProc+1][jProc][kProc]>=0) write_faces_owners(Right)
		if (procIsijk[iProc][jProc-1][kProc]>=0) write_faces_owners(Bottom)
		if (procIsijk[iProc][jProc+1][kProc]>=0) write_faces_owners(Top)
		if (procIsijk[iProc][jProc][kProc-1]>=0) write_faces_owners(Back)
		if (procIsijk[iProc][jProc][kProc+1]>=0) write_faces_owners(Front)


	faces<<")"<<endl;
	faces.close();
	owner<<")"<<endl;
	owner.close();


	neighbour<<nInternalFaces<<endl;
	neighbour<<"("<<endl;
	  for (std::vector<array<int,6> >::iterator ff=faces_Internal.begin();ff<faces_Internal.end();ff++)
		{
			neighbour<<(*ff)[5]<<"\n";
		}

	neighbour<<")"<<endl;
	neighbour.close();

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
	
	cout<<"maxReg  "<<maxReg<<endl;
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

