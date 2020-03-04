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

//~ #define _2D_

#include <sys/stat.h>

#include <fstream>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <iostream>

#include "voxelMesh.h"

using namespace std;


int usage()
{
	cout<<"convert micro-Ct images to OpenFOAM mesh"<<endl;
	cout<<"usage: example:"<<endl;
	cout<<"       rawToFoam imageName.dat_header imageName.dat"<<endl;
	cout<<"       or"<<endl;
	cout<<"       rawToFoam imageName.raw_header imageName.raw"<<endl;
}
int main(int argc, char** argv)
{
  try
    {

    if(argc<3)		return usage();

    //~ char tmpc;
    //~ unsigned int n[3],nStart[3],nEnd[3],nOrig[3],nOrigp2[3],tempi;
    //~ double  xmin[3],xmax[3];



    //~ unsigned int n[3];
    //~ double  xmin[3];
    //~ double dx[3];    
    voxelMesh voxelImage;
	std::string headerName=std::string(argv[1]);
	std::string imageName=std::string(argv[2]);
	cout <<"reading data, header: "<<headerName<<" and image: "<<imageName<<endl;
	readFromHeader(voxelImage, headerName, imageName);
    unsigned int n[3];	voxelImage.getSize(n[0],n[1],n[2]);
    Double3D xmin=voxelImage.X0();
    Double3D dx=voxelImage.dx();

	dx[0]*=1.0e-6; dx[1]*=1.0e-6; dx[2]*=1.0e-6; 
	xmin[0]*=1.0e-6; xmin[1]*=1.0e-6; xmin[2]*=1.0e-6; 


	voxelImage.printInfo();
	
	voxelImage.writeHeader("voxelMesh.dat");
	
	
	
	
	
	cout <<"converting to OpenFOAM format "<<endl;
	
	
	
	
    voxelImage.crop(0,n[0]-1,0,n[1]-1,0,n[2]-1,1,1);	//         XXXXXXXXXXXXXXXXXXXXXXXXXXXX



	
    //~ voxelImage.readBin(std::string(argv[2]),1,nOrig[0],1,nOrig[1],1,nOrig[2]);
    //~ voxelImage.crop( nStart[0]+1,nEnd[0]+1,     nStart[1]+1,nEnd[1]+1,      nStart[2]+1,nEnd[2]+1  , 1, 0  );










    //std::string outdir;
    if (!::mkdir("constant",0777))
    {
        cout  << "Couldn't create 'constant' directory ";
        exit(0);
    //outdir="constant/polyMesh/";
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
"    version     2.0;\n"
"    format      ascii;\n"
"    class       vectorField;\n"
"    location    \"constant/polyMesh\";\n"
"    object      points;\n"
"}\n\n";


//=======================================================================

    voxelField<int> point_mapper(n[0]+1,n[1]+1,n[2]+1,-1);



    int nPoints=0;

{
    std::vector<vector<vector<unsigned char> > >::iterator ddd=voxelImage.begin();//,d1p1=voxelImage.begin()+1;//        int iz=ppp-point_mapper.begin();
    register int iz=0;
    for (;ddd<voxelImage.end()-1;ddd++,iz++)
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
    std::vector<vector<vector<unsigned char> > >::iterator ddd=voxelImage.begin();//,d1p1=voxelImage.begin()+1;//        int iz=ppp-point_mapper.begin();
    register unsigned int iz=0;
    for (;ddd<voxelImage.end()-1;ddd++,iz++)
    {
        register unsigned int iy=0;
        double z=iz*dx[2];
        std::vector<vector<unsigned char> >::iterator dd1=ddd->begin(), dd2=(ddd+1)->begin();
        for (;dd1<ddd->end()-1;dd1++,dd2++,iy++)
        {
            register unsigned int ix=0;
            double y=iy*dx[1];
            std::vector<unsigned char>::iterator d11=dd1->begin(), d12=(dd1+1)->begin();
            std::vector<unsigned char>::iterator d21=dd2->begin(), d22=(dd2+1)->begin();
            for (;d11<dd1->end()-1;d11++,d12++,d21++,d22++,ix++)
            {
                if (*d11==0 || *(d11+1)==0 ||
                    *d12==0 || *(d12+1)==0 ||
                    *d21==0 || *(d21+1)==0 ||
                    *d22==0 || *(d22+1)==0 )
                {
                    double x=ix*dx[0];
                    pointsf<< "("<<x<< ' '<<y<<' '<<z<<")\n";
                    //cout<< "("<<iz<< ' '<<iy<<' '<<ix<<')'<<endl;
                    point_mapper[iz][iy][ix]=iPoints;iPoints++;
                }
            }
        }
    }
    pointsf<<endl<<")"<<endl;

    pointsf.close();

    cout <<"   done"<<endl;    cout.flush();
//=======================================================================



//______________________________________________________________________________________________________________________

    ofstream faces("constant/polyMesh/faces");
    assert(faces);
    faces<<
    "FoamFile\n"
    "{\n"
    "    version     2.0;\n"
    "    format      ascii;\n"
    "    class       faceList;\n"
    "    location    \"constant/polyMesh\";\n"
    "    object      faces;\n"
    "}\n\n";

    ofstream cells("constant/polyMesh/cells");
    assert(cells);
    cells<<
    "FoamFile\n"
    "{\n"
    "    version     2.0;\n"
    "    format      ascii;\n"
    "    class       vectorField;\n"
    "    location    \"constant/polyMesh\";\n"
    "    object      cells;\n"
    "}\n\n";


    ofstream owner("constant/polyMesh/owner");
    assert(owner);
    owner<<
    "FoamFile\n"
    "{\n"
    "    version     2.0;\n"
    "    format      ascii;\n"
    "    class       labelList;\n"
    "    location    \"constant/polyMesh\";\n"
    "    object      owner;\n"
    "}\n\n";

    ofstream neighbour("constant/polyMesh/neighbour");
    assert(neighbour);
    neighbour<<
    "FoamFile\n"
    "{\n"
    "    version     2.0;\n"
    "    format      ascii;\n"
    "    class       labelList;\n"
    "    location    \"constant/polyMesh\";\n"
    "    object      neighbour;\n"
    "}\n\n";


    ofstream boundary("constant/polyMesh/boundary");
    assert(boundary);
    boundary<<
    "FoamFile\n"
    "{\n"
    "    version     2.0;\n"
    "    format      ascii;\n"
    "    class       polyBoundaryMesh;\n"
    "    location    \"constant/polyMesh\";\n"
    "    object      boundary;\n"
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


#define count_cells(ix,iy,iz)            nCells++;                        \
    if (iz!=1)                                                            \
    {                                                                     \
        if (!voxelImage[iz-1][iy][ix])    {}                              \
        else                    mcountFace( Grainwalls)                   \
    }else mcountFace( Back);                                              \
    if (iz!=n[2])                                                          \
    {                                                                      \
        if (!voxelImage[iz+1][iy][ix])    mcountFace( Internal)            \
        else                    mcountFace( Grainwalls)                    \
    }else mcountFace( Front);                                               \
                                                                            \
    if (iy!=1)                                                              \
    {                                                                       \
        if (!voxelImage[iz][iy-1][ix])    {}                               \
        else                    mcountFace( Grainwalls)                    \
    }else mcountFace( Bottom);                                             \
    if (iy!=n[1])                                                           \
    {                                                                       \
        if (!voxelImage[iz][iy+1][ix])    mcountFace( Internal)             \
        else                    mcountFace( Grainwalls)                    \
    }else mcountFace( Top);                                                \
                                                                            \
    if (ix!=1)                                                              \
    {                                                                       \
        if (!voxelImage[iz][iy][ix-1])    {}                               \
        else                    mcountFace( Grainwalls)                    \
    }else mcountFace( Left);                                               \
    if (ix!=n[0])                                                           \
    {                                                                       \
        if (!voxelImage[iz][iy][ix+1])    mcountFace( Internal)            \
        else                    mcountFace( Grainwalls)                    \
    }else mcountFace( Right);



    for (register int iz=1;iz<=n[2];iz++)
    {
        for (register int iy=1;iy<=n[1];iy++)

        {
            for (register int ix=1;ix<=n[0];ix++)
            {
                if (!voxelImage[iz][iy][ix])
                {
                    count_cells(ix,iy,iz);
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
#define write_boundary(name,type)\
boundary<<          \
"    "<<  #name                                         <<endl<<            \
"    {"                                                 <<endl<<            \
"        type            "<<#type<<";"                  <<endl<<            \
"        nFaces          "<<n ##name ##Faces<<';'       <<endl<<            \
"        startFace       "<<i ##name ##StartFace<<';'   <<endl<<            \
"    }"                                                 <<endl;

    boundary<<7  <<endl
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

//std::vector<short> neighbour0_5((n[2]+1)*(n[1]+1)*(n[0]+1)*3, -1);

    //std::vector<vector<vector<int> > > face_mapper1(n[0]+1);
    //std::vector<vector<vector<int> > > face_mapper2(n[0]+1);
//  std::vector<int> owner0((n[2]+1)*(n[1]+1)*(n[0]+1));
//  std::vector<int> owner1((n[2]+1)*(n[1]+1)*(n[0]+1));
//  std::vector<int> owner2((n[2]+1)*(n[1]+1)*(n[0]+1));
//  std::vector<int> neighbour0((n[2]+1)*(n[1]+1)*(n[0]+1));
//  std::vector<int> neighbour1((n[2]+1)*(n[1]+1)*(n[0]+1));
//  std::vector<int> neighbour2((n[2]+1)*(n[1]+1)*(n[0]+1));


//#define init_face_mapper(i)
    //~ std::vector<vector<vector<vector<int> > > > ownerMapper(n[2]+1);
    //~ for (std::vector<vector<vector<vector<int> > > >::iterator ppp=ownerMapper.begin();ppp<ownerMapper.end();ppp++)
    //~ {
        //~ ppp->resize(n[1]+1);
        //~ for (std::vector<vector<vector<int> > >::iterator pp=ppp->begin();pp<ppp->end();pp++)
            //~ pp->resize(n[0]+1, vector<int>(3,-1)) ;
    //~ }
    voxelField<vector<int> > ownerMapper(n[0]+1,n[1]+1,n[2]+1,vector<int>(3,-1));


    cout<<"collecting faces"<<endl;

//~ char aa;
//~ cin >> aa;

//init_face_mapper(0)
//init_face_mapper(1)
//init_face_mapper(2)

//for (register int i=0; i<4;i++)

#define recordFaces_m( l10,l11,l20,l21,l30,l31,dir,ii,jj,kk,type )                   \
  {                                                                     \
    if (ownerMapper[ii][jj][kk][dir]<0)                                         \
    {           ++i##type ##Faces ;                                                                         \
            ownerMapper[ii][jj][kk][dir]=i##type ##Faces + i##type ##StartFace;                             \
            faces_##type[i##type ##Faces][4]=iCells;    /* cell number (for owners) */                                \
            faces_##type[i##type ##Faces][0]=point_mapper[ii][jj][kk];                                  \
            faces_##type[i##type ##Faces][1]=point_mapper[ii+l11][jj+l21][kk+l31];                      \
            faces_##type[i##type ##Faces][2]=point_mapper[ii+l10+l11][jj+l20+l21][kk+l30+l31];          \
            faces_##type[i##type ##Faces][3]=point_mapper[ii+l10][jj+l20][kk+l30];                      \
    }                                                                               \
    else                                                                            \
    {                                                                               \
            faces_##type[ownerMapper[ii][jj][kk][dir]][5]=iCells;        /* cell number (for neighbours) */                            \
    }                                                                               \
    cell[iface++]=i##type ##Faces + i##type ##StartFace;                                                                        \
                                    \
                                    \
}
// cout<<' '<<i##type ##Faces<<'+'<< i##type ##StartFace;
// cout<<'.'<<endl;

#define  iclockwiserecordFaces(type) recordFaces_m( 0,1,1,0,0,0, 2, iz-1,iy-1,ix-1, type)
#define iuclockwiserecordFaces(type) recordFaces_m( 1,0,0,1,0,0, 2, iz-1,iy-1,ix  , type)
#define  jclockwiserecordFaces(type) recordFaces_m( 1,0,0,0,0,1, 1, iz-1,iy-1,ix-1, type)
#define juclockwiserecordFaces(type) recordFaces_m( 0,1,0,0,1,0, 1, iz-1,iy  ,ix-1, type)
#define  kclockwiserecordFaces(type) recordFaces_m( 0,0,0,1,1,0, 0, iz-1,iy-1,ix-1, type)
#define kuclockwiserecordFaces(type) recordFaces_m( 0,0,1,0,0,1, 0, iz  ,iy-1,ix-1, type)


#define write_faces_m       iCells++;                                   \
  {                                                                         \
  int cell[6]; bool internal[]={false,false,false,false,false,false}; int iface=0;      \
    if (ix!=1)                                                              \
    {                                                                       \
        if (!voxelImage[iz][iy][ix-1])    {iclockwiserecordFaces( Internal);internal[iface]=true;}            \
        else                    {iclockwiserecordFaces( Grainwalls)}            \
    }else {iclockwiserecordFaces( Left);}                                       \
    if (ix!=n[0])                                                           \
    {                                                                       \
        if (!voxelImage[iz][iy][ix+1])    {iuclockwiserecordFaces( Internal);internal[iface]=true;}           \
        else                    {iuclockwiserecordFaces( Grainwalls)}           \
    }else {iuclockwiserecordFaces( Right); }                                    \
                                                                            \
    if (iy!=1)                                                              \
    {                                                                       \
        if (!voxelImage[iz][iy-1][ix])    {jclockwiserecordFaces( Internal);internal[iface]=true;}            \
        else                    {jclockwiserecordFaces( Grainwalls)}            \
    }else {jclockwiserecordFaces( Bottom)}                                  \
    if (iy!=n[1])                                                           \
    {                                                                       \
        if (!voxelImage[iz][iy+1][ix])    {juclockwiserecordFaces( Internal);internal[iface]=true;}\
        else                    {juclockwiserecordFaces( Grainwalls)}           \
    }else {juclockwiserecordFaces( Top)}                                    \
                                                                            \
    if (iz!=1)                                                              \
    {                                                                       \
        if (!voxelImage[iz-1][iy][ix])    {kclockwiserecordFaces( Internal);internal[iface]=true;}            \
        else                    {kclockwiserecordFaces( Grainwalls)}            \
    }else {kclockwiserecordFaces( Back)}                                    \
    if (iz!=n[2])                                                           \
    {                                                                       \
        if (!voxelImage[iz+1][iy][ix])    {kuclockwiserecordFaces( Internal);internal[iface]=true;}           \
        else                    {kuclockwiserecordFaces( Grainwalls)}           \
    }else {kuclockwiserecordFaces( Front)}                                      \
    cells<< "(" ;                                                               \
    for (int ii=0;ii<6;ii++)                                                    \
        if  (internal[ii])  cells<< cell[ii]<<' ';                              \
    for (int ii=0;ii<6;ii++)                                                    \
        if  (!internal[ii]) cells<< cell[ii]<<' ';                              \
    cells<< ")\n";                                                              \
    }




    int iCells=-1;

    int iInternalFaces=-1;
    int iGrainwallsFaces=-1;
    int iBackFaces=-1;
    int iFrontFaces=-1;
    int iBottomFaces=-1;
    int iTopFaces=-1;
    int iLeftFaces=-1;
    int iRightFaces=-1;

    for (register unsigned int iz=1;iz<=n[2];iz++)
    {cout<<'.';cout.flush();
        for (register unsigned int iy=1;iy<=n[1];iy++)

        {
            for (register unsigned int ix=1;ix<=n[0];ix++)
            {
                if (!voxelImage[iz][iy][ix])
                {
                    write_faces_m;
                }
            }
        }
    }

    cells<<")"<<endl;

    cells.close();
//______________________________________________

point_mapper.resize(0);



//______________________________________________________________________________________________________________________
#define write_faces_owners(type) \
  for (std::vector<vector<int> >::iterator ff=faces_##type .begin();ff<faces_##type .end();ff++)    \
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














    }
  catch (std::exception &exc)
    {
      std::cerr << endl << endl
                << "----------------------------------------------------"<< endl;
      std::cerr << "Exception on processing: " << endl
                << exc.what() << endl
                << "Aborting!" << endl
                << "----------------------------------------------------"<< endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << endl << endl
                << "----------------------------------------------------"<< endl;
      std::cerr << "Unknown exception!" << endl
                << "Aborting!" << endl
                << "----------------------------------------------------"<< endl;
      return 1;
    }












    return 0;
}
