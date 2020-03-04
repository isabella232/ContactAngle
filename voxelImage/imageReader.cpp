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

#include "voxelMesh.h"
#include <sstream>
//~ #include <streambuf>


namespace MCTProcessing
{
	
	bool ignore( std::stringstream & inputs, voxelMesh & voxelImage)
	{
			return true;
	}
	 
	bool fillHoles( std::stringstream & inputs, voxelMesh & voxelImage)
	{
		unsigned int maxHoleSize;
		inputs>>maxHoleSize;
		
			std::cout<<"fillHoles: eliminating isolated rocks/pores; maxHoleSize:" <<maxHoleSize<<" (default is 2) "<<std::endl;
			voxelImage.fillHoles(maxHoleSize);
			
			voxelImage.FaceMedian(2,4);
			//~ voxelImage.FaceMedian(2,4);
			//~ voxelImage.FaceMedian(2,4);
			return true;
	}
	 
	bool selectPore( std::stringstream & inputs, voxelMesh & voxelImage)
	{
			//~ thresholdImage=true;
			std::cout<<"  converting to binary (0 and 1):"<<std::endl
				 <<"  selecting pore (->0) with values between:";
			unsigned int  thresholdMin=0,thresholdMax=0;
			inputs>>thresholdMin;
			inputs>>thresholdMax;
	 
			std::cout<<" "<<(int)thresholdMin<<"  and "<<(int)thresholdMax<<"  inclusive."<<std::endl;
			voxelImage.threshold101(thresholdMin,thresholdMax);
			return true;
	}
	 
	bool growingThreshold( std::stringstream & inputs, voxelMesh & voxelImage)
	{
			//~ thresholdImage=true;
			std::cout<<"  converting to binary: 0 (pore) and 1 (rock):"<<std::endl
				 <<" start selecting pore with values between:";
			unsigned int  StartThreshold1=0,StartThreshold2=0, finalThreshold1=0,finalThreshold2=0;
			unsigned int nIter=4;
			
			inputs>>StartThreshold1;
			inputs>>StartThreshold2;
			inputs>>finalThreshold1;
			inputs>>finalThreshold2;
			inputs>>nIter;
			std::cout<<" "<<(int)StartThreshold1<<"  and "<<(int)StartThreshold2<<"  inclusive."<<std::endl;
			std::cout<<" growing to voxels  with values between:";
			std::cout<<" "<<(int)finalThreshold1<<"  and "<<(int)finalThreshold2<<", in "<<nIter<<" iterations."<<std::endl;
			voxelImage.segmentPhase
			(
				StartThreshold1,
				StartThreshold2, 
				finalThreshold1, 
				finalThreshold2,
				nIter
			);
			return true;
	}
	 
	bool growPore( std::stringstream & inputs, voxelMesh & voxelImage)
	{
			//~ thresholdImage=true;
			std::cout<<"  growing voxels:"<<std::endl;
			char voxelValueTogrow; inputs>>voxelValueTogrow;
			char growingAlgorithm; inputs>>growingAlgorithm;
			
			while (inputs.good())          // loop while extraction from file is possible
			{
				if (growingAlgorithm!='f')  
				{
					if(voxelValueTogrow=='0')
						voxelImage.growpore();
					else if (voxelValueTogrow=='1')
						voxelImage.shrinkPore();
					else
					{
						std::cerr<<"growing is only implemented for binary images: "<<
						"selected voxel value to grow is "<<voxelValueTogrow << ", which is not acceptable"<<std::endl;
						return false;//error occurred				
					}
				}
				else
				{
					std::cerr<<"selected growing algorithm: "<<growingAlgorithm<<
					" the only implemented algorithm is f which stands for faceGrowing"<<std::endl;
					return false;//error occurred
				}
					
				inputs>>voxelValueTogrow;
				inputs>>growingAlgorithm;
			}
			std::cout<<" done"<<std::endl;
			return true;
	}
	 
	 
	bool resample( std::stringstream & inputs, voxelMesh & voxelImage)
	{
		double nResample=1;
			inputs>>nResample, std::cout<<"  resampling factor: "<<nResample<<std::endl;
			voxelImage.resample(nResample); 
			return true;
	}
	 
	bool rotate( std::stringstream & inputs, voxelMesh & voxelImage)
	{
			char direction;
			inputs>>direction, std::cout<<" swapping x and "<<direction<<" directions"<<std::endl;
			
			voxelImage.rotate(direction);   
			return true;
	}
	 
	bool crop( std::stringstream & inputs, voxelMesh & voxelImage)
	{
		//~ unsigned int nOrig[3];
		unsigned int cropBegin[3], cropEnd[3];
		//~ voxelImage.getSize(nOrig[0],nOrig[1],nOrig[2]);
	 
	 
				//~ std::cout<<"cropping: out of  i=0 to "<<nOrig[0]-1<<",  j=0 to "<<nOrig[1]-1<<", k=0 to"<<nOrig[2]-1<<" "<<std::endl;
				std::cout<<"Crop:   ";
				 for ( unsigned int i=0; i<3;i++)   inputs>>cropBegin[i] >>cropEnd[i],  std::cout<<cropBegin[i]<<' '<<cropEnd[i]<<"    "; 
				  std::cout<<' '<<std::endl;
				//~ std::cout<<"  cropEnd={ ";for ( unsigned int i=0; i<3;i++)  std::cout; std::cout<<'}'<<std::endl;
	 
				for ( unsigned int i=0; i<3;i++)
				{
					//~ n[i]=cropEnd[i]-cropBegin[i];
					//~ double dxi=(xmax[i]-xmin[i])/nOrig[i];
					//~ xmax[i]=xmin[i]+dxi*cropEnd[i];
					//~ xmin[i]=xmin[i]+dx[i]*cropBegin[i];
				}
					//~ std::cout<<"Image size after cropping:"<<std::endl;
					//~ std::cout<<"N:   "<<n[0]<<"    "<<n[1]<<"    "<<n[2]<<std::endl;
					//~ std::cout<<"X0:  "<<  xmin[0]<<"  "<<xmin[1]<<"  "<<xmin[2]<<"  " <<std::endl;
	 
			voxelImage.crop(cropBegin,cropEnd);
			return true;
	}
 
 
 
}
 
 
 
 
bool readHeader
(
	std::ifstream& headerFile,
	unsigned int n[3],
	double  xmin[3],
	double  dx[3]
)
{
    char tmpc;
    for ( unsigned int i=0; i<8;i++)   headerFile>>tmpc, std::cout<<" "<<tmpc;  //ignore the first 8 characters (ascii 3uc)
 
    headerFile>>n[0]>>n[1]>>n[2];                        // number of variables (dimension of
    //~ headerFile>>n[2]>>n[0]>>n[1];                        // original Format, abondoned to avoid confusions it causes
    std::cout<<"Note format conventions changed to Nx(=number of columns) Ny(=number of rows) Nz(=number of layers),  which are as follows in this image"<<std::endl;
    std::cout<<" n          :  "<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl;
    //~ nOrig[0]=n[0],nOrig[1]=n[1],nOrig[2]=n[2];
    headerFile>>    dx[0]>>dx[1]>>dx[2] ;
    headerFile>>    xmin[0]>>xmin[1]>>xmin[2] ;
    std::cout<<" dx  dy  dz :  "<< dx[0]<<"  "<<dx[1]<<"  "<<dx[2]<<std::endl;
    std::cout<<" X0  Y0  Z0 :  "<<  xmin[0]<<"  "<<xmin[1]<<"  "<<xmin[2]<<"  " <<std::endl;;
}


  
  
void	readFromHeader
(

	voxelMesh & voxelImage,
	std::string headerName,
	std::string inputName
)
{
 	unsigned int n[3];
 	double  xmin[3];
	double  dx[3];
	
    // open the file stream
    std::cout<<"Openning header file: "<<headerName<<std::endl;
    std::ifstream headerFile(headerName.c_str());                     
    if(!headerFile)  {std::cout<<"\n\nError: can not open header file, "<<headerName<<std::endl<<std::endl; exit(-1);}



	readHeader(headerFile,n,xmin,dx);


	 voxelImage=voxelMesh(n, dx, xmin,0);		


	if( !inputName.empty() )
	{

      if (inputName.compare(inputName.size()-4,4,".raw") == 0)
      {
        voxelImage.readBin(inputName,0,n[0]-1,0,n[1]-1,0,n[2]-1);
      }
      else
      {
        voxelImage.read(inputName);
      }

	}



	typedef bool(*ProcessP)( std::stringstream&  inputs, voxelMesh & voxelImage);

	std::vector< ProcessP >  Processes;
	std::vector<std::string> processNames;

	processNames.push_back("");  Processes.push_back(& MCTProcessing::ignore);
	processNames.push_back("fillHoles");  Processes.push_back(& MCTProcessing::fillHoles);
	processNames.push_back("pore");  Processes.push_back(& MCTProcessing::selectPore);
	processNames.push_back("resample");  Processes.push_back(& MCTProcessing::resample);
	processNames.push_back("direction");  Processes.push_back(& MCTProcessing::rotate);
	processNames.push_back("crop");  Processes.push_back(& MCTProcessing::crop);
	processNames.push_back("growingThreshold");  Processes.push_back(& MCTProcessing::growingThreshold);
	processNames.push_back("growPore");  Processes.push_back(& MCTProcessing::growPore);

    while (true)
    {
		
		std::string tmpStr;
		headerFile>>tmpStr;	
		std::cout<<tmpStr<<": "<<std::endl;
		bool validKey=false;
		for (unsigned int i=0; i<processNames.size();i++)
		{
			if (headerFile.fail()) break;

			if (tmpStr.compare(processNames[i]) == 0)
			{
				std::stringstream keywordData;
				headerFile.get (*(keywordData.rdbuf()));
				(*Processes[i])(keywordData,voxelImage);
				validKey=true;  break;
			}
		}

		if ( !validKey )
        {
            std::cout<<" skipping the contents of file "<<headerName<<" after entry \""<<tmpStr<<"\" "<<std::endl;
            break;
        }
          
    }
    headerFile.close();

    //~ voxelImage.getSize(n[0],n[1],n[2]);


	//~ if( inputName.empty() ) voxelImage.resize(0); //avoid coding mistakes and save memory
   
      
}
      
