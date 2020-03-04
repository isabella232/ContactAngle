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

    int nOrig[3];
    //trim
    double trimFactorBegin[3];
    double trimFactorEnd[3];
    int nResample=1;
    unsigned char  thresholdMin=0;
    unsigned char  thresholdMax=0;
    bool thresholdImage=false, resampleImage=false,cropImage=false,fillHoles=false;
    
    
    
    
    
    // open the file stream
    char tmpc;
    std::cout<<"Openning header file: "<<headerName<<std::endl;
    std::ifstream headerFile(headerName.c_str());                     // file pointer
    assert(headerFile);

    for ( unsigned int i=0; i<8;i++)   headerFile>>tmpc, std::cout<<" "<<tmpc;  //ignore the first 8 characters (ascii 3uc)

    headerFile>>n[2]>>n[0]>>n[1];                        // number of variables (dimension of
    std::cout<<", "<<n[2]<<", "<<n[1]<<", "<<n[0]<<std::endl;
    nOrig[0]=n[0],nOrig[1]=n[1],nOrig[2]=n[2];
    headerFile>>    xmin[0]>>xmax[0]>>xmin[1]>>xmax[1]>>xmin[2]>>xmax[2] ;
    std::cout<<  xmin[0]<<"-"<<xmax[0]<<", "<<xmin[1]<<"-"<<xmax[1]<<", "<<xmin[2]<<"-"<<xmax[2] <<std::endl;;
    dx[0]=(xmax[0]-xmin[0])/n[0]*1e-6; dx[1]=(xmax[1]-xmin[1])/n[1]*1e-6; dx[2]=(xmax[2]-xmin[2])/n[2]*1e-6;

    while (true)
    {
        std::string tmpStr;
        headerFile>>tmpStr;
        if (headerFile.fail()) break;
        
            std::cout<<tmpStr<<": "<<std::endl;;

        if (tmpStr.compare("crop") == 0)
        {
            cropImage=true;
            std::cout<<"  trimFactorBegin={";for ( unsigned int i=0; i<3;i++)   headerFile>>trimFactorBegin[i], std::cout<<trimFactorBegin[i]<<','; std::cout<<'}'<<std::endl;
            std::cout<<"  trimFactorEnd={";for ( unsigned int i=0; i<3;i++)   headerFile>>trimFactorEnd[i], std::cout<<trimFactorEnd[i]<<','; std::cout<<'}'<<std::endl;

            for ( unsigned int i=0; i<3;i++)
            {
                n[i]=nOrig[i]*trimFactorEnd[i]-nOrig[i]*trimFactorBegin[i];
                double dxi=(xmax[i]-xmin[i]);
                xmax[i]=xmin[i]+dxi*trimFactorEnd[i];
                xmin[i]=xmin[i]+dxi*trimFactorBegin[i];
            }
            std::cout<<"  n: "<<n[0]<<", "<<n[1]<<", "<<n[2]<<std::endl;

            std::cout<<"  x: "<<  xmin[0]<<" - "<<xmax[0]<<", "<<xmin[1]<<" - "<<xmax[1]<<", "<<xmin[2]<<" - "<<xmax[2] <<std::endl;;
        }
        else if (tmpStr.compare("resample") == 0)
        {
            resampleImage=true;
              headerFile>>nResample, std::cout<<"  resampling factor: "<<nResample<<std::endl;
        }
        else if (tmpStr.compare("pore") == 0)
        {
            thresholdImage=true;
            std::cout<<"  making binary inage:"<<std::endl
                     <<"  selecting pore (->0) with values between";
              int thresholdMinInt=0,thresholdMaxInt=0;
              headerFile>>thresholdMinInt;thresholdMin=thresholdMinInt;
              headerFile>>thresholdMaxInt;thresholdMax=thresholdMaxInt;
                
            std::cout<<" "<<(int)thresholdMin<<"  and "<<(int)thresholdMax<<"  inclusive."<<std::endl;
        }
        else if (tmpStr.compare("fillHoles") == 0)
        {
            fillHoles=true;
            std::cout<<"  filling small isolated parts:"<<std::endl;
        }
        else 
        {
            std::cout<<"  skipping the contents of file "<<headerName<<" after word \""<<tmpStr<<"\" "<<std::endl;
            break;
        }
         
        
    }
    headerFile.close();






       voxelImage=voxelMesh(nOrig[0],nOrig[1],nOrig[2],0);

      if (inputName.compare(inputName.size()-4,4,".raw") == 0)
      {
        voxelImage.readBin(inputName,0,nOrig[0]-1,0,nOrig[1]-1,0,nOrig[2]-1);
      }
      else
      {
        voxelImage.read(inputName);
      }



      if (cropImage)
      voxelImage.crop(0,n[0]-1,0,n[1]-1,0,n[2]-1,0,0);
      
      if (resampleImage)
      voxelImage.resample(nResample);
      
      if (thresholdImage)
      {
        voxelImage.threshold101(thresholdMin,thresholdMax);
      }
      
      if (fillHoles)
      {
            std::cout<<"  filling small isolated parts:"<<std::endl;

        voxelMesh dataTmp=voxelImage;
        dataTmp.shrinkSolid(); 
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        voxelImage=dataTmp;
        dataTmp.growSolid();
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        voxelImage=dataTmp;

        dataTmp.shrinkSolid(); 
        dataTmp.shrinkSolid();
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        dataTmp.growSolid();  dataTmp.OR(voxelImage);
        voxelImage=dataTmp;
        dataTmp.growSolid(); 
        dataTmp.growSolid(); 
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        dataTmp.shrinkSolid(); dataTmp.AND(voxelImage);
        voxelImage=dataTmp;
      }
      
      voxelImage.getSize(n[0],n[1],n[2]);
      
      
      
