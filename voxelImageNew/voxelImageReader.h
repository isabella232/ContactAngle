/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.


This file is part of voxelImage library, a small c++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:    a.qaseminejad-raeini09@imperial.ac.uk

\*-------------------------------------------------------------------------*/

#include <sstream>
//~ #include <streambuf>


namespace MCTProcessing
{
template<typename T> bool ignore( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
		return true;
}

template<typename T> bool fillHoles( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
	unsigned int maxHoleSize;
	inputs>>maxHoleSize;

	    std::cout<<"fillHoles: eliminating isolated rocks/pores; maxHoleSize:" <<maxHoleSize<<" (default is 2) "<<std::endl;
        vxlImage.fillHoles(maxHoleSize);

        vxlImage.FaceMedian06(1,5);
        //~ vxlImage.FaceMedian07(2,5);
        //~ vxlImage.FaceMedian07(2,5);
		return true;
}

template<typename T> bool selectPore( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
		std::cout<<"  converting to binary (0 and 1):"<<std::endl
			 <<"  selecting pore (->0) with values between:";
		unsigned int  thresholdMin=0,thresholdMax=0;
		inputs>>thresholdMin;
		inputs>>thresholdMax;

		std::cout<<" "<<int(thresholdMin)<<"  and "<<int(thresholdMax)<<"  inclusive."<<std::endl;
		vxlImage.threshold101(thresholdMin,thresholdMax);
		return true;
}

template<typename T> bool growPore( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
		std::cout<<"  growing voxels:"<<std::endl;
		int voxelValueTogrow; inputs>>voxelValueTogrow;
 		char growingAlgorithm; inputs>>growingAlgorithm;

		while (inputs.good())          // loop while extraction from file is possible
		{
			if (growingAlgorithm!='f')
			{
				if(voxelValueTogrow==0)
					vxlImage.growPore();
				else if (voxelValueTogrow==0)
					vxlImage.shrinkPore();
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


template<typename T> bool resample( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
    double nResample=1;
        inputs>>nResample, std::cout<<"  resampling factor: "<<nResample<<std::endl;
		vxlImage.resample(nResample);
		return true;
}


template<typename T> bool resampleMax( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
    double nResample=1;
        inputs>>nResample, std::cout<<"  resampling factor: "<<nResample<<std::endl;
		vxlImage.resampleMax(nResample);
		return true;
}

template<typename T> bool redirect( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
		char direction;
        inputs>>direction, std::cout<<direction<<", swapping x and "<<direction<<" directions"<<std::endl;

		vxlImage.rotate(direction);
		return true;
}

template<typename T> bool crop( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
    unsigned int cropBegin[3], cropEnd[3];

	std::cout<<"Crop:   ";
	for ( unsigned int i=0; i<3;i++)   inputs>>cropBegin[i] >>cropEnd[i],  std::cout<<cropBegin[i]<<' '<<cropEnd[i]<<"    ";
	std::cout<<' '<<std::endl;

	//cropEnd[0]+=1; cropEnd[1]+=1; cropEnd[2]+=1;
	vxlImage.crop(cropBegin,cropEnd);
	return true;
}


template<typename T> bool cropD( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
	int3 cropBegin{{0,0,0}}, cropEnd=vxlImage.size3();
	int nLayers(0); int value(1);
	std::cout<<"cropD:   ";
	inputs>>cropBegin[0] >>cropBegin[1] >>cropBegin[2];  std::cout<<" "<<cropBegin[0] <<" "<<cropBegin[1] <<" "<<cropBegin[2]<<" --  ";  
	inputs>>cropEnd[0] >>cropEnd[1] >>cropEnd[2];        std::cout<<" "<<cropEnd[0] <<" "<<cropEnd[1] <<" "<<cropEnd[2]<<"  +  ";;  
	inputs >> nLayers >> value;
	std::cout<<nLayers<<" layers with "<<value<<std::endl;
	vxlImage.cropD(cropBegin,cropEnd,nLayers,value);
	return true;
}


}




template<typename T>
void voxelImageT<T>::readFromHeader
(
	std::ifstream& headerFile,
	std::string header,
	int processKeys,
	std::string inputName
)
{
	unsigned int n[3];
	std::string BinaryData="XXX";
	std::cout<<header<<":"<<std::endl;
    if (header.size()>4 && header.compare(header.size()-4,4,".mhd") == 0)
    {
//		BinaryData="True";
		while (true)
		{
			std::string tmpStr;
			std::streampos begLine = headerFile.tellg();
			headerFile>>tmpStr;

			std::cout<<" "<<tmpStr<<": "; std::cout.flush();

			if (headerFile.fail()) break;
			//~ ObjectType = Image
			//~ NDims = 3
			//~ Offset = 0 0 0
			//~ ElementSpacing = 8 8 8
			//~ DimSize = 200 225 153
			//~ ElementType = MET_UCHAR
			//~ ElementDataFile = Ketton100.raw
			std::stringstream keywordData;
			headerFile.get (*(keywordData.rdbuf()));
			std::string tmp;
			if (tmpStr == "ObjectType")
			{
				keywordData >> tmp; keywordData >> tmp;
				if (tmp != "Image") std::cout<<"Warning: ObjectType != Image :="<<tmp;
			}
			else if (tmpStr == "NDims")
			{
				keywordData >> tmp; keywordData >> tmp;
				if (tmp != "3") std::cout<<"Warning: NDims != 3 :="<<tmp;
			}
			else if (tmpStr == "ElementType")
			{
				keywordData >> tmp; keywordData >> tmp;
				if (tmp != "MET_UCHAR") std::cout<<"Warning: ElementType != MET_UCHAR :="<<tmp;
			}
			else if (tmpStr == "Offset")
			{
				keywordData >> tmp; keywordData>>    X0_[0]>>X0_[1]>>X0_[2] ;
				std::cout<<" X0: "<<  X0_[0]<<"  "<<X0_[1]<<"   "<<X0_[2] ;
			}
			else if (tmpStr == "ElementSpacing")
			{
				keywordData >> tmp; keywordData>>    dx_[0]>>dx_[1]>>dx_[2] ;
				std::cout<<" dX: "<< dx_[0]<<"  "<<dx_[1]<<"  "<<dx_[2]<<"   "; std::cout.flush();
			}
			else if (tmpStr == "DimSize")
			{
				keywordData >> tmp; keywordData>>    n[0]>>n[1]>>n[2];
				std::cout<<" Nxyz: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"   "; std::cout.flush();
			}
			else if (tmpStr == "ElementDataFile")
			{
				keywordData >> tmp; if (inputName.empty()) keywordData >> inputName;

				size_t islash=header.find_last_of("\\/");
				if (islash<header.size() && inputName[0]!='/' &&  inputName[1]!=':') inputName=header.substr(0,islash+1)+inputName;
				std::cout<<" ElementDataFile = "<<inputName<<"    ";
			}
			else if (tmpStr == "BinaryData")
			{
				keywordData >> tmp; keywordData >> BinaryData;
				std::cout<<" BinaryData = "<<BinaryData<<"    ";
			}
			else if (tmpStr!="BinaryDataByteOrderMSB" && tmpStr!="CompressedData" &&  tmpStr!="CompressedDataSize" &&  tmpStr!="TransformMatrix" &&
					 tmpStr!="CenterOfRotation" && tmpStr!="AnatomicalOrientation" && tmpStr!="AnatomicalOrientation")
			{
				headerFile.seekg(begLine);
				break;
			}
			std::cout<<std::endl;
		}
	}
	else
	{

		char tmpc;
		for ( unsigned int i=0; i<8;i++)   headerFile>>tmpc, std::cout<<tmpc;  //ignore the first 8 characters (ascii 3uc)

		if (header.size()>7 && header.compare(header.size()-7,7,"_header") == 0)  inputName=header.substr(0,header.size()-7);
		headerFile>>n[0]>>n[1]>>n[2];                        // number of variables (dimension of
		std::cout<<"\n Nxyz: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"   "; std::cout.flush();
		headerFile>>    dx_[0]>>dx_[1]>>dx_[2] ;
		std::cout<<" dX: "<< dx_[0]<<"  "<<dx_[1]<<"  "<<dx_[2]<<"   "; std::cout.flush();
		headerFile>>    X0_[0]>>X0_[1]>>X0_[2] ;
		if (!headerFile)     { std::cout<<"  Incomplete/bad header, aborting"<<std::endl; exit(-1);}
		std::cout<<" X0: "<<  X0_[0]<<"  "<<X0_[1]<<"   "<<X0_[2] << std::endl;
		if (!headerFile)     { std::cout<<"  Incomplete/bad header, continuing anyway"<<std::endl; }

	}

	this->reset(n[0],n[1],n[2],1);
	if( !inputName.empty() && inputName!="NO_READ" && processKeys!=2 )
	{
	  if ((inputName.compare(inputName.size()-4,4,".raw") == 0 && BinaryData!="False") || BinaryData=="True")
	  {
		this->readBin(inputName);
	  }
	  else
	  {
		std::ifstream infile(inputName.c_str());
		assert(infile);
		readAscii(infile);
	  }
	}



	if (processKeys)
	{


		typedef bool(*ProcessP)( std::stringstream&  inputs, voxelImageT<T>& vxlImage);
		std::vector<std::string> processNames;	std::vector< ProcessP >  Processes;

		processNames.push_back("");  			Processes.push_back(& MCTProcessing::ignore);
		processNames.push_back("%");  		Processes.push_back(& MCTProcessing::ignore);
		processNames.push_back(";");  		Processes.push_back(& MCTProcessing::ignore);
		processNames.push_back("fillHoles");  Processes.push_back(& MCTProcessing::fillHoles);
		processNames.push_back("pore");  		Processes.push_back(& MCTProcessing::selectPore);
		processNames.push_back("threshold"); Processes.push_back(& MCTProcessing::selectPore);
		processNames.push_back("resample");	Processes.push_back(& MCTProcessing::resample);
		processNames.push_back("direction");	Processes.push_back(& MCTProcessing::redirect);
		processNames.push_back("crop");		Processes.push_back(& MCTProcessing::crop);
		processNames.push_back("cropD");		Processes.push_back(& MCTProcessing::cropD);
		processNames.push_back("resampleMax");  Processes.push_back(& MCTProcessing::resampleMax);


      while (true)
      {
		std::string tmpStr;
		std::streampos begLine = headerFile.tellg();
		headerFile>>tmpStr;
		(std::cout<<tmpStr<<": ").flush();
		bool validKey=false;
		for (unsigned int i=0; i<processNames.size();i++)
		{
			if (headerFile.fail()) break;

			if (tmpStr.compare(processNames[i])==0 || tmpStr[0]=='%')
			{
				std::stringstream keywordData;
				headerFile.get (*(keywordData.rdbuf()));
				(*Processes[i])(keywordData,*this);
				validKey=true;  break;
			}
		}
		std::cout<<std::endl;

		if ( !validKey )
        {
            std::cout<<" skipping the contents of file "<<header<<" after entry \""<<tmpStr<<"\" "<<std::endl;
            headerFile.seekg(begLine);
            break;
        }
	  }


	}

}










std::unique_ptr<voxelImageTBase> readImage
(
	std::string headerName,
	int processKeys = 1
)
{

	std::cout<<"Openning header file: "<<headerName<<std::endl;
	std::ifstream headerFile(headerName.c_str());
	if(!headerFile)  {std::cout<<"\n\n\nError: can not open header file, "<<headerName<<std::endl<<std::endl; }
	else
	{
    if (headerName.size()>4 && headerName.compare(headerName.size()-4,4,".mhd") == 0)
    {
		while (true)
		{
			std::string tmpStr;
			headerFile>>tmpStr;


			if (headerFile.fail()) break;

			std::stringstream keywordData;
			headerFile.get (*(keywordData.rdbuf()));
			std::string tmp;
			if (tmpStr == "ElementType")
			{
				keywordData >> tmp; keywordData >> tmp;
				headerFile.close();
				if (tmp=="MET_UCHAR")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned char>(headerName, processKeys)); }
				if (tmp=="MET_CHAR")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<char>(headerName, processKeys)); }
				if (tmp=="MET_USHORT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned short>(headerName, processKeys)); }
				if (tmp=="MET_SHORT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<short>(headerName, processKeys)); }
				if (tmp=="MET_UINT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned int>(headerName, processKeys)); }
				if (tmp=="MET_INT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<int>(headerName, processKeys)); }
				if (tmp=="MET_FLOAT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<float>(headerName, processKeys)); }
				if (tmp=="MET_DOUBLE")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<double>(headerName, processKeys)); }
				  
			}

		}
	 }
	}
	
	headerFile.close();
	{ return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned char>(headerName, processKeys)); }

}




template<typename T>
void readConvertFromHeader
(	voxelImageT<T>& vxlImg,
	std::string headerName,
	int processKeys = 1
)
{
	std::unique_ptr<voxelImageTBase> vxlImgTup = readImage(headerName,processKeys);
	voxelImageTBase* vxlImgT = vxlImgTup.get();
	
	bool red = false;
	{auto vxlImage = dynamic_cast<voxelImageT<char>* >(vxlImgT);            if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<" char "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<unsigned char>* >(vxlImgT);   if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<" ucar "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<short>* >(vxlImgT);           if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<" shrt "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<unsigned short>* >(vxlImgT);  if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<" usrt "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<int>* >(vxlImgT);             if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<" intg "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<unsigned int>* >(vxlImgT); 	if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<" uint "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<float>* >(vxlImgT);           if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<" flot "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<double>* >(vxlImgT);          if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<" dobl "; } }
	 
	if(!red) std::cout<<"\n\ncan not convert image\n\n"<<std::endl;
	if(!red) std::cerr<<"\n\ncan not convert image\n\n"<<std::endl;
}
