/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.


This file is part of voxelImage library, a small c++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

For further information please contact us by email:
Ali Q Raeini:    a.qaseminejad-raeini09@imperial.ac.uk

\*-------------------------------------------------------------------------*/


#ifndef voxelImageT_H
#define voxelImageT_H

#include <fstream>
#include <iostream>
#include <vector>
#include <valarray>
#include <cassert>
#include <sstream>
#include <memory>
#include "vec3.h"








template <typename Type> class voxelField
: public std::vector<std::vector<std::vector<Type> > >
{
 protected:
 public:
	voxelField(){};
	voxelField(int3 n, Type value) {  reset(n, value);  };
	voxelField(int n1, int n2, int n3, Type value) {  reset(int3{{n1,n2,n3}}, value);  };
	voxelField(const voxelField<Type>& vm ): std::vector<std::vector<std::vector<Type> > >(vm) {} ;

	void reset(int3 n, Type value);
	void reset( int n1,  int n2,  int n3, Type value) {  reset(int3{{n1,n2,n3}}, value);  };
	void readMicroCT(std::string);
	void readMicroCTHeader(std::ifstream);
	bool readAscii(std::string);
	void readAscii(std::ifstream& in);
	bool readBin(std::string fileName);
	bool readBin(std::string fileName,int iStart,int iEndp1 , int jStart,int jEndp1 , int kStart,int kEndp1 );
	void writeBin(std::string fileName) const;
	void writeBin(std::string fileName,int iStart,int iEndp1 , int jStart,int jEndp1 , int kStart,int kEndp1 ) const;
	void writeAscii(std::string fileName) const;
	void writeAscii(std::string fileName,int iStart,int iEndp1 , int jStart,int jEndp1 , int kStart,int kEndp1) const;
	void writeRotatedXZ(std::ofstream& of) const;
	int3 size3() const;
	void getSize(unsigned int& n1, unsigned int& n2, unsigned int& n3) const;
	void getPrivateSize(unsigned int& n1, unsigned int& n2, unsigned int& n3) const;
	
};




class voxelImageTBase
{
public:
	virtual ~voxelImageTBase() {};
	virtual void write(std::string fileName) = 0;
	virtual void printInfo() const = 0;

};


template <typename T>
class voxelImageT: public voxelField<T>, public voxelImageTBase
{

	vec3	X0_, dx_;

 public:

	voxelImageT():X0_(0.0,0.0,0.0),dx_(1,1,1) {};


	voxelImageT(unsigned int n1, unsigned int n2, unsigned int n3, T value)
	: voxelField<T>( n1,  n2,  n3,  value),  X0_(0.0,0.0,0.0),dx_(1,1,1) {}


	voxelImageT(int3 n, vec3 dx, vec3 xmin, T value)
	: voxelField<T>( n[0],  n[1],  n[2],  value), X0_(xmin),dx_(dx) {}

	voxelImageT(const voxelImageT & vm)
	:  voxelField<T>(vm), X0_(vm.X0_), dx_(vm.dx_){}


	voxelImageT(std::string headerName, int processKeys=1, std::string fileName="") {readFromHeader(headerName, processKeys,fileName);}
	void readFromHeader(std::string headerName, int processKeys=1, std::string fileName="")
	{	if (!headerName.empty())
		{	std::cout<<"Openning header file: "<<headerName<<std::endl;
			std::ifstream headerFile(headerName.c_str());
			if(!headerFile)  {std::cout<<"\n\n\nError: can not open header file, "<<headerName<<std::endl<<std::endl; }
			else
				readFromHeader(headerFile,headerName,processKeys,fileName);
			headerFile.close();
		}
	}

	void readFromHeader( std::ifstream& headerFile,	std::string headerName, int processKeys=1, std::string fileName="");



	bool readAscii(std::string fileName)
	{	///  overwrite as the parent -voxelField<T>- interprets
		///  numerical values as characters not integers

		std::cout<<  " reading "<<fileName<<std::endl;

		//if ( (fileName.compare(fileName.size()-4,4,".dat")==0) || (fileName.compare(fileName.size()-4,4,".txt") == 0) )
		//{
			std::ifstream in(fileName.c_str());
			assert(in);

			char tmpc[8];
			for ( unsigned int i=0; i<8;i++)   in>>tmpc[i];
			if (std::string(tmpc).compare(0,4,"ascii") == 0) //ignore first lines
			{
				int n[3];
				in>>n[2]>>n[0]>>n[1];//ignore first lines
				double  xmin[3],xmax[3];
				in>> xmin[0]>>xmax[0]>>xmin[1]>>xmax[1]>>xmin[2]>>xmax[2] ;
				std::cout<<"Warning: ignoring the header of file "<<fileName<<std::endl;
			}
			else
				in.seekg(0, in.beg);
			readAscii(in);
			in.close();
			return !in.fail();
		//}
		//else
			//this->readBin(fileName);

	}


	void  readAscii(std::ifstream& in)
	{	///  overwrite as the parent -voxelField<T>- interprets
		///  numerical values as characters not integers
		int tmp=0;
		for (typename std::vector<std::vector<std::vector<T> > >::iterator d1=this->begin();  d1<this->end();  d1++)
			for (typename  std::vector<std::vector<T> >::iterator d2=d1->begin();d2<d1->end();d2++)
				for (typename  std::vector<T>::iterator d3=d2->begin();d3<d2->end();d3++)
				{
					in>>tmp;
					*d3=tmp;
				}
	}



	void cropD( int3 cropBegin,  int3 cropEnd,unsigned int emptylayers=0, T emptylayersValue=1) ;
	void crop( unsigned int cropBegin[3],  unsigned int cropEnd[3],unsigned int emptylayers=0, T emptylayersValue=1) ;
	void crop(unsigned int iStart, unsigned int iEnd ,
				 unsigned int jStart, unsigned int jEnd ,
				 unsigned int kStart, unsigned int kEnd ,
				 unsigned int emptylayers=0,T emptylayersValue=1);

	void writeHeader(std::string fileName) const;
	void writeHeader(std::string fileName, int3 iStart, int3 iEnd) const;


	void erodeLayer(int i);
	void resample(double i);
	void resampleMax(double i);
	void rotate(char direction);
	void PointMedian026(unsigned int thereshold0,unsigned int thereshold1);
	void FaceMedian06(unsigned int thereshold0,unsigned int thereshold1);
	void median(short nNeist);

	void AND(const voxelImageT& data2);
	void NOT(const voxelImageT& data2);
	void OR(const voxelImageT& data2);
	void XOR(const voxelImageT& data2);

	void fillHoles(unsigned int maxHoleRadius);

	void shrinkPore();
	void growPore();


   void threshold101(T theresholdMin,T theresholdMax);
	void  segmentPhase
	(
		T StartThreshold1,
		T StartThreshold2,
		T finalThreshold1,
		T finalThreshold2,
		unsigned int nIter
	);

	void writeAConnectedPoreVoxel(std::string fileName) const;
	void setLayer(int k, std::vector<std::vector<T> > &Values);
	void replaceyLayer(int j, int fromj);
	void replacexLayer(int i, int fromi);
	void setBlock(int n1, int n2, int n3, const std::vector<std::vector<std::vector<T> > >&Values);
	template<typename T2> void resetFrom(const voxelImageT<T2>&Values);
	void setFrom(const voxelImageT<T>&Values, int n1, int n2, int n3);
	void growBox(unsigned int nlyr);
	void shrinkBox(unsigned int nlyr)
		{unsigned int beg[3]={nlyr,nlyr,nlyr},
		end[3]={unsigned((*this)[0][0].size()-nlyr-1),unsigned((*this)[0].size()-nlyr-1),unsigned((*this).size()-nlyr-1)};
		crop(beg,end);};

	void write(std::string fileName);
	double volFraction(T vv1,T vv2) const;
	void printInfo() const;

	vec3 X0() const {return X0_;};
	vec3& X0Ch()    {return X0_;};
	vec3 dx() const {return dx_;};
	vec3& dxCh()    {return dx_;};

 };





typedef voxelImageT<unsigned char> voxelImage;

#include "voxelImageI.h"
#include "voxelImageReader.h"


#endif
