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
Ali Q Raeini: a.qaseminejad-raeini09@imperial.ac.uk
\*-------------------------------------------------------------------------*/

#include "voxelImage.h"
#include <map>
#include <limits>   // std::numeric_limits


//using namespace std;

#define forAllNei(r1,r2) \
for (short k_nei_m=r1;k_nei_m<=r2;++k_nei_m) \
for (short j_nei_m=r1;j_nei_m<=r2;++j_nei_m) \
for (short i_nei_m=r1;i_nei_m<=r2;++i_nei_m)

#define nei(dataM,i_M,j_M,k_M) dataM[k_M+k_nei_m][j_M+j_nei_m][i_M+i_nei_m]
#define distSqrNei() (k_nei_m*k_nei_m+j_nei_m*j_nei_m+i_nei_m*i_nei_m)


#define 	forAllkji(datas3s_M)   \
	for ( unsigned int k=0; k<datas3s_M.size() ; ++k )   \
	 for ( unsigned int j=0; j<datas3s_M[k].size() ; ++j )   \
	  for ( unsigned int i=0; i<datas3s_M[k][j].size() ; ++i )


#define forAllNeiInt(ri1,ri2,rj1,rj2,rk1,rk2) \
for (int k_nei_m = rk1;k_nei_m <= rk2;++k_nei_m) \
for (int j_nei_m = rj1;j_nei_m <= rj2;++j_nei_m) \
for (int i_nei_m = ri1;i_nei_m <= ri2;++i_nei_m)

#define forAllNeiU(ri1,ri2,rj1,rj2,rk1,rk2) \
for (unsigned int k_nei_m = rk1;k_nei_m <= rk2;++k_nei_m) \
for (unsigned int j_nei_m = rj1;j_nei_m <= rj2;++j_nei_m) \
for (unsigned int i_nei_m = ri1;i_nei_m <= ri2;++i_nei_m)

#define neib(dataM) dataM[k_nei_m][j_nei_m][i_nei_m]


template<typename  T> inline T maxNei(const voxelImageT<T>& image, int i, int j, int k, int r11, int r22)
{
	T maxx=std::numeric_limits<T>::min();
	forAllNei(r11,r22)
	{
		maxx=std::max(maxx,nei(image,i,j,k));
	}
	return maxx;
}


#define ROCKVV 1





template<typename Type>   void voxelField<Type>::reset(int3 n, Type value)
{
	this->resize(0);
	this->resize(n[2],std::vector<std::vector<Type> >(n[1],std::vector<Type>(n[0],value)));
}



//kji
template<typename Type>   void voxelField<Type>::getSize(unsigned int& n1, unsigned int& n2, unsigned int& n3) const
{
  n3 = (*this).size();
  if (n3>0)
  {
	  n1 = (*this)[0][0].size();
	  n2 = (*this)[0].size();
   }
   else
   {
	  n1 = 0;
	  n2 = 0;
   }
}



template<typename Type>
int3 voxelField<Type>::size3() const
{
  int3 mySize;
  mySize[2] = (*this).size();
  if (mySize[2]>0)
  {
	  mySize[0] = (*this)[0][0].size();
	  mySize[1] = (*this)[0].size();
   }
   return mySize;
}

//kji
template<typename Type>   void voxelField<Type>::getPrivateSize(unsigned int& n1, unsigned int& n2, unsigned int& n3) const
{
  n1 = (*this).size();
  if (n1>0)
  {
	  n2 = (*this)[0].size();
	  n3 = (*this)[0][0].size();
   }
}



//read order sensitive
template<typename Type>   void voxelField<Type>::readAscii(std::ifstream& in)
{
	for (
	typename std::vector<std::vector<std::vector<Type> > >::iterator d1 = this->begin();
	d1<this->end();
	d1++)
	{
		for ( typename std::vector<std::vector<Type> >::iterator d2 = d1->begin();d2<d1->end();d2++)
		{
			for ( typename std::vector<Type>::iterator d3 = d2->begin();d3<d2->end();d3++)
			{
				in>>*d3;
			}
		}
	}
}

//read order sensitive
template<typename Type>   bool voxelField<Type>::readAscii(std::string fileName)
{
	std::cout<<  " reading ascii file "<<fileName<<std::endl;
	std::ifstream in(fileName.c_str());
	if(!in)  {std::cout<<"\n\nError: can not open image file, "<<fileName<<std::endl<<std::endl;
			   return false;}


	readAscii(in);

	in.close();
	return !in.fail();
}




//read order sensitive
template<typename Type>   void voxelField<Type>::readMicroCT(std::string fileName)
{
	std::cout<<  " reading micro-CT file "<<fileName<<std::endl;
	std::ifstream in(fileName.c_str());
	if(!in)  {std::cout<<"\n\nError: can not open image file, "<<fileName<<std::endl<<std::endl; exit(-1);}

	char tmpc;
	for ( unsigned int i = 0; i<8;i++)   in>>tmpc, std::cout<<" "<<tmpc;  //ignore the first 8 characters (ascii 3uc)
	unsigned int n[3];
	double  xmin[3];
	double  xmax[3];

	in>>n[2]>>n[1]>>n[0];						// number of variables (dimension of
	//~ in>>	dx[0]>>dx[1]>>dx[2] ;
	in>>	xmin[0]>>xmin[1]>>xmin[2] ;
	in>>	xmax[0]>>xmax[1]>>xmax[2] ;

	readAscii(in);

	in.close();
}

//kji
template<typename Type>   bool voxelField<Type>::readBin(std::string fileName)
{
		int3 n = size3();
		return readBin(fileName,0,n[0],0,n[1],0,n[2]);
}

//kji
template<typename Type>   bool voxelField<Type>::readBin(std::string fileName,
								int iStart,int iEnd ,
								int jStart,int jEnd ,
								int kStart,int kEnd
)
{
	std::cout<<  " reading binary file "<<fileName;
	std::ifstream in (fileName.c_str(), std::ios::in | std::ios::binary);
	if(!in)  {std::cout<<"\n\n  Error: can not open image file, "<<fileName<<std::endl<<std::endl;
		return false;}
	int k = kStart;
	for ( ;k < kEnd;k++)
	{
		for ( int j = jStart;j < jEnd;j++)
		{
			if(in)  in.read(reinterpret_cast<char*>(&(*this)[k][j][iStart]), (iEnd-iStart)*sizeof(Type) );
		}
		if (!in) break;
	}

	if (!in)
	{
		std::cout<<  "\n\n ***** Error in reading "<<fileName<<" ***** \n"<<"only "<<k<<" layers read"<<std::endl;
		this->resize(k);
	}
	in.close();
	std::cout<<  "."<<std::endl;
	return true;
}


template<typename Type>   void voxelField<Type>::writeBin(std::string fileName) const
{
	int3 imgsize = size3();
	writeBin(fileName,0, imgsize[0] ,0,imgsize[1] ,0, imgsize[2]);
}



//kji
template<typename Type>   void voxelField<Type>::writeBin(std::string fileName,
							   int iStart,int iEnd , int jStart,int jEnd , int kStart,int kEnd ) const
{

	std::cout<<  " writting binary file "<<fileName<<";  i: "<<iStart<<" "<<iEnd<<",  j: "<<jStart<<" "<<jEnd<<",  k: "<<kStart<<" "<<kEnd; std::cout.flush();

	std::ofstream of (fileName.c_str(), std::ios::out | std::ios::binary);
	assert(of);
	for ( int k = kStart;k < kEnd;k++)
	{
		for ( int j = jStart;j < jEnd;j++)
		{
			if (!(*this)[k][j].empty())
				of.write(reinterpret_cast<const char*>(&((*this)[k][j][iStart])), (iEnd-iStart) * sizeof(Type));
			else
				std::cout<<"\n  Error when writting binary file, k:"<<k<<  "j:"<<k<<" "<<std::endl;
		}
	}
	of.flush();
	of.close();
	(std::cout<<  ".  ").flush();

}

template<typename Type>   void voxelField<Type>::writeAscii(std::string fileName,int iStart,int iEnd , int jStart,int jEnd , int kStart,int kEnd) const
{
	std::cout<<  " writting ascii file "<<fileName<<";  "; std::cout.flush();

	std::ofstream of (fileName.c_str());
	assert(of);

	for ( int k = kStart;k < kEnd;k++)
	{
	  for ( int j = jStart;j < jEnd;j++)
	  {
		for ( int i = iStart;i < iEnd;i++)
		{

		   of<<double((*this)[k][j][i])<<' ';
		}
		of<<"\n";
	  }
	}
	of<<std::endl;
	of.close();
	std::cout<<  "."<<std::endl;
}

template<typename Type>   void voxelField<Type>::writeAscii(std::string fileName) const
{
	int3 imgsize = size3();
	writeAscii(fileName,0, imgsize[0] ,0,imgsize[1] ,0, imgsize[2]);
}

template<typename Type>   void voxelField<Type>::writeRotatedXZ(std::ofstream& of) const
{
	for (unsigned int i = 0;i<(*this)[0][0].size();i++) //reversed order with k
	{
		for (unsigned int j = 0;j<(*this)[0].size();j++)
		{
			for (unsigned int k = 0;k<(*this).size();k++) //reversed order with i
			{
				of<<double((*this)[k][j][i])<<' ';
			}
			of<<std::endl;
		}
	}
}




template<typename T>
void voxelImageT<T>::cropD( int3 cropBegin,  int3 cropEnd,unsigned int emptylayers, T emptylayersValue)
{
	crop(cropBegin[0],cropEnd[0]-1,cropBegin[1],cropEnd[1]-1,cropBegin[2],cropEnd[2]-1,emptylayers,emptylayersValue);
}


template<typename T>
void voxelImageT<T>::crop( unsigned int cropBegin[3],  unsigned int cropEnd[3],unsigned int emptylayers, T emptylayersValue)
{
	crop(cropBegin[0],cropEnd[0],cropBegin[1],cropEnd[1],cropBegin[2],cropEnd[2],emptylayers,emptylayersValue);
}

//kji
template<typename T>
void voxelImageT<T>::crop(
							unsigned int iStart, unsigned int iEnd ,
							unsigned int jStart, unsigned int jEnd ,
							unsigned int kStart, unsigned int kEnd  ,
							unsigned int emptylayers, T emptylayersValue
)
{
	(std::cout<<  "  cropping, "<<  "   ["<<iStart<<" "<<iEnd+1 <<  ")  ["<<jStart<<" "<<jEnd+1<< ")  ["<<kStart<<" "<<kEnd+1<<")  ").flush();
	if (emptylayersValue) (std::cout<<  ", adding "<<emptylayers<<" layers with value "<< double(emptylayersValue)<<"  ").flush();

	X0_[0]=X0_[0]+(int(iStart)-int(emptylayers))*dx_[0];   X0_[1]=X0_[1]+(int(jStart)-int(emptylayers))*dx_[1];   X0_[2]=X0_[2]+(int(kStart)-int(emptylayers))*dx_[2];

	voxelImageT<T> tmp=*this;


	this->resize(0);
	this->resize(kEnd+1-kStart+2*emptylayers);
	for ( unsigned int k=0; k<kEnd+1-kStart+2*emptylayers ; k++ )
	{
		(*this)[k].resize(jEnd+1-jStart+2*emptylayers);
		for ( unsigned int j=0; j<jEnd+1-jStart+2*emptylayers; ++j )
		{
			(*this)[k][j].resize(iEnd+1-iStart+2*emptylayers,emptylayersValue);
			for ( unsigned int i=0; i<iEnd-iStart+1+2*emptylayers; ++i )
			{
				(*this)[k][j][i]=emptylayersValue;
			}
		}
	}

	for ( unsigned int k=0; k<kEnd+1-kStart; k++ )
	{
		for ( unsigned int j=0; j<jEnd+1-jStart; ++j )
		{
			for ( unsigned int i=0; i<iEnd+1-iStart; ++i )
			{
				(*this)[k+emptylayers][j+emptylayers][i+emptylayers]=tmp[k+kStart][j+jStart][i+iStart];
			}
		}
	}


}




template<typename T>
void voxelImageT<T>::setLayer(int k, std::vector<std::vector<T> > &Values)
{
		//~ for (unsigned int i=0 ;i<(*this)[k].size(); i++)
		//~ {
			(*this)[k]=Values;
		//~ }
}

template<typename T>
void voxelImageT<T>::replacexLayer(int i, int fromi)
{

	 for ( unsigned int k=0; k<(*this).size() ; k++ )
	{
		for ( unsigned int j=0; j<(*this)[k].size() ; ++j )
		{
			//~ for ( unsigned int i=0; i<Values[k][j].size() ; ++i )
			//~ {
				(*this)[k][j][i]=(*this)[k][j][fromi];
			//~ }
		}
	}

}
template<typename T>
void voxelImageT<T>::replaceyLayer(int j, int fromj)
{

	 for ( unsigned int k=0; k<(*this).size() ; k++ )
	{
		//~ for ( unsigned int j=0; j<Values[k].size() ; ++j )
		//~ {
			for ( unsigned int i=0; i<(*this)[k][j].size() ; ++i )
			{
				(*this)[k][j][i]=(*this)[k][fromj][i];
			}
		//~ }
	}

}

//kji
template<typename T>
void voxelImageT<T>::setBlock(int n1, int n2, int n3, const std::vector<std::vector<std::vector<T> > >&Values)
{
	for ( unsigned int k=0; k<Values.size() ; k++ )
	 for ( unsigned int j=0; j<Values[k].size() ; ++j )
	  for ( unsigned int i=0; i<Values[k][j].size() ; ++i )
			(*this)[k+n3][j+n2][i+n1]=Values[k][j][i];
}

//kji
template<typename T>
template<typename T2>
void voxelImageT<T>::resetFrom(const voxelImageT<T2>&Values)
{
	dx_= Values.dx();
	X0_ = Values.X0();
	this->reset(Values.size3(),0);
	for ( unsigned int k=0; k<(*this).size() ; k++ )
	 for ( unsigned int j=0; j<(*this)[k].size() ; ++j )
	  for ( unsigned int i=0; i<(*this)[k][j].size() ; ++i )
			(*this)[k][j][i]=Values[k][j][i];
}
template<typename T>
void voxelImageT<T>::setFrom(const voxelImageT<T>&Values, int n1, int n2, int n3)
{
	dx_= Values.dx();
	X0_[0] = Values.X0()[0]+n1*dx_[0];
	X0_[1] = Values.X0()[1]+n2*dx_[1];
	X0_[2] = Values.X0()[2]+n3*dx_[2];
	for ( unsigned int k=0; k<(*this).size() ; k++ )
	 for ( unsigned int j=0; j<(*this)[k].size() ; ++j )
	  for ( unsigned int i=0; i<(*this)[k][j].size() ; ++i )
			(*this)[k][j][i]=Values[k+n3][j+n2][i+n1];
}

//kji
template<typename T>
void voxelImageT<T>::growBox(unsigned int nLayers)
{

	int3 n = (*this).size3();
	(*this).crop(0,n[0]-1,0,n[1]-1,0,n[2]-1, nLayers,1);//		 XXXXXXXXXXXXXXXXXXXXXXXXXXXX

	for (unsigned  int i=0; i<nLayers ; i++ )
	{
		(*this).replaceyLayer(n[1]+nLayers+i, n[1]+nLayers-1);
		(*this).replaceyLayer(i, nLayers);
		(*this).replacexLayer(n[0]+nLayers+i, n[0]+nLayers-1);
		(*this).replacexLayer(i, nLayers);
		(*this).setLayer(n[2]+nLayers+i, (*this)[n[2]+nLayers-1]);
		(*this).setLayer(i, (*this)[nLayers]);
	}
}



template<typename T>
void voxelImageT<T>::resample(double nReSampleNotSafe)//  TODO to be tested
{
	if (nReSampleNotSafe < .999)
	{
		int nReSample=1.0/nReSampleNotSafe+0.5;
		voxelImageT<T> tmp=*this;
		//~ reset((*this)[k].size(),int n2,int n3, Type value)
		unsigned int n3=this->size(),n2=(*this)[0].size(),n1=(*this)[0][0].size();
		this->resize(nReSample*n3);
		for ( unsigned int k=0; k<this->size() ; k++ )
		{

			(*this)[k].resize(nReSample*n2);
			for ( unsigned int j=0; j<(*this)[k].size() ; ++j )
			{
				(*this)[k][j].resize(nReSample*n1);
				for ( unsigned int i=0; i<(*this)[k][j].size() ; ++i )
				{
					(*this)[k][j][i]=tmp[(0.5+k)/nReSample][(0.5+j)/nReSample][(0.5+i)/nReSample];
				}
			}
		}
		dx_/=nReSample; //fixed
	}
	else if (nReSampleNotSafe > 1.001)
	{
		int nReSample=nReSampleNotSafe+0.5; /// Warning unsigned doesn't work  wTf
		voxelImageT<T> tmp=*this;
		//~ reset((*this)[k].size(),int n2,int n3, Type value)
		unsigned int n3=this->size(),n2=(*this)[0].size(),n1=(*this)[0][0].size();
		this->resize((n3)/nReSample);
		for ( unsigned int k=0; k<this->size() ; k++ )
		{

			(*this)[k].resize((n2)/nReSample);
			for ( unsigned int j=0; j<(*this)[k].size() ; ++j )
			{
				(*this)[k][j].resize((n1)/nReSample);
				for ( unsigned int i=0; i<(*this)[k][j].size() ; ++i )
				{
					int neiSum=0;
					forAllNei(0,nReSample-1)
					{
						neiSum+=nei(tmp,i*nReSample,j*nReSample,k*nReSample);
					}
					(*this)[k][j][i]=(0.5+double(neiSum)/(nReSample*nReSample*nReSample));

				}
			}
		}
		dx_*=nReSample;
	}

		//~ for (unsigned int i=0 ;i<(*this)[k].size(); i++)
		//~ {
			//~ (*this)[k][i]=Values[i];
		//~ }
}


template<typename T>
void voxelImageT<T>::resampleMax(double nReSampleNotSafe)//  TODO to be tested
{
	if (nReSampleNotSafe < .999)
	{
		int nReSample=1.0/nReSampleNotSafe+0.5;
		voxelImageT<T> tmp=*this;
		//~ reset((*this)[k].size(),int n2,int n3, Type value)
		unsigned int n3=this->size(),n2=(*this)[0].size(),n1=(*this)[0][0].size();
		this->resize(nReSample*n3);
		for ( unsigned int k=0; k<this->size() ; k++ )
		{

			(*this)[k].resize(nReSample*n2);
			for ( unsigned int j=0; j<(*this)[k].size() ; ++j )
			{
				(*this)[k][j].resize(nReSample*n1);
				for ( unsigned int i=0; i<(*this)[k][j].size() ; ++i )
				{
					(*this)[k][j][i]=tmp[(0.5+k)/nReSample][(0.5+j)/nReSample][(0.5+i)/nReSample];
				}
			}
		}
		dx_/=nReSample; //fixed
	}
	else if (nReSampleNotSafe > 1.001)
	{
		int nReSample=nReSampleNotSafe+0.5; /// Warning unsigned doesn't work  wTf
		voxelImageT<T> tmp=*this;
		//~ reset((*this)[k].size(),int n2,int n3, Type value)
		unsigned int n3=this->size(),n2=(*this)[0].size(),n1=(*this)[0][0].size();
		this->resize((n3)/nReSample);
		for ( unsigned int k=0; k<this->size() ; k++ )
		{

			(*this)[k].resize((n2)/nReSample);
			for ( unsigned int j=0; j<(*this)[k].size() ; ++j )
			{
				(*this)[k][j].resize((n1)/nReSample);
				for ( unsigned int i=0; i<(*this)[k][j].size() ; ++i )
				{
					T neiSum=0;
					forAllNei(0,nReSample-1)
					{
						neiSum=std::max(neiSum, nei(tmp,i*nReSample,j*nReSample,k*nReSample));
					}
					(*this)[k][j][i]=neiSum;//(0.5+double(neiSum)/(nReSample*nReSample*nReSample));

				}
			}
		}
		dx_*=nReSample;
	}
}



template<typename T>
void voxelImageT<T>::rotate(char direction)
{// wrong X0
	unsigned int n1,n2,n3;

	voxelField<T>::getSize(n1,n2,n3);
	if (direction=='z')
	{
		//~ int nMinTmp=nMin_[0];
		//~ nMin_[0]=nMin_[2];
		//~ nMin_[2]=nMinTmp;
		double X0Tmp=X0_[0];
		X0_[0]=X0_[2];
		X0_[2]=X0Tmp;
		double dxTmp=dx_[0];
		dx_[0]=dx_[2];
		dx_[2]=dxTmp;

		voxelImageT<T> tmp=*this;
		this->reset(n3,n2,n1,0);
		for ( unsigned int k=0; k<n3 ; k++ )
		{
			for ( unsigned int j=0; j<n2 ; ++j )
			{
				for ( unsigned int i=0; i<n1 ; ++i )
				{
					(*this)[i][j][k]=tmp[k][j][i];
				}
			}
		}
	}
	else if (direction=='y')
	{
		//~ int nMinTmp=nMin_[0];
		//~ nMin_[0]=nMin_[1];
		//~ nMin_[1]=nMinTmp;
		double X0Tmp=X0_[0];
		X0_[0]=X0_[1];
		X0_[1]=X0Tmp;
		double dxTmp=dx_[0];
		dx_[0]=dx_[1];
		dx_[1]=dxTmp;

		voxelImageT<T> tmp=*this;
		this->reset(n2,n1,n3,0);
		for ( unsigned int k=0; k<n3 ; k++ )
			for ( unsigned int j=0; j<n2 ; ++j )
				for ( unsigned int i=0; i<n1 ; ++i )
					(*this)[k][i][j]=tmp[k][j][i];
	}
	else if (direction=='-')
	{
		std::cout<<" -> flipping image,  x origin will be invalid "<<std::endl;
		voxelImageT<T> tmp=*this;
		for ( unsigned int k=0; k<n3 ; k++ )
			for ( unsigned int j=0; j<n2 ; ++j )
				for ( unsigned int i=0; i<n1 ; ++i )
					(*this)[k][j][n1-1-i]=tmp[k][j][i];
	}
	else
	{
		std::cout<<"\n\nSwapping "<<direction<<" and x directions(!?!), sorry can't do that >-( "<<std::endl;
		std::cerr<<"Swapping "<<direction<<" and x directions(!?!), sorry can't do that >-( \n\n"<<std::endl;
	}

}


template<typename T>
void voxelImageT<T>::PointMedian026(unsigned int thereshold0,unsigned int thereshold1)
{
	unsigned long nChanged(0);
	vec3 doubletmp(0,0,0);
	int3 n = voxelField<T>::size3();
	for (unsigned int i=0;i<3;i++) n[i]=n[i]+2;
	
	voxelImageT<T> voxls(n,doubletmp,doubletmp,1);
		voxls.setBlock(0, 0, 0, (*this));
		voxls.setBlock(2, 2, 2, (*this));
		voxls.setBlock(1, 1, 1, (*this));


	for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
	{
		for ( unsigned int j=1; j<voxls[k].size()-1 ; ++j )
		{
			for ( unsigned int i=1; i<voxls[k][j].size()-1 ; ++i )
			{
				unsigned int neiSum=0;
				forAllNei(-1,1)
				{
					neiSum+=nei(voxls,i,j,k);
				}
				neiSum-=voxls[k][j][i];
				if (neiSum <= thereshold0  && (*this)[k-1][j-1][i-1])
				{
					(*this)[k-1][j-1][i-1]=0;
					++nChanged;
				}
				else if (neiSum >= thereshold1  && !((*this)[k-1][j-1][i-1]))
				{
					(*this)[k-1][j-1][i-1]=1;
					++nChanged;
				}
		   }
		}
	}
	std::cout<<"PointMedian026  changed: "<<nChanged<<std::endl;

}





//~
//~ #define forAllFaceNei \/
//~ for (int k_nei_m=-1;k_nei_m<2;k_nei_m++) \/
//~ for (int j_nei_m=-1;j_nei_m<2;j_nei_m++) \/
//~ for (int i_nei_m=-1;i_nei_m<2;i_nei_m++)
//~
//~ #define nei(i,j,k) (*this)[k+k_nei_m][j+j_nei_m][i+i_nei_m]


template<typename T>
void voxelImageT<T>::FaceMedian06(unsigned int thereshold0,unsigned int thereshold1)
{
	unsigned long nChanged(0);
	vec3 doubletmp(0,0,0);
	int3 n = voxelField<T>::size3();
	for (unsigned int i=0;i<3;i++) n[i]=n[i]+2;

	voxelImageT<T> voxls(n,doubletmp,doubletmp,1);
		voxls.setBlock(0, 0, 0, (*this));
		voxls.setBlock(2, 2, 2, (*this));
		voxls.setBlock(1, 1, 1, (*this));

	for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
	{
		for ( unsigned int j=1; j<voxls[k].size()-1 ; ++j )
		{
			for ( unsigned int i=1; i<voxls[k][j].size()-1 ; ++i )
			{
				 unsigned int neiSum = //voxls[k][j][i]
									  voxls[k][j][i-1]+voxls[k][j][i+1]
									  +voxls[k][j-1][i]+voxls[k][j+1][i]
									  +voxls[k+1][j][i]+voxls[k-1][j][i];

				if (neiSum <= thereshold0 && (*this)[k-1][j-1][i-1])
				{
					(*this)[k-1][j-1][i-1]=0;
					++nChanged;
				}
				else if (neiSum >= thereshold1 && !((*this)[k-1][j-1][i-1]))
				{
					(*this)[k-1][j-1][i-1]=1;
					++nChanged;
				}
			}
		}
	}
	std::cout<<"FaceMedian06  changed: "<<nChanged<<std::endl;
}


template<typename T>
void voxelImageT<T>::shrinkPore()
{
	voxelImageT<T> voxls=*this;


	for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
	{
		for ( unsigned int j=1; j<voxls[k].size()-1 ; ++j )
		{
			for ( unsigned int i=1; i<voxls[k][j].size()-1 ; ++i )
			{

				if (voxls[k][j][i]==0 && ( voxls[k][j][i-1]  || voxls[k][j][i+1] ||
									  voxls[k][j-1][i] || voxls[k][j+1][i] || voxls[k-1][j][i] || voxls[k+1][j][i] ) )
					(*this)[k][j][i]=1;
		   }
		}
	}


	for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
	{
		for ( unsigned int j=1; j<voxls[k].size()-1 ; ++j )
		{
			//~ for ( unsigned int i=0; i<voxls[k][j].size()-1 ; ++i )
			{	unsigned int i=0;
				if (voxls[k][j][i]==0 && ( voxls[k][j+1][i] || voxls[k][j-1][i] || voxls[k+1][j][i] || voxls[k-1][j][i] ) )
					(*this)[k][j][i]=1;

				i=voxls[k][j].size()-1;
				if (voxls[k][j][i]==0 && ( voxls[k][j+1][i] || voxls[k][j-1][i] || voxls[k+1][j][i] || voxls[k-1][j][i] ) )
					(*this)[k][j][i]=1;
		   }
		}
	}


	for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
	{
		//~ for ( unsigned int j=0; j<voxls[k].size()-1 ; ++j )
		{
			for ( unsigned int i=1; i<voxls[k][0].size()-1 ; ++i )
			{
				unsigned int j=0;
				if (voxls[k][j][i]==0 && ( voxls[k][j][i-1] || voxls[k][j][i+1] || voxls[k-1][j][i] || voxls[k+1][j][i] ) )
					(*this)[k][j][i]=1;

				j=voxls[k].size()-1;
				if (voxls[k][j][i]==0 && ( voxls[k][j][i-1] || voxls[k][j][i+1] || voxls[k-1][j][i] || voxls[k+1][j][i] ) )
					(*this)[k][j][i]=1;
		   }
		}
	}

	//~ for ( unsigned int k=0; k<voxls.size()-1 ; k++ )
	{
		for ( unsigned int j=1; j<voxls[0].size()-1 ; ++j )
		{
			for ( unsigned int i=1; i<voxls[0][0].size()-1 ; ++i )
			{
				unsigned int k=0;
				if (voxls[k][j][i]==0 && ( voxls[k][j][i-1] || voxls[k][j][i+1] || voxls[k][j-1][i] || voxls[k][j+1][i] ) )
					(*this)[k][j][i]=1;

				k=voxls.size()-1;
				if (voxls[k][j][i]==0 && ( voxls[k][j][i-1] || voxls[k][j][i+1] || voxls[k][j-1][i] || voxls[k][j+1][i] ) )
					(*this)[k][j][i]=1;
		   }
		}
	}


}

template<typename T>
class mapComparer  {  public: bool operator() (std::pair<const T,short>& i1, std::pair<const T,short> i2) {return i1.second<i2.second;}  };

template<typename T>
void voxelImageT<T>::median(short nNeist)
{
	voxelImageT<T> voxls=*this;
	long long nChanges = 0;
	for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
	for ( unsigned int j=1; j<voxls[k].size()-1 ; j++ )
	for ( unsigned int i=1; i<voxls[k][j].size()-1 ; i++ )
	{
		register T pID = voxls[k][j][i];

		short nSames(0);
		std::map<T,short> neis;///.  ID-counter

		register T
		neiPID = voxls[k][j][i-1];
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls[k][j][i+1];
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls[k][j-1][i];
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls[k][j+1][i];
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls[k-1][j][i];
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls[k+1][j][i];
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;

		if(nSames<nNeist)
		{
			typename std::map<T,short>::iterator neitr = max_element(neis.begin(), neis.end(), mapComparer<T>());
			if (neitr->second>nSames)
			{
				++nChanges;
				(*this)[k][j][i] = neitr->first;
				//~ std::cout<<"  * "<<int(pID)<<" "<<int(neitr->first)<<" "<<int(nSames)<<" "<<(neitr->second)<<" "<<neis[0]<<" "<<neis[1]<<" "<<neis[2]<<"  * ";
			}
			else if ( pID!=ROCKVV && nSames==neis[1])
			{
				++nChanges;
				(*this)[k][j][i] = 1;
				std::cout<<"  * "<<int(pID)<<":"<<int(nSames)<<" "<<int(neitr->first)<<":"<<(neitr->second)<<" "<<neis[0]<<" "<<neis[1]<<" "<<neis[2]<<"  * ";
			}
			else if ( pID!=1 )
			{ std::cout<<"  X "<<int(pID)<<":"<<int(nSames)<<" "<<int(neitr->first)<<":"<<(neitr->second)<<" "<<neis[0]<<" "<<neis[1]<<" "<<neis[2]<<"  X ";
			}
		 }
	  }

	( std::cout<<"  nMedian: "<< std::left<<nChanges<<"  \n").flush();

}





template<typename T>
void voxelImageT<T>::growPore() // optimized function, should be further optimized as it is frequently used
{

	voxelImageT<T> voxls=*this;


	for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
	{
		for ( unsigned int j=1; j<voxls[k].size()-1 ; ++j )
		{
			for ( unsigned int i=1; i<voxls[k][j].size()-1 ; ++i )
			{

				if (voxls[k][j][i] && ( !voxls[k][j][i-1] || !voxls[k][j][i+1] ||
									  !voxls[k][j-1][i] || !voxls[k+1][j][i] || !voxls[k][j+1][i] || !voxls[k-1][j][i] ) )
					(*this)[k][j][i]=0;
		   }
		}
	}


	for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
	{
		for ( unsigned int j=1; j<voxls[k].size()-1 ; ++j )
		{
			// for ( unsigned int i=0; i<voxls[k][j].size()-1 ; ++i )
			{	unsigned int i=0;
				if (voxls[k][j][i] && ( !voxls[k][j][i+1] || !voxls[k][j-1][i]  || !voxls[k+1][j][i] || !voxls[k][j+1][i] || !voxls[k-1][j][i] ) )
					(*this)[k][j][i]=0;

				i=voxls[k][j].size()-1;
				if (voxls[k][j][i] && ( !voxls[k][j][i-1] || !voxls[k][j-1][i]  || !voxls[k+1][j][i] || !voxls[k][j+1][i] || !voxls[k-1][j][i] ) )
					(*this)[k][j][i]=0;
		   }
		}
	}


	for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
	{
		// for ( unsigned int j=0; j<voxls[k].size()-1 ; ++j )
		{
			for ( unsigned int i=1; i<voxls[k][0].size()-1 ; ++i )
			{	unsigned int j=0;

				if (voxls[k][j][i] && ( !voxls[k][j][i-1] || !voxls[k][j][i+1]  || !voxls[k+1][j][i] || !voxls[k][j+1][i] || !voxls[k-1][j][i] ) )
					(*this)[k][j][i]=0;

				j=voxls[k].size()-1;
				if (voxls[k][j][i] && ( !voxls[k][j][i-1] || !voxls[k][j][i+1] || !voxls[k][j-1][i]  || !voxls[k+1][j][i] || !voxls[k-1][j][i] ) )
					(*this)[k][j][i]=0;
		   }
		}
	}

	// for ( unsigned int k=0; k<voxls.size()-1 ; k++ )
	{
		for ( unsigned int j=1; j<voxls[0].size()-1 ; ++j )
		{
			for ( unsigned int i=1; i<voxls[0][0].size()-1 ; ++i )
			{	unsigned int k=0;

				if (voxls[k][j][i] && ( !voxls[k][j][i-1] || !voxls[k][j][i+1] || !voxls[k+1][j][i] || !voxls[k][j+1][i] || !voxls[k][j-1][i] ) )
					(*this)[k][j][i]=0;

				k=voxls.size()-1;
				if (voxls[k][j][i] && ( !voxls[k][j][i-1] || !voxls[k][j][i+1] || !voxls[k][j-1][i] || !voxls[k][j+1][i] || !voxls[k-1][j][i] ) )
					(*this)[k][j][i]=0;
		   }
		}
	}

}








template<typename T>
void voxelImageT<T>::NOT(const voxelImageT& data2)
{
	forAllkji((*this))	(*this)[k][j][i]= (*this)[k][j][i] && !data2[k][j][i];
}
template<typename T>
void voxelImageT<T>::AND(const voxelImageT& data2)
{
	forAllkji((*this))
				(*this)[k][j][i]= (*this)[k][j][i] && data2[k][j][i];
}
template<typename T>
void voxelImageT<T>::OR(const voxelImageT& data2)
{
	forAllkji((*this))
				(*this)[k][j][i]= (*this)[k][j][i] || data2[k][j][i];
}

template<typename T>
void voxelImageT<T>::XOR(const voxelImageT& data2)
{
	forAllkji((*this))
		(*this)[k][j][i]= (*this)[k][j][i] != data2[k][j][i];
}


template<typename T>
void voxelImageT<T>::threshold101(T theresholdMin,T  theresholdMax)
{
	forAllkji((*this))
		(*this)[k][j][i]=	( (*this)[k][j][i] < theresholdMin )  || ( (*this)[k][j][i] > theresholdMax  );

}

template<typename T>
void voxelImageT<T>::segmentPhase
(
	T StartThreshold1,
	T StartThreshold2,
	T finalThreshold1,
	T finalThreshold2,
	unsigned int nIter
)
{
	voxelImageT<T> tmpFinal=*this;
	this->threshold101(StartThreshold1, StartThreshold2);
	this->fillHoles(2);
	voxelImageT<T> tmpStart=*this;

	tmpFinal.threshold101(finalThreshold1, finalThreshold2);
	tmpFinal.FaceMedian06(2,4);
	tmpFinal.FaceMedian06(2,4);
	tmpFinal.FaceMedian06(2,4);
	for ( unsigned int i=0 ; i < nIter ; ++i )
	{ this->growPore();  this->OR(tmpFinal); }

	this->AND(tmpStart);
}




template<typename T>
void voxelImageT<T>::fillHoles(unsigned int maxHoleRadius)
{
	std::cout<<"  filling small isolated parts: "<<std::flush;
	voxelImageT<T> dataTmp=*this;
		std::cout<<"-"<<std::flush;

	dataTmp.shrinkPore(); std::cout<<".";std::cout.flush();
	for ( unsigned int i=0 ; i < 6 ; ++i )
		{ dataTmp.growPore();  dataTmp.OR(*this); std::cout<<".";std::cout.flush();}
	*this=dataTmp;
	std::cout<<"-"<<std::flush;

	dataTmp.growPore();
	for ( unsigned int i=0 ; i < 4 ; ++i )
		{ dataTmp.shrinkPore(); dataTmp.AND(*this); std::cout<<".";std::cout.flush();}
	*this=dataTmp;
	std::cout<<"-"<<std::flush;

	if ( maxHoleRadius > 1)
	{
		for ( unsigned int i=0 ; i < maxHoleRadius ; ++i )
			{ dataTmp.shrinkPore(); std::cout<<".";std::cout.flush();}
		for ( unsigned int i=0 ; i < maxHoleRadius*6 ; ++i )
			{ dataTmp.growPore();  dataTmp.OR(*this); std::cout<<".";std::cout.flush();}
		*this=dataTmp;
		std::cout<<"-"<<std::flush;

		for ( unsigned int i=0 ; i < maxHoleRadius ; ++i )
			{ dataTmp.growPore(); std::cout<<".";std::cout.flush();}
		for ( unsigned int i=0 ; i < maxHoleRadius*4 ; ++i )
			{ dataTmp.shrinkPore(); dataTmp.AND(*this); std::cout<<".";std::cout.flush();}
		*this=dataTmp;
		std::cout<<"-"<<std::flush;
	}
	std::cout<<"."<<std::endl;

}

template<typename T>
void voxelImageT<T>::writeAConnectedPoreVoxel(std::string fileName) const
{

	std::cout<<" finding a connected pore voxel:";
	//~ bool foundAConnectedPoreVoxel=false;
	for ( int nShrink=4; nShrink>=0  ; nShrink-- )
	{
		voxelImageT<T> dataTmp=*this;
		for ( int iShrink=1; iShrink<=nShrink ; iShrink++ )
		{
			dataTmp.shrinkPore();
		}

		std::cout<<" nShrink: "<<nShrink<<std::endl;;
		dataTmp.printInfo();

		//~ dataTmp.fillHoles();
		//~ dataTmp.AND(*this);

		for (unsigned int k=dataTmp.size()*3/4-2; k>dataTmp.size()*1/8+1  ; k-- )
		{
			for (unsigned int j=dataTmp[k].size()*3/4-2; j>dataTmp[k].size()*1/8+1  ; j-- )
			{
				for (unsigned int i=dataTmp[k][j].size()*3/4-2; i>dataTmp[k][j].size()*1/8+1  ; i-- )
				{

					if (dataTmp[k][j][i]==0 && !(															  dataTmp[k+1][j+1][i+1]
							|| dataTmp[k-1][j-1][i-1] || dataTmp[k-1][j-1][i+1] ||  dataTmp[k-1][j+1][i-1] || dataTmp[k-1][j+1][i+1]
							|| dataTmp[k+1][j-1][i-1] || dataTmp[k+1][j-1][i+1] ||  dataTmp[k+1][j+1][i-1]
							|| dataTmp[k][j][i-1] || dataTmp[k][j][i+1]
							|| dataTmp[k][j-1][i] || dataTmp[k][j+1][i]
							|| dataTmp[k-1][j][i] || dataTmp[k+1][j][i]
								  )
						)
					{
						std::ofstream of(fileName.c_str());
						assert(of);
						of<<(i+0.5)*dx_[0]+X0_[0]<<" "<<(j+0.5)*dx_[1]+X0_[1]<<" "<<(k+0.5)*dx_[2]+X0_[2]<<std::endl;
						of.close();
						std::cout<<" found  ("<<i<<" "<<j<<" "<<k<<") -> "<<double(dataTmp[k][j][i])<<std::endl;
						std::cout<<" found  xyz:  "<<i*dx_[0]+X0_[0]<<" "<<j*dx_[1]+X0_[1]<<" "<<k*dx_[2]+X0_[2]<<std::endl;
						return ;
					}
				}
			}
		}

	}

	std::cout<<" \n			  ----  ERRORR   -----	  \n"<<std::endl;

	std::cout<<" \n\n ----  didn't find  a connected pore voxel   -----  \n\n"<<std::endl;

}



template<typename T>
void voxelImageT<T>::printInfo() const
{
	unsigned long long nPores=0;
	int3 n=this->size3();
	std::cout<<"calculating image porosity:"<<std::endl;
	forAllkji((*this))
				if ((*this)[k][j][i]==0 )
					nPores++;

	std::cout << " total porosity: " << nPores<<"/"<<"("<<n[0]<<"*"<<n[1]<<"*"<<n[2]<<") = "<< double(nPores)/(n[0]*n[1]*n[2]) << std::endl;

}

template<typename T>
double voxelImageT<T>::volFraction(T vv1,T vv2) const
{
	unsigned long long nPores=0;
	int3 n=this->size3();
	for ( unsigned int k=0; k<this->size() ; k++ )
		for ( unsigned int j=0; j<(*this)[k].size() ; ++j )
			for ( unsigned int i=0; i<(*this)[k][j].size() ; ++i )
						if ( vv1<=(*this)[k][j][i] && (*this)[k][j][i]<=vv2 )
							nPores++;

	return double(nPores)/(n[0]*n[1]*n[2]);

}


template<typename T>
void voxelImageT<T>::writeHeader(std::string outputName) const
{
	this->voxelImageT<T>::writeHeader(outputName, int3{{0,0,0}}, voxelField<T>::size3());
}

template<typename T>
void voxelImageT<T>::writeHeader(std::string outName, int3 iStart, int3 iEnd) const
{
	vec3 xmin(X0_+(iStart*dx_));
	int3	n(iEnd-iStart);
	if (outName.size()<7 || outName.compare(outName.size()-7,7,"_header")!=0)
	{
		int islash=outName.find_last_of("\\/"); if (islash>=int(outName.size())) islash=-1;
		std::string title=outName.substr(islash+1);
		if (outName.size()>=4 && outName.compare(outName.size()-4,4,".mhd")==0)
			title=title.substr(0,title.size()-4)+".raw";
		else
			outName=outName.substr(0,outName.find_last_of("."))+".mhd";

		std::string typeNmeVTK="MET_UCHAR";
		if (typeid(T)==typeid(char)) typeNmeVTK="MET_CHAR";
		else if (typeid(T)==typeid(short)) typeNmeVTK="MET_SHORT";
		else if (typeid(T)==typeid(unsigned short)) typeNmeVTK="MET_USHORT";
		else if (typeid(T)==typeid(int)) typeNmeVTK="MET_INT";
		else if (typeid(T)==typeid(unsigned int)) typeNmeVTK="MET_UINT";
		else if (typeid(T)==typeid(float)) typeNmeVTK="MET_FLOAT";
		else if (typeid(T)==typeid(double)) typeNmeVTK="MET_DOUBLE";

		std::ofstream outputHeaderFile(outName.c_str());
		assert(outputHeaderFile);
		outputHeaderFile
			 <<"ObjectType =  Image"<<std::endl
			 <<"NDims =	   3"<<std::endl
			 <<"ElementType = "<<typeNmeVTK<<std::endl <<std::endl
			 <<"DimSize =		"<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl
			 <<"ElementSpacing = "<<dx_[0]<<" "<<"  " <<dx_[1]<<" "<<"  " <<dx_[2]<<std::endl
			 <<"Offset =		 "<<X0_[0]<<" "<<"  " <<X0_[1]<<" "<<"  " <<X0_[2]<<std::endl <<std::endl
			 <<"ElementDataFile = "<<title<<std::endl <<std::endl
			 <<std::endl;
	}
	else
	{
	  std::ofstream outputHeaderFile(outName.c_str());					 // file pointer
	  assert(outputHeaderFile);
	  outputHeaderFile
		 <<"Nxyz"<<std::endl
		 <<"dxX0"<<std::endl
		 <<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl
		 <<dx_[0]<<" "<<"  " <<dx_[1]<<" "<<"  " <<dx_[2]<<std::endl
		 <<xmin[0]<<" "<<"  " <<xmin[1]<<" "<<"  " <<xmin[2]<<std::endl
		 <<"\n\nComments:"<<std::endl
		 <<" first 9 entries above are:"<<std::endl
		 <<"	Nx Ny Nz"<<std::endl
		 <<"	dx dy dz"<<std::endl
		 <<"	Xo Yo Zo"<<std::endl
		 <<" Nx, Ny and Nz  count for the number of columns, rows and layers respectively as written in the file"<<std::endl
		 <<" Optional keywords (move above Comments to activate):"<<std::endl
		 <<"	crop		0  299   0  299   0  299 "<<std::endl
		 <<"	pore 		0 0 "<<std::endl
		 <<"	resample	1"<<std::endl
		 <<"	direction	z"<<std::endl
		 <<"	..... "<<std::endl
		 <<std::endl;
	}
}


template<typename T>
void voxelImageT<T>::write(std::string outName)
{
	if (outName.compare(outName.size()-4,4,".mhd") == 0 || outName.compare(outName.size()-4,4,".raw") == 0)
	{
		this->writeBin(outName.substr(0,outName.size()-4)+".raw");
		writeHeader(outName);
	}
	else
	{
		this->writeAscii(outName);
		writeHeader(outName+"_header");
	}

}

