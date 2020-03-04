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


#ifndef voxelMesh_H
#define voxelMesh_H

    #include <fstream>
    #include <iostream>
    #include <vector>
	#include <valarray>
	
    #include <cassert>
    
template<typename Type>
class threeD: public std::valarray<Type>
{
  public:
	threeD(): std::valarray<Type>(3) {};
	threeD(Type x, Type y, Type z): std::valarray<Type>(3) {(*this)[0] = x;(*this)[1] = y;(*this)[2] = z;};
	threeD(Type x[3]): std::valarray<Type>(3) {(*this)[0] = x[0];(*this)[1] = x[1];(*this)[2] = x[2];};
	threeD(std::valarray<Type> x): std::valarray<Type>(x) {};
	void operator = (std::valarray<Type> x) {std::valarray<Type>::operator = (x); };
};

 
typedef   threeD<long long> Long3D;//with fancy constructors
typedef   threeD<int> Int3D;//with fancy constructors
typedef   threeD<double> Double3D;//with fancy constructors

//after initialization, switch back to vallarray
#define int3D  std::valarray<int> 
#define double3D  std::valarray<double> 

inline Double3D operator*( int3D n, Double3D x){ Double3D res; for (int i=0;i<3;i++){res[i]=n[i]*x[i];}; return res;} 

template<typename T1>
inline Double3D operator*( T1 n, double3D x){ Double3D res; for (int i = 0;i<3;i++){res[i] = n[i]*x[i];}; return res;}

template<typename T2>
inline Double3D operator/( Long3D n, T2 x){ Double3D res; for (int i = 0;i<3;i++){res[i] = n[i]/x;}; return res;}



template <typename Type> class voxelField
: public std::vector<std::vector<std::vector<Type> > >
{
  protected:
  public:
    voxelField(){};
    voxelField(unsigned int n1, unsigned int n2, unsigned int n3, Type value) ;
    voxelField(const voxelField<Type>& vm ): std::vector<std::vector<std::vector<Type> > >(vm) {} ;
    void reset(unsigned int n1,unsigned int n2,unsigned int n3, Type value);
    void readMicroCT(std::string); 
    void readMicroCTHeader(std::ifstream); 
    bool read(std::string);
    void read(std::ifstream& in);
    bool readBin(std::string fileName);
    bool readBin(std::string fileName,int iStart,int iEnd , int jStart,int jEnd , int kStart,int kEnd );
    void writeBin(std::string fileName) const;
    void writeBin(std::string fileName,int iStart,int iEnd , int jStart,int jEnd , int kStart,int kEnd ) const;
    void writeData(std::ofstream& of) const;
    void writeRotatedXZ(std::ofstream& of) const; 
    void write(std::string) const;
    Int3D size3D() const;
    void getSize(unsigned int& n1, unsigned int& n2, unsigned int& n3) const;
    void getPrivateSize(unsigned int& n1, unsigned int& n2, unsigned int& n3) const;

    //~ void operator = (voxelField<Type>);
};





template<typename Type>
voxelField<Type>::voxelField(unsigned int n1, unsigned int n2, unsigned int n3, Type value)
{
    reset(n1,n2,n3, value);
}



//kji-----check needed
template<typename Type>   void voxelField<Type>::reset(unsigned int n1,unsigned int n2,unsigned int n3, Type value)
{
    //~ this->resize(n3,std::vector<std::vector<Type> >(n2,std::vector<Type>(n1,value)));
    this->resize(0);
    this->resize(n3,std::vector<std::vector<Type> >(n2,std::vector<Type>(n1,value)));
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
Int3D voxelField<Type>::size3D() const
{
  Int3D mySize;
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
template<typename Type>   void voxelField<Type>::read(std::ifstream& in)
{
	//~ register double tmp = 0;
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
                //~ *d3 = tmp;
            }
        }
    }
}

//read order sensitive
template<typename Type>   bool voxelField<Type>::read(std::string fileName)
{
    std::cout<<  " reading "<<fileName<<std::endl;
    std::ifstream in(fileName.c_str());
    if(!in)  
		{std::cout<<"\n\nError: can not open image file, "<<fileName<<std::endl<<std::endl; 
		return false;}

    
    read(in);
    
    in.close();
    return true;
}




//read order sensitive
template<typename Type>   void voxelField<Type>::readMicroCT(std::string fileName)
{
    std::cout<<  " reading "<<fileName<<std::endl;
    std::ifstream in(fileName.c_str());
    if(!in)  {std::cout<<"\n\nError: can not open image file, "<<fileName<<std::endl<<std::endl; exit(-1);}

    char tmpc;
    for ( unsigned int i = 0; i<8;i++)   in>>tmpc, std::cout<<" "<<tmpc;  //ignore the first 8 characters (ascii 3uc)
	unsigned int n[3];
	double  xmin[3];
	double  xmax[3];
	
    in>>n[2]>>n[1]>>n[0];                        // number of variables (dimension of
    //~ in>>    dx[0]>>dx[1]>>dx[2] ;
    in>>    xmin[0]>>xmin[1]>>xmin[2] ;
    in>>    xmax[0]>>xmax[1]>>xmax[2] ;

    read(in);
    
    in.close();
}

//kji
template<typename Type>   bool voxelField<Type>::readBin(std::string fileName)
{
		int3D n = size3D();
        return readBin(fileName,0,n[0]-1,0,n[1]-1,0,n[2]-1);
}

//kji
template<typename Type>   bool voxelField<Type>::readBin(std::string fileName,
                                int iStart,int iEnd ,
                                int jStart,int jEnd ,
                                int kStart,int kEnd
)
{
    std::cout<<  " reading "<<fileName;
    std::ifstream in (fileName.c_str(), std::ios::in | std::ios::binary);
    if(!in)  {std::cout<<"\n\n  Error: can not open image file, "<<fileName<<std::endl<<std::endl;
		return false;}

    for ( int k = kStart;k <= kEnd;k++)
    {
        for ( int j = jStart;j <= jEnd;j++)
        {

            if(!
                in.read((char*)(&(*this)[k][j][iStart]),
                (iEnd-iStart+1) * sizeof(Type) )
              )
            {
                    std::cout<<  "  error in reading "<<fileName<<std::endl; std::cout.flush();
            }
        }
    }
    in.close();
    std::cout<<  "."<<std::endl;
	return true;
}


template<typename Type>   void voxelField<Type>::writeBin(std::string fileName) const
{
	int3D imgsize = size3D();
	writeBin(fileName,0, imgsize[0]-1 ,0,imgsize[1]-1 ,0, imgsize[2]-1);
}



//kji
template<typename Type>   void voxelField<Type>::writeBin(std::string fileName,
                               int iStart,int iEnd , int jStart,int jEnd , int kStart,int kEnd ) const
{
    //~ std::cout<<  "writting "<<' '<<(*this).size()<<' '<<(*this)[400].size()<<' '<<(*this)[0][0].size()<<' '<<std::endl;

    std::cout<<  " writting binary file "<<fileName<<";  ";
    std::cout<<  "  i: "<<iStart<<" "<<iEnd<<",  ";
    std::cout<<  "  j: "<<jStart<<" "<<jEnd<<",  ";
    std::cout<<  "  k: "<<kStart<<" "<<kEnd; std::cout.flush();
    //~ std::cout<<  "  size: "<<kEnd-kStart+1<<" layers of "<<iEnd-iStart+1<<" by "<<jEnd-jStart+1<<" "<<std::endl;

    std::ofstream of (fileName.c_str(), std::ios::out | std::ios::binary);
    assert(of);


    for ( int k = kStart;k <= kEnd;k++)
    {
        for ( int j = jStart;j <= jEnd;j++)
        {
            if (!(*this)[k][j].empty())
                of.write((char*)(&((*this)[k][j][iStart])),
                (iEnd-iStart+1) * sizeof(Type));
            else
                std::cout<<"\n  Error when writting binary file, k:"<<k<<  "j:"<<k<<" "<<std::endl;


        }

    }
    of.flush();
    of.close();
    std::cout<<  "."<<std::endl;

}


template<typename Type>   void voxelField<Type>::writeData(std::ofstream& of) const
{
    for (
    typename std::vector<std::vector<std::vector<Type> > >::const_iterator d1 = this->begin();
    d1<this->end();
    d1++)
    {
        for ( typename std::vector<std::vector<Type> >::const_iterator d2 = d1->begin();d2<d1->end();d2++)
        {
            for ( typename std::vector<Type>::const_iterator d3 = d2->begin();d3<d2->end();d3++)
            {
                of<<double(*d3)<<' ';
            }
            of<<"\n";
        }
    }
    of<<std::endl;
}

template<typename Type>   void voxelField<Type>::writeRotatedXZ(std::ofstream& of) const
{
   
	for (unsigned int i = 0;i<(*this)[0][0].size();i++) //reversed order with k
    {
        for (unsigned int j = 0;j<(*this)[0].size();j++)
        {
            for (unsigned int k = 0;k<(*this).size();k++) //reversed order with i
			{
                of<<(double)(*this)[k][j][i]<<' ';
			}
            of<<std::endl;
        }

    }
}

template<typename Type>   void voxelField<Type>::write(std::string fileName) const
{

    std::ofstream of(fileName.c_str());
    assert(of);
    writeData(of);
    of.close();


}


class voxelMesh: public voxelField<unsigned char>
{

    void Gauss(unsigned int  k);
    //~ Int3D  nMin_;
    Double3D    X0_, dx_;

  public:

    voxelMesh(): X0_(0.0,0.0,0.0),dx_(1,1,1) {};
    //~ voxelMesh(unsigned int n1, unsigned int n2, unsigned int n3, unsigned char value) ;
    //~ voxelMesh(unsigned int n[3], double dx[3], double xmin[3], unsigned char value) ;
    //~ voxelMesh(const voxelMesh & ) ;




    voxelMesh(unsigned int n1, unsigned int n2, unsigned int n3, unsigned char value) 
: voxelField<unsigned char>( n1,  n2,  n3,  value),  X0_(0.0,0.0,0.0),dx_(1,1,1) {}


    voxelMesh(unsigned int n[3], double dx[3], double xmin[3], unsigned char value)
    : voxelField<unsigned char>( n[0],  n[1],  n[2],  value), X0_(xmin),dx_(dx) {}

    voxelMesh(const voxelMesh & vm) 
:  voxelField<unsigned char>(vm), X0_(vm.X0_),dx_(vm.dx_){}




void     read(std::ifstream& in)
{
	int tmp=0;
    for (
     std::vector<std::vector<std::vector<unsigned char> > >::iterator d1=this->begin();
    d1<this->end();
    d1++)
    {
        for (  std::vector<std::vector<unsigned char> >::iterator d2=d1->begin();d2<d1->end();d2++)
        {
            for (  std::vector<unsigned char>::iterator d3=d2->begin();d3<d2->end();d3++)
            {
                in>>tmp;
                *d3=tmp;
            }
        }
    }
}

void     read(std::string fileName)
{
//  overwrite as the parent -voxelField<unsigned char>- interprets
//  numerical values as characters not integers

    std::cout<<  " reading "<<fileName<<std::endl;


	if ( (fileName.compare(fileName.size()-4,4,".dat")==0) || (fileName.compare(fileName.size()-4,4,".txt") == 0) )     
	{
		std::ifstream in;
		in.open(fileName.c_str());
		assert(in);
		
		char tmpc[8];
		for ( unsigned int i=0; i<8;i++)   in>>tmpc[i];  
		if (std::string(tmpc).compare(0,4,"ascii") == 0) //ignore first lines
		{
			int n[3];
			in>>n[2]>>n[0]>>n[1];//ignore first lines
			double  xmin[3],xmax[3];
			in>>    xmin[0]>>xmax[0]>>xmin[1]>>xmax[1]>>xmin[2]>>xmax[2] ;
			std::cout<<"Warning: ignoring the header of file "<<fileName<<std::endl;
		}
		else
		{
			in.seekg(0, std::ios::beg);
		}
		read(in);
		in.close();	
    }
    else
    {
        readBin(fileName);
	}

}








    void crop( unsigned    int cropBegin[3],  unsigned    int cropEnd[3]) ;
    void crop(unsigned int iStart, unsigned int iEnd ,
                unsigned int jStart, unsigned int jEnd ,
                unsigned int kStart, unsigned int kEnd ,
                unsigned int emptylayers,unsigned char emptylayersValue);

    void writeHeader(std::string fileName) const;
	void writeHeader(std::string fileName, int3D iStart, int3D iEnd) const;


    void erodeLayer(int i);
    void resample(double i);
    void rotate(char direction);
    void PointMedian(unsigned int thereshold0,unsigned int thereshold1);
    void FaceMedian(unsigned int thereshold0,unsigned int thereshold1);
       
    void AND(const voxelMesh& data2);
    void OR(const voxelMesh& data2);
    void XOR(const voxelMesh& data2);
    
    void threshold101(unsigned char theresholdMin,unsigned char theresholdMax);

    void Gauss();

    void fillHoles(unsigned int maxHoleRadius);

    void shrinkPore();
    void growpore();
	//void ReplaceVoxels(int nVoxels);
	//template<typename Type> 	void replaceValue(voxelMesh& vImage, Type v, Type newv);
	//void ParticleLabel();
	
	void median(short nNeist);
	void  segmentPhase
	(
		unsigned char StartThreshold1,
		unsigned char StartThreshold2, 
		unsigned char finalThreshold1, 
		unsigned char finalThreshold2,
		unsigned int nIter
	);
    //~ void Segment3Phase
    //~ (
      //~ char StartThresholdW1,char StartThresholdW2, char finalThresholdW1 , char finalThresholdW2 
      //~ char StartThresholdO1,char StartThresholdO2, char finalThresholdO1 , char finalThresholdO2 
      //~ char StartThresholdG1,char StartThresholdG2, char finalThresholdG1 , char finalThresholdG2 
      //~ char StartThresholdR1,char StartThresholdR2, char finalThresholdR1 , char finalThresholdR2 
    //~ );

    //~ void addX(int k, int distance);
    //~ void writeSTLBINARY(std::string fileName) const;
    void writeAConnectedPoreVoxel(std::string fileName) const;
    void setLayer(int k, std::vector<std::vector<unsigned char> > &Values);
    void replaceyLayer(int j, int fromj);
    void replacexLayer(int i, int fromi);
    void setBlock(int n1, int n2, int n3, std::vector<std::vector<std::vector<unsigned char> > >&Values);
	void growBox(unsigned int nLayers);

	//~ void read(std::ifstream& in);
    //~ void read(std::string fileName);
    void write(std::string fileName);
    void printInfo() const;
    
    Double3D X0() const {return X0_;};
    Double3D dx() const  {return dx_;};
    //~ Int3D nMin() const {return nMin_;};

 };





void  readFromHeader
(voxelMesh& voxelImage, std::string headerName, std::string inputName = "" );



#endif
