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
#include <map>
#include <algorithm>


using namespace std;

//~ template<typename Type>
//~ void voxelField<Type>::operator=(voxelField<Type> VF)
//~ {
    //~ this->vector<vector<vector<Type> > >::operator=(VF);
//~ }

//read order sensitive

void voxelMesh::crop( unsigned int cropBegin[3],  unsigned int cropEnd[3])
{
	crop(cropBegin[0],cropEnd[0],cropBegin[1],cropEnd[1],cropBegin[2],cropEnd[2],0,0);
}
//kji
void voxelMesh::crop(
                            unsigned int iStart, unsigned int iEnd ,
                            unsigned int jStart, unsigned int jEnd ,
                            unsigned int kStart, unsigned int kEnd  ,
                            unsigned int emptylayers, unsigned char emptylayersValue
)
{
    cout<<  " Cropping, ";
    cout<<  "  i: "<<iStart<<" "<<iEnd <<  "  j: "<<jStart<<" "<<jEnd<<  "  k: "<<kStart<<" "<<kEnd<<endl;
    cout<<  "  adding "<<emptylayers<<" layers with value "<< (double) emptylayersValue <<endl;
    
	//~ nMin_[0]=nMin_[0]+iStart-emptylayers;   nMin_[1]=nMin_[1]+jStart-emptylayers;   nMin_[2]=nMin_[2]+kStart-emptylayers;
	X0_[0]=X0_[0]+(int(iStart)-int(emptylayers))*dx_[0];   X0_[1]=X0_[1]+(int(jStart)-int(emptylayers))*dx_[1];   X0_[2]=X0_[2]+(int(kStart)-int(emptylayers))*dx_[2];
	
    //cout<<  " -> ";
    //cout<<  "  nmin: "<<nMin_[0]<<" "<<nMin_[1] <<  " "<<nMin_[2]<<endl;
   // cout<<  "  dx: "<<dx_[0]<<" "<<dx_[1] <<  " "<<dx_[2]<<endl;
    //cout<<  "  X0: "<<X0_[0]<<" "<<X0_[1] <<  " "<<X0_[2]<<endl;
    
    voxelMesh tmp=*this;
 


    //~ reset((*this)[k].size(),int n2,int n3, Type value)
    this->resize(kEnd+1-kStart+2*emptylayers);
    for ( unsigned int k=0; k<kEnd+1-kStart+2*emptylayers ; k++ )
    {
        (*this)[k].resize(jEnd+1-jStart+2*emptylayers);
        for ( unsigned int j=0; j<jEnd+1-jStart+2*emptylayers; j++ )
        {
            (*this)[k][j].resize(iEnd+1-iStart+2*emptylayers,emptylayersValue);
            for ( unsigned int i=0; i<iEnd-iStart+1+2*emptylayers; i++ )
            {
                (*this)[k][j][i]=emptylayersValue;
            }
        }
    }

    for ( unsigned int k=0; k<kEnd+1-kStart; k++ )
    {
        for ( unsigned int j=0; j<jEnd+1-jStart; j++ )
        {
            for ( unsigned int i=0; i<iEnd+1-iStart; i++ )
            {
                (*this)[k+emptylayers][j+emptylayers][i+emptylayers]=tmp[k+kStart][j+jStart][i+iStart];
            }
        }
    }

        //~ for (unsigned int i=0 ;i<(*this)[k].size(); i++)
        //~ {
            //~ (*this)[k][i]=Values[i];
        //~ }
        
        
}




void voxelMesh::setLayer(int k, vector<vector<unsigned char> > &Values)
{
        //~ for (unsigned int i=0 ;i<(*this)[k].size(); i++)
        //~ {
            (*this)[k]=Values;
        //~ }
}

void voxelMesh::replacexLayer(int i, int fromi)
{

     for ( unsigned int k=0; k<(*this).size() ; k++ )
    {
        for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
        {
            //~ for ( unsigned int i=0; i<Values[k][j].size() ; i++ )
            //~ {
                (*this)[k][j][i]=(*this)[k][j][fromi];
            //~ }
        }
    }

}
void voxelMesh::replaceyLayer(int j, int fromj)
{

     for ( unsigned int k=0; k<(*this).size() ; k++ )
    {
        //~ for ( unsigned int j=0; j<Values[k].size() ; j++ )
        //~ {
            for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
            {
                (*this)[k][j][i]=(*this)[k][fromj][i];
            }
        //~ }
    }

}
//kji
void voxelMesh::growBox(unsigned int nLayers)
{

	
	Int3D n = (*this).size3D();
    (*this).crop(0,n[0]-1,0,n[1]-1,0,n[2]-1, nLayers,1);//         XXXXXXXXXXXXXXXXXXXXXXXXXXXX

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

//kji
void voxelMesh::setBlock(int n1, int n2, int n3, vector<vector<vector<unsigned char> > >&Values)
{

     for ( unsigned int k=0; k<Values.size() ; k++ )
    {
        for ( unsigned int j=0; j<Values[k].size() ; j++ )
        {
            for ( unsigned int i=0; i<Values[k][j].size() ; i++ )
            {
                (*this)[k+n3][j+n2][i+n1]=Values[k][j][i];
            }
        }
    }

}



#define forAllNei(r1,r2) \
for (int k_nei_m=r1;k_nei_m<=r2;++k_nei_m) \
for (int j_nei_m=r1;j_nei_m<=r2;++j_nei_m) \
for (int i_nei_m=r1;i_nei_m<=r2;++i_nei_m)

#define nei(dataM,i_M,j_M,k_M) dataM[k_M+k_nei_m][j_M+j_nei_m][i_M+i_nei_m]




void voxelMesh::resample(double nReSampleNotSafe)//  TODO to be tested
{
    if (nReSampleNotSafe < .999)
    {
        int nReSample=1.0/nReSampleNotSafe+0.5;
        voxelMesh tmp=*this;
        //~ reset((*this)[k].size(),int n2,int n3, Type value)
        unsigned int n3=this->size(),n2=(*this)[0].size(),n1=(*this)[0][0].size();
        this->resize(nReSample*n3);
        for ( unsigned int k=0; k<this->size() ; k++ )
        {

            (*this)[k].resize(nReSample*n2);
            for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
            {
                (*this)[k][j].resize(nReSample*n1);
                for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
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
        voxelMesh tmp=*this;
        //~ reset((*this)[k].size(),int n2,int n3, Type value)
        unsigned int n3=this->size(),n2=(*this)[0].size(),n1=(*this)[0][0].size();
        this->resize((n3)/nReSample);
        for ( unsigned int k=0; k<this->size() ; k++ )
        {

            (*this)[k].resize((n2)/nReSample);
            for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
            {
                (*this)[k][j].resize((n1)/nReSample);
                for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
                {
                    //~ int neiSum=0;
                    //~ forAllNei(0,nReSample-1)
                    //~ {
                        //~ neiSum+=nei(tmp,i*nReSample,j*nReSample,k*nReSample);
                    //~ } 
                    //~ (*this)[k][j][i]=(0.5+double(neiSum)/(nReSample*nReSample*nReSample));
                       (*this)[k][j][i]=tmp[0.5+k*nReSample][0.5+j*nReSample][0.5+i*nReSample];


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



void voxelMesh::rotate(char direction)
{// wrong X0	
	unsigned int n1,n2,n3;

	getSize(n1,n2,n3);
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

		voxelMesh tmp=*this;
		this->reset(n3,n2,n1,0);
		for ( unsigned int k=0; k<n3 ; k++ )
		{
			for ( unsigned int j=0; j<n2 ; j++ )
			{
				for ( unsigned int i=0; i<n1 ; i++ )
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
		
		voxelMesh tmp=*this;
		this->reset(n2,n1,n3,0);
		for ( unsigned int k=0; k<n3 ; k++ )
		{
			for ( unsigned int j=0; j<n2 ; j++ )
			{
				for ( unsigned int i=0; i<n1 ; i++ )
				{
					(*this)[k][i][j]=tmp[k][j][i];
				}
			}
		}
	}
	else if (direction=='-') 
	{
		cout<<" -> flipping image,  x origin will be invalid "<<endl;
		voxelMesh tmp=*this;
		for ( unsigned int k=0; k<n3 ; k++ )
		{
			for ( unsigned int j=0; j<n2 ; j++ )
			{
				for ( unsigned int i=0; i<n1 ; i++ )
				{
					(*this)[k][j][n1-1-i]=tmp[k][j][i];
				}
			}
		}
		
	}
	else
	{
		cout<<"\n\nSwapping "<<direction<<" and x directions(!?!), sorry can't do that >-( "<<endl;
		cerr<<"Swapping "<<direction<<" and x directions(!?!), sorry can't do that >-( \n\n"<<endl;
	}

}


void voxelMesh::PointMedian(unsigned int thereshold0,unsigned int thereshold1)
{
    voxelMesh voxls=(*this);
    //~ // reset((*this)[k].size(),int n2,int n3, Type value)
    for ( unsigned int k=2; k<this->size()-2 ; k++ )
    {
        for ( unsigned int j=2; j<(*this)[k].size()-2 ; j++ )
        {
            for ( unsigned int i=2; i<(*this)[k][j].size()-2 ; i++ )
            {
                unsigned int neiSum=0;
                forAllNei(-1,1)
                {
                    neiSum+=nei(voxls,i,j,k);
                }
                if (neiSum <= thereshold0 )
                    (*this)[k][j][i]=0;
                if   (neiSum >= thereshold1)
                    (*this)[k][j][i]=1;
           }
        }
    }
}





//~
//~ #define forAllFaceNei \/
//~ for (int k_nei_m=-1;k_nei_m<2;k_nei_m++) \/
//~ for (int j_nei_m=-1;j_nei_m<2;j_nei_m++) \/
//~ for (int i_nei_m=-1;i_nei_m<2;i_nei_m++)
//~
//~ #define nei(i,j,k) (*this)[k+k_nei_m][j+j_nei_m][i+i_nei_m]


void voxelMesh::FaceMedian(unsigned int thereshold0,unsigned int thereshold1)
{
    unsigned int n[3]={0,0,0};
    double doubletmp[3]={0,0,0};
    getSize(n[0],n[1],n[2]);
	for (unsigned int i=0;i<3;i++) n[i]=n[i]+2;
	
    voxelMesh voxls(n,doubletmp,doubletmp,1);
        voxls.setBlock(0, 0, 0, (*this));
        voxls.setBlock(2, 2, 2, (*this));
        voxls.setBlock(1, 1, 1, (*this));

    for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
    {
        for ( unsigned int j=1; j<voxls[k].size()-1 ; j++ )
        {
            for ( unsigned int i=1; i<voxls[k][j].size()-1 ; i++ )
            {
                 unsigned int neiSum = voxls[k][j][i-1]+voxls[k][j][i+1]
									  +voxls[k][j-1][i]+voxls[k][j+1][i]
                                      +voxls[k+1][j][i]+voxls[k-1][j][i];

                if (neiSum <= thereshold0 )
                    (*this)[k-1][j-1][i-1]=0;
                if (neiSum >= thereshold1)
                    (*this)[k-1][j-1][i-1]=1;
            }
        }
    }
}


void voxelMesh::shrinkPore() 
{
         //~ {  FaceMedian(0,1);  }
//~ return;
    voxelMesh voxls=*this;


    for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
    {
        for ( unsigned int j=1; j<voxls[k].size()-1 ; j++ )
        {
            for ( unsigned int i=1; i<voxls[k][j].size()-1 ; i++ )
            {

                if (voxls[k][j][i]==0 && ( voxls[k][j][i-1]  || voxls[k][j][i+1] || 
                                      voxls[k][j-1][i] || voxls[k][j+1][i] || voxls[k-1][j][i] || voxls[k+1][j][i] ) )
                    (*this)[k][j][i]=1;
           }
        }
    }


    for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
    {
        for ( unsigned int j=1; j<voxls[k].size()-1 ; j++ )
        {
            //~ for ( unsigned int i=0; i<voxls[k][j].size()-1 ; i++ )
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
        //~ for ( unsigned int j=0; j<voxls[k].size()-1 ; j++ )
        {
            for ( unsigned int i=1; i<voxls[k][0].size()-1 ; i++ )
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
        for ( unsigned int j=1; j<voxls[0].size()-1 ; j++ )
        {    	
            for ( unsigned int i=1; i<voxls[0][0].size()-1 ; i++ )
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

void voxelMesh::growpore() // optimized function, should be further optimized as it is frequently used
{
    
    voxelMesh voxls=*this;


    for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
    {
        for ( unsigned int j=1; j<voxls[k].size()-1 ; j++ )
        {
            for ( unsigned int i=1; i<voxls[k][j].size()-1 ; i++ )
            {

                if (voxls[k][j][i]==1 && ( voxls[k][j][i-1]+voxls[k][j][i+1] < 2 || 
                                      voxls[k][j-1][i]+voxls[k+1][j][i]+voxls[k][j+1][i]+voxls[k-1][j][i] < 4 ) )
                    (*this)[k][j][i]=0;
           }
        }
    }


    for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
    {
        for ( unsigned int j=1; j<voxls[k].size()-1 ; j++ )
        {
            // for ( unsigned int i=0; i<voxls[k][j].size()-1 ; i++ )
            {	unsigned int i=0;
                if (voxls[k][j][i]==1 && (voxls[k][j][i+1]+voxls[k][j-1][i] +voxls[k+1][j][i]+voxls[k][j+1][i]+voxls[k-1][j][i] < 5 ) )
                    (*this)[k][j][i]=0;
           
				i=voxls[k][j].size()-1;
                if (voxls[k][j][i]==1 && (voxls[k][j][i-1]+voxls[k][j-1][i] +voxls[k+1][j][i]+voxls[k][j+1][i]+voxls[k-1][j][i] < 5 ) )
                    (*this)[k][j][i]=0;
           }
        }
    }
    

    for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
    {
        // for ( unsigned int j=0; j<voxls[k].size()-1 ; j++ )
        {
            for ( unsigned int i=1; i<voxls[k][0].size()-1 ; i++ )
            {	unsigned int j=0;

                if (voxls[k][j][i]==1 && (voxls[k][j][i-1]+voxls[k][j][i+1] +voxls[k+1][j][i]+voxls[k][j+1][i]+voxls[k-1][j][i] < 5 ) )
                    (*this)[k][j][i]=0;
           
				j=voxls[k].size()-1;
                if (voxls[k][j][i]==1 && (voxls[k][j][i-1]+voxls[k][j][i+1]+voxls[k][j-1][i] +voxls[k+1][j][i]+voxls[k-1][j][i] < 5 ) )
                    (*this)[k][j][i]=0;
           }
        }
    }
    
    // for ( unsigned int k=0; k<voxls.size()-1 ; k++ )
    {
        for ( unsigned int j=1; j<voxls[0].size()-1 ; j++ )
        {
            for ( unsigned int i=1; i<voxls[0][0].size()-1 ; i++ )
            {	unsigned int k=0;

                if (voxls[k][j][i]==1 && (voxls[k][j][i-1]+voxls[k][j][i+1]
                                      +voxls[k+1][j][i]+voxls[k][j+1][i]+voxls[k][j-1][i] < 5 ) )
                    (*this)[k][j][i]=0;
           
				k=voxls.size()-1;
                if (voxls[k][j][i]==1 && (voxls[k][j][i-1]+voxls[k][j][i+1]+voxls[k][j-1][i]
                                      +voxls[k][j+1][i]+voxls[k-1][j][i] < 5 ) )
                    (*this)[k][j][i]=0;
           }
        }
    }
    
    
}



void voxelMesh::AND(const voxelMesh& data2)
{
    //~ // reset((*this)[k].size(),int n2,int n3, Type value)
    for ( unsigned int k=0; k<this->size() ; k++ )
    {
        for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
        {
            for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
            {
                (*this)[k][j][i]= (*this)[k][j][i] && data2[k][j][i];
            }
        }
    }
}
void voxelMesh::OR(const voxelMesh& data2)
{
    //~ // reset((*this)[k].size(),int n2,int n3, Type value)
    for ( unsigned int k=0; k<this->size() ; k++ )
    {
        for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
        {
            for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
            {
                (*this)[k][j][i]= (*this)[k][j][i] || data2[k][j][i];
            }
        }
    }
}

void voxelMesh::XOR(const voxelMesh& data2)
{
    //~ // reset((*this)[k].size(),int n2,int n3, Type value)
    for ( unsigned int k=0; k<this->size() ; k++ )
    {
        for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
        {
            for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
            {
                (*this)[k][j][i]= (*this)[k][j][i] != data2[k][j][i];
            }
        }
    }
}


void voxelMesh::threshold101(unsigned char theresholdMin,unsigned char  theresholdMax)
{
    for ( unsigned int k=0; k<this->size() ; k++ )
    {
        for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
        {
            for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
            {
                    (*this)[k][j][i]=
                    ( (*this)[k][j][i] < theresholdMin )  || ( (*this)[k][j][i] > theresholdMax  );
            }
        }
    }
}

void  voxelMesh::segmentPhase
(
	unsigned char StartThreshold1,
	unsigned char StartThreshold2, 
	unsigned char finalThreshold1, 
	unsigned char finalThreshold2,
	unsigned int nIter
)
{
	voxelMesh tmpFinal=*this;
	this->threshold101(StartThreshold1, StartThreshold2);
	this->fillHoles(2);
	voxelMesh tmpStart=*this;
	
	tmpFinal.threshold101(finalThreshold1, finalThreshold2);
	tmpFinal.FaceMedian(2,4);
	tmpFinal.FaceMedian(2,4);
	tmpFinal.FaceMedian(2,4);
	for ( unsigned int i=0 ; i < nIter ; i++ )
	{ this->growpore();  this->OR(tmpFinal); }
	
	this->AND(tmpStart); 
}






void voxelMesh::Gauss()
{
    //~ // reset((*this)[k].size(),int n2,int n3, Type value)
    for ( unsigned int k=2; k<this->size()-2 ; k++ )
    {
        Gauss(k);
    }
}
void voxelMesh::Gauss(unsigned int k)
{



        for ( unsigned int j=2; j<(*this)[k].size()-2 ; j++ )
        {
            for ( unsigned int i=2; i<(*this)[k][j].size()-2 ; i++ )
            {
                int neiSum=0;
                forAllNei(-1,1)
                {
                    neiSum+=nei((*this),i,j,k);
                }
                neiSum=(0.5+neiSum/27.0);
                (*this)[k][j][i]=neiSum;
            }
        }
}


void voxelMesh::fillHoles(unsigned int maxHoleRadius)
{
	cout<<"  filling small isolated parts: "<<flush;
	voxelMesh dataTmp=*this;
		cout<<"-"<<flush;
 
	dataTmp.shrinkPore(); cout<<".";cout.flush();
	for ( unsigned int i=0 ; i < 6 ; i++ )
		{ dataTmp.growpore();  dataTmp.OR(*this); cout<<".";cout.flush();}
	*this=dataTmp;
	cout<<"-"<<flush;

	dataTmp.growpore();
	for ( unsigned int i=0 ; i < 4 ; i++ )
		{ dataTmp.shrinkPore(); dataTmp.AND(*this); cout<<".";cout.flush();}
	*this=dataTmp;
	cout<<"-"<<flush;	
	
	if ( maxHoleRadius > 1)
	{	
		for ( unsigned int i=0 ; i < maxHoleRadius ; i++ )
			{ dataTmp.shrinkPore(); cout<<".";cout.flush();}
		for ( unsigned int i=0 ; i < maxHoleRadius*6 ; i++ )
			{ dataTmp.growpore();  dataTmp.OR(*this); cout<<".";cout.flush();}
		*this=dataTmp;
		cout<<"-"<<flush;

		for ( unsigned int i=0 ; i < maxHoleRadius ; i++ )
			{ dataTmp.growpore(); cout<<".";cout.flush();}
		for ( unsigned int i=0 ; i < maxHoleRadius*4 ; i++ )
			{ dataTmp.shrinkPore(); dataTmp.AND(*this); cout<<".";cout.flush();}
		*this=dataTmp;
		cout<<"-"<<flush;
	}
	cout<<"."<<endl;

}

void voxelMesh::writeAConnectedPoreVoxel(string fileName) const
{

        cout<<" finding a connected pore voxel:";
        //~ bool foundAConnectedPoreVoxel=false;
        for ( int nShrink=4; nShrink>=0  ; nShrink-- )
        {
            voxelMesh dataTmp=*this;
            for ( int iShrink=1; iShrink<=nShrink ; iShrink++ )
            {            
                dataTmp.shrinkPore(); 
            }

			cout<<" nShrink: "<<nShrink<<endl;;
			dataTmp.printInfo();
			
            //~ dataTmp.fillHoles();
            //~ dataTmp.AND(*this);
            
            for (unsigned int k=dataTmp.size()*3/4-2; k>dataTmp.size()*1/8+1  ; k-- )
            {
                for (unsigned int j=dataTmp[k].size()*3/4-2; j>dataTmp[k].size()*1/8+1  ; j-- )
                {
                    for (unsigned int i=dataTmp[k][j].size()*3/4-2; i>dataTmp[k][j].size()*1/8+1  ; i-- )
                    {

                        if (dataTmp[k][j][i]==0 && !(                                                              dataTmp[k+1][j+1][i+1] 
								|| dataTmp[k-1][j-1][i-1] || dataTmp[k-1][j-1][i+1] ||  dataTmp[k-1][j+1][i-1] || dataTmp[k-1][j+1][i+1] 
                                || dataTmp[k+1][j-1][i-1] || dataTmp[k+1][j-1][i+1] ||  dataTmp[k+1][j+1][i-1]
								|| dataTmp[k][j][i-1] || dataTmp[k][j][i+1] 
								|| dataTmp[k][j-1][i] || dataTmp[k][j+1][i] 
								|| dataTmp[k-1][j][i] || dataTmp[k+1][j][i] 
                                      )
                                      )
                        {
                            ofstream of(fileName.c_str());
                            assert(of);
                            of<<(i+0.5)*dx_[0]+X0_[0]<<" "<<(j+0.5)*dx_[1]+X0_[1]<<" "<<(k+0.5)*dx_[2]+X0_[2]<<endl;   
                            of.close();    
                            cout<<" found  ("<<i<<" "<<j<<" "<<k<<") -> "<<(double)dataTmp[k][j][i]<<endl;
                            cout<<" found  xyz:  "<<i*dx_[0]+X0_[0]<<" "<<j*dx_[1]+X0_[1]<<" "<<k*dx_[2]+X0_[2]<<endl;
                            return ;              
                        }
                    }
                }
            }
                
        }
        
        cout<<" \n              ----  ERRORR   -----      \n"<<endl;
        
        cout<<" \n\n ----  didn't find  a connected pore voxel   -----  \n\n"<<endl;

}


 
void voxelMesh::printInfo() const
{
	unsigned int nPores=0;
	unsigned int n[3]={0,0,0};
    getSize(n[0],n[1],n[2]);
    cout<<"calculating image porosity:"<<endl;
	for ( unsigned int k=0; k<this->size() ; k++ )
    {
        for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
        {
            for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
            {

                        if ((*this)[k][j][i]==(unsigned char)0 )
                        {
                            nPores++;              
                        }
                        if ((*this)[k][j][i]==(unsigned char)2 )
                        {
                            nPores++;              
                        }
            }
        }
    }
    cout << " total porosity: " << nPores<<"/"<<"("<<n[0]<<"*"<<n[1]<<"*"<<n[2]<<") = "<< ((double)nPores)/(n[0]*n[1]*n[2]) << endl;
	///Porosity_x.txt
	ofstream  ofPores("Porosity.txt");
	ofPores<<"Total porosity: "<<nPores<<"/"<<"("<<n[0]<<"*"<<n[1]<<"*"<<n[2]<<") = "<<((double)nPores)/(n[0]*n[1]*n[2])<<endl;

}


void voxelMesh::writeHeader(string outputName) const
{
	this->writeHeader(outputName, Int3D(0, 0, 0 ), size3D());
}

void voxelMesh::writeHeader(string outputName, int3D iStart, int3D iEnd) const
{

	//~ cout <<"1 "<<dx_[0]<<" "<<X0_[0]<<" "<<iStart[0]<<endl;
	//~ cout <<"1 "<<iEnd[0]<<"-"<<iStart[0]<<"="<<(iEnd-iStart)[0]<<endl;

	//~ Double3D dx=dx_*1.0e+6;
	Double3D xmin(X0_+(iStart*dx_));
	//~ cout <<"1.5 "<<endl;
	Int3D    n(iEnd-iStart);
	//~ cout <<"2 "<<endl;
    string headerName;
    if (outputName.compare(outputName.size()-7,7,"_header") == 0)
		headerName=outputName;
	else
		headerName=outputName+"_header";

    ofstream outputHeaderFile(headerName.c_str());                     // file pointer
    assert(outputHeaderFile);
    outputHeaderFile
         <<"Nxyz"<<endl  
         <<"dxX0"<<endl  
         <<n[0]<<" "<<n[1]<<" "<<n[2]<<endl
         //~ <<"0.0 "<<xmax[0]-xmin[0]<<"   0.0  "<<xmax[1]-xmin[1]<<"    0.0 "<<xmax[2]-xmin[2]<<endl            
         <<dx_[0]<<" "<<"  " <<dx_[1]<<" "<<"  " <<dx_[2]<<endl                    
         <<xmin[0]<<" "<<"  " <<xmin[1]<<" "<<"  " <<xmin[2]<<endl                    
         <<"\n\nComments:"<<endl  
         <<" first 9 entries above are:"<<endl  
         <<"    Nx Ny Nz"<<endl  
         <<"    dx dy dz"<<endl  
         <<"    Xo Yo Zo"<<endl  
         <<" Nx, Ny and Nz  count for the number of columns, rows and layers respectively as written in the file"<<endl         
         <<" Optional keywords (move above Comments to activate):"<<endl  
         <<"    crop	    0  299   0  299   0  299 "<<endl  
         <<"    pore 	    0 0 "<<endl  
         <<"    resample	1"<<endl         
         <<"    ..... "<<endl         
         //~ <<"\n image absolute bounds:"<<endl         
         //~ <<xmin[0]<<" "<<xmax[0]<<"  " <<xmin[1]<<" "<<xmax[1]<<"  " <<xmin[2]<<" "<<xmax[2]<<endl            
         <<endl;
}


void voxelMesh::write(string outputName)
{

	unsigned int n[3];
	getSize(n[0],n[1],n[2]);

	if (outputName.compare(outputName.size()-4,4,".raw") == 0)
	{
		writeBin(outputName,0,n[0]-1,0,n[1]-1,0,n[2]-1);
		writeHeader(outputName+"_header");

	}
	else if (outputName.compare(outputName.size()-9,9,"_newf.dat") == 0)
	{		///for Branko

		ofstream rockFile((outputName).c_str());                   
		assert(rockFile);
		writeRotatedXZ(rockFile);
		rockFile.close();

		ofstream outputHeaderFile((outputName+"_header").c_str());                     // file pointer
		assert(outputHeaderFile);
		outputHeaderFile
		 <<"Nxyz"<<endl  
		 <<"dxX0"<<endl  
		 <<n[2]<<" "<<n[1]<<" "<<n[0]<<endl
		 //~ <<"0.0 "<<xmax[0]-xmin[0]<<"   0.0  "<<xmax[1]-xmin[1]<<"    0.0 "<<xmax[2]-xmin[2]<<endl            
		 <<dx_[2]<<" "<<"  " <<dx_[1]<<" "<<"  " <<dx_[0]<<endl                    
		 <<X0_[2]<<" "<<"  " <<X0_[1]<<" "<<"  " <<X0_[0]<<endl                    
		 <<"\n\nNote: image is rotated"<<endl  
		 <<"\n\nComments:"<<endl  
		 <<" first 9 entries above are:"<<endl  
		 <<"    Nx Ny Nz"<<endl  
		 <<"    dx dy dz"<<endl  
		 <<"    Xo Yo Zo"<<endl  
		 <<" Nx, Ny and Nz  count for the number of columns, rows and layers respectively"<<endl         
		 <<" Optional keywords (move above Comments to activate):"<<endl  
		 <<"    crop	    0  300   0  300   0  300 "<<endl  
		 <<"    pore 	    0 0 "<<endl  
		 <<"    resample	1"<<endl         
		 <<endl;
		outputHeaderFile.close();
	}
	else
	{
		voxelField::write(outputName);
		writeHeader(outputName+"_header");
	}

}


void voxelMesh::median(short nNeist) 
{
	const voxelMesh voxls(*this);
    long long nChanges = 0;
    for ( unsigned int k=1; k<voxls.size()-1 ; k++ )
		for ( unsigned int j=1; j<voxls[k].size()-1 ; j++ )
			for ( unsigned int i=1; i<voxls[k][j].size()-1 ; i++ )
			{
			   unsigned char pID = voxls[k][j][i];

				 short nSames(0);
				 map<unsigned char,short> neis;///.  ID-counter

				  unsigned char 
				 neiPID = voxls[k][j][i-1];
				 if (neiPID != pID  ) 	 ++(neis.insert(pair<unsigned char,short>(neiPID,0)).first->second); else ++nSames;
				 neiPID = voxls[k][j][i+1];
				 if (neiPID != pID  ) 	 ++(neis.insert(pair<unsigned char,short>(neiPID,0)).first->second); else ++nSames;
				 neiPID = voxls[k][j-1][i];
				 if (neiPID != pID  ) 	 ++(neis.insert(pair<unsigned char,short>(neiPID,0)).first->second); else ++nSames;
				 neiPID = voxls[k][j+1][i];
				 if (neiPID != pID  ) 	 ++(neis.insert(pair<unsigned char,short>(neiPID,0)).first->second); else ++nSames;
				 neiPID = voxls[k-1][j][i];
				 if (neiPID != pID  ) 	 ++(neis.insert(pair<unsigned char,short>(neiPID,0)).first->second); else ++nSames;
				 neiPID = voxls[k+1][j][i];
				 if (neiPID != pID  ) 	 ++(neis.insert(pair<unsigned char,short>(neiPID,0)).first->second); else ++nSames;

				 if(nSames<3)
				 {
					 pair<unsigned char,short> neiBest(pID,0);
					 map<unsigned char,short>::iterator  neitr=neis.begin();
					 for( ;  neitr!=neis.end(); ++neitr) if (neitr->second>neiBest.second) neiBest=*neitr;
					 //~ pair<unsigned char,short>* neitr=&neiBest;
					 //~ map<unsigned char,short>::iterator neitr = max_element(neis.begin(), neis.end(), neis.value_comp());
					 if ( (neiBest.second>=nSames))
					 {
						 if ( neiBest.second>nNeist)
						 {
							++nChanges; 
							(*this)[k][j][i] = neiBest.first;
							//~ cout<<"  "<<neis.size()<<"    "<<int(pID)<<" "<<nSames<<"    :   "<<int(neiBest.first)<<" "<<neiBest.second<<">="<<nNeist<<"    "<<endl;
						 }
						 else if ( neiBest.second==nNeist && pID!=1 )  //
						 {
							++nChanges; 
							(*this)[k][j][i] = neiBest.first;
							//~ cout<<"  "<<neis.size()<<"    "<<int(pID)<<" "<<nSames<<"    :   "<<int(neiBest.first)<<" "<<neiBest.second<<">="<<nNeist<<"    "<<endl;
						 }
					 }
				 }
			  }

	( cout<<"  nMedian"<<nNeist<<": "<< left<<nChanges<<"  \n").flush();

}



//////////////////////////////////////////////////////////////////////
/*
template<typename T>
void replaceValue(voxelImageT<T>& vImage, T v, T newv)
{
	 (cout<<int(v)<<"->"<<int(newv)<<"    ").flush();
    forAllkji(vImage)  if (vImage[k][j][i]==v)  vImage[k][j][i]=newv;
}
void voxelMesh::ParticleLabel()
{
	long long nGanglia=0;
	//unsigned int nLabels = 0;
	 voxelMesh voxls(*this);
	//const voxelMesh nLabels(*this) == (unsigned char)0;
	unsigned int n[3]={0,0,0};
    getSize(n[0],n[1],n[2]);
    cout<<"calculating ganglia size (voxels):"<<endl;
    long long nChanges = 1;
    while(nChanges)
    { nChanges = 0;
		
		for ( unsigned int k=0; k<this->size() ; k++ )
			for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
				for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
					
					if (voxls[k][j][i]==(unsigned char)2 )
					{
						++nGanglia;
					 
						unsigned char minv = voxls[k][j][i];
						minv=min(minv,voxls[k][j][i-1]);
						minv=min(minv,voxls[k][j][i+1]);
						minv=min(minv,voxls[k][j-1][i]);
						minv=min(minv,voxls[k][j+1][i]);
						minv=min(minv,voxls[k-1][j][i]);
						minv=min(minv,voxls[k+1][j][i]);
						if(minv<voxls[k][j][i]) {voxls[k][j][i]=minv; ++nChanges;}

					}
                
	 cout<<": "<<nChanges<<" ============  "; cout.flush();

	}
	
	///write the file
	ofstream out;            // output
	out.open("Labeled.raw",ios::binary);
	assert(out);
	for ( unsigned int k=0; k<this->size() ; k++ )
		for ( unsigned int j=0; j<(*this)[k].size() ; j++ )
			for ( unsigned int i=0; i<(*this)[k][j].size() ; i++ )
				
				out << voxls[k][j][i];' ';
				
	out.close();
}
*/
