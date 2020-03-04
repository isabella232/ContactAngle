

#include <sys/stat.h>
#include <math.h>
    #include <fstream>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <iostream>
using namespace std; 



int n[3];

class shape  
{	
public:
	virtual unsigned char value(double , double , double )=0;
};


class sphere : public shape
{
	double a, b, c,r;
public:
	sphere(double aa, double bb, double cc, double rr)
	{
		a=aa;		b=bb;		c=cc;;		r=rr;
		
		cout <<"sphere: a="<<a<<"   b="<<b<<"    c="<<c<<"    r="<<r<<endl;
	};
		
	unsigned char value(double xx, double yy, double zz)
	{
		if  ( (xx-a)*(xx-a) + (yy-b)*(yy-b) + (zz-c)*(zz-c)<=r*r
			)
			return 2; //Oil
		else 
			return 0; //Brine
	}
};




int main(int argc, char** argv)
{
  try
    {	
	char* fname = argv[1];  // file name
	ifstream in;            // file pointer
	cout<<"Openning data file: "<<fname<<endl;
	in.open(fname);
	assert(in);

	char tmpc;
	for ( unsigned int i=0; i<8;i++)   in>>tmpc, cout<<" "<<tmpc;//, out<<tmpc<<' '; out<<endl; //ignore the first 8 characters (ascii 3uc)

	in>>n[0]>>n[1]>>n[2];   //out<<n[0]<<' '<<n[1]<<' '<<n[2]<<endl;   	                 // number of variables (dimension of
	cout<<", "<<n[0]<<", "<<n[1]<<", "<<n[2]<<endl;  

	double cccccc;
	double xmin[]={0, 0, 0};
	double xmax[]={n[0]+0.0, n[1]+0.0, n[2]+0.0};  
	double dx[]={1.0, 1.0, 1.0};
	
	in>>    cccccc>>cccccc>>cccccc>>cccccc>>cccccc>>cccccc ;


shape * voxelizedshape[1];
for (int i=0;i<1;i++)
{
	in>>tmpc;
	if (tmpc=='s')
	{
		double rr;
		in>>rr;  
		voxelizedshape[i]=new sphere(n[0]/2,n[1]/2,n[2]/2,rr);
		cout <<"\nsphere: "<<"  r="<<rr<<"\n"<<endl;
	}	
}
	in.close();

///write the file
	ofstream out;            // output
	out.open("output.raw",ios::binary);
	assert(out);
	for (register double x2=xmin[2]+dx[2]/2; x2< xmax[2];x2+=dx[2])
	{
		for (register double x1=xmin[1]+dx[1]/2;x1< xmax[1];x1+=dx[1])
		{
			for (register double x0=xmin[0]+dx[0]/2;x0< xmax[0];x0+=dx[0])
			{
				out<<//voxelizedshape[1]->value(x0,x1,x2)*
					 voxelizedshape[0]->value(x0,x1,x2);//  <<' ';
			}
		}
	}	
	out.close();






	
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


