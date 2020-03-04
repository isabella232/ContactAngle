

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

class paralellPlate : public shape
{
	double a1,b1,a2,b2,h,c1,c2;
public:
	paralellPlate(double a, double b, double c)
	{
		a1=a;
		b1=0;
		c1=c;
		a2=a;
		b2=b;
		c2=c;
			
		cout <<"paralellPlate: a1="<<a1<<"   a2="<<a2<<"    b1="<<b1<<"    b2="<<b2<<endl;
	};
		
	unsigned char value(double x, double y, double z)
	{
		if ( y > a1*x+c1*z+b1 && y < a2*x+c2*z+b2 )
		{
			return 0; //pore
			//~ cout<<' '<<0;
		}
		else 
		{
			return 1; //rock
			//~ cout<<' '<<1;		
		}
	}
};


	double mag(double a, double b, double c)
	{
		return sqrt(a*a+b*b+c*c);
	}
		
	double magCrossProduct(double a1, double a2, double a3, double b1, double b2, double b3)
	{
		double aa=a2*b3 - a3*b2;
		double bb=a3*b1 - a1*b3;
		double cc=a1*b2 - a2*b1;
		return mag(aa,bb,cc);
	}

	
	
class voxelCylinder : public shape
{
	double a1,b1,a2,b2,h,c1,c2,d1,r;
	double x1,y1,z1,x2,y2,z2,magx1_x2;
public:
	voxelCylinder(double a, double b, double c)
	{//http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		//y=a1*x+b/2, z=c1*x : two lines describing centre of cylinder
		a1=a;		b1=b;		c1=c;		a2=a;		b2=b;		c2=c;
		//~ r=sqrt(pow(b/2,2)/(1+pow(a,2)+pow(c,2)));
		 r=b/2;
		x1=0; y1=2*r; z1=r; 
		x2=x1+1; y2=y1+a1; z2=z1+c1; 
		magx1_x2=mag(x2-x1,y2-y1,z2-z1);
		 		//~ d1=1;//sqrt(1+pow(a,2)+pow(c,2));

		cout <<"voxelCylinder: a1="<<a1<<"    b1="<<b1<<"   c1="<<c1<<"    r="<<r<<endl;
	};

	
	unsigned char value(double x, double y, double z)
	{
		//~ cout<<"y "<<y <<"x"<<
		if ( magCrossProduct(x-x1,y-y1,z-z1,x-x2,y-y2,z-z2)/magx1_x2 <= r )
		{
			 //~ cout<<"x:  " <<x<<' '<<y<<' '<<z<<endl;
			return 0; //pore
		}
		else 
		{
						 //~ cout<<"x:  " <<x<<' '<<y<<' '<<z<<endl;
			return 1; //rock

			//~ cout<<' '<<1;		
		}
	}
};

class plate : public shape
{
	double a1,b1,a2,b2,h,c1,c2,xo_,zo_;
public:
	plate(double a, double b, double c, double xo, double zo)
	{
		a1=a;		b1=0;		c1=c;
		a2=a;		b2=b;		c2=c;
		xo_=xo;		zo_=zo;
		cout <<"plate: a1="<<a1<<"   a2="<<a2<<"    r="<<b1<<"    b2="<<b2<< " x0="<<xo<<"   zo="<<zo<<endl;
	};
		
	unsigned char value(double x, double y, double z)
	{
		if  ( y > a1*x+c1*z+b1 //&& y < a2*x+c2*z+b2 
			  && mag(x - xo_ ,   y - (a1*xo_+c1*zo_+b1)  ,   z - zo_) <= b2 
			)
			return 0; //pore
		else 
			return 1; //rock
	}
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
	
	  //~ xmin[3],xmax[3],dx[3];

	double cccccc;
	double xmin[]={0, 0, 0};
	double xmax[]={n[0]+0.0, n[1]+0.0, n[2]+0.0};  double dx[]={1.0, 1.0, 1.0};
	
	in>>    cccccc>>cccccc>>cccccc>>cccccc>>cccccc>>cccccc ;
	//~ out<<   xmin[0]<<' '<<xmax[0]<<' '<<xmin[1]<<' '<<xmax[1]<<' '<<xmin[2]<<' '<<xmax[2]<<endl ;
	//~ cout<<  xmin[0]<<", "<<xmax[0]<<", "<<xmin[1]<<", "<<xmax[1]<<", "<<xmin[2]<<", "<<xmax[2] <<endl;;
//initialize
	//~ dx[0]=(xmax[0]-xmin[0])/n[0];dx[1]=(xmax[1]-xmin[1])/n[1];dx[2]=(xmax[2]-xmin[2])/n[2];
	//~ cout<<  "dx"<<": "<<dx[0]<<", "<<dx[1]<<", "<<dx[2]<<endl;



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
	else if (tmpc=='p')
	{
		double a,b,c;
		in>>a >> b>>c;  
		voxelizedshape[i]=new plate(a,b,c,(xmax[0]-xmin[0])/2,(xmax[2]-xmin[2])/2);
		cout <<"\nplate: slope="<<a<<"  r="<<b<<"  rr="<<b/sqrt(1+a*a)<<"\n"<<endl;
	}	
	else if (tmpc=='f')
	{
		double a,b,c;
		in>>a >> b>>c;  
		voxelizedshape[i]=new paralellPlate(a,b,c);
		cout <<"\nslope="<<a<<"  seperation="<<b<<"  h="<<b/sqrt(1+a*a)<<"\n"<<endl;
	}
	else if (tmpc=='c')
	{
		double a,b,c;
		in>>a >> b>>c;  
		voxelizedshape[i]=new voxelCylinder(a,b,c);
		cout <<"\nslope="<<a<<"  seperation="<<b<<"  rr="<<b/sqrt(1+a*a)<<"\n"<<endl;
	}

	else
	{
		cout <<"unsupported shape type: "<<tmpc<<"in file "<<fname<<", only 'p' supported, standing for paralell plate"<<endl;
	}
}
	in.close();
	//~ paralellPlate   pp(a,b,c);
//write the file
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
			//~ out<<endl;
		}
		//~ out<<endl;
	}	
	//~ out<<endl;
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


