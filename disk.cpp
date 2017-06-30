#define PI 3.14159265358979323846
#define DISK 1
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex>
#include<iostream>
#include<cstdlib>
#include<time.h>
using namespace std;  



double Fract(double n)					//Fractorial, invalid for n>30, for it leaks out.
{
	double s=1;
	if((double)(n)<0.01&&(double)(n)>-0.01)
		return 1.;
	if(n<0.0)
		return 0.;
	for(double i=n;i>0;i=i-1.0)
	{
		s=s*i;
	}
	
	return s;
}
double C(double up, double low)			//Number of combinations
{
	if(low<0)
		return 0;
	if(up<low)
		return 0;
	if(up==low&&up==0)
		return 1;
	double multi=1;
	for(double i=up;i>up-low;i--)
	{
		multi=multi*i;
	}
	multi=multi/Fract(low);
	return multi;
	
}


complex<double> LY(int n, complex<double> *z, int p)	//Laughlin wave-function, with particle number n and their coordinates z, CF capturing 2p vortixes.
{
	int i,j;
	complex<double> multi;
	complex<double> sum;
	
	multi=polar(1.,0.);
	sum=polar(0.,0.);
	for(i=0;i<n;i++)
	{
		for(j=i+1;j<n;j++)
		{
			multi=multi*pow((z[i]-z[j]),double(2*p+1));
		}
	}
	for(i=0;i<n;i++)
	{
		sum=sum+z[i]*conj(z[i]);
	}
	return multi*exp(-0.25*sum)/double(DISK);
}
double EllipticE(double x)	//Elliptic function, for the use of the evaluation for some wave functions.
{
	return PI/2.-PI/8.*x-3.*PI/128.*x*x-5.*PI/512.*pow(x,3)-175.*PI/32768.*pow(x,4)-441.*PI/131072.*pow(x,5)-4851*PI/2097152.*pow(x,6)-14157.*PI/8388608.*pow(x,7)-2760615.*PI/2147483648.*pow(x,8)-8690825.*PI/8589934592*pow(x,9)-112285459.*PI/137438953472.*pow(x,10);

}
/***********************************/
/**Three potentials: e-e, e-b, b-b**/
/***********************************/
double Vee(int n, complex<double> *z)
{
	double V=0;
	for(int i=0;i<n;i++)
		for(int j=i+1;j<n;j++)
		{
			V=V+1./abs(z[i]-z[j]);
		}
	return V;
}
double Vbb(int n, double niu)
{
	return n*8./3./PI*sqrt(niu*n/2);
}
double Vbe(int n, double niu, double RN, complex<double> *z)
{
	double sum=0.0;
	for(int i=0;i<n;i++)
	{
		sum=sum+EllipticE((abs(z[i])/RN)*(abs(z[i])/RN));
	}
	return -sqrt(2.*niu*n)*sum*2/PI;
}
/*****************************************************************************/
/* For general disk functions.
/*****************************************************************************/

complex<double> Jas(int n_p, complex<double> *z)	//Jastrow factor without oder 2p, also written as (phi_1)^2p (LLL function). n_p=particle number
{
	complex<double> product;
	product=polar(1.0,0.0);
	for(int i=0;i<n_p;i++)
	{
		for(int j=i+1;j<n_p;j++)
		{
			product=product*(z[i]-z[j]);
		}
	}
	return product;
}
double JasRatio(int n_p, complex<double> *z, complex<double> zn, int nw)
{
	double sum1=0.,sum2=0.;
	for(int j=0;j<n_p;j++)
	{
		if(j!=nw)
		{
			sum1=sum1+log(norm(z[j]-zn))-log(norm(z[j]-z[nw]));
		}
	}
	return exp(sum1);//*exp(-0.5*(norm(zn)-norm(z[nw])));
}

double Nm(int n, int m)	//Normalization factor
{
	return pow(-1,n)*sqrt(Fract(n)/(2.*PI*pow(2,m)*Fract(n+m)));
}
complex<double> ynm(int n, int m, complex<double> z)	//Disk single particle wave-function without Gaussian Factor
{
	complex<double> sum;
	complex<double> subsm;
	sum=polar(0.,0.);
	
	for(int k=0;k<=n;k++)
	{
		subsm=polar(0.,0.);
		for(int l=0;l<=k;l++)
		{
			subsm=subsm+Fract(k)/Fract(l)*C(k+m,k-l)*pow(z,m+l);
		}
		sum=sum+subsm*pow(-1,k)*C(n+m,n-k)/Fract(k);
	}
	return sum*Nm(n,m);
}
complex<double> CFY(int n, int n_p, complex<double> *z, int p)	//CFY, filling factor=niu=n/(2pn+1)
{
	int i=0,j;
	double RN;
	double gaus=0.;
	double niu=(double)n/(2.*(double)p*(double)n+1.);
	
	RN=sqrt(2.*n_p/niu);		//magnetic length=l not l*, so we use niu not n.
	//cout<<(int)sqrt(RN*RN/2.)<<endl;
	
	complex<double> product;
	complex<double> sum;
	
	
	sum=polar(0.,0.);
	for(j=0;j<n_p;j++)
	{
		i=j;
		product=polar(1.,0.);
		for(int nl=0;nl<n;nl++)
		{
			for(int m=-nl;m<(int)(RN*RN/2.)&&i<n_p;m++)
			{
				product=product*ynm(nl,m,z[i%n_p])/RN;
				product=product*pow(-1,i);
				i=i+1;
			}
		}
		sum=sum+product;
	}
	for(int k=0;k<n_p;k++)
	{
		gaus=gaus+norm(z[k]);
	}
	return sum;//*exp(-0.25*gaus);//*Jas(n_p, p, z);
}

double MonCa(int max, int n, int ln, int p/*, double R*/)		//max means the max steps of computation, n=particle number, R=radis of the plate, CF capturing 2p vortixes.
{
	int i;
	
	double RN;
	double niu;
	niu=1./(2*p+1);
	RN=sqrt(2*n/niu);
	
	complex<double> z[n];
	complex<double> r[n];
	double testp;
	double sum=0;
	for(i=0;i<n;i++)
	{
		z[i]=polar((double)(double(rand())/double(RAND_MAX))*RN,(double)(double(rand())/double(RAND_MAX))*2.*PI);
		//cout<<z[i]<<endl;
	}
	srand(time(NULL));
	for(i=0;i<max;i=i+1)
	{
		for(int j=0;j<n;j++)
		{
			r[j]=z[j];
			
			
		}
		r[i%n]=r[i%n]+polar((double(rand())/double(RAND_MAX))*0.2*RN,(double(rand())/double(RAND_MAX))*2.*PI);	//Slow here, change is small.
		//cout<<r[i%n]<<"         "<<z[i%n]<<endl;
		
		//testp=norm(LY(n,r,p))/norm(LY(n,z,p));
		testp=norm(CFY(ln, n, r, p))/norm(CFY(ln, n, z, p))*exp(-0.5*(norm(r[i%n])-norm(z[i%n])))*pow(JasRatio(n,z,r[i%n],i%n),2*p);
		//cout<<"NORM="<<norm(CFY(ln, n, r, p))<<"  "<<norm(CFY(ln, n, z, p))<<endl;
		if(testp>(double(rand())/double(RAND_MAX)))
		{				
			//cout<<"WES"<<endl;
			for(int m=0;m<n;m++)
			{
				z[m]=r[m];		//Metropolis acceptance.
			}
		}
		sum=sum+Vee(n,z)+Vbe(n,niu,RN,z)+Vbb(n,niu);	//cout<<"Number="<<i<<"  "<<Coulomb(n,z)<<endl;
		
	}
	cout<<"Filling Factor is :"<<ln/(2*p*ln+1.)<<endl;
	return sum/double(max)/n;
	
}

int main()
{
	int n;
	complex<double> c1;
	c1=polar(1.,PI/3);
	cout<<pow(c1,3)<<endl<<RAND_MAX<<endl;	
	for(int i=0;i<10;i++)
	{
		//cout<<rand()%100<<" "<<rand()%100<<endl;
	}
	
	for(;;)
	{
		cin>>n;
		cout<<MonCa(500000,n , 1 ,2)<<endl;
	}
	//cout<<(MonCa(1000,8,2,5)-MonCa(1000000,8,0,5))/10.<<endl;
	//cout<<MonCa(10000,10,0,10)<<endl;
}
