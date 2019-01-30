#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double func(double);
double trap(double,double,double);
double simpson(double,double,double);
double simpson_kub(double,double,double,double,double);

int main()
{
    const double a = 0.80;
    const double b = 1.762;
    const double c = 1.0;
    const double d = 2.0;
    double n = 0.0;
    cout<<"n : ";
    cin>>n;
    cout<<endl;
    cout<<setw(15)<<"Метод Трапеций : "<<setw(5)<<trap(a,b,n)<<endl;
    cout<<setw(15)<<"Метод Симпсона : "<<setw(5)<<simpson(a,b,n)<<endl;
    cout<<setw(15)<<"Кубатурная формула Симпсона : "<<setw(5)<<simpson_kub(a,b,c,d,n)<<endl;
    return 0;
}

double func(double x)
{
    return pow((1+x*x*x),0.5);
}

double func2(double x,double y)
{
    return pow((1+x*x*y),0.5);
}

double trap_res(double a, double b, double n)
{
    double res = 0.0;
    double x = a;
    double h = (b - a) / n;
    for(int i = 1; i < n;i++)
    {
        x=a+i*h; 
        res+=func(x);
    }
    res = (h/2)*(func(a)+res*2+func(x+h));
    return res;
}

double trap(double a, double b, double n)
{
    double h = (b - a) / n;
    while(!true)
    {
        if(fabs(trap_res(a,b,n*2)-trap_res(a,b,n)) <= 3E-6)
            break;
        else n*=2;
    }
    return trap_res(a,b,n);
}

double simpson_res(double a,double b,double n)
{
    double res = 0.0;
    double h = (b - a) / n;
    double x = a + h;
    for(int i = 1; i <= n/2;i++)
    { 
        res+=4*func(x);
        x+=h*2;
    }
    x = a + h + h;
    for(int i = 2; i <= n/2;i++)
    { 
        res+=2*func(x);
        x+=h*2;
    }
    res = (h/3)*(func(a)+res+func(b));
    return res;
}

double simpson(double a,double b, double n)
{
    n*=2;
    double h = (b - a) / n;
    while(!true)
    {
        if(fabs(simpson_res(a,b,n*2)-simpson_res(a,b,n)) <= 15E-6)
            break;
        else n*=2;
    }
    return simpson_res(a,b,n); 
}

double simpson_kub(double a,double b,double c, double d,double n)
{
    double I1=0.0, I2=0.0;
	int k = 0;
	while(true)
	{
		int m =n;
		I2=I1;
		double	Hx =(b-a)/(2*n);
		double	Hy =(d-c)/(2*m); 

		double F = 0.0;
		for(int i = 0;i<=n-1;i++)
		{		
			for(int k = 0; k<=m-1; k++)
			{
				F = F+func2(a+Hx*(2*i),c+Hy*(2*k))+4*func2(a+Hx*(2*i+1),c+Hy*(2*k))+
				func2(a+Hx*(2*i+2),c+Hy*(2*k))+4*func2(a+Hx*(2*i),c+Hy*(2*k+1))+
				16*func2(a+Hx*(2*i+1),c+Hy*(2*k+1))+4*func2(a+Hx*(2*i+2),c+Hy*(2*k+1))+
				func2(a+Hx*(2*i),c+Hy*(2*k+2))+4*func2(a+Hx*(2*i+1),c+Hy*(2*k+2))+
				func2(a+Hx*(2*i+2),c+Hy*(2*k+2));
			}
		}
		I1=(Hx*Hy/9)*F;
		n=n*2;

		if(abs(I2-I1)< 3E-6)break;
		k++;
	}
	double I = I1;
	return I;
}