#include "gauss.cpp"

const float eps1 = 1E-9;
const float eps2 = 1E-9;
const float NIT = 15;

void jacobicalc(double**,double*);
void changedelta(double&,double&,double*);

double f1(double x1, double x2)
{
    double result = 0.0;
    //result = log(1+(x1+x2)/5)-sin(x2/3)-x1+1.1;
    result = x1 + 2*x2*x2-9;
    return result;
}

double f2(double x1, double x2)
{
    double result = 0.0;
    //result = cos(x1*x2/6)-x2+0.5;
    result = x1 + x2 - 3;
    return result;
}

int main()
{
    int n = 2,k = 0;
    double *x,*vectornev,**jacobi,d1,d2;
    x = new double[n];
    vectornev = new double[n];
    jacobi = new double*[n];
    for(int i = 0; i < n; i++)
    {
        jacobi[i] = new double[n];
        cout<<"x["<<i+1<<"] = ";
        cin>>x[i];
    }
    cout<<setw(4)<<'k'<<setw(20)<<"d1"<<setw(20)<<"d2"<<endl;
    do{
        vectornev[0] = f1(x[0],x[1]);
        vectornev[1] = f2(x[0],x[1]);
        /*for(int i = 0; i < n; i++)
            cout<<vectornev[i]<<endl;*/
        jacobicalc(jacobi,x);
        ofstream out("jacobi.txt");
        out<<n<<endl;
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
                out<<jacobi[i][j]<<"  ";
            out<<endl;
        }
        for(int i = 0; i < n; i++)
            out<<vectornev[i]<<endl;
        out.close();
        gauss();
        ifstream in("gaussres.txt");
        for(int i = 0; i < n; i++)
        {
            double a;
            in>>a;
            x[i]+=a;
            //cout<<x[i]<<endl;
        }
        in.close();
        changedelta(d1,d2,x);;
        cout<<setw(4)<< k+1 <<setw(20)<< d1 <<setw(20)<< d2 <<endl;
        k++;
        if(k>NIT) 
        {
            cout<<"IER = 2"<<endl;
            return -1;
        }
    }while(d1>eps1 || d2>eps2);
    cout<<"Ответ : "<<endl;
    for(int i = 0; i < n; i++)
        cout<<setw(4)<<'x'<<i<<" = "<<x[i]<<endl;
    return 0;
}

void jacobicalc(double** jacobi,double* x)
{
    jacobi[0][0] = (f1(x[0]-eps1,x[1])-f1(x[0],x[1]))/eps1;
    jacobi[0][1] = (f1(x[0],x[1]-eps1)-f1(x[0],x[1]))/eps1;
    jacobi[1][0] = (f2(x[0]-eps1,x[1])-f2(x[0],x[1]))/eps1;
    jacobi[1][1] = (f2(x[0],x[1]-eps1)-f2(x[0],x[1]))/eps1;
    return;
}

void changedelta(double& d1,double& d2,double* x)
{
    d1 = abs(f1(x[0],x[1]));
    if(d1<abs(f2(x[0],x[1]))) d1 = abs(f2(x[0],x[1]));
    ifstream in("gaussres.txt");
    double a[2],max;
    in>>a[0];
    in>>a[1];
    in.close();
    max = abs(a[0]);
    int k = 0;
    if(max<abs(a[1])) 
    {
        max = abs(a[1]);
        k = 1;
    }
    if(x[k]<1) d2 = max;
    else d2 = max/x[k];
    return;
}
