#include "gauss.cpp"

void fill_XY(double*,double*, int&);
void print_XY(double*,double*, int);
void make_SLAY(double*,double*,double*,double*,double**,int,int);
void fill_matrix_for_gauss(double**,double*,int,int);
void print_result(double*,double*,int,int);

int main()
{
    int m = 2;
    int n = 0;
    ifstream infile("table.txt");
    infile>>n;
    infile.close();
    double *x,*y,*powerx,**sumx,*praw;
    x = new double[n];
    y = new double[n];
    powerx = new double[m*2-2];
    praw = new double[m];
    sumx = new double*[m];
    for (int i = 0; i < m; i++)
        sumx[i]=new double[m];
    fill_XY(x,y,n);
    print_XY(x,y,n);
    make_SLAY(x,y,powerx,praw,sumx,n,m);
    delete [] powerx;
    delete [] praw;
    for(int i = 0; i < m; i++)
        delete [] sumx[i];
    delete [] sumx;
    gauss();
    print_result(x,y,n,m);
    delete [] x;
    delete [] y;
    return 0;
}

void fill_XY(double* x,double* y, int& n)
{
    ifstream infile("table.txt");
    infile>>n;
    for (int i = 0; i < n ; i++)
        infile>>x[i];
    for (int i = 0; i < n ; i++)
        infile>>y[i];
    infile.close();
    return;
}

void print_XY(double* x,double* y, int n)
{
    cout<<"n = "<<n<<'\n';
    for (int i = 0; i < n ; i++)
        cout<<setw(13)<<x[i];
    cout<<endl;
    for (int i = 0; i < n ; i++)
        cout<<setw(13)<<y[i];
    cout<<endl;
    return;
}

void make_SLAY(double* x, double* y, double* powerx, double* praw, double** sumx, int n, int m)
{
    double *x_res;
    x_res = new double[n];
    powerx[0] = 0.0;
    for(int i = 0; i < n; i++)
    {
        x_res[i] = x[i]; 
        powerx[0] += x_res[i]; 
    }
    for(int i = 1; i < m*2-2; i++)
    {
        powerx[i] = 0.0;
        for(int j = 0; j < n; j++)
        {
            x_res[j]*=x[j];
            powerx[i]+= x_res[j];
        }
    }
    delete [] x_res;
    for(int i = 0; i < m; i++)
        for(int j = 0; j < m; j++)
            sumx[i][j] = powerx[i+j-1];
    sumx[0][0] = n;
    for(int i = 0; i < m; i++)
    {
        praw[i] = 0.0;
        for(int j = 0; j < n; j++){
            //praw[i] += y[j]*pow(x[j],i);
            praw[i] += log10(y[j])*pow(x[j],i);
        }
    }
    fill_matrix_for_gauss(sumx,praw,n,m);
    return;
}

void fill_matrix_for_gauss(double** sumx, double* praw,int n, int m)
{
    ofstream file("matrixforgauss.txt");
    file<<m<<endl<<endl;
    for(int i = 0 ; i < m ; i++)
    {
        for(int j = 0; j < m; j++)
            file<<sumx[i][j]<<"  ";
        file<<endl;
    }
    file<<endl;
    for(int i = 0; i<m; i++)
        file<<praw[i]<<endl;
    return;
}

void print_result(double* x,double* y,int n,int m)
{
    ifstream file("gaussres.txt");
    double* a = new double[m];
    for(int i = 0 ; i < m ; i++)
        file>>a[i];
    file.close();
    a[0] = pow(10,a[0]);
    for(int i = 0; i < n; i++)
        cout<<setw(13)<<a[0]*pow(10.0,a[1]*x[i]);
    cout<<endl;
    cout<<"Ответ : y = "<<a[0]<<"*10^("<<a[1]<<"*x)";
    double d_res=0.0;
    for(int i = 0; i < n; i++)
    {
        double q = y[i];
        q-=a[0]*pow(10.0,a[1]*x[i]);
        d_res +=q*q;
    }
    cout<<endl;
    d_res *= 1.0/(n-m-1);
    cout<<"d = "<<sqrt(d_res)<<endl;
    return;
}