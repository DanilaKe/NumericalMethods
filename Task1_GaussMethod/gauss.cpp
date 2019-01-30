//Метод Гаусса
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

void filla(int,float**,float*,ifstream&);
void gaussmethod(int,float**,float*);
void swap(float&,float&);
void printgauss(int,float**,float*);
void vectornev(int,float**,float*,float*);

//int main()  //if gauss 
int gauss()
{
    int n = 0;
    ifstream infile("matrixforgauss.txt");  //if gauss 
    //ifstream infile("jacobi.txt");
    infile>>n;
    infile.close();
    //cout<<"Порядок матрицы : "<<n<<endl; //if gauss 
    float **array,*vector;
    array=new float*[n];   //Матрица коэффициентов A
    for (int i = 0; i < n; i++)
        array[i]=new float[n];
    vector = new float[n];  //Вектор b
    ifstream in("matrixforgauss.txt");  //if gauss 
    //ifstream in("jacobi.txt");
    filla(n,array,vector,in);
    in.close();
    //printgauss(n,array,vector); //if gauss 
    gaussmethod(n,array,vector);
    delete [] vector;
    for(int i = 0; i < n; i++)
        delete [] array[i];
    delete [] array;
    return 0;
}

void printgauss(int n,float** array,float* vector) // Вывод на экран матрицу А и вектор b
{
    cout<<"Матрица A : "<<endl;
    for(int i = 0;i < n; i++)
    {
        for(int j = 0;j < n;j++)
            cout<<setw(4)<<array[i][j]<<"  ";
        cout<<endl;
    }
    cout<<"Вектор b : "<<endl;
    for(int i = 0;i < n; i++)
    {
        cout<<setw(4)<<vector[i]<<endl;
    }
}

void filla(int n,float **array,float *vector,ifstream &infile)  //заполнение массива и вектора из файла
{
    if(!infile) 
    {
        cout<<"Нету файла!"<<endl;
        return;
    }
    infile>>n;
    for(int i = 0;i < n; i++)
        for(int j = 0;j < n;j++)
            infile>>array[i][j];
    for(int i = 0;i < n; i++)
        infile>>vector[i];
}


void gaussmethod(int n,float** array,float* vector) 
{
    float max = 0.0,*result;
    result = new float[n];
    int max_index = 0;
    for(int k = 0;k < n; k++)
    {
        max_index = k;
        for(int i = k; i < n; i++)            
            if(abs(array[i][k]) > max)
            {
                max = abs(array[i][k]);
                max_index = i;
            }
        if(max==0)
        {
            cout<<"Матрица вырожденна!"<<endl;
            return;
        }
        if(k != max_index)
        {
            swap(vector[k],vector[max_index]);
            for(int i = k; i < n; i++)
                swap(array[max_index][i],array[k][i]);
        }
        if(abs(array[k][k])<=0.00001)
        {
            cout<<"Матрица вырожденна!"<<endl;
            return;
        }
        result[k] = 0.0;
        vector[k] /= array[k][k];
        for(int i = k+1; i < n; i++)
        {
            vector[i] -= array[i][k]*vector[k];
        }
        for(int j = n-1; j >= k; j--)
        {
            array[k][j]/=array[k][k];
            for(int i = 1+k; i < n;i++)
            {   
                array[i][j] -= array[k][j]*array[i][k];
            }
        }
    }
    result[n-1] = vector[n-1];
    for(int i = n-2; i >= 0; i--)
    {
        result[i] = vector[i];
        for (int j = i+1; j < n; j++) result[i]-=array[i][j]*result[j];
    }
    //cout<<"Ответ : "<<endl; //if gauss 
    ofstream out("gaussres.txt");
    for(int i = 0; i < n; i++)
    {
        //cout<<"x"<<i+1<<" = "<<result[i]<<endl;//if gauss 
        out<<result[i]<<endl;
    }
    out.close();
    /*ifstream infile("file1.txt");  //if gauss 
    filla(n,array,vector,infile);
    infile.close();
    vectornev(n,array,vector,result);*/
    delete [] result;
}

void vectornev(int n,float** array,float* vector,float* result)  //вычисленние вектора невязки
{
    float *vectornev = new float[n];
    cout<<"Вектор невязки :"<<endl;
    fstream out("gaussres.txt",ios::in | ios::out);
    for(int i = 0; i < n ; i++)
    {
        float a;
        out>>a;
    }
    out<<'\n';
    for(int i = 0; i < n; i++)
    {
        vectornev[i] = vector[i];
        for(int j = 0; j < n; j++)
        {
            vectornev[i] -= array[i][j]*result[j];
        }
        cout<<"("<<vectornev[i]<<')'<<'e'<<i+1<<endl;
        out<<vectornev[i]<<endl;
    }
    out.close();
    float a = 0.0;
    for(int i = 0; i < n; i++)
    {
        if(a<abs(vectornev[i])) a = vectornev[i];
    }
    cout<<"Норма вектора невязки : "<<a<<endl;
    delete [] vectornev;
}

void swap(float &a,float &b)
{
    float q = 0.0;
    q = a;
    a = b;
    b = q;
}