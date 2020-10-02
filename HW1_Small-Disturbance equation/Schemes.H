#ifndef SCHEMES_H
#define SCHEMES_H
#include <iostream>
#include <vector>
#include "Matrix.H"
#include <cmath>
#include <map>
#include "Solve.H"
#include "Field.H"


void pointJacobi(Field& phi)
{
    for(int i=1;i<phi.mesh.ny-1;i++)
    {
        for(int j=1;j<phi.mesh.nx-1;j++)
        {
            {   
                double dxp = phi.mesh.x[j+1]-phi.mesh.x[j];
                double dxm = phi.mesh.x[j]-phi.mesh.x[j-1];
                double dyp = phi.mesh.y[i+1]-phi.mesh.y[i];
                double dym = phi.mesh.y[i]-phi.mesh.y[i-1];
                double dxc = phi.mesh.x[j+1]-phi.mesh.x[j-1];
                double dyc = phi.mesh.y[i+1]-phi.mesh.y[i-1];

                double aii= 0.5 /(A/dxp/dxm+1/dyp/dym);
                phi.value1[i][j] = aii * (2.0*A/dxp/dxc * phi.value[i][j+1] + 2.0*A/dxm/dxc * phi.value[i][j-1] + 2.0/dyp/dyc * phi.value[i+1][j]+ 2.0/dym/dyc * phi.value[i-1][j]);
            }
        } 
    }
}

void pointGaussSeidel(Field& phi)
{
    for(int i=1;i<phi.mesh.ny-1;i++)
    {
        for(int j=1;j<phi.mesh.nx-1;j++)
        {
            {   
                double dxp = phi.mesh.x[j+1]-phi.mesh.x[j];
                double dxm = phi.mesh.x[j]-phi.mesh.x[j-1];
                double dyp = phi.mesh.y[i+1]-phi.mesh.y[i];
                double dym = phi.mesh.y[i]-phi.mesh.y[i-1];
                double dxc = phi.mesh.x[j+1]-phi.mesh.x[j-1];
                double dyc = phi.mesh.y[i+1]-phi.mesh.y[i-1];

                double aii= 0.5 /(A/dxp/dxm+1/dyp/dym);
                phi.value1[i][j] = aii * (2.0*A/dxp/dxc * phi.value[i][j+1] + 2.0*A/dxm/dxc * phi.value1[i][j-1] + 2.0/dyp/dyc * phi.value[i+1][j]+ 2.0/dym/dyc * phi.value1[i-1][j]);                                                                                    
            }
        }
    }
}

void lineJacobi(Field& phi)
{
    int n = phi.mesh.ny;
    std::vector<double> a(n);
    std::vector<double> soln(n) ,b(n-1),c(n-1),d(n);
    int nstart =0;
    int nend = n-1;
    for(int j=1;j<phi.mesh.nx-1;j++)
    {
    //  a.clear
        nstart = 0;
        nend = n-1;
        for(int i=0;i<phi.mesh.ny;i++)
        {   
            
            if(phi.mesh.onBoundary[i][j] == 0)
            {   
                double dxp = phi.mesh.x[j+1]-phi.mesh.x[j];
                double dxm = phi.mesh.x[j]-phi.mesh.x[j-1];
                double dyp = phi.mesh.y[i+1]-phi.mesh.y[i];
                double dym = phi.mesh.y[i]-phi.mesh.y[i-1];
                double dxc = phi.mesh.x[j+1]-phi.mesh.x[j-1];
                double dyc = phi.mesh.y[i+1]-phi.mesh.y[i-1];
                if(i==n-2)
                {
                    d[i]=( 2.0*A/dxp/dxc * phi.value[i][j+1] + 2.0*A/dxm/dxc * phi.value[i][j-1]);
                    a[i] = (2.0/dyp/dym + 2.0 * A/dxp/dxm);
                    d[i]+=(2.0/dyp/dyc)*phi.value[i+1][j];
                    c[i-1]=(-2.0/dym/dyc);
                }
                else
                {
                    d[i]=( 2.0*A/dxp/dxc * phi.value[i][j+1] + 2.0*A/dxm/dxc * phi.value[i][j-1]);
                    a[i] = (2.0/dyp/dym + 2.0 * A/dxp/dxm);
                    b[i]=(-2.0/dyp/dyc);
                    c[i-1]=(-2.0/dym/dyc);
                }
                
                
            }
            else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 1)
            {
             //   a[i] = 1;
                soln[i] = phi.mesh.bcValueMap[i][j];
                if(i==0)
                {
                    nstart++;
                }
                else if(i==n-1)
                {
                    nend--;
                }
            }
            else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 2)  
            {
                a[i] = 1;
                b[i] = -1;
                d[i] = phi.mesh.bcValueMap[i][j];
            }
            else
            {
                std::cout<<"wrong BC type"<<std::endl;
                exit;
            }
        }
        nchase(a,b,c,d,soln,nstart,nend);
        for(int k=0;k<phi.mesh.nx;k++)
        {
            phi.value1[k][j]=soln[k];
        }
    }
}

void gaussSeidelLineRelaxation(Field& phi)
{
    int n = phi.mesh.ny;
    std::vector<double> a(n);
    std::vector<double> soln(n) ,b(n-1),c(n-1),d(n);
    int nstart =0;
    int nend = n-1;
    for(int j=1;j<phi.mesh.nx-1;j++)
    {
    //  a.clear
        nstart = 0;
        nend = n-1;
        for(int i=0;i<phi.mesh.ny;i++)
        {   
            
            if(phi.mesh.onBoundary[i][j] == 0)
            {   
                double dxp = phi.mesh.x[j+1]-phi.mesh.x[j];
                double dxm = phi.mesh.x[j]-phi.mesh.x[j-1];
                double dyp = phi.mesh.y[i+1]-phi.mesh.y[i];
                double dym = phi.mesh.y[i]-phi.mesh.y[i-1];
                double dxc = phi.mesh.x[j+1]-phi.mesh.x[j-1];
                double dyc = phi.mesh.y[i+1]-phi.mesh.y[i-1];
                if(i==n-2)
                {
                    d[i]=( 2.0*A/dxp/dxc * phi.value[i][j+1] + 2.0*A/dxm/dxc * phi.value1[i][j-1]);
                    a[i] = (2.0/dyp/dym + 2.0 * A/dxp/dxm);
                    d[i]+=(2.0/dyp/dyc)*phi.value[i+1][j];
                    c[i-1]=(-2.0/dym/dyc);
                }
                else
                {
                    d[i]=( 2.0*A/dxp/dxc * phi.value[i][j+1] + 2.0*A/dxm/dxc * phi.value1[i][j-1]);
                    a[i] = (2.0/dyp/dym + 2.0 * A/dxp/dxm);
                    b[i]=(-2.0/dyp/dyc);
                    c[i-1]=(-2.0/dym/dyc);
                }
                
                
            }
            else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 1)
            {
             //   a[i] = 1;
                soln[i] = phi.mesh.bcValueMap[i][j];
                if(i==0)
                {
                    nstart++;
                }
                else if(i==n-1)
                {
                    nend--;
                }
            }
            else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 2)  
            {
                a[i] = 1;
                b[i] = -1;
                d[i] = phi.mesh.bcValueMap[i][j];
            }
            else
            {
                std::cout<<"wrong BC type"<<std::endl;
                exit;
            }
        }
        nchase(a,b,c,d,soln,nstart,nend);
        for(int k=0;k<phi.mesh.nx;k++)
        {
            phi.value1[k][j]=soln[k];
        }
            
    }
}

void ADITwoStep(Field& phi)
{
    gaussSeidelLineRelaxation(phi);
    phi.value.swap(phi.value1);
    int n = phi.mesh.nx;
    std::vector<double> a(n);
    std::vector<double> soln(n) ,b(n-1),c(n-1),d(n);
    int nstart =0;
    int nend = n-1;
    for(int i=1;i<phi.mesh.nx-1;i++)
    {
    //  a.clear
        nstart = 0;
        nend = n-1;
        for(int j=0;j<phi.mesh.ny;j++)
        {   
            
            if(phi.mesh.onBoundary[i][j] == 0)
            {   
                double dxp = phi.mesh.x[j+1]-phi.mesh.x[j];
                double dxm = phi.mesh.x[j]-phi.mesh.x[j-1];
                double dyp = phi.mesh.y[i+1]-phi.mesh.y[i];
                double dym = phi.mesh.y[i]-phi.mesh.y[i-1];
                double dxc = phi.mesh.x[j+1]-phi.mesh.x[j-1];
                double dyc = phi.mesh.y[i+1]-phi.mesh.y[i-1];
                if(i==n-2)
                {
                    d[j]=( 2.0/dyp/dyc * phi.value[i+1][j] + 2.0/dym/dyc * phi.value[i-1][j]);
                    a[j] = (2.0/dyp/dym + 2.0 * A/dxp/dxm);
                    d[j]+=(2.0*A/dxp/dxc)*phi.value[i][j+1];
                    c[j-1]=(-2.0*A/dxm/dxc);
                }
                else if(i ==1)
                {
                    d[j] =( 2.0/dyp/dyc * phi.value[i+1][j] + 2.0/dym/dyc * phi.value[i-1][j]);
                    a[j] = (2.0/dyp/dym + 2.0 * A/dxp/dxm );
                    b[j] = (-2.0*A/dxp/dxc);
                  //  d[j]+=(2.0/dxm/dxc)*phi.value[i][j-1];
                    c[j-1]=(-2.0*A/dxm/dxc);
                }
                else
                {
                    d[j] = ( 2.0/dyp/dyc * phi.value[i+1][j] + 2.0/dym/dyc * phi.value[i-1][j]);
                    a[j] = (2.0/dyp/dym + 2.0 * A/dxp/dxm);
                    b[j] = (-2.0*A/dxp/dxc);
                    c[j-1]=(-2.0*A/dxm/dxc);
                }
                
                
            }
            else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 1)
            {
             //  a[i] = 1;
                soln[j] = phi.mesh.bcValueMap[i][j];
                if(j==0)
                {
                    nstart++;
                }
                else if(j==n-1)
                {
                    nend--;
                }
            }
            else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 2)  
            {
                a[j] = 1;
                b[j] = -1;
                d[j] = phi.mesh.bcValueMap[i][j];
            }
            else
            {
                std::cout<<"wrong BC type"<<std::endl;
                exit;
            }
        }
        nchase(a,b,c,d,soln,nstart,nend);
        for(int k=0;k<phi.mesh.nx;k++)
        {
            phi.value1[i][k]=soln[k];
        }
            
    }
}
#endif