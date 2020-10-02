#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
const double M_inf = 0.735;
const double gamma = 1.4, rho_inf = 1, p_inf = 1/gamma/pow(M_inf,2); 
const double v_inf = 1;
const double a_inf = v_inf/M_inf;
#include "Field.H"
using namespace std;
//#include "Schemes.H"
/*int main(void)
{
    std::vector<double> a,b,c,d;
    for(int i=0;i<5;i++)
    {
        a.push_back(5);
    } 
    for(int i=0;i<4;i++)
    {
        b.push_back(6);
        c.push_back(1);
    }
    d.push_back(1);
    d.push_back(0);
    d.push_back(0);
    d.push_back(0);
    d.push_back(1);
    std::vector<double> x(a.size());
  //  chase(a,b,c,d,x);
    Jacobi(a,b,c,d,x);
    std::cout<<"x="<<std::endl;
    for(int i=0;i<x.size();i++)
    {
        std::cout<<x[i]<<"         ";
    }
    return 0;
}*/
double newtonsolve(double xmin,double xtotal,int n)
{
    double k = 2;
    double d,f,ff;
    for(;;)
    {
        f  = xtotal * (k - 1) - xmin * (pow(k,n) - 1);
        ff = xtotal - xmin*n*pow(k,n-1);
        d = f/ff;
        k  = k - f/ff;
        if(std::abs(f/ff)<1e-10)
        {
            return k;
        }
    }
}

int main()
{
    int Nx=51,Ny=51;
    double D = 50;
    int N1=15,N2=20;
    double c = 1, th = 0.06;
    double ymin = th/10;
    std::vector<double> xlist(Nx);
    std::vector<double> ylist(Ny);
    
        double k = newtonsolve(c/N2,D,N1);
        xlist[N1] = 0;
        for(int i=N1+1;i<Nx-N1;i++)
        {
            xlist[i] = xlist[i-1] + c/N2;
        }
        for(int i=Nx-N1;i<Nx;i++)
        {
            xlist[i] = xlist[Nx-1-N1] + D*(pow(k,i-Nx+1+N1) - 1)/(pow(k,N1)-1);
            xlist[Nx-1-i] = xlist[N1] - D*(pow(k,i-Nx+1+N1) - 1)/(pow(k,N1)-1);
        }
        ylist[0]=0;
        k = newtonsolve(ymin,D,Ny-1);
        for(int i=1;i<Ny;i++)
        {
            ylist[i] = ylist[0] + D * ( pow(k,i) -1 )/(pow(k,Ny-1) -1);
        }
    std::vector<std::vector<int> >bcTypemap;
    std::vector<std::vector<bool> >boundarymap;
    std::vector<std::vector<double> >bcvaluemap;

    for(int i=0;i<Ny;i++)
    {
        boundarymap.push_back(std::vector<bool>());
        bcTypemap.push_back(std::vector<int>());
        bcvaluemap.push_back(std::vector<double>());
        for(int j=0;j<Nx;j++)
        {
            bcTypemap[i].push_back(0);
            boundarymap[i].push_back(false);
            bcvaluemap[i].push_back(0);
        }
    }
    //initialize boundarypatches
    std::vector<boundaryPatch> bd;
    //top
    {
        std::vector<std::vector<int> > index,neighbour;
        std::vector<double> value;
        std::vector<double> dist;
        int type = 1;
        for(int j=0;j<Nx;j++)
        {
            index.push_back(std::vector<int>(2));
            neighbour.push_back(std::vector<int>(2));
            index[j][0] = Ny-1;
            index[j][1] = j;
            neighbour[j][0] = Ny-2;
            neighbour[j][1] = j;
            dist.push_back(std::sqrt(pow((xlist[index[j][1]]-xlist[neighbour[j][1]]),2)+pow((ylist[index[j][0]]-ylist[neighbour[j][0]]),2)));
            value.push_back(0);
            boundarymap[index[j][0]][index[j][1]] = true;
            bcTypemap[index[j][0]][index[j][1]] = 1;
            bcvaluemap[index[j][0]][index[j][1]] = value[j];
        }
        boundaryPatch bp(index,neighbour,value,dist,type); 
        bd.push_back(bp);
    }
    //bottom
    {
        std::vector<std::vector<int> > index,neighbour;
        std::vector<double> value;
        std::vector<double> dist;
        double a;
        int type = 2;
        for(int j=0;j<Nx;j++)
        { 
            index.push_back(std::vector<int>(2));
            neighbour.push_back(std::vector<int>(2));
            index[j][0] = 0;
            index[j][1] = j;
            neighbour[j][0] = 1; 
            neighbour[j][1] = j;
            
            dist.push_back(std::sqrt(pow((xlist[index[j][1]]-xlist[neighbour[j][1]]),2)+pow((ylist[index[j][0]]-ylist[neighbour[j][0]]),2)));
            if(j >= N1 && j<Nx-N1)
            {
                value.push_back(v_inf * (xlist[index[j][1]] - 0.5 * c)/std::sqrt( pow(25.09/6*c,2)-pow((xlist[index[j][1]]-0.5*c),2) )); 
            }
            else
            {
                value.push_back(0);
            }
            boundarymap[index[j][0]][index[j][1]] = true;
            bcTypemap[index[j][0]][index[j][1]] = 2;
            bcvaluemap[index[j][0]][index[j][1]] = value[j]*dist[j];
        }
        boundaryPatch bp(index,neighbour,value,dist,type);
        bd.push_back(bp);
    }
    //left
    {
        std::vector<std::vector<int> > index,neighbour;
        std::vector<double> value;
        std::vector<double> dist;
        int type = 1;
        for(int j=0;j<Ny;j++)
        {
            index.push_back(std::vector<int>(2));
            neighbour.push_back(std::vector<int>(2));
            index[j][0] = j;
            index[j][1] = 0;
            neighbour[j][0] = j;
            neighbour[j][1] = 1;
            dist.push_back(std::sqrt(pow((xlist[index[j][1]]-xlist[neighbour[j][1]]),2)+pow((ylist[index[j][0]]-ylist[neighbour[j][0]]),2)));
            value.push_back(0);
            boundarymap[index[j][0]][index[j][1]] = true;
            bcTypemap[index[j][0]][index[j][1]] = 1;
            bcvaluemap[index[j][0]][index[j][1]] = value[j];
        }
        boundaryPatch bp(index,neighbour,value,dist,type);
        bd.push_back(bp);
    }
    //right
    {
        std::vector<std::vector<int> > index,neighbour;
        std::vector<double> value;
        std::vector<double> dist;
        int type = 1;
        for(int j=0;j<Ny;j++)
        {
            index.push_back(std::vector<int>(2));
            neighbour.push_back(std::vector<int>(2));
            index[j][0] = j;
            index[j][1] = Nx-1;
            neighbour[j][0] = j;
            neighbour[j][1] = Nx-2;
            dist.push_back(std::sqrt(pow((xlist[index[j][1] ]-xlist[neighbour[j][1]]),2)+pow((ylist[index[j][0]]-ylist[neighbour[j][0]]),2)));
            value.push_back(0);
            boundarymap[index[j][0]][index[j][1]] = true; 
            bcTypemap[index[j][0]][index[j][1]] = 1;
            bcvaluemap[index[j][0]][index[j][1]] = value[j];
        }
        boundaryPatch bp(index,neighbour,value,dist,type);
        bd.push_back(bp);
    }
    Mesh mesh(xlist,ylist,Nx,Ny,bd,boundarymap,bcTypemap,bcvaluemap);  

    
    
    Field phi(mesh),A(mesh),mu(mesh);
    A.initialize(0);
    phi.initialize(0.000001);
    mu.initialize(0);
   // std::cout<<"x="<<mesh.nx<<std::endl;
    
    phi.bcUpdate();
    for(int i=0;i<mesh.ny;i++)
    {
        A.value[i][0] = 1-pow(M_inf,2);
        for(int j=1;j<mesh.nx-1;j++)
        {
            double dphidx = (phi.value[i][j+1] - phi.value[i][j-1] )/(xlist[j+1]-xlist[j-1]);
            A.value[i][j] = 1-pow(M_inf,2)-(gamma+1)*pow(M_inf,2)*dphidx/v_inf;
        }
        A.value[i][mesh.nx-1] = 1-pow(M_inf,2);
    }
    double eps = 1e-6;
    for(int i=0;i<mesh.ny;i++)
    {
        for(int j=0;j<mesh.nx-1;j++)
        {
            if(A.value[i][j]<0 && std::abs(A.value[i][j])>eps)
            {
                mu.value[i][j] = 0;
            }
            else if(A.value[i][j]<0 && std::abs(A.value[i][j])>eps)
            {
                mu.value[i][j] = 1;
            }
            else
            {
                mu.value[i][j]=eps;
            }
             
        }
    }

    for(int i=0;i<400;i++)
    {
         phi.traniter(A);
        if((i+1)%10==0)
        {
   //        std::cout/*<<"i="<<i<<"   "<<"Res = "*/<<std::abs(phi.res)<<std::endl;
          // std::cout<<"i="<<phi.ii<<"                    "<<"j="<<phi.jj<<std::endl;
        }
    }
    std::vector<double> u(mesh.nx),v(mesh.nx),M(mesh.nx),a(mesh.nx);
    std::vector<std::vector<double> > fu(mesh.ny),fv(mesh.ny),fa(mesh.ny),fm(mesh.ny);
   
    std::cout<<" =============================="<<std::endl;

    std::vector<double> cp(mesh.nx); 
    for(int j=1;j<mesh.nx-1;j++)
    {
        double u = (phi.value[1][j+1]-phi.value[1][j-1])/(phi.mesh.x[j+1]-phi.mesh.x[j-1])+v_inf;
        double v = (phi.value[2][j]-phi.value[0][j])/(phi.mesh.y[2]-phi.mesh.y[0]);
        double p =  p_inf * pow( (  1- (gamma-1)/2 * pow(M_inf,2)*((pow(u,2)+pow(v,2))/pow(v_inf,2)-1)  ), (gamma/(gamma-1))); 
         cp[j] = ( p - p_inf )/0.5/rho_inf/pow(v_inf,2);
    }
    for(int i=0;i<mesh.ny;i++)
    {
        fu[i]=(std::vector<double>(mesh.nx));
        fv[i]=(std::vector<double>(mesh.nx));
        fm[i]=(std::vector<double>(mesh.nx));
        fa[i]=(std::vector<double>(mesh.nx));
        fu[i][0] = v_inf;
        fv[i][0] = 0;
        fa[i][0] = a_inf;
        fm[i][0] = M_inf;
        
        for(int j=0;j<mesh.nx;j++)
        {
            if(j==0 || j == mesh.nx-1 || i == mesh.ny-1)
            {
                fu[i][mesh.nx-1] = v_inf;
                fv[i][mesh.nx-1] = 0;
                fa[i][mesh.nx-1] = a_inf;
                fm[i][mesh.nx-1] = M_inf;    
            }
            else
            {
                if(i==0)
                {
                    fv[i][j]==0;    
                }
                else
                {
                    fv[i][j]=( phi.value[i+1][j]-phi.value[i-1][j] )/(ylist[i+1]-ylist[i-1]);
                } 
                fu[i][j]=v_inf + ( phi.value[i][j+1]-phi.value[i][j-1] )/(xlist[j+1]-xlist[j-1]);

                double p =  p_inf * pow( (  1- (gamma-1)/2 * pow(M_inf,2)*((pow(fu[i][j],2)+pow(fv[i][j],2))/pow(v_inf,2)-1)  ), (gamma/(gamma-1)));
                double rho =  rho_inf * pow( (  1- (gamma-1)/2 * pow(M_inf,2)*((pow(fu[i][j],2)+pow(fv[i][j],2))/pow(v_inf,2)-1)  ), (1/(gamma-1)));
                fa[i][j]=sqrt(p/rho*gamma);
                fm[i][j]=fu[i][j]/fa[i][j];
            } 
            
        }
        
    }
    std::ofstream foutfull;
    
    foutfull.open("C:\\Users\\Zousf\\Desktop\\AEM8251\\Codes\\data.csv",ios::out);

    for (int i = 0;i < Ny;i++)
    {
        for (int j = 0;j < Nx;j++)
        {
            foutfull << xlist[j] << "," << ylist[i] << "," << fu[i][j] << "," << fv[i][j] <<","<<fm[i][j] <<std::endl;
        }
    }
    foutfull.close();

    cp[0]=cp[1];
     
    cp[phi.mesh.nx]=cp[phi.mesh.nx-1];

    for(int i = N1;i<mesh.nx-N1;i++)
    {
        std::cout<</*"x="<<xlist[i]<<"       "<< */cp[i]<<std::endl;
    }
    
    getchar();
    //Field res(mesh);
    return 0;
}