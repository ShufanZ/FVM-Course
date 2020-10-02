#include <iostream>
#include <cmath>
#include <vector>
const double M_inf = 0.735;
const double A0 = 1 - pow(M_inf,2);
const double gamma = 1.4, rho_inf = 1, p_inf = 1; 
const double a_inf = std::sqrt(gamma * p_inf/rho_inf);
const double v_inf = a_inf * M_inf;
#include "Field.H"
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
        f  = xtotal * (k - 1 ) - xmin * (pow(k,n) - 1);
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
            xlist[i] = xlist[i-1] + 0.05;
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
            value.push_back(v_inf * xlist[index[j][1]]);
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
            value.push_back(v_inf * xlist[index[j][1]]);
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
            value.push_back(v_inf * xlist[index[j][1] ]);
            boundarymap[index[j][0]][index[j][1]] = true; 
            bcTypemap[index[j][0]][index[j][1]] = 1;
            bcvaluemap[index[j][0]][index[j][1]] = value[j];
        }
        boundaryPatch bp(index,neighbour,value,dist,type);
        bd.push_back(bp);
    }
    Mesh mesh(xlist,ylist,Nx,Ny,bd,boundarymap,bcTypemap,bcvaluemap);  

    
    
    Field phi(mesh),A(mesh);
    phi.initialize(0);
    A.initialize(A0);
   // std::cout<<"x="<<mesh.nx<<std::endl;
    for(int i=0;i<mesh.ny;i++)
    {
        for(int j=0;j<mesh.nx;j++)
        {
            phi.value[i][j] = v_inf * xlist[j];
            phi.value1[i][j] = phi.value[i][j];
        }
    }
    

    
    phi.bcUpdate();
    std::string scheme="gaussSeidelLineRelaxation";
    
    for(int i=0;i<1000;i++)
    {
        phi.iter(scheme,A);
        if((i+1)%10==0)
        {
            std::cout<</*"i="<<(i+1)<<"   "<<"Res = "<<*/std::abs(phi.res)<<std::endl;
        }
    }
  std::cout<<"============================================="<<std::endl;
    std::vector<double> cp(mesh.nx); 
    for(int j=1;j<mesh.nx-1;j++)
    {
        double u = (phi.value[1][j+1]-phi.value[1][j-1])/(phi.mesh.x[j+1]-phi.mesh.x[j-1]);
        double v = (phi.value[2][j]-phi.value[0][j])/(phi.mesh.y[2]-phi.mesh.y[0]);
        double p =  p_inf * pow( (  1- (gamma-1)/2 * pow(M_inf,2)*((pow(u,2)+pow(v,2))/pow(v_inf,2)-1)  ), (gamma/(gamma-1))); 
         cp[j] = ( p - p_inf )/0.5/rho_inf/pow(v_inf,2);
    }
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