#ifndef FIELD_H
#define FIELD_H
#include <iostream>
#include <vector>
#include "Matrix.H"
#include <cmath>
#include <map>
#include "Solve.H"
#include "eigen/Eigen/Core"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Sparse"

enum SCHEMETYPE
    {
        unknown          = 0,
        pJ               = 1, 
        pGS              = 2,
        lJ               = 3,
        gSLR             = 4,
        ADI              = 5
    };
static std::map<std::string,SCHEMETYPE> typeList = 
{
    {"pointJacobi",pJ},
    {"pointGaussSeidel",pGS},
    {"lineJacobi",lJ},
    {"gaussSeidelLineRelaxation",gSLR},
    {"ADITwoStep",ADI}
};


class Field
{
    public:
    std::vector<std::vector<double> > value;
    std::vector<std::vector<double> > value1;
    const Mesh& mesh;
    double res;
    int ii,jj;
    int scheme;
    // constructors 
    Field(const Mesh& mesh)
    :mesh(mesh)
    {}

    Field(const Field& Field_1)
    :mesh(Field_1.mesh)
    {
        for(int i=0;i<mesh.nx;++i)
        {
            value.push_back(std::vector<double>() );
            for(int j=0;j<mesh.ny;++j)
            {
                value[i].push_back(Field_1.value[i][j]);
            }
        }  
        value1.assign(value.begin(),value.end() );  
    }

    ~Field()
    {}    
    
    void initialize(double uniformValue)
    {
        for(int i=0;i<mesh.nx;++i)
        {
            value.push_back(std::vector<double>());
            for(int j=0;j<mesh.ny;++j)
            {
                value[i].push_back(uniformValue);
            }
        }
        value1.assign(value.begin(),value.end() );   
    }
    
    void iter(const std::string scheme,Field& A);
    void traniter(Field& A);
    void bcUpdate()
    {
        for(int i=0;i<mesh.boundary.size();i++)
        {
            if(mesh.boundary[i].type == 1)
            {
                for(int j=0;j<mesh.boundary[i].index.size();j++)
                {
                    value[mesh.boundary[i].index[j][0]][mesh.boundary[i].index[j][1]] = mesh.boundary[i].bcValue[j];
                    value1[mesh.boundary[i].index[j][0]][mesh.boundary[i].index[j][1]] = value[mesh.boundary[i].index[j][0]][mesh.boundary[i].index[j][1]];
                }
            }
            else if(mesh.boundary[i].type == 2)
            { 
                for(int j=0;j<mesh.boundary[i].index.size();j++)
                {
                    value[mesh.boundary[i].index[j][0]][mesh.boundary[i].index[j][1]] =  value[mesh.boundary[i].neighbour[j][0]][mesh.boundary[i].neighbour[j][1]] +  mesh.boundary[i].dist[j] * mesh.boundary[i].bcValue[j];
                    value1[mesh.boundary[i].index[j][0]][mesh.boundary[i].index[j][1]] =  value[mesh.boundary[i].index[j][0]][mesh.boundary[i].index[j][1]];
                }
            }
            else
            {
                std::cout<<"wrong boundary type"<<std::endl;
                exit;
            }
                        
        }
    }
    void subsonicRes(Field& A)
    {
        res = 0;
        for(int i=1;i<mesh.ny-1;i++)
        {
            for(int j=1;j<mesh.nx-1;j++)
            {
                double dxp = mesh.x[j+1]-mesh.x[j];
                double dxm = mesh.x[j]-mesh.x[j-1];
                double dyp = mesh.y[i+1]-mesh.y[i];
                double dym = mesh.y[i]-mesh.y[i-1];
                double dxc = mesh.x[j+1]-mesh.x[j-1];
                double dyc = mesh.y[i+1]-mesh.y[i-1];

                double R = 2.0 * A.value[i][j]*(value[i][j+1]/dxp/dxc-value[i][j]/dxp/dxm+value[i][j-1]/dxc/dxm)
                         + 2.0 * (value[i+1][j]/dyp/dyc-value[i][j]/dyp/dym+value[i-1][j]/dyc/dym);

                if(std::abs(res)<std::abs(R))
                {
                    res = R;
                }
            }
        }
    }
    void transonicRes(Field& A)
    {
        res = 0;
        ii = 0;
        jj = 0;
        for(int i=1;i<mesh.ny-1;i++)
        {
            for(int j=1;j<mesh.nx-1;j++)
            {
                double dxp = mesh.x[j+1]-mesh.x[j];
                double dxm = mesh.x[j]-mesh.x[j-1];
                double dxcm = mesh.x[j]-mesh.x[j-2];
                double dxmm = mesh.x[j-1] - mesh.x[j-2];
                double dyp = mesh.y[i+1]-mesh.y[i];
                double dym = mesh.y[i]-mesh.y[i-1];
                double dxc = mesh.x[j+1]-mesh.x[j-1];
                double dyc = mesh.y[i+1]-mesh.y[i-1];
                
                double R=0;
                if(A.value[i][j]>0)
                {
                    R += 2.0 * A.value[i][j]*(value[i][j+1]/dxp/dxc-value[i][j]/dxp/dxm+value[i][j-1]/dxc/dxm);
                }
                if(j>1 && A.value[i][j-1]<0)
                {
                    R += 2.0 * A.value[i][j-1]*(value[i][j]/dxm/dxcm-value[i][j-1]/dxm/dxmm+value[i][j-2]/dxcm/dxmm);
                }
                R += 2.0 * (value[i+1][j]/dyp/dyc-value[i][j]/dyp/dym+value[i-1][j]/dyc/dym);

                if(std::abs(res)<std::abs(R))
                {
                    res = R;
                    ii = i;
                    jj = j;
                    
                }
                
            }
        }
    }
};





void pointJacobi(Field& phi,Field& A)
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

                double aii= 0.5 /(A.value[i][j]/dxp/dxm+1/dyp/dym);
                phi.value1[i][j] = aii * (2.0*A.value[i][j]/dxp/dxc * phi.value[i][j+1] + 2.0*A.value[i][j]/dxm/dxc * phi.value[i][j-1] + 2.0/dyp/dyc * phi.value[i+1][j]+ 2.0/dym/dyc * phi.value[i-1][j]);
            }
        } 
    }
}

void pointGaussSeidel(Field& phi,Field& A)
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

                double aii= 0.5 /(A.value[i][j]/dxp/dxm+1/dyp/dym);
                phi.value1[i][j] = aii * (2.0*A.value[i][j]/dxp/dxc * phi.value[i][j+1] + 2.0*A.value[i][j]/dxm/dxc * phi.value1[i][j-1] + 2.0/dyp/dyc * phi.value[i+1][j]+ 2.0/dym/dyc * phi.value1[i-1][j]);                                                                                    
            }
        }
    }
}

void lineJacobi(Field& phi,Field& A)
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
                    d[i]=( 2.0*A.value[i][j]/dxp/dxc * phi.value[i][j+1] + 2.0*A.value[i][j]/dxm/dxc * phi.value[i][j-1]);
                    a[i] = (2.0/dyp/dym + 2.0 * A.value[i][j]/dxp/dxm);
                    d[i]+=(2.0/dyp/dyc)*phi.value[i+1][j];
                    c[i-1]=(-2.0/dym/dyc);
                }
                else
                {
                    d[i]=( 2.0*A.value[i][j]/dxp/dxc * phi.value[i][j+1] + 2.0*A.value[i][j]/dxm/dxc * phi.value[i][j-1]);
                    a[i] = (2.0/dyp/dym + 2.0 * A.value[i][j]/dxp/dxm);
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

void gaussSeidelLineRelaxation(Field& phi,Field& A)
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
                    d[i]=( 2.0*A.value[i][j]/dxp/dxc * phi.value[i][j+1] + 2.0*A.value[i][j]/dxm/dxc * phi.value1[i][j-1]);
                    a[i] = (2.0/dyp/dym + 2.0 * A.value[i][j]/dxp/dxm);
                    d[i]+=(2.0/dyp/dyc)*phi.value[i+1][j];
                    c[i-1]=(-2.0/dym/dyc);
                }
                else
                {
                    d[i]=( 2.0*A.value[i][j]/dxp/dxc * phi.value[i][j+1] + 2.0*A.value[i][j]/dxm/dxc * phi.value1[i][j-1]);
                    a[i] = (2.0/dyp/dym + 2.0 * A.value[i][j]/dxp/dxm);
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

void ADITwoStep(Field& phi,Field& A)
{
    gaussSeidelLineRelaxation(phi,A);
    phi.value.swap(phi.value1);
    int n = phi.mesh.nx;
    using namespace Eigen;
    Eigen::VectorXd soln(n), matB(n);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    for(int i=1;i<phi.mesh.ny-1;i++)
    {
    //  a.clear
     //   nstart = 0;
      //  nend = n-1;
        int j=0;  
        tripletList.clear();
        if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 1)
        {
             //   a[i] = 1;
            tripletList.push_back(T(j, j, 1.0));
            matB(j) = phi.mesh.bcValueMap[i][j];
        }
        else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 2)  
        {
            tripletList.push_back(T(j, j, 1.0));
            tripletList.push_back(T(j, j+1, -1.0));
            matB(j) = phi.mesh.bcValueMap[i][j];
        }
        else
        {
            std::cout<<"wrong BC type"<<std::endl;
            exit;
        }
        for(j++;j<phi.mesh.nx-1;j++)
        {  
            double dxp = phi.mesh.x[j+1]-phi.mesh.x[j];
            double dxm = phi.mesh.x[j]-phi.mesh.x[j-1];
            double dxcm = phi.mesh.x[j]-phi.mesh.x[j-2];
            double dxmm = phi.mesh.x[j-1] - phi.mesh.x[j-2];
            double dyp = phi.mesh.y[i+1]-phi.mesh.y[i];
            double dym = phi.mesh.y[i]-phi.mesh.y[i-1];
            double dxc = phi.mesh.x[j+1]-phi.mesh.x[j-1];
            double dyc = phi.mesh.y[i+1]-phi.mesh.y[i-1];
           
            tripletList.push_back(T(j, j, (A.value[i][j]*2.0/dxp/dxm + 2.0/dyp/dym)));
            tripletList.push_back(T(j, j - 1, -2.0*A.value[i][j]/dxm/dxc));
            tripletList.push_back(T(j, j + 1, -2.0*A.value[i][j]/dxp/dxc));
            
            matB(j) = 2.0/dyp/dyc * phi.value[i+1][j] + 2.0/dym/dyc * phi.value[i-1][j];
        }
        if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 1)
        {
             //   a[i] = 1;
            tripletList.push_back(T(j, j, 1.0));
            matB(j) = phi.mesh.bcValueMap[i][j];
        }
        else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 2)  
        {
            tripletList.push_back(T(j, j, 1.0));
            tripletList.push_back(T(j, j - 1, -1.0));
            matB(j) = phi.mesh.bcValueMap[i][j];
        }
        else
        {
            std::cout<<"wrong BC type"<<std::endl;
            exit;
        }

        SparseMatrix<double> matA(n,n);
        matA.setFromTriplets(tripletList.begin(), tripletList.end());

        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
        solver.compute(matA);
        soln = solver.solve(matB);
        std::vector<double> x(n),b(n);
        for(int k=0;k<phi.mesh.nx;k++)
        {
            x[k] = soln(k);
            b[k] = matB(k);
        }


       // nchase(a,b,c,d,soln,nstart,nend);
        for(int k=0;k<phi.mesh.nx;k++)
        {
            phi.value[i][k]=soln(k);
        }
              
    }
}


void transonicGSLR(Field& phi, const Field& A)
{
    int n = phi.mesh.ny;
    using namespace Eigen;
  //  int nstart =0;
  //  int nend = n-1;
    Eigen::VectorXd soln(n), matB(n);
    typedef Eigen::Triplet<double> T;
    std::vector<T > tripletList;

    for(int j=1;j<phi.mesh.nx-1;j++)
    {
    //  a.clear
     //   nstart = 0;
      //  nend = n-1;
        int i=0;  
        tripletList.clear();
        if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 1)
        {
             //   a[i] = 1;
            tripletList.push_back(T(i, i, 1.0));
            matB(i) = phi.mesh.bcValueMap[i][j];
        }
        else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 2)  
        {
            tripletList.push_back(T(i, i, 1.0));
            tripletList.push_back(T(i, i+1, -1.0));
            matB(i) = phi.mesh.bcValueMap[i][j];
        }
        else
        {
            std::cout<<"wrong BC type"<<std::endl;
            exit;
        }
        for(i++;i<phi.mesh.ny-1;i++)
        {  
            double dxp = phi.mesh.x[j+1]-phi.mesh.x[j];
            double dxm = phi.mesh.x[j]-phi.mesh.x[j-1];
            double dxcm = phi.mesh.x[j]-phi.mesh.x[j-2];
            double dxmm = phi.mesh.x[j-1] - phi.mesh.x[j-2];
            double dyp = phi.mesh.y[i+1]-phi.mesh.y[i];
            double dym = phi.mesh.y[i]-phi.mesh.y[i-1];
            double dxc = phi.mesh.x[j+1]-phi.mesh.x[j-1];
            double dyc = phi.mesh.y[i+1]-phi.mesh.y[i-1];

            double aa = 0, dd = 0;
                
            if(A.value[i][j]>0)
            {
                aa +=  ( A.value[i][j]*2.0/dxp/dxm);
                dd += A.value[i][j]*( 2.0/dxp/dxc * phi.value[i][j+1] + 2.0/dxm/dxc * phi.value1[i][j-1]);
            }
            if(j>1 && A.value[i][j-1]<0)
            {
                aa += -2.0*A.value[i][j-1]/dxm/dxcm; 
                dd += 2.0*A.value[i][j-1] * (-1.0*phi.value1[i][j-1]/dxm/dxmm+phi.value1[i][j-2]/dxcm/dxmm); 
                
            }
            tripletList.push_back(T(i, i, aa+2.0/dyp/dym));
            tripletList.push_back(T(i, i - 1, -2.0/dym/dyc ));
            tripletList.push_back(T(i, i + 1, -2.0/dyp/dyc));
            
            matB(i) = dd;
        }
        if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 1)
        {
             //   a[i] = 1;
            tripletList.push_back(T(i, i, 1.0));
            matB(i) = phi.mesh.bcValueMap[i][j];
        }
        else if(phi.mesh.onBoundary[i][j] == 1 && phi.mesh.bcType[i][j] == 2)  
        {
            tripletList.push_back(T(i, i, 1.0));
            tripletList.push_back(T(i, i - 1, -1.0));
            matB(i) = phi.mesh.bcValueMap[i][j];
        }
        else
        {
            std::cout<<"wrong BC type"<<std::endl;
            exit;
        }

        SparseMatrix<double> matA(n,n);
        matA.setFromTriplets(tripletList.begin(), tripletList.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
        solver.compute(matA);
        soln = solver.solve(matB);
        std::vector<double> x(n),b(n);
        for(int k=0;k<phi.mesh.nx;k++)
        {
            x[k] = soln(k);
            b[k] = matB(k);
        }


       // nchase(a,b,c,d,soln,nstart,nend);
        for(int k=0;k<phi.mesh.nx;k++)
        {
            phi.value1[k][j]=soln(k);
        }      
    }
}
void Field::iter(const std::string scheme, Field& A)
{
    switch(typeList[scheme])
    {
        case pJ:
            pointJacobi(*this,A);
            value.swap(value1);
            this->bcUpdate();
            break;
        case pGS:
            pointGaussSeidel(*this,A);
            value.swap(value1);
            this->bcUpdate();
            break;
        case lJ:
            lineJacobi(*this,A);
            value.swap(value1);
            break;
        case gSLR:
            gaussSeidelLineRelaxation(*this,A);
            value.swap(value1);
            break;
        case ADI:
            ADITwoStep(*this,A);
            this->bcUpdate();
          //  value.swap(value1);
            break;
        default:
            std::cout<<"Unknown scheme"<<std::endl;
            exit;
    }
    subsonicRes(A);
//    value.swap(value1);
//    this->bcUpdate();
    //value.swap(value1);
}
void Field::traniter(Field& A)
{
    transonicGSLR(*this,A);
    value.swap(value1);
    for(int i=0;i<mesh.ny-1;i++)
    {
        for(int j=1;j<mesh.nx-1;j++)
        {
        //    if(A.value[i][j]>0)
          //  {
                A.value[i][j] = 1-pow(M_inf,2)-(gamma+1)*pow(M_inf,2)/v_inf*(this->value[i][j+1] - this->value[i][j-1] )/(this->mesh.x[j+1]-this->mesh.x[j-1]);    
                if(A.value[i][j]<0)
                {
                    double u = (this->value[i][j+1] - this->value[i][j-1] )/(this->mesh.x[j+1]-this->mesh.x[j-1]); 
                 //   std::cout<<"i="<<i<<"    j="<<j<<std::endl;
                  //  std::cout<<"velocity="<<u<<std::endl;
                }
          //  }
          //  else if(A.value[i][j]<0 && j>1)
          //  {
          //      A.value[i][j] = 1 - pow(M_inf,2)-(gamma+1)*pow(M_inf,2)/v_inf*(this->value[i][j] - this->value[i][j-2] )/(this->mesh.x[j]-this->mesh.x[j-2]);
          //  }
            
        }
 //   A.value[i][mesh.nx-1] = 1-pow(M_inf,2);
    }
    transonicRes(A);
}



#endif
