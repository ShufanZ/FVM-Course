#include <cmath>
#include <vector>

void chase(const std::vector<double>& a, const std::vector<double>& b,const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x)
{
    {
        std::vector<double> sol;
        int n = a.size();
        std::vector<double> y(n);
        y[0]=d[0];
        std::vector<double> L(n-1),u(n);
        u[0] = a[0];
        for(int i=1;i<n;i++)
        {
            L[i-1] = c[i-1]/u[i-1];
            u[i] = a[i]-L[i-1]*b[i-1];
            y[i]=d[i]-L[i-1] * y[i-1];
        }
        x[n-1] = y[n-1]/u[n-1];
        for(int i=n-2;i>=0;i--)
        {
            x[i] = (y[i]-b[i]*x[i+1])/u[i];
        }
    }
}

void nchase(const std::vector<double>& a, const std::vector<double>& b,const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x,int start, int end)
{
        int n = end-start+1;
        std::vector<double> sol;
        std::vector<double> y(n);
        y[start]=d[start];
        std::vector<double> L(n-1),u(n);
        u[start] = a[start];
        for(int i=start+1;i<n;i++)
        {
            L[i-1] = c[i-1]/u[i-1];
            u[i] = a[i]-L[i-1]*b[i-1];
            y[i]=d[i]-L[i-1] * y[i-1];
        }
        x[n-1] = y[n-1]/u[n-1];
        for(int i=n-2;i>=start;i--)
        {
            x[i] = (y[i]-b[i]*x[i+1])/u[i];
        }
}

void Jacobi(const std::vector<double>& a, const std::vector<double>& b,const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x)
{
    {
        int n = a.size();
        std::vector<double> sol(n);
        x[0]= 1;
        for(int i=1;i<n;i++)
        {
            x[i] = 0;
        }
        double res;
        for(;;)
        {
            res = 0;
            sol[0] = 1/a[0]*(-b[0]*x[1]+d[0]);
            for(int j=1;j<n-1;j++)
            {
                sol[j] =1/a[j]*(-c[j-1]*x[j-1]-b[j]*x[j+1]+d[j]);
            }
            sol[n-1] = 1/a[n-1] *(-c[n-2]*x[n-2]+d[n-1]);
            for(int j=0;j<n;j++)
            {
                res += pow( (x[j]-sol[j])/(x[j]),2);
                x[j]=sol[j];
            }
            if(res<0.000001)
            {
                break;
            }
        }
    }
}