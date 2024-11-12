#include <bits/stdc++.h>

class LinearEquation
{
    //This class is built for solve linear equations
    public:
        int n;  //scale
        std::vector<std::vector<double>> A;
        std::vector<double> b;  //Ax=b
        /*
        member functions:
            MultiA: calculate z=Ax
            TransF: D=A.T
            Hilbert: initialize A as a Hilbert matrix
            displayA: print A; displayVec: print any vector with dimension n
            NormL1: calculate the L1-norm of two vectors with common dimension n
        methods:
            DooliteLU,CholeskyLU: LU decomposition
            UpperSolve,LowerSolve: unsafe, be only used in LU solvers
            Jacobi iteration, Gauss-Seider iteration, common gradient descent
        */
        LinearEquation(int m)
        {
            n=m;
            b.resize(n,0);
            A.resize(n);
            for(auto it=A.begin();it!=A.end();++it)
                (*it).resize(n,0);
        }

        void MultiA(std::vector<double> x, std::vector<double>& z)
        {
            //calculate z=Ax
            for(int i=0;i<n;i++)
            {
                z[i]=0;
                for(int j=0;j<n;j++)
                    z[i]+=A[i][j]*x[j];
            }
        }
        void TransF(LinearEquation& D)
        {
            //D=A'
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<n;j++)
                    (D.A)[i][j]=A[j][i];
            }
        }
        void Hilbert()
        {
            //Sommon A as Hilbert(n)
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<n;j++)
                    A[i][j]=1/double(i+j+1);
            }
        }
        void displayA()
        {
            //Display A
            for(auto i=A.begin();i!=A.end();++i)
            {
                for(int j=0;j<n;j++)
                    printf("%1.4f ", (*i)[j]);
                std::cout<<std::endl;
            }
        }
        void displayVec(std::vector<double> x)
        {
            //Display any vector with dimension n
            for(auto i=x.begin();i!=x.end();++i)
                printf("%.12f ",*i);
            std::cout<<'\n';
        }
        double NormL1(std::vector<double> x,std::vector<double> y)
        {
            //L^1 distance of x and y
            double s=0;
            for(int i=0;i<n;i++)
                s+=abs(x[i]-y[i]);
            return s;
        }
        double InnProd(std::vector<double> x,std::vector<double> y)
        {
            //<x,y>=y'x
            double sum=0;
            for(int i=0;i<n;i++)
                sum+=x[i]*y[i];
            return sum;
        }

    //Doolite
        bool DooliteLU(std::vector<double>& x,LinearEquation& L,LinearEquation& U);
    //Cholesky
        bool CholeskyLU(std::vector<double>& x,LinearEquation& L,LinearEquation& U);

    //Jacobi iteration
        int JacobiIter(std::vector<double>& x,std::vector<double> xStart={0});
    //Gauss-Seider iteration
        int GSIter(std::vector<double>& x,std::vector<double> xStart={0});
    //Gradient Descent
        int GD(std::vector<double>& x,std::vector<double> xStart={0});

    private:
        const int MaxIt=5000;
        const double error=1e-5;
        void DiagSolve(std::vector<double>& x)
        {
            //Solve while A is of diagram form
            for(int i=0;i<n;i++)
                x[i]=b[i]/A[i][i];
        }
        void UpperSolve(std::vector<double>& x)
        {
            //Solve while A has only upper non-zero eles
            for(int i=n-1;i>=0;i--)
            {
                double s=0;
                for(int j=n-1;j>i;j--)
                    s+=A[i][j]*x[j];
                x[i]=(b[i]-s)/A[i][i];
            }
        }
        void LowerSolve(std::vector<double>& x)
        {
            //Solve while A has only lower non-zero eles
            for(int i=0;i<n;i++)
            {
                double s=0;
                for(int j=0;j<i;j++)
                    s+=A[i][j]*x[j];
                x[i]=(b[i]-s)/A[i][i];
            }
        }
};

bool LinearEquation::DooliteLU(std::vector<double>& x,LinearEquation& L,LinearEquation& U)
{
    for(int i=0;i<n;i++)
    {
        (L.A)[i][i]=1;
        for(int j=i;j<n;j++)
        {
            (U.A)[i][j]=A[i][j];
            for(int l=0;l<i;l++)
                (U.A)[i][j]-=(L.A)[i][l]*(U.A)[l][j];
        }
        for(int j=i+1;j<n;j++)
        {
            (L.A)[j][i]=A[j][i];
            for(int l=0;l<i;l++)
                (L.A)[j][i]-=(L.A)[j][l]*(U.A)[l][i];
            (L.A)[j][i]/=(U.A)[i][i];
        }
    }
    std::vector<double> y(n);
    L.b=b;
    L.LowerSolve(y);
    U.b=y;
    U.UpperSolve(x);
    return true;
}

bool LinearEquation::CholeskyLU(std::vector<double>& x,LinearEquation& L,LinearEquation& U)
{
    for(int i=0;i<n;i++)
    {
        (L.A)[i][i]=A[i][i];
        for(int j=0;j<i;j++)
            (L.A)[i][i]-=pow((L.A)[i][j],2);
        (L.A)[i][i]=sqrt((L.A)[i][i]);
        for(int j=i+1;j<n;j++)
        {
            (L.A)[j][i]=A[j][i];
            for(int k=0;k<i;k++)
                (L.A)[j][i]-=(L.A)[j][k]*(L.A)[i][k];
            (L.A)[j][i]/=(L.A)[i][i];
        }
    }
    L.TransF(U);
    std::vector<double> y(n);
    L.b=b;
    L.LowerSolve(y);
    U.b=y;
    U.UpperSolve(x);
    return true;
}

int LinearEquation::JacobiIter(std::vector<double>& x,std::vector<double> xStart)
{
    std::vector<double> y=xStart;
    if(y.size()!=n)
        y.resize(n,0);
    x=y;
    LinearEquation M(n),N(n);
    int i,j,k=0;
    for(i=0;i<n;i++)
    {
        (M.A)[i][i]=A[i][i];
        for(j=0;j<n;j++)
        {
            if(j!=i)
                (N.A)[i][j]=-A[i][j];
        }
    }
    while(k<=MaxIt)
    {
        N.MultiA(y,M.b);
        for(i=0;i<n;i++)
            (M.b)[i]+=b[i];
        M.DiagSolve(x);
        if(NormL1(x,y)<error)
            return k;
        y=x;
        k++;
    }
    return 0;
}

int LinearEquation::GSIter(std::vector<double>& x,std::vector<double> xStart)
{
    std::vector<double> y=xStart;
    if(y.size()!=n)
        y.resize(n,0);
    x=y;
    LinearEquation M(n),N(n);
    int i,j,k=0;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(j<=i)
                (M.A)[i][j]=A[i][j];
            if(j!=i)
                (N.A)[i][j]=-A[i][j];
        }
    }
    while(k<=MaxIt)
    {
        N.MultiA(y,M.b);
        for(i=0;i<n;i++)
            (M.b)[i]+=b[i];
        M.LowerSolve(x);
        if(NormL1(x,y)<error)
            return k;
        y=x;
        k++;
    }
    return 0;
}

int LinearEquation::GD(std::vector<double>& x,std::vector<double> xStart)
{
    std::vector<double> y=xStart;
    if(y.size()!=n)
        y.resize(n,0);
    x=y;
    std::vector<double> r(n,0);
    int i,k=0;
    while(k<=MaxIt)
    {
        MultiA(y,r);
        for(i=0;i<n;i++)
            r[i]=b[i]-r[i];
        std::vector<double> s(n,0);
        MultiA(r,s);
        double al=InnProd(r,r)/InnProd(s,r);
        for(i=0;i<n;i++)
            x[i]=y[i]+al*r[i];
        if(NormL1(x,y)<error)
            return k;
        y=x;
        k++;
    }
    return 0;
}