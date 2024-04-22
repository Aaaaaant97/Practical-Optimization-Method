#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip> 
#define theta 10000.0
#define PI 3.1415926535

using namespace std;
using namespace Eigen;

double f(Vector4d x){
        double x1=x(0),x2=x(1),x3=x(2),x4=x(3),ans;
        double a=1/2*(x1*x1+x2*x2+x3*x3+x4*x4);
        double b=x1*(5*x1+x2+1/2*x4)+(x1+4*x2+1/2*x3)*x2+(1/2*x2+3*x3)*x3+(1/2*x1+2*x4)*x4;
        ans=a+theta/4*b;
        return ans;
}

Vector4d grad(Vector4d x){
        Vector4d ans;
        double x1=x(0),x2=x(1),x3=x(2),x4=x(3);
        ans(0)=x1+theta/4*(5*x1+x2+1/2*x4+5*x1+x2+1/2*x4);
        ans(1)=x2+theta/4*(x1+4*x2+x1+4*x2+1/2*x3+1/2*x3);
        ans(2)=x3+theta/4*(1/2*x2+3*x3+1/2*x2+3*x3);
        ans(3)=x4+theta/4*(1/2*x1+2*x4+1/2*x1+2*x4);
        return ans;
}

Matrix4d hessian(Vector4d x){
        double x1=x(0),x2=x(1),x3=x(2),x4=x(3);
        Matrix4d hess;
        hess(0,0)=1+theta/2*5;
        hess(0,1)=theta/2;
        hess(0,2)=0;
        hess(0,3)=theta/4;
        hess(1,0)=theta/2;
        hess(1,1)=1+2*theta;
        hess(1,2)=theta/4;
        hess(1,3)=0;
        hess(2,0)=0;
        hess(2,1)=theta/4;
        hess(2,2)=1+3/2*theta;
        hess(2,3)=0;
        hess(3,0)=theta/4;
        hess(3,1)=0;
        hess(3,2)=0;
        hess(3,3)=1+theta;
        return hess;
}


void Newton(Vector4d x,double e){
        int iter=0;
        while(iter<=100000){
                Vector4d p=-hessian(x).inverse()*grad(x);
                double alpha=1;
                x+=alpha*p;
                if(grad(x).norm()<e){
                        break;
                }
                cout<<f(x)<<endl;
                iter++;
        }
        /*
        cout<<"The number of iterations is "<<iter<<endl;
        cout <<"Minimum found at x1 = "<<fixed<<setprecision(8)<<x(0) << ", x2 = " << x(1)  <<  ", x3 = " << x(2)  << ", x4 = " << x(3)  <<endl;
        cout<<"The value of the function is "<<f(x)<<endl;
        cout<<"hess:"<<endl;
                cout<<hessian(x)(0,0)<<"  "<<hessian(x)(0,1)<<"  "<<hessian(x)(0,2)<<"  "<<hessian(x)(0,3)<<endl;
                cout<<hessian(x)(1,0)<<"  "<<hessian(x)(1,1)<<"  "<<hessian(x)(1,2)<<"  "<<hessian(x)(1,3)<<endl;
                cout<<hessian(x)(2,0)<<"  "<<hessian(x)(2,1)<<"  "<<hessian(x)(2,2)<<"  "<<hessian(x)(2,3)<<endl;
                cout<<hessian(x)(3,0)<<"  "<<hessian(x)(3,1)<<"  "<<hessian(x)(3,2)<<"  "<<hessian(x)(3,3)<<endl;
                */
}

int main(){
        double e=1e-6;
        double ang=50.0*PI/180.0;
        double a=cos(ang),b=sin(ang);
        Vector4d x(a,b,a,b);
        Newton(x,e);
        return 0;
}