#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip> 

using namespace std;
using namespace Eigen;

double f(Vector2d x){
        double x1=x(0),x2=x(1);
        double a=1.5-x1+x1*x2,b=2.25-x1+x1*pow(x2,2),c=2.625-x1+x1*pow(x2,3),ans;
        ans=a*a+b*b+c*c;
        return ans;
}

Vector2d grad(Vector2d x){
        Vector2d ans;
        double x1=x(0),x2=x(1);
        double a=1.5-x1+x1*x2,b=2.25-x1+x1*pow(x2,2),c=2.625-x1+x1*pow(x2,3);
        ans(0)=2*a*(x2-1)+2*b*(pow(x2,2)-1)+2*c*(pow(x2,3)-1);
        ans(1)=2*a*x1+4*x1*x2*b+6*x1*x2*x2*c;
        return ans;
}

Matrix2d hessian(Vector2d x){
        double x1 = x(0),x2 = x(1);
        Matrix2d hess;
        hess(0, 0)=2*pow(x2-1,2)+2*pow(x2*x2-1,2)+2*pow(x2*x2*x2-1,2);
        hess(0, 1)=2*(1.5-x1+x1*x2)+2*(x2-1)*x1+4*(2.25-x1+x1*x2*x2)*x2+4*(x2*x2-1)*x1*x2+6*(2.2625-x1+x1*x2*x2*x2)*x2*x2+6*(x2*x2*x2-1)*x1*x2*x2;
        hess(1, 0)=2*(1.5-x1+x1*x2)+2*(x2-1)*x1+4*(2.25-x1+x1*x2*x2)*x2+4*(x2*x2-1)*x1*x2+6*(2.2625-x1+x1*x2*x2*x2)*x2*x2+6*(x2*x2*x2-1)*x1*x2*x2;
        hess(1, 1)=2*x1*x1+4*(2.25-x1+x1*x2*x2)*x1+8*pow(x1*x2,2)+12*(2.2625-x1+x1*x2*x2*x2)*x1*x2+18*pow(x1*x2*x2,2);
        return hess;
}

double LineSearch(Vector2d x,Vector2d p){
        double alpha=1,c=1e-4;
        Vector2d y=x+alpha*p;
        while(f(y)>f(x)+c*alpha*(grad(x)(0)*p(0)+grad(x)(1)*p(1))){
                alpha*=0.5;
                y=x+alpha*p;
        }
        return alpha;
}

void Newton(Vector2d x,double e){
        int iter=0;
        while(iter<=100000){
                Vector2d p=-hessian(x).inverse()*grad(x);
                double alpha=LineSearch(x,p);
                x+=alpha*p;
                if(grad(x).norm()<e){
                        break;
                }
                cout<<f(x)<<" "<<x(0)<<" "<<x(1)<<endl;
                iter++;
        }
        /*
        cout<<"The number of iterations is "<<iter<<endl;
        cout <<"Minimum found at x1 = "<<fixed<<setprecision(10)<<x(0) << ", x2 = " << x(1)  << endl;
        cout<<"The value of the function is "<<f(x)<<endl;
        cout<<"hess:"<<endl;
                cout<<hessian(x)(0,0)<<"  "<<hessian(x)(0,1)<<endl;
                cout<<hessian(x)(1,0)<<"  "<<hessian(x)(1,1)<<endl;
                cout<<endl;
                */
}

int main(){
        double e=1e-6;
        Vector2d x(-1.2,1); 
        Newton(x,e);
        return 0;
}