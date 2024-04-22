#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip> 

using namespace std;
using namespace Eigen;

double f(Vector2d x){
        double a=x(1)-x(0)*x(0),b=1-x(0),ans;
        ans=100*pow(a,2)+pow(b,2);
        return ans;
}

Vector2d grad(Vector2d x){
        Vector2d ans;
        ans(0)=-400*(x(1)-x(0)*x(0))*x(0)-2*(1-x(0));
        ans(1)=200*(x(1)-x(0)*x(0));
        return ans;
}


Matrix2d hessian(Vector2d x){
        Matrix2d hess;
        hess(0,0)=-400*(x(1)-3*x(0)*x(0))+2;
        hess(0,1)=-400*x(0);
        hess(1,0)=-400*x(0);
        hess(1,1)=200;
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
        while(true){
                Vector2d p=-hessian(x).inverse()*grad(x);
                double alpha=LineSearch(x,p);
                x+=alpha*p;
                if(grad(x).norm()<e){
                        break;
                }
                cout<<f(x)<<endl;
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