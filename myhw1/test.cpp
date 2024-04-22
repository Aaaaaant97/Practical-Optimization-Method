#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip> 

using namespace std;
using namespace Eigen;

double f(Vector2d x){
        double x1=x(0),x2=x(1);
        double a=1.5-x1+x1*x2,b=2.25-x1+x1*x2*x2,c=2.625-x1+x1*x2*x2*x2,ans;
        ans=a*a+b*b+c*c;
        return ans;
}

Vector2d grad(Vector2d x){
        Vector2d ans;
        double x1 = x(0);
        double x2 = x(1);
        ans(0) = -2 * (1.5 - x1 + x1 * x2) * (x2 - 1) - 2 * (2.25 - x1 + x1 * pow(x2, 2)) * (x2 * x2 - 1) - 2 * (2.625 - x1 + x1 * pow(x2, 3)) * (x2 * x2 * x2 - 1);
        ans(1) = -2 * (1.5 - x1 + x1 * x2) * x1 - 4 * (2.25 - x1 + x1 * pow(x2, 2)) * x1 * x2 - 6 * (2.625 - x1 + x1 * pow(x2, 3)) * x1 * pow(x2, 2);
        return ans;
}

Matrix2d hessian(Vector2d x){
        double x1 = x(0);
        double x2 = x(1);
        Matrix2d hess;
        hess(0, 0) = 2 * pow(x2 - 1, 2) + 2 * pow(x2 * x2 - 1, 2) + 2 * pow(x2 * x2 * x2 - 1, 2);
        hess(0, 1) = -2 * (1.5 - x1 + x1 * x2) + 4 * x1 * (2.25 - x1 + x1 * pow(x2, 2)) + 6 * x1 * pow(x2, 2) * (2.625 - x1 + x1 * pow(x2, 3));
        hess(1, 0) = -2 * (1.5 - x1 + x1 * x2) + 4 * x1 * (2.25 - x1 + x1 * pow(x2, 2)) + 6 * x1 * pow(x2, 2) * (2.625 - x1 + x1 * pow(x2, 3));
        hess(1, 1) = 2 * pow(x1, 2) + 4 * pow(x1 * x2, 2) + 6 * pow(x1 * pow(x2, 2), 2);
        return hess;
}

double LineSearch(Vector2d x,Vector2d p){
        double alpha=1,c=1e-4;
        Vector2d y=x+alpha*p;
        while(f(y)>f(x)-c*alpha*(p(0)*p(0)+p(1)*p(1))){
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
                iter++;
                cout<<iter<<":"<<x(0)<<","<<x(1)<<endl;
        }
        cout << "Minimum found at x1 = "<<fixed<<setprecision(10)<< x(0) << ", x2 = " << x(1) << " after " << iter << " iterations." << endl;
}

int main(){
        double e=1e-6;
        Vector2d x(-1.2,1); 
        Newton(x,e);
        return 0;
}