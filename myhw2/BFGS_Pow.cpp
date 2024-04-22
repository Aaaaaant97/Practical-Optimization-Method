#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

double f(VectorXd x){
    double ans=0,a,b,c,d;
    a=x(0)+10*x(1);
    b=x(2)-x(3);
    c=x(1)-2*x(2);
    d=x(0)-x(3);
    ans=a*a+5*b*b+c*c*c*c+10*d*d*d*d;
    return ans;
}

VectorXd gradient(VectorXd x){
    VectorXd ans(x.size());
    double a,b,c,d;
    a=x(0)+10*x(1);
    b=x(2)-x(3);
    c=x(1)-2*x(2);
    d=x(0)-x(3);
    ans(0)=2*a+40*d*d*d;
    ans(1)=20*a+4*c*c*c;
    ans(2)=10*b-8*c*c*c;
    ans(3)=-10*b-40*d*d*d;
    return ans;
}

double LineSearch(VectorXd x,VectorXd p){
    double alpha=1,c1=1e-4,c2=0.9;
    double alpha_min=0,alpha_max=1;
    for(int i=0;i<10000;i++){
        VectorXd y=x+alpha*p;
        VectorXd grad=gradient(x);
        if(f(y)>f(x)+c1*alpha*grad.dot(p)){
            alpha*=0.5;
            alpha_max=alpha;
        }else if(gradient(y).dot(p)<c2*grad.dot(p)){
            alpha_min=alpha;
            alpha=(alpha_min+alpha_max)/2.0;
        }else{
            return alpha;
        }
    }
    return alpha;
}

void bfgs(VectorXd x,double e){
    int iter=0;
    MatrixXd H=MatrixXd::Identity(x.size(),x.size());
    MatrixXd In=MatrixXd::Identity(x.size(),x.size());
    while(iter<1000000){
        cout<<f(x)<<endl;
        VectorXd p=-1*H*gradient(x);
        double alpha=LineSearch(x,p);
        VectorXd x_pre=x;
        x+=alpha*p;
        if(gradient(x).norm()<e){
            break;
        }
        VectorXd s=x-x_pre;
        VectorXd y=gradient(x)-gradient(x_pre);
        double rho=s.dot(y);
        MatrixXd A=In-s*y.transpose()/rho;
        MatrixXd B=In-y*s.transpose()/rho;
        MatrixXd C=s*s.transpose()/rho;
        H=A*H*B+C;
        iter++;
    }
    cout<<iter<<endl;
    cout<<"Solution x:\n"<<x<<endl;
    cout<<"Function Value:\n"<<f(x)<<endl;
}

int main(){
    int size=4;
    double e=1e-5;
    VectorXd x(size);
    x(0)=3;x(1)=-1;x(2)=0;x(3)=1;
    bfgs(x,e);
}