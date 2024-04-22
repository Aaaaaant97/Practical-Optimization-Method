#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

double f(VectorXd x){
    double ans=0,a,b;
    for(int i=0;i<x.size()-1;i++){
        a=x(i+1)-x(i)*x(i);
        b=1-x(i);
        ans+=100*a*a+b*b;
    }
    return ans;
}

VectorXd gradient(VectorXd x){
    VectorXd ans(x.size());
    for(int i=0;i<x.size();i++){
        if(i==0){
            ans(i)=-400*(x(i+1)-x(i)*x(i))*x(i)-2+2*x(i);
        }else if(i<x.size()-1){
            ans(i)=200*(x(i)-x(i-1)*x(i-1))-400*(x(i+1)-x(i)*x(i))*x(i)-2+2*x(i);
        }else{
            ans(i)=200*(x(i)-x(i-1)*x(i-1));
        }
    }
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
    int size=6;
    double e=1e-5;
    VectorXd x(size);
    for(int i=0;i<size;i++){
        if(i%2==0){
            x(i)=-1.2;
        }else{
            x(i)=1;
        }
    }
    bfgs(x,e);
}