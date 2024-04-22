#include <iostream>
#include <cmath>
#include <iomanip>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#define miu 0.001
using namespace std;
using namespace Eigen;

double AbsoluteValue(double a){
    if(a>=0){
        return a;
    }else{
        return -1*a;
    }
}

double L(VectorXd x){
    double delta=0.001*miu;
    double ret=0;
    for(int i=0;i<1024;i++){
        double value=x(i);
        if(AbsoluteValue(value)<delta){
            ret+=1.0/(2*delta)*value*value;
        }else{
            ret+=AbsoluteValue(value)-delta/2.0;
        }
    }
    return ret;
}

double f(MatrixXd A,VectorXd x,VectorXd b){
        VectorXd c(b.size());
        c=A*x-b;
        double M=c.norm();
        double ans=M*M/2.0+miu*L(x);
        return ans;
}

VectorXd g(MatrixXd A,VectorXd x,VectorXd b){
    VectorXd ans(x.size());
    for(int i=0;i<1024;i++){
        ans(i)=0;
    }
    VectorXd Ax=A*x-b;
    for(int i=0;i<1024;i++){
        for(int j=0;j<512;j++){
            ans(i)+=A(j,i)*Ax(j);
        }
    }
    return ans;
}

double Linesearch(MatrixXd A,VectorXd x,VectorXd b,VectorXd p){
        double alpha=1,c=1e-4;
        VectorXd y=x-alpha*p;
        while(f(A,y,b)>f(A,x,b)-c*alpha*p.norm()*p.norm()){
            alpha*=0.618;
            y=x-alpha*p;
        }
        return alpha;
}

void SteepestDescent(MatrixXd A,VectorXd x,VectorXd b,double e){
    int iter=0;
    while(iter<=500){
        VectorXd p(1024);
        if(iter==0){
            p=g(A,x,b);
        }
        double alpha=Linesearch(A,x,b,p);
        x=x-alpha*p;
        p=g(A,x,b);
        if(p.norm()<e){
            break;
        }
        iter++;
        if(iter%10==0){
            //cout<<"The number of iterations is "<<iter<<endl;
            //cout<<"The value of the function is "<<f(A,x,b)<<endl;
            cout<<f(A,x,b)<<endl;
        }
    }
    cout<<f(A,x,b)<<endl;
    cout<<iter<<endl;
    //cout<<"The number of iterations is "<<iter<<endl;
    //cout<<"The value of the function is "<<f(A,x,b)<<endl;
}

int main(){
    double e=1e-6;
    int m=512,n=1024;
    MatrixXd A=MatrixXd::Random(m,n);

    SparseVector<double> y(n);
    for (int i = 0; i < n; i++) {
        if (rand() / (double)RAND_MAX < 0.1) {
            y.insert(i) = rand() / (double)RAND_MAX; 
        }
    }
    VectorXd x(n);
    x=y.toDense();

    VectorXd b=VectorXd::Random(m);
    SteepestDescent(A,x,b,e);
    return 0;
}