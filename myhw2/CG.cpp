#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void CG(MatrixXd A, VectorXd b, VectorXd& x, double e){
    VectorXd r=A*x-b,r_pre;
    VectorXd p=-1.0*r;
    int iter=0;
    while(iter<=100000&&r.norm()>=e){
        double alpha=-1.0*(r.dot(p))/(p.dot(A*p));
        x=x+alpha*p;
        r_pre=r;
        r+=alpha*A*p;
        cout<<r.norm()<<endl;
        double beta=r.dot(r)/r_pre.dot(r_pre);
        p=-1*r+beta*p;
        iter++;
    }
    //cout<<"total iteration times:"<<iter<<endl;
    //cout<<r.norm()<<endl;
}

int main() {
    int size=20;
    double e=1e-6;
    MatrixXd A(size,size);
    VectorXd b(size);
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            A(i,j)=1.0/(i+j+1);
        }
        b(i)=1;
    }

    VectorXd x = VectorXd::Zero(size);

    CG(A,b,x,e);
    //std::cout << "Solution x:\n" << x << std::endl;
    return 0;
}
