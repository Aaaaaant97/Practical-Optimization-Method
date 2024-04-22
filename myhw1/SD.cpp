#include <iostream>
#include <cmath>
#include <iomanip> 
using namespace std;

double f(double x1,double x2){
        double a=x2-x1*x1,b=1-x1,ans;
        ans=100*a*a+b*b;
        return ans;
}

double g1(double x1,double x2){
        double a=x2-x1*x1,b=1-x1;
        double ans=-400*a*x1-2*b;
        return ans;
}

double g2(double x1,double x2){
        double ans=200*(x2-x1*x1);
        return ans;
}

double norm(double a,double b){
        double ans=sqrt(a*a+b*b);
        return ans;
}

double LineSearch(double x1,double x2,double p1,double p2){
        double alpha=1,c=1e-4;
        double y1=x1-alpha*p1,y2=x2-alpha*p2;
        while(f(y1,y2)>f(x1,x2)-c*alpha*(p1*p1+p2*p2)){
                alpha*=0.5;
                y1 = x1 - alpha * p1;
                y2 = x2 - alpha * p2;
        }
        return alpha;
}

void SteepestDescent(double x1,double x2,double e){
        int iter=0;
        while(iter<=100000){
                double p1,p2;
                if(iter==0){
                        p1=g1(x1,x2);
                        p2=g2(x1,x2);
                }
                double alpha=LineSearch(x1,x2,p1,p2);
                x1=x1-alpha*p1;
                x2=x2-alpha*p2;
                p1=g1(x1,x2);
                p2=g2(x1,x2);
                if(norm(p1,p2)<e){
                        break;
                }
                iter++;
                if(iter%100==0){
                        cout<<f(x1,x2)<<endl;
                }
        }
        //cout<<"The number of iterations is "<<iter<<endl;
        //cout <<"Minimum found at x1 = "<<fixed<<setprecision(10)<<x1 << ", x2 = " << x2  << endl;
        //cout<<"The value of the function is "<<f(x1,x2)<<endl;
}

int main(){
        double e=1e-6;
        double x1=-1.2,x2=1;
        SteepestDescent(x1,x2,e);
        return 0;
}