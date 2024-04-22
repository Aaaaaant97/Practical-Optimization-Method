#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main() {
    MatrixXd A(3,3);
    A << 1, 2, 3,
        4, 5, 6,
        7, 8, 10;

    std::cout << "Matrix A:\n" << A << std::endl;

    MatrixXd A_inv = A.inverse();

    if(A_inv.determinant() != 0) { // 检查矩阵是否可逆
        std::cout << "Inverse of A:\n" << A_inv << std::endl;
    }
    else {
        std::cout << "Matrix A is not invertible." << std::endl;
    }

    return 0;
}
