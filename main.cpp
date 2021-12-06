#include <iostream>
#include "Cholesky.cpp"
#include "IO.cpp"

int main() {
    { //пример сравнения нормы матрицы результата метода Холецкого и вычисления СЛАУ с ней методом Гаусса
        const Lobaev::Math::Matrix<double> matrix_a({
            {1, 2, 4},
            {2, 13, 23},
            {4, 23, 77}
        });
        const Lobaev::Math::Matrix<double> matrix_l = Lobaev::Math::cholesky(matrix_a);
        const Lobaev::Math::Matrix<double> matrix_llt = matrix_l * matrix_l.transpose();

        std::cout << "Matrix A:" << std::endl;
        Lobaev::IO::print_matrix(std::cout, matrix_a);
        std::cout << std::endl;

        std::cout << "Matrix L:" << std::endl;
        Lobaev::IO::print_matrix(std::cout, matrix_l);
        std::cout << std::endl;

        std::cout << "Matrix (L * L^T) ~= A:" << std::endl;
        Lobaev::IO::print_matrix(std::cout, matrix_llt);
        std::cout << std::endl;

        std::cout << "Matrix A euclidean norm: " << matrix_a.norm_euclidean<double>() << std::endl;
        std::cout << "Matrix (L * L^T) euclidean norm: " << matrix_llt.norm_euclidean<double>() << std::endl;
    }
    return 0;
}
