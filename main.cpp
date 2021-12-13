#include <iostream>
#include "Cholesky.cpp"
#include "Gauss.cpp"
#include "IO.cpp"

int main() {
    { //пример сравнения нормы матрицы результата метода Холецкого
        std::cout << "--- Task #1 ---" << std::endl;
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
        std::cout << std::endl;
    }
    { //пример сравнения норм векторов результата метода Холецкого и вычисления СЛАУ с ним методом Гаусса
        std::cout << "--- Task #2 ---" << std::endl;
        const Lobaev::Math::Matrix<double> matrix_a({
            {1, 2, 4},
            {2, 13, 23},
            {4, 23, 77}
        });
        const Lobaev::Math::Vector<double> vector_x({
            10,
            11,
            2
        });
        const Lobaev::Math::Vector<double> vector_b = matrix_a * vector_x;
        const Lobaev::Math::Matrix<double> matrix_l = Lobaev::Math::cholesky(matrix_a);
        const Lobaev::Math::Matrix<double> matrix_l_transposed = matrix_l.transpose();
        const Lobaev::Math::Vector<double> vector_y = matrix_l_transposed * vector_x;

        const Lobaev::Math::Vector<double> gauss_result_1 = Lobaev::Math::Gauss::gauss(matrix_l, vector_b);
        const Lobaev::Math::Vector<double> gauss_result_2 = Lobaev::Math::Gauss::gauss(matrix_l_transposed, vector_y);

        std::cout << "Vector Y euclidean norm: " << vector_y.norm_euclidean<double>() << std::endl;
        std::cout << "Gauss vector Y euclidean norm: " << gauss_result_1.norm_euclidean<double>() << std::endl;
        std::cout << std::endl;
        std::cout << "Vector X euclidean norm: " << vector_x.norm_euclidean<double>() << std::endl;
        std::cout << "Gauss vector X euclidean norm: " << gauss_result_2.norm_euclidean<double>() << std::endl;
        std::cout << std::endl;
    }
    return 0;
}
