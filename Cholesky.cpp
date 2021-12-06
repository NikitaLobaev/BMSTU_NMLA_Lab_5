#include "Matrix.cpp"

namespace Lobaev::Math {

    template <class T>
    Matrix<T> cholesky(const Matrix<T> matrix_a) {
        if (!matrix_a.is_square()) {
            throw "cholesky: only square matrix A is supported.";
        }

        Matrix<T> matrix_l(matrix_a.rows_count(), matrix_a.columns_count());

        matrix_l(0, 0) = (T) std::pow(matrix_a(0, 0), (T) std::pow((T) 2, -1));

        for (size_t j = 1; j < matrix_l.rows_count(); j++) {
            matrix_l(j, 0) = matrix_a(j, 0) / matrix_l(0, 0);
        }

        for (size_t i = 1; i < matrix_l.rows_count(); i++) {
            T sum = (T) 0;
            for (size_t p = 0; p < i; p++) {
                sum += matrix_l(i, p) * matrix_l(i, p);
            }
            matrix_l(i, i) = (T) std::pow(matrix_a(i, i) - sum, (T) std::pow((T) 2, -1));

            for (size_t j = i + 1; j < matrix_l.rows_count(); j++) {
                sum = (T) 0;
                for (size_t p = 0; p < i; p++) {
                    sum += matrix_l(i, p) * matrix_l(j, p);
                }
                matrix_l(j, i) = (matrix_a(j, i) - sum) / matrix_l(i, i);
            }
        }

        return matrix_l;
    }

}