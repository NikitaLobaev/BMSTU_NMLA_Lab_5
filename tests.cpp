#include "gtest/gtest.h"
#include "Cholesky.cpp"

using namespace Lobaev::Math;

TEST(cholesky, test_1_ok) {
    const Matrix<double> matrix_a({
        {1, 3},
        {3, 13}
    });
    EXPECT_EQ(matrix_a, matrix_a.transpose());

    const Matrix<double> matrix_l = cholesky(matrix_a);

    const Matrix<double> matrix_llt = matrix_l * matrix_l.transpose();

    ASSERT_EQ(matrix_a, matrix_llt);
}

TEST(cholesky, test_2_ok) {
    const Matrix<double> matrix_a({
        {1, 2, 4},
        {2, 13, 23},
        {4, 23, 77}
    });
    EXPECT_EQ(matrix_a, matrix_a.transpose());

    const Matrix<double> matrix_l = cholesky(matrix_a);

    const Matrix<double> matrix_llt = matrix_l * matrix_l.transpose();

    ASSERT_EQ(matrix_a, matrix_llt);
}
