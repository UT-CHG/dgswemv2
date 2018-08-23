#ifndef SERIALIZATION_EIGEN_MATRIX_HPP
#define SERIALIZATION_EIGEN_MATRIX_HPP

namespace hpx {
namespace serialization {

// Vector serialization
template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
std::enable_if_t<(Rows == 1) || (Cols == 1)>
    serialize(output_archive& ar, const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& vec, unsigned) {
    ar << vec.size();

    for (uint i = 0; i < vec.size(); ++i) {
        ar << vec(i);
    }
}

// Dense matrix serialization
template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
std::enable_if_t<!((Rows == 1) || (Cols == 1))>
    serialize(output_archive& ar, const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& mat, unsigned) {
    ar << mat.rows() << mat.cols();

    // Eigen defaults to column major
    for (uint j = 0; j < mat.cols(); ++j) {
        for (uint i = 0; i < mat.rows(); ++i) {
            ar << mat(i, j);
        }
    }
}

template <typename Scalar, int Options, typename StorageIndex>
void serialize(output_archive& ar, const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& mat, unsigned) {
    using Iterator = typename Eigen::SparseMatrix<Scalar, Options, StorageIndex>::InnerIterator;

    ar << mat.rows() << mat.cols() << mat.nonZeros();

    for (uint k = 0; k < mat.outerSize(); ++k) {
        for (Iterator it(mat, k); it; ++it) {
            ar << it.row() << it.col() << it.value();
        }
    }
}

// Vector serialization
template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
std::enable_if_t<(Rows == 1) || (Cols == 1)>
    serialize(input_archive& ar, Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& vec, unsigned) {
    size_t size;
    ar >> size;
    vec.resize(size);

    Scalar value{};

    for (uint i = 0; i < size; ++i) {
        ar >> value;
        vec(i) = value;
    }
}

// Dense matrix serialization
template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
std::enable_if_t<!((Rows == 1) || (Cols == 1))>
    serialize(input_archive& ar, Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& mat, unsigned) {
    size_t n_rows, n_cols;
    ar >> n_rows >> n_cols;
    mat.resize(n_rows, n_cols);

    Scalar value{};

    for (size_t j = 0; j < n_cols; ++j) {
        for (size_t i = 0; i < n_rows; ++i) {
            ar >> value;
            mat(i, j) = value;
        }
    }
}

template <typename Scalar, int Options, typename StorageIndex>
void serialize(input_archive& ar, Eigen::SparseMatrix<Scalar, Options, StorageIndex>& mat, unsigned) {
    size_t rows, cols, n_non_zeros;
    ar >> rows >> cols >> n_non_zeros;
    mat.resize(rows, cols);
    mat.reserve(n_non_zeros);

    size_t i, j;
    Scalar value{};

    for (size_t id = 0UL; id < n_non_zeros; ++id) {
        ar >> i >> j >> value;
        mat.insert(i, j) = value;
    }
}
}
}
#endif