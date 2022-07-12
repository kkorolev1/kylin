#ifndef KYLIN_MATRIX_HPP
#define KYLIN_MATRIX_HPP

#include <array>
#include <iostream>
#include <vector>
#include <type_traits>

namespace kylin {

    template <typename Container>
    class AbstractStorage {
    public:
        using T = typename Container::value_type;

        T& operator[](size_t i) {
            return data_[i];
        }

        const T& operator[](size_t i) const {
            return data_[i];
        }

        void swap(AbstractStorage& other) {
            data_.swap(other.data_);
        }
    protected:
        Container data_;
    };

    template <typename T, int Rows, int Cols>
    class StaticStorage : public AbstractStorage<std::array<T, Rows * Cols>> {
    public:
        using ValueType = T;

        constexpr static inline int kStaticStorageMaxSize = 4096;

        explicit StaticStorage(size_t, const T& value = T()) {
            this->data_.fill(value);
        }
    };

    template <typename T>
    class DynamicStorage : public AbstractStorage<std::vector<T>> {
    public:
        using ValueType = T;

        explicit DynamicStorage(size_t size, const T& value = T()) {
            this->data_.resize(size, value);
        }
    };

    template <typename T, int Rows, int Cols>
    using OptimizedStaticStorage = std::conditional_t<Rows * Cols <= StaticStorage<T, Rows, Cols>::kStaticStorageMaxSize,
            StaticStorage<T, Rows, Cols>, DynamicStorage<T>>;

    template <typename MatrixStorage>
    class OrderedMatrix : protected MatrixStorage {
    public:
        using T = typename MatrixStorage::ValueType;
        OrderedMatrix(int nrows, int ncols, const T& value = T()) : MatrixStorage(nrows * ncols, value), nrows_(nrows), ncols_(ncols) {}

        int nrows() const {
            return nrows_;
        }

        int ncols() const {
            return ncols_;
        }
    protected:
        int nrows_;
        int ncols_;
    };

    template <typename MatrixStorage>
    class RowMatrix : public OrderedMatrix<MatrixStorage> {
    public:
        using T = typename MatrixStorage::ValueType;
        using OrderedMatrix<MatrixStorage>::OrderedMatrix;

        RowMatrix(const std::initializer_list<std::initializer_list<T>>& list) :
            OrderedMatrix<MatrixStorage>(list.size(), list.begin()->size()) {
            if (list.size() != this->nrows_) {
                throw std::invalid_argument("Rows dimensions differs");
            }

            size_t index = 0;
            for (auto& row : list) {
                if (row.size() != this->ncols_) {
                    throw std::invalid_argument("Columns dimensions differs");
                }
                for (auto& value : row) {
                    this->data_[index++] = value;
                }
            }
        }

        T& operator()(size_t i, size_t j) {
            return at(i, j);
        }

        const T& operator()(size_t i, size_t j) const {
            return at(i, j);
        }

        T& at(size_t i, size_t j) {
            return this->data_[i * this->ncols_ + j];
        }

        const T& at(size_t i, size_t j) const {
            return this->data_[i * this->ncols_ + j];
        }
    };

    template <typename MatrixStorage>
    class ColumnMatrix : public OrderedMatrix<MatrixStorage> {
    public:
        using T = typename MatrixStorage::ValueType;
        using OrderedMatrix<MatrixStorage>::OrderedMatrix;

        ColumnMatrix(const std::initializer_list<std::initializer_list<T>>& list) :
            OrderedMatrix<MatrixStorage>(list.size(), list.begin()->size()) {
            if (list.size() != this->nrows_) {
                throw std::invalid_argument("Rows dimensions differs");
            }

            size_t row_index = 0;
            for (auto& row : list) {
                if (row.size() != this->ncols_) {
                    throw std::invalid_argument("Columns dimensions differs");
                }
                size_t col_index = 0;
                for (auto& value : row) {
                    (*this)(row_index, col_index) = value;
                    ++col_index;
                }
                ++row_index;
            }
        }

        T& operator()(size_t i, size_t j) {
            return at(i, j);
        }

        const T& operator()(size_t i, size_t j) const {
            return at(i, j);
        }

        T& at(size_t i, size_t j) {
            return this->data_[i + j * this->nrows_];
        }

        const T& at(size_t i, size_t j) const {
            return this->data_[i + j * this->nrows_];
        }
    };

    template <typename MatrixClass>
    class MatrixOp {
    public:

        template <typename MatrixClass2>
        auto& operator+=(const MatrixClass2& other) {
            auto& self = static_cast<MatrixClass&>(*this);
            for (int i = 0; i < self.nrows(); ++i) {
                for (int j = 0; j < self.ncols(); ++j) {
                    self(i, j) += other(i, j);
                }
            }
            return *this;
        }

        template <typename MatrixClass2>
        auto& operator-=(const MatrixClass2& other) {
            auto& self = static_cast<MatrixClass&>(*this);
            for (int i = 0; i < self.nrows(); ++i) {
                for (int j = 0; j < self.ncols(); ++j) {
                    self(i, j) -= other(i, j);
                }
            }
            return *this;
        }

        template <typename MatrixClass2>
        auto operator+(const MatrixClass2& other) {
            auto& self = static_cast<MatrixClass&>(*this);
            MatrixClass res(self);
            res += other;
            return res;
        }

        template <typename MatrixClass2>
        auto operator-(const MatrixClass2& other) {
            auto& self = static_cast<MatrixClass&>(*this);
            MatrixClass res(self);
            res -= other;
            return res;
        }

        template <typename MatrixClass2>
        auto operator*(const MatrixClass2& other) {
            auto& self = static_cast<MatrixClass&>(*this);
            MatrixClass res(self.nrows(), other.ncols());

            for (int i = 0; i < self.nrows(); ++i) {
                for (int j = 0; j < self.ncols(); ++j) {
                    for (int k = 0; k < self.ncols(); ++k) {
                        res(i, j) += self(i, k) * other(k, j);
                    }
                }
            }

            return res;
        }

        template <typename MatrixClass2>
        auto& operator*=(const MatrixClass& other) {
            auto& self = static_cast<MatrixClass&>(*this);
            auto res = self * other;
            self = std::move(res);
            return self;
        }

        auto transpose() {
            auto& self = static_cast<MatrixClass&>(*this);
            MatrixClass res(self.ncols(), self.nrows());

            for (int i = 0; i < self.nrows(); ++i) {
                for (int j = 0; j < self.ncols(); ++j) {
                    res(j, i) = self(i, j);
                }
            }

            return res;
        }

        friend std::ostream& operator<<(std::ostream& os, const MatrixClass& m) {
            for (int i = 0; i < m.nrows(); ++i) {
                for (int j = 0; j < m.ncols(); ++j) {
                    os << m(i, j) << " ";
                }
                os << "\n";
            }
            return os;
        }
    };

    template <typename T, int Rows, int Cols, template <class> class MatrixOrder = RowMatrix>
    class StaticMatrix : public MatrixOp<StaticMatrix<T, Rows, Cols, MatrixOrder>>, public MatrixOrder<OptimizedStaticStorage<T, Rows, Cols>> {
    public:
        using Base = MatrixOrder<OptimizedStaticStorage<T, Rows, Cols>>;
        explicit StaticMatrix(const T& value = T()) : Base(Rows, Cols, value) {}
        StaticMatrix(const std::initializer_list<std::initializer_list<T>>& elems) : Base(elems) {}
    };

    template <typename T, template <class> class MatrixOrder = RowMatrix>
    class DynamicMatrix : public MatrixOp<DynamicMatrix<T, MatrixOrder>>, public MatrixOrder<DynamicStorage<T>> {
    public:
        using Base = MatrixOrder<DynamicStorage<T>>;
        DynamicMatrix(int nrows, int ncols, const T& value = T()) : Base(nrows, ncols, value) {}
        DynamicMatrix(const std::initializer_list<std::initializer_list<T>>& elems) : Base(elems) {}
    };

    template <typename T>
    auto diag(const std::vector<T>& elems) {
        DynamicMatrix<T, RowMatrix> res(elems.size(), elems.size());

        for (size_t i = 0; i < elems.size(); ++i) {
            res(i, i) = elems[i];
        }

        return res;
    }

    template <typename T>
    auto eye(size_t size, const T& value) {
        DynamicMatrix<T, RowMatrix> res(size, size);

        for (size_t i = 0; i < size; ++i) {
            res(i, i) = value;
        }

        return res;
    }

    template <typename T>
    class MatrixView {

    };
}

#endif