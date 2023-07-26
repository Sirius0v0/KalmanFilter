#pragma once

#include <cmath>
#include <algorithm>
#include <functional>

#include <vector>
#include <string>
#include <stack>
#include <tuple>
#include <initializer_list>

#include <sstream>
#include <istream>
#include <ostream>
#include <iomanip>

#include <KFExceptions/KFExceptions.hpp>
#include <matrix/IOFormat.hpp>

namespace kalmans
{

    template<typename T>
    using matrixdata = std::vector<T>;

    template<typename T>
    class Matrix
    {
    public:
        Matrix();
        Matrix(int row, int col);
        Matrix(int row, int col, T e);
        Matrix(int row, int col, const std::initializer_list<T>& init_lst);
        Matrix(const Matrix& mat);
        Matrix(Matrix&& mat) noexcept;
        virtual ~Matrix() = default;

        void write(std::ostream& os) const;
        void read(std::istream& is);

        Matrix& operator=(const T& e);
        Matrix& operator=(const Matrix& mat);
        Matrix operator-() const;
        T& operator()(int i, int j);
        T& operator()(int i);
        T operator()(int i, int j) const;
        T operator()(int i) const;

        int get_row() const { return row_; }
        int get_col() const { return col_; }
        int get_size() const { return row_ * col_; }

        void swap(Matrix& mat) noexcept;
        int assign(std::initializer_list<T> data_lst);
        int assign(int num, T d);
        void resize(int new_row, int new_col);
        T& at(int i, int j);
        T& at(int i);
        T at(int i, int j) const;
        T at(int i) const;

        void row_swap(int i, int j);
        void row_swap(int i, int j, int start, int end);
        void col_swap(int i, int j);
        void col_swap(int i, int j, int start, int end);

        Matrix transpose() const;
        Matrix operator+(const Matrix& mat) const;
        Matrix operator-(const Matrix& mat) const;

        Matrix operator*(const Matrix& mat) const;

        std::tuple<Matrix<double>, Matrix<double>, Matrix<double>> lup_decomposition() const;
        Matrix<double> inverse() const;

        static Matrix Eye(int dim);

    private:

        int row_;
        int col_;
        matrixdata<T> data_;

        void init(int row = 0, int col = 0);
    };

    template<typename T>
    std::ostream& operator<< (std::ostream& os, const Matrix<T>& mat)
    {
        mat.write(os);
        return os;
    }

    template<typename T>
    std::istream& operator>> (std::istream& is, const Matrix<T>& mat)
    {
        return is;
    }

    template<typename T>
    Matrix<T>& operator<< (Matrix<T>& mat, std::initializer_list<T>&& init_lst)
    {
        mat.assign(init_lst);
        return mat;
    }

    template<typename T>
    void Matrix<T>::init(int row, int col)
    {
        this->row_ = row;
        this->col_ = col;
        this->data_.resize(row * col, 0);
    }

    template <typename T>
    Matrix<T>::Matrix()
    {
        init();
    }

    template <typename T>
    Matrix<T>::Matrix(int row, int col)
    {
        init(row, col);
    }

    template <typename T>
    Matrix<T>::Matrix(int row, int col, T e) :Matrix(row, col)
    {
        for (int i = 0; i < this->get_size(); ++i)
        {
            this->data_.at(i) = e;
        }
    }

    template <typename T>
    Matrix<T>::Matrix(int row, int col, const std::initializer_list<T>& init_lst)
        :Matrix(row, col)
    {
        this->assign(init_lst);
    }

    template <typename T>
    Matrix<T>::Matrix(const Matrix& mat)
    {
        init(mat.get_row(), mat.get_col());
        for (int i = 0; i < mat.get_size(); ++i)
        {
            this->data_.at(i) = mat(i);
        }
    }

    template <typename T>
    Matrix<T>::Matrix(Matrix&& mat) noexcept :Matrix()
    {
        this->swap(mat);
    }

    template <typename T>
    void Matrix<T>::write(std::ostream& os) const
    {
        if (this->get_row() == 0 || this->get_col() == 0)
            return;

        os.setf(std::ios::fixed, std::ios::floatfield);
        os.precision(IOFormat::precision);

        os << IOFormat::start_delimiter;
        size_t num = 0;
        for (const auto& elem : this->data_)
        {
            os << IOFormat::elem_delimiter << elem;
            num++;
            if (num == this->get_size())
                os << IOFormat::elem_delimiter << IOFormat::end_delimiter;
            if (num % this->get_col() == 0)
                os << IOFormat::row_delimiter;
        }
    }

    template <typename T>
    void Matrix<T>::read(std::istream& is)
    {
        // 更新当前读取格式
        IOFormat::update_cache_flag();

        std::string buf;
        int row = 0;
        int col = 0;
        matrixdata<T> data;

        if (IOFormat::cache_start_delim)
            is >> buf;

        if (IOFormat::cache_row_delim)
        {
            if (IOFormat::cache_elem_delim)
            {

            }
            else
            {

            }
        }
        else
        {
            if (IOFormat::cache_elem_delim)
            {

            }
            else
            {

            }
        }

        if (IOFormat::cache_end_delim)
            is >> buf;

        this->data_.swap(data);
        this->row_ = row;
        this->col_ = col;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator=(const T& e)
    {
        this->assign(this->get_size(), e);
        return *this;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator=(const Matrix& mat)
    {
        this->swap(Matrix(mat));
        return *this;
    }

    template <typename T>
    Matrix<T> Matrix<T>::operator-() const
    {
        Matrix tmp(*this);
        std::for_each(tmp.data_.begin(), tmp.data_.end(), [](T& x) {x *= -1; });
        return tmp;
    }

    template <typename T>
    T& Matrix<T>::operator()(int i, int j)
    {
        if (i < 0 || j < 0 || i > this->row_ - 1 || j > this->col_ - 1)
        {
            std::ostringstream oss;
            oss << "尝试访问 (" << i << ", " << j << ") 位置的值，但该索引超出范围："
                << "允许的最大索引为 (" << this->get_row() - 1 << ", " << this->get_col() - 1 << ")\n";
            throw kalmans::OutOfRangeError(oss.str());
        }
        return this->data_.at(this->col_ * i + j);
    }

    template <typename T>
    T& Matrix<T>::operator()(int i)
    {
        if (i < 0 || i > this->get_size() - 1)
        {
            std::ostringstream oss;
            oss << "尝试访问 (" << i << ") 位置的值，但该索引超出范围："
                << "允许的最大索引为 (" << this->get_size() - 1 << ")\n";
            throw kalmans::OutOfRangeError(oss.str());
        }
        return this->data_.at(i);
    }

    template <typename T>
    T Matrix<T>::operator()(int i, int j) const
    {
        if (i < 0 || j < 0 || i > this->row_ - 1 || j > this->col_ - 1)
        {
            std::ostringstream oss;
            oss << "尝试访问 (" << i << ", " << j << ") 位置的值，但该索引超出范围："
                << "允许的最大索引为 (" << this->get_row() - 1 << ", " << this->get_col() - 1 << ")\n";
            throw kalmans::OutOfRangeError(oss.str());
        }
        return this->data_.at(this->col_ * i + j);
    }

    template <typename T>
    T Matrix<T>::operator()(int i) const
    {
        if (i < 0 || i > this->get_size() - 1)
        {
            std::ostringstream oss;
            oss << "尝试访问 (" << i << ") 位置的值，但该索引超出范围："
                << "允许的最大索引为 (" << this->get_size() - 1 << ")\n";
            throw kalmans::OutOfRangeError(oss.str());
        }
        return this->data_.at(i);
    }

    template <typename T>
    void Matrix<T>::swap(Matrix& mat) noexcept
    {
        std::swap(this->row_, mat.row_);
        std::swap(this->col_, mat.col_);
        this->data_.swap(mat.data_);
    }

    template <typename T>
    int Matrix<T>::assign(std::initializer_list<T> data_lst)
    {
        if (data_lst.size() > this->get_size())
        {
            std::ostringstream oss;
            oss << "尝试给大小为" << this->get_size() << "的矩阵赋值个数为"
                << data_lst.size() << "的值，数据大小超过容量\n";
            throw kalmans::LengthError(oss.str());
        }

        int idx = 0;
        for (auto&& elem : data_lst)
        {
            this->data_.at(idx) = std::move(elem);
            idx++;
        }
        return idx;
    }

    template <typename T>
    int Matrix<T>::assign(int num, T d)
    {
        if (num > this->get_size())
        {
            std::ostringstream oss;
            oss << "尝试给大小为" << this->get_size() << "的矩阵赋值个数为"
                << num << "的值，数据大小超过容量\n";
            throw kalmans::LengthError(oss.str());
        }

        int idx = 0;
        for (int i = 0; i < num; ++i)
        {
            this->data_.at(idx) = d;
            idx++;
        }
        return idx;
    }

    template <typename T>
    void Matrix<T>::resize(int new_row, int new_col)
    {
        this->row_ = new_row;
        this->col_ = new_col;
        this->data_.resize(new_row * new_col, 0);
        this->data_.shrink_to_fit();
    }


    template <typename T>
    T& Matrix<T>::at(int i, int j)
    {
        return (*this)(i, j);
    }

    template <typename T>
    T& Matrix<T>::at(int i)
    {
        return (*this)(i);
    }

    template <typename T>
    T Matrix<T>::at(int i, int j) const
    {
        return (*this)(i, j);
    }

    template <typename T>
    T Matrix<T>::at(int i) const
    {
        return (*this)(i);
    }

    template <typename T>
    void Matrix<T>::row_swap(int i, int j, int start, int end)
    {
        for (int k = start; k < end; ++k)
        {
            T tmp = (*this)(i, k);
            (*this)(i, k) = (*this)(j, k);
            (*this)(j, k) = tmp;
        }
    }

    template <typename T>
    void Matrix<T>::row_swap(int i, int j)
    {
        this->row_swap(i, j, 0, this->get_col());
    }

    template <typename T>
    void Matrix<T>::col_swap(int i, int j, int start, int end)
    {
        for (int k = start; k < end; ++k)
        {
            T tmp = (*this)(k, i);
            (*this)(k, i) = (*this)(k, j);
            (*this)(k, j) = tmp;
        }
    }

    template <typename T>
    void Matrix<T>::col_swap(int i, int j)
    {
        this->col_swap(i, j, 0, this->get_row());
    }

    template <typename T>
    Matrix<T> Matrix<T>::transpose() const
    {
        Matrix<T> tmp(this->get_col(), this->get_row());
        for (int i = 0; i < tmp.get_row(); ++i)
            for (int j = 0; j < tmp.get_col(); ++j)
            {
                tmp(i, j) = (*this)(j, i);
            }
        return tmp;
    }

    template <typename T>
    Matrix<T> Matrix<T>::operator+(const Matrix& mat) const
    {
        Matrix tmp(*this);
        if (tmp.get_col() != mat.get_col() || tmp.get_row() != mat.get_row())
        {
            std::ostringstream oss;
            oss << "矩阵形状 (" << tmp.get_row() << ", " << tmp.get_col() << ") 与 ("
                << mat.get_row() << ", " << mat.get_col() << ") 不匹配\n";
            throw kalmans::ValueError(oss.str());
        }
        for (int i = 0; i < mat.get_size(); ++i)
        {
            tmp(i) += mat(i);
        }
        return tmp;
    }

    template <typename T>
    Matrix<T> Matrix<T>::operator-(const Matrix& mat) const
    {
        Matrix tmp(*this);
        return tmp + (-mat);
    }

    template <typename T>
    Matrix<T> Matrix<T>::operator*(const Matrix<T>& mat) const
    {
        if (this->get_col() != mat.get_row())
        {
            std::ostringstream oss;
            oss << "矩阵形状 (" << this->get_row() << ", " << this->get_col() << ") 与 ("
                << mat.get_row() << ", " << mat.get_col() << ") 不匹配："
                << this->get_col() << "!=" << mat.get_row() << '\n';
            throw kalmans::ValueError(oss.str());
        }
        Matrix<T> tmp(this->get_row(), mat.get_col());
        for (int i = 0; i < tmp.get_row(); ++i)
            for (int j = 0; j < tmp.get_col(); ++j)
                for (int k = 0; k < mat.get_row(); ++k)
                {
                    tmp(i, j) += (*this)(i, k) * mat(k, j);
                }
        return tmp;
    }

    template <typename T>
    std::tuple<Matrix<double>, Matrix<double>, Matrix<double>> Matrix<T>::lup_decomposition() const
    {
        if (this->get_row() != this->get_col())
        {
            std::ostringstream oss;
            oss << "矩阵形状 (" << this->get_row() << ", " << this->get_col() << ") 非方阵，LUP分解需要矩阵为方阵\n";
            throw kalmans::ValueError(oss.str());
        }
        Matrix<double> U(*this);
        Matrix<double> L = Matrix<double>::Eye(this->get_row());
        Matrix<double> P = Matrix<double>::Eye(this->get_row());

        for (int j = 0; j < this->get_col() - 1; ++j)
        {
            // 选择列主元
            int principal = j;
            T max = fabs(U(principal, j));
            // 只遍历当前主元下方的元素
            for (int i = j + 1; i < this->get_row(); ++i)
            {
                if (fabs(U(i, j)) > max)
                {
                    principal = i;
                    max = fabs(U(i, j));
                }
            }
            if (j != principal) // 选择主元交换
            {
                // 交换 U(j,j:n) <-> U(i,j:n)
                U.row_swap(j, principal, j, this->get_col());
                // 交换 L(j,1:j-1) <-> L(i,1:j-1)
                L.row_swap(j, principal, 0, j - 1);
                // 交换 P(j,1:n) <-> P(i,1:n)
                P.row_swap(j, principal);
            }

            // 列主元如果为 0，则跳过该列的消元
            if (fabs(U(j, j)) < 1e-7)
                continue;
            // 计算 L 和 U
            for (int i = j + 1; i < this->get_row(); ++i)
            {
                L(i, j) = U(i, j) / U(j, j) / 1.0;
                for (int k = j; k < this->get_row(); ++k)
                {
                    U(i, k) -= L(i, j) * U(j, k);
                }
            }
        }

        return std::make_tuple(L, U, P);
    }

    template <typename T>
    Matrix<double> Matrix<T>::inverse() const
    {
        auto [L_, U_, P_] = this->lup_decomposition();

        // L 求逆
        for (int i = 0; i < L_.get_row(); ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                if (i == j)
                {
                    L_(i, j) = 1.0 / L_(i, j);
                }
                else
                {
                    double sum = 0.0;
                    for (int k = j; k < i; ++k)
                    {
                        sum += L_(i, k) * L_(k, j);
                    }
                    L_(i, j) = -1.0 * sum / L_(i, i);
                }
            }
        }

        // U 求逆
        for (int j = U_.get_col() - 1; j >= 0; --j)
        {
            for (int i = j; i >= 0; --i)
            {
                if (i == j)
                {
                    U_(i, j) = 1.0 / U_(i, j);
                }
                else
                {
                    double sum = 0.0;
                    for (int k = j; k >= i + 1; --k)
                    {
                        sum += U_(i, k) * U_(k, j);
                    }
                    U_(i, j) = -1.0 * sum / U_(i, i);
                }
            }
        }

        return Matrix<double>(U_ * L_ * P_);
    }

    template <typename T>
    Matrix<T> Matrix<T>::Eye(int dim)
    {
        Matrix<T> eye(dim, dim);

        for (int i = 0; i < dim; ++i)
        {
            eye(i, i) = 1;
        }

        return eye;
    }

}
