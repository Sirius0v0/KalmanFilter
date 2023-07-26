#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <stack>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <initializer_list>
#include <functional>
#include <algorithm>

#include <KFExceptions/KFExceptions.hpp>

namespace kalmans
{

    template<typename T>
    using matrixdata = std::vector<T>;

    template<typename T, int Row = 0, int Col = 0>
    class Matrix
    {
        friend std::ostream& operator<< (std::ostream& out, const Matrix<T, Row, Col>& mat)
        {
            std::cout << "[";
            int num = 0;
            for (const auto& elem : mat.data_)
            {
                std::cout << std::fixed << std::setprecision(4) << "\t" << (fabs(elem) < 1e-8 ? 0 : elem);
                num++;
                if (num == mat.get_size())
                    std::cout << "\t]";
                if (num % mat.get_col() == 0)
                    std::cout << '\n';
            }
            return out;
        }
    public:
        Matrix();
        Matrix(const std::initializer_list<T>& init_lst);
        Matrix(const Matrix& mat);
        Matrix(Matrix&& mat) noexcept;
        Matrix(T e);
        virtual ~Matrix() = default;

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
        std::string get_row_name() const { return row_name_; }
        std::string get_col_name() const { return col_name_; }
        std::string get_data_name() const { return data_name_; }

        void swap(Matrix& mat) noexcept;
        int assign(std::initializer_list<T> data_lst);
        int assign(int num, T d);
        void resize(int new_row, int new_col);
        T& at(int i, int j);
        T& at(int i);
        T at(int i, int j) const;
        T at(int i) const;

        Matrix<T, Col, Row> transpose() const;
        Matrix operator+(const Matrix& mat) const;
        Matrix operator-(const Matrix& mat) const;

        template<int NewCol>
        Matrix<T, Row, NewCol> operator*(const Matrix<T, Col, NewCol>& mat) const;

        std::tuple<Matrix<double, Row, Col>, Matrix<double, Row, Col>, Matrix<double, Row, Col>> lup_decomposition() const;
        Matrix<double, Row, Col> inverse() const;

        static Matrix Eye();

    private:

        int row_;
        std::string row_name_;
        int col_;
        std::string col_name_;
        matrixdata<T> data_;
        std::string data_name_;

        void init(const std::string& data_str, const std::string& row_str,
            const std::string& col_str, int row = Row, int col = Col);

        void row_swap(int i, int j);
        void row_swap(int i, int j, int start, int end);
        void col_swap(int i, int j);
        void col_swap(int i, int j, int start, int end);
    };

    template<typename T, int Row, int Col>
    Matrix<T, Row, Col>& operator<< (Matrix<T, Row, Col>& mat, std::initializer_list<T>&& init_lst)
    {
        mat.assign(init_lst);
        return mat;
    }

    template <typename T, int Row, int Col>
    void Matrix<T, Row, Col>::swap(Matrix& mat) noexcept
    {
        std::swap(this->row_, mat.row_);
        std::swap(this->col_, mat.col_);
        this->row_name_.swap(mat.row_name_);
        this->col_name_.swap(mat.col_name_);
        this->data_.swap(mat.data_);
        this->data_name_.swap(mat.data_name_);
    }

    template <typename T, int Row, int Col>
    int Matrix<T, Row, Col>::assign(std::initializer_list<T> data_lst)
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

    template <typename T, int Row, int Col>
    int Matrix<T, Row, Col>::assign(int num, T d)
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

    template <typename T, int Row, int Col>
    void Matrix<T, Row, Col>::resize(int new_row, int new_col)
    {
        this->row_ = new_row;
        this->col_ = new_col;
        this->data_.resize(new_row * new_col, 0);
    }

    template<typename T, int Row, int Col>
    void Matrix<T, Row, Col>::init(const std::string& data_str, const std::string& row_str,
        const std::string& col_str, int row, int col)
    {
        this->row_ = row;
        this->row_name_ = row_str;
        this->col_ = col;
        this->col_name_ = col_str;
        this->data_.resize(row * col, 0);
        this->data_name_ = data_str;
    }

    template <typename T, int Row, int Col>
    void Matrix<T, Row, Col>::row_swap(int i, int j, int start, int end)
    {
        for (int k = start; k < end; ++k)
        {
            T tmp = (*this)(i, k);
            (*this)(i, k) = (*this)(j, k);
            (*this)(j, k) = tmp;
        }
    }

    template <typename T, int Row, int Col>
    void Matrix<T, Row, Col>::row_swap(int i, int j)
    {
        this->row_swap(i, j, 0, this->get_col());
    }

    template <typename T, int Row, int Col>
    void Matrix<T, Row, Col>::col_swap(int i, int j, int start, int end)
    {
        for (int k = start; k < end; ++k)
        {
            T tmp = (*this)(k, i);
            (*this)(k, i) = (*this)(k, j);
            (*this)(k, j) = tmp;
        }
    }

    template <typename T, int Row, int Col>
    void Matrix<T, Row, Col>::col_swap(int i, int j)
    {
        this->col_swap(i, j, 0, this->get_row());
    }

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col>::Matrix()
    {
        init("__data__", "__row__", "__col__");
    }

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col>::Matrix(const std::initializer_list<T>& init_lst) :Matrix()
    {
        this->assign(init_lst);
    }

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col>::Matrix(const Matrix& mat)
    {
        init(mat.get_data_name(),
            mat.get_row_name(),
            mat.get_col_name(),
            mat.get_row(),
            mat.get_col());
        for (int i = 0; i < mat.get_size(); ++i)
        {
            this->data_.at(i) = mat(i);
        }
    }

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col>::Matrix(Matrix&& mat) noexcept :Matrix()
    {
        this->swap(mat);
    }

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col>::Matrix(T e) :Matrix()
    {
        for (int i = 0; i < this->get_size(); ++i)
        {
            this->data_.at(i) = e;
        }
    }

    template <typename T, int Row, int Col>
    T& Matrix<T, Row, Col>::operator()(int i, int j)
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

    template <typename T, int Row, int Col>
    T& Matrix<T, Row, Col>::operator()(int i)
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

    template <typename T, int Row, int Col>
    T Matrix<T, Row, Col>::operator()(int i, int j) const
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

    template <typename T, int Row, int Col>
    T Matrix<T, Row, Col>::operator()(int i) const
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

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col>& Matrix<T, Row, Col>::operator=(const T& e) : Matrix()
    {
        this->assign(this->get_size(), e);
        return *this;
    }

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col>& Matrix<T, Row, Col>::operator=(const Matrix& mat)
    {
        this->swap(Matrix(mat));
        return *this;
    }

    template <typename T, int Row, int Col>
    T& Matrix<T, Row, Col>::at(int i, int j)
    {
        return (*this)(i, j);
    }

    template <typename T, int Row, int Col>
    T& Matrix<T, Row, Col>::at(int i)
    {
        return (*this)(i);
    }

    template <typename T, int Row, int Col>
    T Matrix<T, Row, Col>::at(int i, int j) const
    {
        return (*this)(i, j);
    }

    template <typename T, int Row, int Col>
    T Matrix<T, Row, Col>::at(int i) const
    {
        return (*this)(i);
    }

    template <typename T, int Row, int Col>
    Matrix<T, Col, Row> Matrix<T, Row, Col>::transpose() const
    {
        Matrix<T, Col, Row> tmp;
        for (int i = 0; i < tmp.get_row(); ++i)
            for (int j = 0; j < tmp.get_col(); ++j)
            {
                tmp(i, j) = (*this)(j, i);
            }
        return tmp;
    }

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col> Matrix<T, Row, Col>::operator+(const Matrix& mat) const
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

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col> Matrix<T, Row, Col>::operator-() const
    {
        Matrix tmp(*this);
        std::for_each(tmp.data_.begin(), tmp.data_.end(), [](T& x) {x *= -1; });
        return tmp;
    }

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col> Matrix<T, Row, Col>::operator-(const Matrix& mat) const
    {
        Matrix tmp(*this);
        return tmp + (-mat);
    }

    template <typename T, int Row, int Col>
    template <int NewCol>
    Matrix<T, Row, NewCol> Matrix<T, Row, Col>::operator*(const Matrix<T, Col, NewCol>& mat) const
    {
        if (this->get_col() != mat.get_row())
        {
            std::ostringstream oss;
            oss << "矩阵形状 (" << this->get_row() << ", " << this->get_col() << ") 与 ("
                << mat.get_row() << ", " << mat.get_col() << ") 不匹配："
                << this->get_col() << "!=" << mat.get_row() << '\n';
            throw kalmans::ValueError(oss.str());
        }
        Matrix<T, Row, NewCol> tmp;
        for (int i = 0; i < tmp.get_row(); ++i)
            for (int j = 0; j < tmp.get_col(); ++j)
                for (int k = 0; k < mat.get_row(); ++k)
                {
                    tmp(i, j) += (*this)(i, k) * mat(k, j);
                }
        return tmp;
    }

    template <typename T, int Row, int Col>
    std::tuple<Matrix<double, Row, Col>, Matrix<double, Row, Col>, Matrix<double, Row, Col>> Matrix<T, Row, Col>::lup_decomposition() const
    {
        if (this->get_row() != this->get_col())
        {
            std::ostringstream oss;
            oss << "矩阵形状 (" << this->get_row() << ", " << this->get_col() << ") 非方阵，LUP分解需要矩阵为方阵\n";
            throw kalmans::ValueError(oss.str());
        }
        Matrix<double, Row, Col> U(*this);
        auto L = Matrix<double, Row, Col>::Eye();
        auto P = Matrix<double, Row, Col>::Eye();

        for (int j = 0; j < Col - 1; ++j)
        {
            // 选择列主元
            int principal = j;
            T max = fabs(U(principal, j));
            // 只遍历当前主元下方的元素
            for (int i = j + 1; i < Row; ++i)
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
                U.row_swap(j, principal, j, Col);
                // 交换 L(j,1:j-1) <-> L(i,1:j-1)
                L.row_swap(j, principal, 0, j - 1);
                // 交换 P(j,1:n) <-> P(i,1:n)
                P.row_swap(j, principal);
            }

            // 列主元如果为 0，则跳过该列的消元
            if (fabs(U(j, j)) < 1e-7)
                continue;
            // 计算 L 和 U
            for (int i = j + 1; i < Row; ++i)
            {
                L(i, j) = U(i, j) / U(j, j) / 1.0;
                for (int k = j; k < Row; ++k)
                {
                    U(i, k) -= L(i, j) * U(j, k);
                }
            }
        }

        return std::make_tuple(L, U, P);
    }

    template <typename T, int Row, int Col>
    Matrix<double, Row, Col> Matrix<T, Row, Col>::inverse() const
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

        return Matrix<double, Row, Col>(U_ * L_ * P_);
    }

    template <typename T, int Row, int Col>
    Matrix<T, Row, Col> Matrix<T, Row, Col>::Eye()
    {
        Matrix<T, Row, Col> eye;

        for (int i = 0; i < Row; ++i)
        {
            eye(i, i) = 1;
        }

        return eye;
    }

}
