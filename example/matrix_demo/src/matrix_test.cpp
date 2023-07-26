#include <iostream>
#include <Matrix/Matrix.hpp>
#include <vector>

int main()
{
    std::cout << "Hello My Matrix~\n";
    std::cout << "=========== 构造一个矩阵 ===========\n";
    // 1.1. 默认构造
    kalmans::Matrix<int> my_mat1;
    kalmans::Matrix<int, 3, 3> my_mat2;
    kalmans::Matrix<int, 3, 3> my_mat3(6);
    std::cout << my_mat3;
    // 1.2. 拷贝构造
    kalmans::Matrix my_mat4(my_mat3);
    // 1.3. 移动构造
    kalmans::Matrix my_mat5(std::move(my_mat3));
    std::cout << my_mat5;
    std::cout << my_mat3;

    std::cout << "=========== 为矩阵赋值 ===========\n";
    // 2.1. 赋值 6 个 6
    my_mat2.assign(6, 6);
    std::cout << my_mat2;
    // 2.2. 使用 << 运算符配合初始化列表初始化
    my_mat2 << std::initializer_list({ 6, 5, 4, 3, 2 });
    std::cout << my_mat2;
    // 2.3. 初始化列表赋值
    my_mat3.assign({ 1,2,3,4,5,6,7,8,9 });
    std::cout << my_mat3;
    // 2.4. 所有元素赋值 8
    kalmans::Matrix<int, 2, 3> my_mat6 = 8;
    std::cout << my_mat6;
    // 2.5. 等号运算符赋值
    kalmans::Matrix my_mat7 = my_mat6;
    std::cout << my_mat7;

    std::cout << "=========== 索引与修改 ===========\n";
    std::cout << "  + 修改前mat3 = \n" << my_mat3 << '\n';
    // 3.1. 使用 () 运算符访问与修改
    my_mat3(2, 2) = 66;
    my_mat3(6) = 99;
    // 3.2. 使用 at 函数访问与修改
    my_mat3.at(0, 0) = 100;
    my_mat3.at(2) = 12;
    std::cout << "  + 修改后mat3 = \n" << my_mat3 << '\n';

    std::cout << "=========== 交换矩阵 ===========\n";
    std::cout << "交换前mat2矩阵：\n" << my_mat2 << '\n';
    std::cout << "交换前mat3矩阵：\n" << my_mat3 << '\n';
    // 4.1. swap 函数交换矩阵
    my_mat2.swap(my_mat3);
    std::cout << "交换后mat2矩阵：\n" << my_mat2 << '\n';
    std::cout << "交换后mat3矩阵：\n" << my_mat3 << '\n';

    std::cout << "=========== 矩阵运算 ===========\n";
    // 5.1. 转置运算
    std::cout << "mat3^T = \n" << my_mat3.transpose() << '\n';
    // 5.2. 负号运算
    std::cout << "-mat3 = \n" << -my_mat3 << '\n';
    // 5.3. 加法运算
    std::cout << "mat2 + mat3 = \n" << my_mat2 + my_mat3 << '\n';
    // 5.4. 减法运算
    std::cout << "mat2 - mat3 = \n" << my_mat2 - my_mat3 << '\n';
    std::cout << "-mat3 + mat2 = \n" << -my_mat3 + my_mat2 << '\n';
    // 5.5. 乘法运算
    kalmans::Matrix<int, 2, 3> A{2, 3, 4, 5, 6, 7};
    kalmans::Matrix B = A;
    std::cout << "A * B^T = \n" << A * B.transpose() << '\n';
    // 5.6. LU/LUP分解
    kalmans::Matrix<double, 3, 3> lu_mat1{2, 1, 2, 4, 5, 4, 6, -3, 5};
    kalmans::Matrix<double, 2, 2> lu_mat2{0,1,0,1};
    auto [L1, U1, P1] = lu_mat1.lup_decomposition();
    std::cout << "LUP分解示例1：\n";
    std::cout << lu_mat1 << "=\n";
    std::cout << L1 << "*\n";
    std::cout << U1 << "*\n";
    std::cout << P1;
    auto [L2, U2, P2] = lu_mat2.lup_decomposition();
    std::cout << "LUP分解示例2：\n";
    std::cout << lu_mat2 << "=\n";
    std::cout << L2 << "*\n";
    std::cout << U2 << "*\n";
    std::cout << P2;
    // 5.7. 求逆运算
    auto lu_mat_inv = lu_mat1.inverse();
    std::cout << "验证求逆：\nlu_mat * lu_mat_inv = \n";
    std::cout << lu_mat1 * lu_mat_inv;

    return EXIT_SUCCESS;
}