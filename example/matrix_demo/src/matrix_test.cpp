#include <fstream>
#include <iostream>
#include <Matrix/Matrix.hpp>
#include <filesystem>

int main()
{
    std::cout << "Hello My Matrix~\n";
    kalmans::IOFormat::set_add_head(false);
    std::cout << "=========== 构造一个矩阵 ===========\n";
    // 1.1. 默认构造
    kalmans::Matrix<int> my_mat1;
    // 1.2. 其他构造
    kalmans::Matrix<int> my_mat2(3, 3);
    kalmans::Matrix<int> my_mat3(3, 3, 6);
    kalmans::Matrix<double> my_mat4(2, 3, { 1,2,3,4,5,6 });
    std::cout << my_mat4;
    // 1.3. 拷贝构造
    kalmans::Matrix my_mat5(my_mat3);
    // 1.4. 移动构造
    kalmans::Matrix my_mat6(std::move(my_mat3));
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
    my_mat3.resize(3, 3);
    my_mat3.assign({ 1,2,3,4,5,6,7,8,9 });
    std::cout << my_mat3;
    // 2.4. 所有元素赋值 8
    kalmans::Matrix<int> my_mat7(3, 4);
    my_mat7 = 8;
    my_mat1 = my_mat7;
    std::cout << my_mat7;
    std::cout << my_mat1;
    // 2.5. 等号运算符赋值
    kalmans::Matrix my_mat8 = my_mat7;
    std::cout << my_mat8;
    // 2.6. 利用矩阵为为矩阵赋值
    auto submat = kalmans::Matrix<double>::Random(1, 3);
    kalmans::Matrix<double> mat9(3, 3);
    mat9.assign(submat, 3);
    std::cout << mat9;

    std::cout << "=========== 索引与修改 ===========\n";
    std::cout << "  + 修改前mat3 = \n" << my_mat3 << '\n';
    // 3.1. 使用 () 运算符访问与修改
    my_mat3(2, 2) = 66;
    my_mat3(6) = 99;
    // 3.2. 使用 at 函数访问与修改
    my_mat3.at(0, 0) = 100;
    my_mat3.at(2) = 12;
    std::cout << "  + 修改后mat3 = \n" << my_mat3 << '\n';

    std::cout << "=========== 矩阵交换 ===========\n";
    std::cout << "交换前mat2矩阵：\n" << my_mat2 << '\n';
    std::cout << "交换前mat3矩阵：\n" << my_mat3 << '\n';
    // 4.1. swap 函数交换矩阵
    my_mat2.swap(my_mat3);
    std::cout << "交换后mat2矩阵：\n" << my_mat2 << '\n';
    std::cout << "交换后mat3矩阵：\n" << my_mat3 << '\n';
    //4.2. row_swap col_swap对整行/列交换
    my_mat3.row_swap(1, 2);
    std::cout << "交换后mat3矩阵：\n" << my_mat3 << '\n';
    //4.3. row_swap col_swap对部分行/列交换
    my_mat3.row_swap(0, 2, 1, 2);
    std::cout << "交换后mat3矩阵：\n" << my_mat3 << '\n';

    std::cout << "=========== 矩阵运算 ===========\n";
    // 5.1. 转置运算
    std::cout << "mat3^T = \n" << my_mat3.transpose() << '\n';
    // 5.2. 负号运算
    std::cout << "-mat3 = \n" << -my_mat3 << '\n';
    // 5.3. 加法运算
    std::cout << "mat2 + mat3 = \n" << my_mat2 + my_mat3 << '\n';
    my_mat2 += my_mat3;
    std::cout << "mat2矩阵：\n" << my_mat2 << '\n';
    // 5.4. 减法运算
    my_mat2 -= my_mat3;
    std::cout << "mat2矩阵：\n" << my_mat2 << '\n';
    std::cout << "mat2 - mat3 = \n" << my_mat2 - my_mat3 << '\n';
    std::cout << "-mat3 + mat2 = \n" << -my_mat3 + my_mat2 << '\n';
    // 5.5. 乘法运算
    kalmans::Matrix<int> A(2, 3, { 2, 3, 4, 5, 6, 7 });
    kalmans::Matrix B = A;
    std::cout << "A * B^T = \n" << A * B.transpose() << '\n';
    auto valuemul = kalmans::Matrix<double>::Random(3, 4);
    std::cout << "随机矩阵：\n" << valuemul;
    valuemul *= 3;
    std::cout << "随机矩阵：\n" << valuemul;
    // 5.6. LU/LUP分解
    kalmans::Matrix<double> lu_mat1(3, 3, { 2, 1, 2, 4, 5, 4, 6, -3, 5 });
    kalmans::Matrix<double> lu_mat2(2, 2, { 0, 1, 0, 1 });
    auto [L1, U1, P1] = lu_mat1.lup_decomposition();
    std::cout << "LUP分解示例1：\n" << lu_mat1 << "=\n" << L1 << "*\n" << U1 << "*\n" << P1;
    auto [L2, U2, P2] = lu_mat2.lup_decomposition();
    std::cout << "LUP分解示例2：\n" << lu_mat2 << "=\n" << L2 << "*\n" << U2 << "*\n" << P2;
    // 5.7. 求逆运算
    auto lu_mat_inv = lu_mat1.inverse();
    std::cout << "验证求逆：\nlu_mat * lu_mat_inv = \n";
    std::cout << kalmans::IOFormat::set_precision(4) << kalmans::IOFormat::set_start_delim("[") <<
        kalmans::IOFormat::set_end_delim("]") << lu_mat1 * lu_mat_inv;

    std::cout << "=========== 特殊矩阵构造 ===========\n";
    // 6.1. 单位阵
    auto eye1 = kalmans::Matrix<int>::Eye(3);
    std::cout << eye1;
    // 6.2. 对角阵
    auto diag1 = kalmans::Matrix<double>::Diag({ 1,2,3,4,5,6 });
    std::cout << diag1;
    // 6.3. 随机阵
    unsigned int seed = 10;
    auto rand1 = kalmans::Matrix<double>::Random(4, 3);
    std::cout << rand1;
    auto rand2 = kalmans::Matrix<double>::Random(3, 3, seed);
    std::cout << rand2;
    auto rand3 = kalmans::Matrix<int>::Random(2, 2, 0, 10);
    std::cout << rand3;
    // 6.4. 对角阵块
    auto blkdiag = kalmans::Matrix<double>::BlkDiag({ rand1, rand2 });
    std::cout << blkdiag;
    auto blkdiag2 = kalmans::Matrix<int>::BlkDiag(kalmans::Matrix{rand3}, 3);
    std::cout << blkdiag2;
    // 6.5. 垂直拼接矩阵
    auto ver_cat_mat = kalmans::Matrix<double>::ver_mat_cat({ rand1, rand2 });
    std::cout << ver_cat_mat;

    std::cout << "=========== 数据读写 ===========\n";
    kalmans::IOFormat::set_add_head(true);
    namespace fs = std::filesystem;
    std::string filename = "./text_mat.txt";

    if (!fs::exists(filename)) {
        // 文件不存在，打开文件
        std::ofstream file(filename);
        if (file.is_open()) {
            // 文件成功打开，进行操作
            file << kalmans::IOFormat::set_precision(8)
                << kalmans::IOFormat::set_elem_delim(" ")
                << kalmans::IOFormat::set_row_delim("\n")
                << kalmans::IOFormat::set_start_delim("[")
                << kalmans::IOFormat::set_end_delim("]")
                << lu_mat_inv;
            file.close();
        }
        else {
            std::cerr << "文件打开失败\n";
        }
    }

    kalmans::Matrix<double> read_mat;
    if (fs::exists(filename)) {
        std::ifstream ifs(filename);
        if (ifs.is_open()) {
            // 文件成功打开，进行操作
            ifs >> kalmans::IOFormat::set_precision(6)
                >> kalmans::IOFormat::set_elem_delim(" ")
                >> kalmans::IOFormat::set_row_delim("\n")
                >> kalmans::IOFormat::set_start_delim("[")
                >> kalmans::IOFormat::set_end_delim("]")
                >> read_mat;
            ifs.close();
        }
        else {
            std::cerr << "文件打开失败\n";
        }
    }
    std::cout << kalmans::IOFormat::set_add_head(false)
        << kalmans::IOFormat::set_elem_delim(" ")
        << kalmans::IOFormat::set_row_delim("\n")
        << "读取到的矩阵：\n" << read_mat;

    return EXIT_SUCCESS;
}