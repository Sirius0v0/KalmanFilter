#include <fstream>
#include <iostream>
#include <filesystem>
#include <chrono>

#include <Matrix/Matrix.hpp>
#include <filter_demo/ekf_demo.h>
namespace fs = std::filesystem;

int main()
{
    std::cout << "Hello My Filter~\n";
    constexpr int STATE_DIM = 9;
    constexpr int U_DIM = 3;
    constexpr int MEASURE_DIM = 4;
    MyEKFilter filter(STATE_DIM, U_DIM, MEASURE_DIM);

    auto dir = fs::path(__FILE__).parent_path();
    const fs::path measure_filename = fs::absolute(dir / "measure.txt");

    kalmans::Matrix<double> measure_data;
    if (fs::exists(measure_filename)) {
        std::ifstream ifs(measure_filename);
        if (ifs.is_open()) {
            // 文件成功打开，进行操作
            ifs >> kalmans::IOFormat::set_precision(6)
                >> kalmans::IOFormat::set_elem_delim(" ")
                >> kalmans::IOFormat::set_row_delim("\n")
                >> kalmans::IOFormat::set_start_delim("")
                >> kalmans::IOFormat::set_end_delim("")
                >> measure_data;
            ifs.close();
        }
        else {
            std::cerr << "文件打开失败\n";
        }
    }

    auto P0 = kalmans::Matrix<double>::Diag({
        1, 0.1, 0.01,
        1, 0.1, 0.01,
        1, 0.1, 0.01
        });
    kalmans::Matrix<double> x(STATE_DIM, 1);
    x(0) = 1500e3;
    x(3) = 10e3;
    x(7) = -250;

    kalmans::IOFormat::set_add_head(false);
    kalmans::IOFormat::set_precision(2);
    int i = 0;
    // 滤波器初始化
    filter.init(x, P0);
    std::cout << "[Time=" << 2 * i << "]" << filter.get_state().transpose() << '\n';

    kalmans::Matrix<double> z(MEASURE_DIM, 1);
    kalmans::Matrix<double> u;
    kalmans::Matrix<double> filter_value(measure_data.get_row() + 1, STATE_DIM);
    filter_value.assign(filter.get_state(), STATE_DIM * i);

    // 计时开始
    auto begin = std::chrono::high_resolution_clock::now();

    for (; i < measure_data.get_row(); ++i)
    {

        for (int j = 0; j < MEASURE_DIM; ++j)
            z(j) = measure_data(i, j);

        filter.step(u, z);
        filter_value.assign(filter.get_state(), STATE_DIM * (i + 1));
        std::cout << "[Time=" << 2 * (i + 1) << "]" << filter.get_state().transpose() << '\n';
    }

    // 计时结束
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    auto time_stamp = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::system_clock::now().time_since_epoch()
    );
    const std::string filter_filename = ("filter_res_" + std::to_string(time_stamp.count()) + ".txt");

    if (!fs::exists(filter_filename)) {
        // 文件不存在，打开文件
        std::ofstream file(filter_filename);
        if (file.is_open()) {
            // 文件成功打开，进行操作
            file << kalmans::IOFormat::set_precision(6)
                << kalmans::IOFormat::set_elem_delim("\t")
                << kalmans::IOFormat::set_row_delim("\n")
                << kalmans::IOFormat::set_start_delim("")
                << kalmans::IOFormat::set_end_delim("")
                << filter_value;
            std::cout << "滤波完成，用时" << elapsed.count() * 1e-9 << "秒, 数据保存在"
                << filter_filename << "中\n";
            file.close();
        }
        else {
            std::cerr << "文件打开失败\n";
        }
    }
    return EXIT_SUCCESS;
}