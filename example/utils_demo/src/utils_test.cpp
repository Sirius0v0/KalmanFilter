#include <iostream>
#include <iomanip>
#include <map>
#include <Utils/random.hpp>

int main()
{
    std::cout << "Hello My Utils~\n";
    std::cout << "=========== 均匀分布随机数生成 ===========\n";
    // 1.1. 默认 [0,1) double 均匀分布
    std::cout << kalmans::Random<>::uniform() << '\n';
    // 1.2. 指定 seed 的 [0,1) 均匀分布
    unsigned int seed = 1;
    std::cout << kalmans::Random<>::uniform(seed) << '\n';
    // 1.3. 指定 [min, max) 均匀分布
    std::cout << kalmans::Random<int>::uniform(1,3) << '\n';
    // 1.4. 指定 seed 的 [min, max) 均匀分布
    std::cout << kalmans::Random<int>::uniform(1,3,seed) << '\n';

    std::cout << "=========== 正态分布随机数生成 ===========\n";
    // 1.5. 正态分布的生成方法与均匀分布类似
    std::map<int, int> hist{};
    for (int i = 0; i < 1e5; ++i)
    {
        ++hist[kalmans::Random<int>::normal(5, 2)];
    }
    for (auto p : hist)
    {
        std::cout << std::setw(2)
            << p.first << " " << std::string(p.second / 2000, '*') << '\n';
    }

    return EXIT_SUCCESS;
}