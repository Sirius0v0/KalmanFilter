#pragma once
#include <random>

namespace kalmans
{
    template<typename T = double>
    class Random
    {
    public:
        static T uniform(unsigned int seed = std::random_device()());
        static T uniform(T min, T max, unsigned int seed = std::random_device()());

        static T normal(unsigned int seed = std::random_device()());
        static T normal(T mean, T stddev, unsigned int seed = std::random_device()());
    };

    template <typename T = double>
    T Random<T>::uniform(unsigned seed)
    {
        T min = 0;
        T max = 1;
        static std::mt19937 engine(seed);
        static unsigned int  last_uniform_seed = 0;
        if (last_uniform_seed != seed)
        {
            engine.seed(seed);
            last_uniform_seed = seed;
        }
        std::uniform_real_distribution<> dist(min, max);
        return dist(engine);
    }

    template <typename T = double>
    T Random<T>::uniform(T min, T max, unsigned seed)
    {
        static std::mt19937 engine(seed);
        static unsigned int last_uniform_seed = 0;
        if (last_uniform_seed != seed)
        {
            engine.seed(seed);
            last_uniform_seed = seed;
        }
        std::uniform_real_distribution<> dist(min, max);
        return dist(engine);
    }

    template <typename T = double>
    T Random<T>::normal(unsigned seed)
    {
        T min = 0;
        T max = 1;
        static std::mt19937 engine(seed);
        static unsigned int last_normal_seed = 0;
        if (last_normal_seed != seed)
        {
            engine.seed(seed);
            last_normal_seed = seed;
        }
        std::normal_distribution<> dist(min, max);
        return dist(engine);
    }

    template <typename T = double>
    T Random<T>::normal(T mean, T stddev, unsigned seed)
    {
        static std::mt19937 engine(seed);
        static unsigned int last_normal_seed = 0;
        if (last_normal_seed != seed)
        {
            engine.seed(seed);
            last_normal_seed = seed;
        }
        std::normal_distribution<> dist(mean, stddev);
        return dist(engine);
    }
}
