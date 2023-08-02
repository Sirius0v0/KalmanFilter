#pragma once

#include <string>
#include <iostream>

namespace kalmans
{
    class IOFormat
    {
    public:
        static std::string elem_delimiter;
        static std::string row_delimiter;
        static std::string start_delimiter;
        static std::string end_delimiter;
        static unsigned precision;
        static bool add_head;
        // 如果分隔符为空白，读取时忽略
        static bool cache_elem_delim;
        static bool cache_row_delim;
        static bool cache_start_delim;
        static bool cache_end_delim;

        static void update_cache_flag();

        static IOFormat set_elem_delim(std::string elem_delim);
        static IOFormat set_row_delim(std::string row_delim);
        static IOFormat set_start_delim(std::string strat_delim);
        static IOFormat set_end_delim(std::string end_delim);
        static IOFormat set_precision(const unsigned pre);
        static IOFormat set_add_head(const bool add);

        friend std::ostream& operator<< (std::ostream& os, IOFormat);
        friend std::istream& operator>> (std::istream& is, IOFormat);
    };
}
