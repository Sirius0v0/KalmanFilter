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
        static IOFormat set_precision(unsigned pre);
        static IOFormat set_add_head(bool add);
    };

    std::string IOFormat::elem_delimiter = " ";
    std::string IOFormat::row_delimiter = "\n";
    std::string IOFormat::start_delimiter = std::string();
    std::string IOFormat::end_delimiter = std::string();
    unsigned IOFormat::precision = 8;
    bool IOFormat::add_head = true;
    bool IOFormat::cache_elem_delim
        = (elem_delimiter.find_first_not_of(std::string(" \t\n")) != std::string::npos);
    bool IOFormat::cache_row_delim
        = (row_delimiter.find_first_not_of(std::string(" \t\n")) != std::string::npos);
    bool IOFormat::cache_start_delim
        = (start_delimiter.find_first_not_of(std::string(" \t\n")) != std::string::npos);
    bool IOFormat::cache_end_delim
        = (end_delimiter.find_first_not_of(std::string(" \t\n")) != std::string::npos);

    inline void IOFormat::update_cache_flag()
    {
        const std::string white_str(" \t\n");
        IOFormat::cache_elem_delim =
            (elem_delimiter.find_first_not_of(white_str) != std::string::npos);
        IOFormat::cache_row_delim =
            (row_delimiter.find_first_not_of(white_str) != std::string::npos);
        IOFormat::cache_start_delim =
            (start_delimiter.find_first_not_of(white_str) != std::string::npos);
        IOFormat::cache_end_delim =
            (end_delimiter.find_first_not_of(white_str) != std::string::npos);
    }

    inline IOFormat IOFormat::set_elem_delim(std::string elem_delim)
    {
        IOFormat::elem_delimiter = elem_delim;
        return IOFormat();
    }

    inline IOFormat IOFormat::set_row_delim(std::string row_delim)
    {
        IOFormat::row_delimiter = row_delim;
        return IOFormat();
    }

    inline IOFormat IOFormat::set_start_delim(std::string start_delim)
    {
        IOFormat::start_delimiter = start_delim;
        return IOFormat();
    }

    inline IOFormat IOFormat::set_end_delim(std::string end_delim)
    {
        IOFormat::end_delimiter = end_delim;
        return IOFormat();
    }

    inline IOFormat IOFormat::set_precision(unsigned pre)
    {
        IOFormat::precision = pre;
        return IOFormat();
    }

    inline IOFormat IOFormat::set_add_head(bool add)
    {
        if (add)
            IOFormat::add_head = true;
        else
            IOFormat::add_head = false;
        return IOFormat();
    }

    inline std::ostream& operator<< (std::ostream& os, IOFormat)
    {
        return os;
    }

    inline std::istream& operator>> (std::istream& is, IOFormat)
    {
        return is;
    }
    
}
