#include <Matrix/IOFormat.h>
#include <string>
#include <iostream>
#include <utility>

namespace kalmans
{
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

    void IOFormat::update_cache_flag()
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

    IOFormat IOFormat::set_elem_delim(std::string elem_delim)
    {
        IOFormat::elem_delimiter = std::move(elem_delim);
        return {};
    }

    IOFormat IOFormat::set_row_delim(std::string row_delim)
    {
        IOFormat::row_delimiter = std::move(row_delim);
        return {};
    }

    IOFormat IOFormat::set_start_delim(std::string start_delim)
    {
        IOFormat::start_delimiter = std::move(start_delim);
        return {};
    }

    IOFormat IOFormat::set_end_delim(std::string end_delim)
    {
        IOFormat::end_delimiter = std::move(end_delim);
        return {};
    }

    IOFormat IOFormat::set_precision(const unsigned pre)
    {
        IOFormat::precision = pre;
        return {};
    }

    IOFormat IOFormat::set_add_head(const bool add)
    {
        if (add)
            IOFormat::add_head = true;
        else
            IOFormat::add_head = false;
        return {};
    }

    std::ostream& operator<< (std::ostream& os, IOFormat)
    {
        return os;
    }

    std::istream& operator>> (std::istream& is, IOFormat)
    {
        return is;
    }

}
