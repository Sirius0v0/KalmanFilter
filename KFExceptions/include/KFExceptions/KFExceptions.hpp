#pragma once

#include <string>
#include <iostream>

namespace kalmans
{
    struct KFException : public std::logic_error
    {
        explicit KFException(const std::string& message)
            :logic_error(message){}
    };

    struct LengthException : public KFException
    {
        explicit LengthException(const std::string& message)
            :KFException("[LengthException] "+message) {}
    };

    struct OutOfRangeException : public KFException
    {
        explicit OutOfRangeException(const std::string& message)
            :KFException("[OutOfRangeException] "+message){}
    };
}
