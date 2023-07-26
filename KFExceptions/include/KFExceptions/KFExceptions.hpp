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

    struct LengthError : public KFException
    {
        explicit LengthError(const std::string& message)
            :KFException("[LengthError] "+message) {}
    };

    struct OutOfRangeError : public KFException
    {
        explicit OutOfRangeError(const std::string& message)
            :KFException("[OutOfRangeError] "+message){}
    };

    struct ValueError : public KFException
    {
        explicit ValueError(const std::string& message)
            :KFException("[ValueError] "+message){}
    };

    struct ZeroDivisionError : public KFException
    {
        explicit ZeroDivisionError(const std::string& message)
            :KFException("[ZeroDivisionError] "+message){}
    };
}
