#ifndef __MACROS_HPP
#define __MACROS_HPP


#ifdef __MP
#include <gmpxx.h>
    typedef mpf_class value_t;
    #define ZERO( fl) (!cmp( (fl), 0))
    #define EQUAL( fl, val) (!cmp( fl, val))
#else
    #define EPSILON ( 1e-5)
    #define INT( val) static_cast<int>( val)
    #define ABS( val) ((1 - 2*((val) < 0)) * (val))
    #define ZERO( fl) (ABS(fl) <= EPSILON)
    #define EQUAL( al, bl) ZERO( (al) - (bl))
    typedef double value_t;
#endif


using namespace std::string_literals;
using result_t = std::unordered_map<std::string, 
            std::pair<value_t, bool>>;

#endif
