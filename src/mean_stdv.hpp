#ifndef __MEAN_STDV_HPP
#define __MEAN_STDV_HPP

#include <cmath>
#include <tuple>

#include "algorithm_extra.hpp"

/*
 * Compute mean and stdv of a sequence of samples.
 */
template < typename Float_Type, class Input_Iterator, class Unary_Function = detail::Identity >
std::pair< Float_Type, Float_Type >
get_mean_stdv(Input_Iterator it_begin, Input_Iterator it_end,
              Unary_Function&& f = Unary_Function())
{
    Float_Type s = 0.0;
    Float_Type s2 = 0.0;
    long unsigned n = 0;
    for (Input_Iterator it = it_begin; it != it_end; ++it)
    {
        s += f(*it);
        s2 += f(*it) * f(*it);
        ++n;
    }
    Float_Type mean = n > 0? s / n : (Float_Type)0;
    Float_Type stdv = n > 1? std::sqrt((s2 - s * mean * 2.0 + mean * mean * (Float_Type)n)/(n - 1)) : (Float_Type)0;
    return std::make_pair(mean, stdv);
} // get_mean_stdv

template < typename Float_Type, class Input_Range, class Unary_Function = detail::Identity  >
std::pair< Float_Type, Float_Type >
get_mean_stdv(const Input_Range& rg, Unary_Function&& f = Unary_Function())
{
    return get_mean_stdv< Float_Type >(rg.begin(), rg.end(), f);
}

#endif
