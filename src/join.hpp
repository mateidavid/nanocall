#ifndef __JOIN_HPP
#define __JOIN_HPP

#include <iostream>
#include <string>

namespace join_ns
{

namespace detail
{

struct Identity
{
    template < typename U >
    U operator () (const U& v) const { return v; }
}; // struct Identity

/**
 * Temporary object that will cause the given range to be printed to an ostream.
 * Create such objects implicitly using the join() function below.
 */
template < typename Input_Iterator, typename Unary_Function >
struct join_object
{
    Input_Iterator it_begin;
    Input_Iterator it_end;
    const std::string& sep;
    Unary_Function& f;

    join_object(Input_Iterator _it_begin, Input_Iterator _it_end, const std::string& _sep, Unary_Function&& _f)
        : it_begin(_it_begin), it_end(_it_end), sep(_sep), f(_f) {}

    friend std::ostream& operator << (std::ostream& os, const join_object& o)
    {
        bool first = true;
        for (auto it = o.it_begin; it != o.it_end; ++it)
        {
            if (not first)
            {
                os << o.sep;
            }
            first = false;
            os << o.f(*it);
        }
        return os;
    }
}; // struct join_object

} // namespace detail

/**
 * Create join_object given a pair of begin/end iterators.
 */
template < typename Input_Iterator, typename Unary_Function = detail::Identity >
detail::join_object< Input_Iterator, Unary_Function >
join(Input_Iterator it_begin, Input_Iterator it_end, const std::string& sep,
     Unary_Function&& f = Unary_Function())
{
    return detail::join_object< Input_Iterator, Unary_Function >(it_begin, it_end, sep, std::move(f));
} // join()

/**
 * Create join_object given a range.
 */
template < typename Input_Range, typename Unary_Function = detail::Identity >
auto
join(Input_Range&& rg, const std::string& sep,
     Unary_Function&& f = Unary_Function())
    -> detail::join_object< decltype(rg.end()), Unary_Function >
{
    return detail::join_object< decltype(rg.end()), Unary_Function >(rg.begin(), rg.end(), sep, std::move(f));
} // join()

} // namespace join

#endif

#ifdef SAMPLE_JOIN

/*

Compile with:

g++ -std=c++11 -Wall -Wextra -pedantic -D SAMPLE_JOIN -x c++ join.hpp -o sample-join

*/

#include <vector>

using namespace std;
using namespace join_ns;

int main()
{
    vector< unsigned > v{ 42, 1, 15 };
    cout << "v: " << join(v.begin(), v.end(), ", ") << endl;
    cout << "v[0]: " << join(v.begin(), v.begin() + 1, ", ") << endl;
    cout << "v+1: " << join(v.begin(), v.end(), ", ", [] (unsigned i) { return i + 1; }) << endl;
    cout << "v_rg: " << join(v, ",") << endl;
    cout << "v+2_rg: " << join(v, ",", [] (unsigned i) { return i + 2; }) << endl;
}

#endif
