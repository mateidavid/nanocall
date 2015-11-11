#ifndef __KMER_HPP
#define __KMER_HPP

#include <array>
#include <string>

template < unsigned Kmer_Size >
class Kmer
{
public:
    static size_t to_int(const std::string& s)
    {
        static std::array< int8_t, 256 > base_to_int;
        static bool table_initialized = false;
        if (not table_initialized)
        {
            for (unsigned i = 0; i < 256; ++i)
            {
                base_to_int[i] = -1;
            }
            base_to_int['A'] = 0;
            base_to_int['C'] = 1;
            base_to_int['G'] = 2;
            base_to_int['T'] = 3;
            table_initialized = true;
        }
        size_t res = 0;
        for (size_t i = 0; i < s.size(); ++i)
        {
            res <<= 2;
            res += base_to_int[static_cast< unsigned >(s[i])];
        }
        return res;
    }
    static std::string to_string(size_t k)
    {
        static const std::string int_to_base("ACGT");
        std::string res;
        for (size_t j = 0; j < Kmer_Size; ++j)
        {
            res += int_to_base[(k >> (2 * (Kmer_Size - j - 1))) & 0x3];
        }
        return res;
    }
    static unsigned min_skip(unsigned k1, unsigned k2)
    {
        if (k1 == k2)
        {
            return 0;
        }
        else
        {
            for (unsigned k = Kmer_Size - 1; k > 0; --k)
            {
                if ((k1 & ((1u << (2 * k)) - 1)) == (k2 >> (2 * (Kmer_Size - k))))
                {
                    return Kmer_Size - k;
                }
            }
            return Kmer_Size;
        }
    }
}; // class Kmer

#endif
