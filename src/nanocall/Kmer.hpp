#ifndef __KMER_HPP
#define __KMER_HPP

#include <array>
#include <string>
#include <mutex>

template < unsigned Kmer_Size >
class Kmer
{
public:
    static const unsigned n_states = (1u << (2 * Kmer_Size));
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
    static size_t to_int(const std::array< char, Kmer_Size >& a)
    {
        return to_int(std::string(a.begin(), a.end()));
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
    static unsigned prefix(unsigned i, unsigned k)
    {
        return i >> (2 * (Kmer_Size - k));
    }
    static unsigned suffix(unsigned i, unsigned k)
    {
        return i & ((1u << (2 * k)) - 1);
    }

    /*
     * Precompute, for every kmer i, the maximum k such that suffix(i, k) == prefix(i, k)
     */
    static unsigned max_self_overlap(unsigned i)
    {
        assert(i < n_states);
        std::array< unsigned, n_states > _max_self_overlap;
        bool _inited = false;
        if (not _inited)
        {
            static std::mutex _mutex;
            {
                std::lock_guard< std::mutex > _lock(_mutex);
                if (not _inited) // recheck
                {
                    for (unsigned i = 0; i < n_states; ++i)
                    {
                        _max_self_overlap[i] = 0;
                        for (unsigned k = Kmer_Size - 1; k >= 1; --k)
                        {
                            if (suffix(i, k) == prefix(i, k))
                            {
                                _max_self_overlap[i] = k;
                                break;
                            }
                        }
                    }
                    _inited = true;
                }
            }
        }
        return _max_self_overlap[i];
    }

    /*
     * Precompute neighbours at distance 1 and 2.
     */
    static const std::vector< unsigned >& neighbour_list(unsigned i, unsigned d)
    {
        assert(i < n_states);
        assert(d == 1 or d == 2);
        static std::array< std::array< std::vector< unsigned >, 2 >, 4096 > _neighbour_list;
        static bool _inited = false;
        if (not _inited)
        {
            static std::mutex _mutex;
            {
                std::lock_guard< std::mutex > _lock(_mutex);
                if (not _inited) // recheck
                {
                    for (unsigned i = 0; i < n_states; ++i)
                    {
                        _neighbour_list[i][0].clear();
                        _neighbour_list[i][1].clear();
                        for (unsigned b1 = 0; b1 < 4; ++b1)
                        {
                            unsigned i1 = (suffix(i, Kmer_Size - 1) << 2) + b1;
                            _neighbour_list[i][0].push_back(i1);
                            for (unsigned b2 = 0; b2 < 4; ++b2)
                            {
                                unsigned i2 = (suffix(i1, Kmer_Size - 1) << 2) + b2;
                                _neighbour_list[i][1].push_back(i2);
                            }
                        }
                    }
                    _inited = true;
                }
            }
        }
        return _neighbour_list[i][d - 1];
    }
}; // class Kmer

#endif
