//  filter base classes
//  Copyright (C) 2010 Tim Blechmann
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; see the file COPYING.  If not, write to
//  the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
//  Boston, MA 02111-1307, USA.

#ifndef NOVA_DSP_DETAIL_FILTER_BASE_HPP
#define NOVA_DSP_DETAIL_FILTER_BASE_HPP

#include <cmath>

#include "../cache_aligned_array.hpp"
#include "../branch_hints.hpp"

namespace nova
{
namespace detail
{

/** array of samples, usable as argument by muladd classes
 *
 * */
template <typename sample_type, unsigned int n>
struct sample_block:
    public nova::aligned_array<sample_type, n>
{
    sample_block operator+=(sample_type val)
    {
        for (unsigned int i = 0; i != n; ++i)
            aligned_array<sample_type, n>::operator[](i) += val;

        return *this;
    }

    sample_block operator*=(sample_type val)
    {
        for (unsigned int i = 0; i != n; ++i)
            aligned_array<sample_type, n>::operator[](i) *= val;

        return *this;
    }
};

template <typename internal_type,
          bool checked>
class filter_base
{
protected:
    static inline internal_type check(internal_type arg)
    {
        if (checked)
            if (likely(std::isnormal(arg)))
                return arg;
            else
                return 0;
        else
            return arg;
    }
};

} /* namespace detail */
} /* namespace nova */

#endif /* NOVA_DSP_DETAIL_FILTER_BASE_HPP */
