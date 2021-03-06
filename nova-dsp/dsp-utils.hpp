//  utility functions
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

#ifndef NOVA_DSP_UTILS_HPP
#define NOVA_DSP_UTILS_HPP

#include <algorithm>

#include "branch_hints.hpp"

namespace nova
{

template <typename I>
inline I wrap_optimistic(I val, const I mod)
{
    if (unlikely(val >= mod))
        do
            val -= mod;
        while (unlikely(val >= mod));
    else
        while (unlikely(val < I(0)))
            val += mod;

    return val;
}

template <typename T>
inline T clip(T t, T low, T hi)
{
    return std::min(hi, std::max(t, low));
}

} /* namespace nova */

#endif /* NOVA_DSP_UTILS_HPP */
