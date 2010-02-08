//  some utility functions
//  Copyright (C) 2006, 2007 Tim Blechmann
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

#ifndef NOVA_DSP_UTILITIES_HPP
#define NOVA_DSP_UTILITIES_HPP

#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>

#include <cassert>
#include <cmath>

namespace nova
{

namespace detail
{


/** \brief compute a decay factor from for a feedback delay line */
template <typename T>
inline T decay_factor(T delay, T decay)
{
    BOOST_STATIC_ASSERT(boost::is_floating_point<T>::value);
    const T log_threshold = log(0.001); // db2amp(-60)
    return exp(log_threshold * delay / decay);
}

/* decay factor for a delay of 1 */
template <typename T>
inline T decay_factor(T decay)
{
    BOOST_STATIC_ASSERT(boost::is_floating_point<T>::value);
    const T log_threshold = log(0.001); // db2amp(-60)
    return exp(log_threshold / decay);
}


template <typename T>
inline T phasor_increment(T frequency, T samplerate)
{
    BOOST_STATIC_ASSERT(boost::is_floating_point<T>::value);
    return frequency / samplerate;
}

} /* namespace detail */
} /* namespace nova */

#endif /* _UTILITIES_HPP */
