//  sample extractor function
//  Copyright (C) 2007, 2008 Tim Blechmann
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

#ifndef NOVA_DSP_SAMPLE_EXTRACTOR_HPP
#define NOVA_DSP_SAMPLE_EXTRACTOR_HPP

#include <boost/array.hpp>

namespace nova
{

/* @{ */
inline float * get_samples (float * arg)
{
    return arg;
}

inline const float * get_samples (const float * arg)
{
    return arg;
}

inline double * get_samples (double * arg)
{
    return arg;
}

inline const double * get_samples (const double * arg)
{
    return arg;
}

template <typename sample_type, std::size_t n>
#ifdef __GCC__
inline typename boost::array<sample_type, n>::value_type * __attribute__((aligned(64))) * get_samples(boost::array<sample_type, n> & arg)
#else
inline typename boost::array<sample_type, n>::value_type * get_samples(boost::array<sample_type, n> & arg)
#endif
{
    return arg.begin();
}

template <typename sample_type, std::size_t n>
#ifdef __GCC__
inline const  __attribute__((aligned(64))) * get_samples(boost::array<sample_type, n> const & arg)
#else
inline const sample_type * get_samples(boost::array<sample_type, n> const & arg)
#endif
{
    return arg.begin();
}

inline float * get_samples (float ** arg)
{
    return arg[0];
}

inline const float * get_samples (const float ** arg)
{
    return arg[0];
}

inline double * get_samples (double ** arg)
{
    return arg[0];
}

inline const double * get_samples (const double ** arg)
{
    return arg[0];
}
/* @} */

}

#endif /* NOVA_DSP_SAMPLE_EXTRACTOR_HPP */
