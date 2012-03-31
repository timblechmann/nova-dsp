//  mitra-regalia filter classes
//  Copyright (C) 2005, 2007, 2010 Tim Blechmann
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

#ifndef NOVA_DSP_MITRA_REGALIA_FILTERS_HPP
#define NOVA_DSP_MITRA_REGALIA_FILTERS_HPP

#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <algorithm>

#include <boost/math/constants/constants.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>


#include "muladd_helpers.hpp"
#include "sample_extractor.hpp"
#include "parameter.hpp"

#include "branch_hints.hpp"
#include "dsp-utils.hpp"

#include "detail/filter_base.hpp"

namespace nova
{

/* the mathematical constants */
#define pi boost::math::constants::pi<internal_type>()

namespace detail
{

template <typename internal_type>
class mitra_regalia_eq
{
protected:
    static inline internal_type compute_sample(internal_type in, internal_type allpass_sample, internal_type K)
    {
        return ( (in + allpass_sample) +
                 (in - allpass_sample) * K) * 0.5;
    }
};

/** \brief low/high-shelf filter based on:
 *
 *  Tunable digital frequency response equalization filters
 *  Regalia, P.; Mitra, S.
 *  Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on
 *  Volume 35, Issue 1, Jan 1987 Page(s): 118 - 120
 *
 * */
template <bool low,
          typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating>
class mitra_regalia_shelf_base:
    protected detail::mitra_regalia_eq<internal_type>,
    protected detail::filter_base<internal_type, checked>
{
protected:
    typedef typename boost::mpl::if_c<interpolating,
                                      detail::exponential_interpolating_parameter<internal_type, 1, 2205>, /* 50ms @ 44100 */
                                      detail::noninterpolating_parameter<internal_type>
                                      >::type gain_parameter_type;

    typedef typename boost::mpl::if_c<interpolating,
                                      detail::linear_interpolating_parameter<internal_type, 2205>, /* 50ms @ 44100 */
                                      detail::noninterpolating_parameter<internal_type>
                                      >::type filter_parameter_type;

    mitra_regalia_shelf_base(void):
        a_(0), K(1), z_(0)
    {}

public:
    template <typename input_buffer_type, typename output_buffer_type>
    inline void perform(input_buffer_type const & in, output_buffer_type & out, uint n)
    {
        muladd_helper_nop<sample_type> ma;
        perform(in, out, n, ma);
    }

    template <typename input_buffer_type, typename output_buffer_type, typename muladd_helper>
    inline void perform(input_buffer_type const & in, output_buffer_type & out, uint n, muladd_helper & ma)
    {
        const sample_type * in_sample = get_samples(in);
        sample_type * out_sample = get_samples(out);

        internal_type z = z_;

        gain_parameter_type k (K);
        filter_parameter_type a (a_);

        assert(n);
        do
            *out_sample++ = ma(compute_sample(*in_sample++, z, a, k));
        while(--n);

        z_ = z;
        K.update(k);
        a_.update(a);
    }

    void set_factor(internal_type k)
    {
        K = k;
    }

    void reset_factor(internal_type k)
    {
        K.reset(k);
    }

    void set_frequency(internal_type freq)
    {
        a_ = compute_a(freq * 2 * pi);
    }

    void reset_frequency(internal_type freq)
    {
        a_.reset(compute_a(freq * 2 * pi));
    }


protected:
    static inline sample_type compute_sample(sample_type in, internal_type & z, const internal_type a,
                                             const internal_type K)
    {
        internal_type allpassed = compute_allpass(in, a, z);
        return detail::mitra_regalia_eq<internal_type>::compute_sample(in, allpassed, K);
    }

    static inline internal_type compute_allpass(internal_type in, const internal_type a, internal_type & z)
    {
        internal_type new_z = in - z*a;
        internal_type ret = z + new_z*a;

        z = detail::filter_base<internal_type, checked>::check(new_z);
        if (low)
            return -ret;
        else
            return ret;
    }

    static inline internal_type compute_a(internal_type omega)
    {
        internal_type tan_omega_2 = std::tan(omega * 0.5);
        return (tan_omega_2 - 1) / (tan_omega_2 + 1);
    }

protected:
    filter_parameter_type a_;
    gain_parameter_type K;
    internal_type z_;
};

template <bool low,
          typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating>
class mitra_regalia_first_order_filter:
    protected detail::mitra_regalia_eq<internal_type>,
    protected detail::filter_base<internal_type, checked>
{
protected:
    typedef typename boost::mpl::if_c<interpolating,
                                      detail::exponential_interpolating_parameter<internal_type, 1, 2205>, /* 50ms @ 44100 */
                                      detail::noninterpolating_parameter<internal_type>
                                      >::type gain_parameter_type;

    typedef typename boost::mpl::if_c<interpolating,
                                      detail::linear_interpolating_parameter<internal_type, 2205>, /* 50ms @ 44100 */
                                      detail::noninterpolating_parameter<internal_type>
                                      >::type filter_parameter_type;

    mitra_regalia_first_order_filter(void):
        a_(0), z_(0)
    {}

public:
    template <typename input_buffer_type, typename output_buffer_type>
    inline void perform(input_buffer_type const & in, output_buffer_type & out, uint n)
    {
        muladd_helper_nop<sample_type> ma;
        perform(in, out, n, ma);
    }

    template <typename input_buffer_type, typename output_buffer_type, typename muladd_helper>
    inline void perform(input_buffer_type const & in, output_buffer_type & out, uint n, muladd_helper & ma)
    {
        const sample_type * in_sample = get_samples(in);
        sample_type * out_sample = get_samples(out);

        internal_type z = z_;

        filter_parameter_type a (a_);

        for (uint i = 0; i != n; ++i)
            out_sample[i] = ma(compute_sample(in_sample[i], z, a));

        z_ = z;
        a_.update(a);
    }

    void set_frequency(internal_type freq)
    {
        a_ = compute_a(freq * 2 * pi);
    }

    void reset_frequency(internal_type freq)
    {
        a_.reset(compute_a(freq * 2 * pi));
    }


protected:
    static inline sample_type compute_sample(sample_type in, internal_type & z, const internal_type a)
    {
        in = in * 0.5;
        internal_type allpassed = compute_allpass(in, a, z);

        if (low)
            return in - allpassed;
        else
            return in + allpassed;
    }

    static inline internal_type compute_allpass(internal_type in, const internal_type a, internal_type & z)
    {
        internal_type new_z = in - z*a;
        internal_type ret = z + new_z*a;

        z = detail::filter_base<internal_type, checked>::check(new_z);
        return -ret;
    }

    static inline internal_type compute_a(internal_type omega)
    {
        internal_type tan_omega_2 = std::tan(omega * 0.5);
        return (tan_omega_2 - 1) / (tan_omega_2 + 1);
    }

protected:
    filter_parameter_type a_;
    internal_type z_;
};

template <bool pass,
          typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating=true>
class mitra_regalia_2nd_order:
    detail::filter_base<internal_type, checked>
{
    typedef boost::array<internal_type, 2> z_type;

    typedef typename boost::mpl::if_c<interpolating,
                                      detail::exponential_interpolating_parameter<internal_type, 1, 2205>, /* 50ms @ 44100 */
                                      detail::noninterpolating_parameter<internal_type>
                                      >::type gain_parameter_type;

    typedef typename boost::mpl::if_c<interpolating,
                                      detail::linear_interpolating_parameter<internal_type, 2205>, /* 50ms @ 44100 */
                                      detail::noninterpolating_parameter<internal_type>
                                      >::type filter_parameter_type;

public:
    mitra_regalia_2nd_order(void):
        a_(0), b_(0)
    {
        z_.assign(0);
    }

    template <typename input_buffer_type, typename output_buffer_type>
    inline void perform(input_buffer_type const & in, output_buffer_type & out, uint n)
    {
        muladd_helper_nop<sample_type> ma;
        perform(in, out, n, ma);
    }

    template <typename input_buffer_type, typename output_buffer_type, typename muladd_helper>
    inline void perform(input_buffer_type const & in, output_buffer_type & out, uint n, muladd_helper & ma)
    {
        const sample_type * in_sample = get_samples(in);
        sample_type * out_sample = get_samples(out);

        z_type z = z_;
        filter_parameter_type a (a_);
        filter_parameter_type b (b_);

        for (uint i = 0; i != n; ++i)
            out_sample[i] = ma(compute_sample(in_sample[i], z, a, b));

        z_ = z;
        a_.update(a);
        b_.update(b);
    }

protected:
    static inline internal_type compute_allpass(const internal_type in, const internal_type a, const internal_type b,
                                                z_type & z)
    {
        internal_type step1 = in - z[1] * a;
        internal_type step2 = step1 - z[0] * b;

        internal_type back1 = z[0] + b * step2;
        internal_type back2 = z[1] + a * step1;

        z[0] = detail::filter_base<internal_type, checked>::check(step2);
        z[1] = detail::filter_base<internal_type, checked>::check(back1);

        return back2;
    }

    static inline sample_type compute_sample(sample_type in, z_type & z, const internal_type a,
                                             const internal_type b)
    {
        in = in * 0.5;
        internal_type allpassed = compute_allpass(in, a, b, z);

        if (pass)
            return in - allpassed;
        else
            return in + allpassed;
    }


    static inline internal_type compute_a(internal_type omega)
    {
        internal_type tan_omega_2 = std::tan(omega * 0.5);

        return (1 - tan_omega_2) / (1 + tan_omega_2);
    }

    static inline internal_type compute_b(internal_type omega_0)
    {
        return -std::cos(omega_0);
    }

public:
    void set_frequency(internal_type omega_0)
    {
        b_ = compute_b(omega_0 * 2 * pi);
    }

    void set_bandwidth(internal_type omega)
    {
        a_ = compute_a(omega * 2 * pi);
    }

    void reset_frequency(internal_type omega_0)
    {
        b_.reset(compute_b(omega_0 * 2 * pi));
    }

    void reset_bandwidth(internal_type omega)
    {
        a_.reset(compute_a(omega * 2 * pi));
    }

protected:
    filter_parameter_type a_;
    filter_parameter_type b_;
    boost::array<internal_type, 2> z_;
};


}

/** \brief low-shelf filter
 * */
template <typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating>
class mitra_regalia_low_shelf:
    public detail::mitra_regalia_shelf_base<true, sample_type, internal_type, checked, interpolating>
{};

/** \brief high-shelf filter
 * */
template <typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating>
class mitra_regalia_high_shelf:
    public detail::mitra_regalia_shelf_base<false, sample_type, internal_type, checked, interpolating>
{};


/** \brief low-pass filter
 * */
template <typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating>
class mitra_regalia_low_pass:
    public detail::mitra_regalia_first_order_filter<true, sample_type, internal_type, checked, interpolating>
{};

/** \brief high-pass filter
 * */
template <typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating>
class mitra_regalia_high_pass:
    public detail::mitra_regalia_first_order_filter<false, sample_type, internal_type, checked, interpolating>
{};

/** \brief band-pass filter
 * */
template <typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating>
class mitra_regalia_band_pass:
    public detail::mitra_regalia_2nd_order<true, sample_type, internal_type, checked, interpolating>
{};

/** \brief high-pass filter
 * */
template <typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating>
class mitra_regalia_band_reject:
    public detail::mitra_regalia_2nd_order<false, sample_type, internal_type, checked, interpolating>
{};


/** \brief bandpass/reject filter based on:
 *
 *  Tunable digital frequency response equalization filters
 *  Regalia, P.; Mitra, S.
 *  Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on
 *  Volume 35, Issue 1, Jan 1987 Page(s): 118 - 120
 *
 *  using Gray and Markel allpass lattice filter structure
 *
 * */
template <typename sample_type,
          typename internal_type,
          bool checked,
          bool interpolating=true>
class mitra_regalia_eq:
    detail::mitra_regalia_eq<internal_type>,
    detail::filter_base<internal_type, checked>
{
    typedef boost::array<internal_type, 2> z_type;

    typedef typename boost::mpl::if_c<interpolating,
                                      detail::exponential_interpolating_parameter<internal_type, 1, 2205>, /* 50ms @ 44100 */
                                      detail::noninterpolating_parameter<internal_type>
                                      >::type gain_parameter_type;

    typedef typename boost::mpl::if_c<interpolating,
                                      detail::linear_interpolating_parameter<internal_type, 2205>, /* 50ms @ 44100 */
                                      detail::noninterpolating_parameter<internal_type>
                                      >::type filter_parameter_type;

public:
    mitra_regalia_eq(void):
        a_(0), b_(0), K(1)
    {
        z_.assign(0);
    }

    template <typename input_buffer_type, typename output_buffer_type>
    inline void perform(input_buffer_type const & in, output_buffer_type & out, uint n)
    {
        detail::muladd_helper_nop<sample_type> ma;
        perform(in, out, n, ma);
    }

    template <typename input_buffer_type, typename output_buffer_type, typename muladd_helper>
    inline void perform(input_buffer_type const & in, output_buffer_type & out, uint n, muladd_helper & ma)
    {
        const sample_type * in_sample = get_samples(in);
        sample_type * out_sample = get_samples(out);

        z_type z = z_;

        gain_parameter_type k (K);
        filter_parameter_type a (a_);
        filter_parameter_type b (b_);

        assert(n);
        do
            *out_sample++ = ma(compute_sample(*in_sample++, z, a, b, k));
        while(--n);

        z_ = z;
        K.update(k);
        a_.update(a);
        b_.update(b);
    }


protected:
    static inline internal_type compute_allpass(internal_type in, const internal_type a, const internal_type b,
                                                z_type & z)
    {
        internal_type step1 = in - z[1] * a;
        internal_type step2 = step1 - z[0] * b;

        internal_type back1 = z[0] + b * step2;
        internal_type back2 = z[1] + a * step1;

        z[0] = detail::filter_base<internal_type, checked>::check(step2);
        z[1] = detail::filter_base<internal_type, checked>::check(back1);

        return back2;
    }


    static inline sample_type compute_sample(sample_type in, z_type & z, const internal_type a,
                                             const internal_type b, const internal_type K)
    {
        internal_type allpassed = compute_allpass(in, a, b, z);
        return detail::mitra_regalia_eq<internal_type>::compute_sample(in, allpassed, K);
    }

    static inline internal_type compute_a(internal_type omega)
    {
        internal_type tan_omega_2 = std::tan(omega * (internal_type)0.5);

        return (1 - tan_omega_2) / (1 + tan_omega_2);
    }

    static inline internal_type compute_b(internal_type omega_0)
    {
        return -std::cos(omega_0);
    }

public:
    void set_frequency(internal_type omega_0)
    {
        b_ = compute_b(omega_0 * 2 * pi);
    }

    void set_bandwidth(internal_type omega)
    {
        a_ = compute_a(omega * 2 * pi);
    }

    void set_factor(internal_type k)
    {
        K = k;
    }

    void reset_frequency(internal_type omega_0)
    {
        b_.reset(compute_b(omega_0 * 2 * pi));
    }

    void reset_bandwidth(internal_type omega)
    {
        a_.reset(compute_a(omega * 2 * pi));
    }

    void reset_factor(internal_type k)
    {
        K.reset(k);
    }

protected:
    filter_parameter_type a_;
    filter_parameter_type b_;
    gain_parameter_type K;
    boost::array<internal_type, 2> z_;
};

#undef pi

} /* namespace nova */

#endif /* NOVA_DSP_MITRA_REGALIA_FILTERS_HPP */
