#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../nova-dsp/filter.hpp"
#include "../nova-dsp/sample_extractor.hpp"

#include "../nova-dsp/median_filter.hpp"

#include <cmath>

using namespace nova;
using namespace std;

BOOST_AUTO_TEST_CASE( biquad_test )
{
    biquad<float> filter;
    filter.set(-1.95279, 0.956632, 0.00096087,  0.00192174, 0.00096087);

    /* step response */
    for (int i = 0; i != 1024; ++i)
    {
        using namespace std;
        float in = 1.f;
        float out;

        float * ip = &in;
        float * op = &out;

        filter.perform(ip, op, 1);

        BOOST_REQUIRE (abs(out) < 10); /* a reasonable bound */
    }
}


BOOST_AUTO_TEST_CASE( svf_test )
{
    for (float q = 0.5; q < 1000; q*=2)
    {
        state_variable_filter<float> filter;
        filter.set_frequency(0.01);
        filter.set_q(q);

        /* step response */
        for (int i = 0; i != 1024; ++i)
        {
            using namespace std;
            float in = 1.f;
            float low, high, band, notch;

            float * ip = &in;
            float * lp = &low;
            float * hp = &high;
            float * bp = &band;
            float * np = &notch;

            filter.perform(ip, lp, hp, bp, np, 1);

            BOOST_REQUIRE (abs(low) < 10); /* a reasonable bound */
            BOOST_REQUIRE (abs(high) < 10); /* a reasonable bound */
            BOOST_REQUIRE (abs(band) < 10); /* a reasonable bound */
            BOOST_REQUIRE (abs(notch) < 10); /* a reasonable bound */
        }
    }
}



BOOST_AUTO_TEST_CASE( fractional_median_test )
{
    nova::fractional_median_filter<float, 10> filter_;

    filter_.resize(4.5);

    for (int i = 0; i != 1024; ++i)
        filter_.step(i);
}

BOOST_AUTO_TEST_CASE( median_test )
{
    {
        nova::median_filter<float, 10> filter_(4);

        for (int i = 0; i != 6; ++i)
            filter_.step(1.f/i);

        filter_.resize(8);
    }

    {
        nova::median_filter<float, 10> filter_(8);

        for (int i = 0; i != 10; ++i)
            filter_.step(1.f/i);

        filter_.resize(4);
    }
}

BOOST_AUTO_TEST_CASE( static_median_test )
{
    nova::static_median_filter<float, 10> filter_;

    for (int i = 0; i != 1024; ++i)
        filter_.step(1.f/i);
}
