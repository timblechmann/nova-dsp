#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../nova-dsp/mitra_regalia_filters.hpp"
#include "../nova-dsp/sample_extractor.hpp"

using namespace nova;

namespace
{
template <typename shelf_type>
void test_shelf(void)
{
    shelf_type filter;

    filter.set_factor(2);
    filter.set_frequency(0.1);

    /* step response */
    for (int i = 0; i != 10000; ++i)
    {
        using namespace std;
        float in = 1.f;
        float out;

        float * ip = &in;
        float * op = &out;

        filter.perform(ip, op, 1);

        BOOST_REQUIRE (fabs(out) < 10);  /* a reasonable bound */
    }
}
}


BOOST_AUTO_TEST_CASE( mitra_shelf_tests )
{
    test_shelf<mitra_regalia_low_shelf<float, float, true, true> >();
    test_shelf<mitra_regalia_low_shelf<float, float, true, false> >();
    test_shelf<mitra_regalia_low_shelf<float, float, false, true> >();
    test_shelf<mitra_regalia_low_shelf<float, float, false, false> >();

    test_shelf<mitra_regalia_high_shelf<float, float, true, true> >();
    test_shelf<mitra_regalia_high_shelf<float, float, true, false> >();
    test_shelf<mitra_regalia_high_shelf<float, float, false, true> >();
    test_shelf<mitra_regalia_high_shelf<float, float, false, false> >();
}


BOOST_AUTO_TEST_CASE( mitra_eq_test )
{
    mitra_regalia_eq<float, float, true> filter3;

    filter3.set_frequency(0.01);
    filter3.set_bandwidth(0.1);
    filter3.set_factor(2);

    /* step response */
    for (int i = 0; i != 1024; ++i)
    {
        using namespace std;
        float in = 1.f;
        float band;

        float * ip = &in;
        float * bp = &band;

        filter3.perform(ip, bp, 1);

        BOOST_REQUIRE (fabs(band) < 10); /* a reasonable bound */
    }
}
