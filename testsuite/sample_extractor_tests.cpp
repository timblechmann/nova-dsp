#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../nova-dsp/sample_extractor.hpp"
#include "../nova-dsp/cache_aligned_array.hpp"

using namespace nova;

BOOST_AUTO_TEST_CASE( ptr_tests )
{
    float * f;
    BOOST_REQUIRE_EQUAL(get_samples(f), f);

    const float * cf;
    BOOST_REQUIRE_EQUAL(get_samples(cf), cf);
}

BOOST_AUTO_TEST_CASE( aligned_array_tests )
{
    aligned_array<float, 64> array;
    BOOST_REQUIRE_EQUAL(get_samples(array), array.begin());
}
