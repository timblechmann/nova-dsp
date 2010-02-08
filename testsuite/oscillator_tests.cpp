#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../nova-dsp/sample_extractor.hpp"
#include "../nova-dsp/oscillator.hpp"
#include "../nova-dsp/sine_oscillator.hpp"

using namespace nova;

BOOST_AUTO_TEST_CASE( table_oscillators )
{
    table_lookup_oscillator<> osc;
    fixed_table_lookup_oscillator<> osc2;
}

BOOST_AUTO_TEST_CASE( sine_oscillators )
{
    sine_oscillator<> sinosc;
}
