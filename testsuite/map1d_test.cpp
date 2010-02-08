#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#define BOOST_DISABLE_ASSERTS 1


const int size = 1024;


#include "../nova-dsp/chaos/map1d.hpp"


namespace
{
float dummy;
}


using namespace nova::chaos;

BOOST_AUTO_TEST_CASE( map1d_test )
{
    logistic_map<float, double> lm;
    select_interpolator<0, logistic_map<float, double> >::type interp;

    boost::array<float, 1024> out;

    lm.perform(out, 1024);
    lm.perform(out.begin(), 1024);

    lm.perform(interp, out, 1024);
    lm.perform(interp, out.begin(), 1024);

    dummy = out.back();
}

BOOST_AUTO_TEST_CASE( map1d_test_sample_hold )
{
    logistic_map<float, double> lm;
    select_interpolator<1, logistic_map<float, double> >::type interp;

    boost::array<float, 1024> out;

    lm.perform(interp, out, 1024);
    lm.perform(interp, out.begin(), 1024);

    dummy = out.back();
}

BOOST_AUTO_TEST_CASE( map1d_test_linear )
{
    logistic_map<float, double> lm;
    select_interpolator<2, logistic_map<float, double> >::type interp;

    boost::array<float, 1024> out;

    lm.perform(interp, out, 1024);
    lm.perform(interp, out.begin(), 1024);

    dummy = out.back();
}

BOOST_AUTO_TEST_CASE( map1d_test_quadratic )
{
    logistic_map<float, double> lm;
    select_interpolator<3, logistic_map<float, double> >::type interp;

    boost::array<float, 1024> out;

    lm.perform(interp, out, 1024);
    lm.perform(interp, out.begin(), 1024);

    dummy = out.back();
}

BOOST_AUTO_TEST_CASE( delayed_logistic_test )
{
    typedef delayed_logistic_map<float, double> map_type;
    map_type lm;
    select_interpolator<3, map_type >::type interp;

    boost::array<float, 1024> out;

    lm.perform(out.begin(), 1024);
    lm.perform(interp, out.begin(), 1024);

    dummy = out.back();
}
