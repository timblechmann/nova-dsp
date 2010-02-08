#include "benchmark_helpers.hpp"
#include "../nova-dsp/median_filter.hpp"

using namespace nova;

aligned_array<float, 64> in, out;
fractional_median_filter<float, 100> filter_;

void __noinline__ bench_1(unsigned int n)
{
    filter_.perform(in, out, 64);
}


int main()
{
    out.assign(0);
    fill_container(in);

    const unsigned int iterations = 500000;
    filter_.resize(45);

    run_bench(boost::bind(bench_1, 64), iterations);
}
