set(benchmarks
    median_filter_benchmarks.cpp
)

# build benchmarks
foreach(benchmark ${benchmarks})
  string(REPLACE .cpp "" benchmark_name ${benchmark} )
  add_executable(${benchmark_name} ${benchmark})
endforeach(benchmark)
