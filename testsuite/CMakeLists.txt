set(tests
    filter.cpp
    interpolation.cpp
    map1d_test.cpp
    mitra_regalia_filters.cpp
    oscillator_tests.cpp
    parameter.cpp
    sample_extractor_tests.cpp
)

# build benchmarks
foreach(test ${tests})
  string(REPLACE .cpp "" test_name ${test} )
  add_executable(${test_name} ${test})
  add_test(${test_name}_run ${EXECUTABLE_OUTPUT_PATH}/${test_name})
endforeach(test)
