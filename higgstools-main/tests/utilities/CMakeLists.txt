add_executable(testUtilities Algorithm.cpp ArithmeticArray.cpp Json.cpp
                             LinearInterpolator.cpp Logging.cpp)
target_link_libraries(testUtilities Catch2::Catch2WithMain HiggsUtilities
                      prettyprint)
include(CTest)
include(Catch)
catch_discover_tests(testUtilities TEST_PREFIX "Utilities::" PROPERTIES LABELS
                     "Utilities")
