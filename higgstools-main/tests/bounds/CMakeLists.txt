add_executable(
  testBounds
  Acceptance.cpp
  AppliedLimit.cpp
  BasicLimit.cpp
  Bounds.cpp
  Likelihood.cpp
  LikelihoodLimit.cpp
  Limit.cpp
  ProcessLimits.cpp
  WidthLimits.cpp)
target_link_libraries(testBounds Catch2::Catch2WithMain HiggsUtilities
                      HiggsPredictions HiggsBounds prettyprint)
include(CTest)
include(Catch)
catch_discover_tests(testBounds TEST_PREFIX "Bounds::" PROPERTIES LABELS
                     "Bounds")
target_compile_definitions(
  testBounds PRIVATE TESTLIMITS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/testlimits/")
