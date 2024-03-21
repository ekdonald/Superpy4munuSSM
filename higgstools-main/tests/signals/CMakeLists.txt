add_executable(testSignals Evaluation.cpp JsonSupport.cpp Measurement.cpp
                           RateMeasurement.cpp Signals.cpp Uncertainties.cpp)
target_link_libraries(
  testSignals
  Catch2::Catch2WithMain
  HiggsUtilities
  HiggsSignals
  HiggsPredictions
  prettyprint
  Eigen3::Eigen)
target_compile_definitions(
  testSignals
  PRIVATE TESTMEASUREMENTS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/testmeasurements/")

include(CTest)
include(Catch)
catch_discover_tests(testSignals TEST_PREFIX "Signals::" PROPERTIES LABELS
                     "Signals")
