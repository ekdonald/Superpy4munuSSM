add_library(
  HiggsSignals OBJECT
  JsonSupport.cpp
  Evaluation.cpp
  Measurement.cpp
  Uncertainties.cpp
  MeasurementImpl.cpp
  Signals.cpp
  Measurements/RateMeasurement.cpp)
target_include_directories(
  HiggsSignals PUBLIC ${HiggsBounds_SOURCE_DIR}/include
                      ${HiggsBounds_BINARY_DIR}/include)
set_target_properties(HiggsSignals PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(
  HiggsSignals
  PRIVATE HiggsUtilities Eigen3::Eigen
  PUBLIC HiggsPredictions)
target_compile_options(
  HiggsSignals
  PRIVATE -Wall
          -Wextra
          -pedantic
          -pedantic-errors
          -Wmissing-declarations
          ${HiggsTools_EXTRA_COMPILE_ARGS}
          -Werror)
