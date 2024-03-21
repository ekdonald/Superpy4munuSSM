add_library(
  HiggsBounds OBJECT
  Acceptance.cpp
  Bounds.cpp
  Likelihood.cpp
  Limit.cpp
  limits/BasicLimit.cpp
  limits/LikelihoodLimit.cpp
  limits/ProcessLimits.cpp
  limits/WidthLimits.cpp)
target_include_directories(HiggsBounds PUBLIC ${PROJECT_SOURCE_DIR}/include
                                              ${PROJECT_BINARY_DIR}/include)
set_target_properties(HiggsBounds PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(
  HiggsBounds
  PRIVATE HiggsUtilities
  PUBLIC HiggsPredictions)
target_compile_options(
  HiggsBounds
  PRIVATE -Wall
          -Wextra
          -pedantic
          -pedantic-errors
          -Wmissing-declarations
          ${HiggsTools_EXTRA_COMPILE_ARGS}
          -Werror)
