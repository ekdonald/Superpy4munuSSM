add_library(
  HiggsPredictions OBJECT
  CalcHgamgam.cpp
  Channels.cpp
  Clustering.cpp
  Constraints.cpp
  EffectiveCouplings.cpp
  Basics.cpp
  FormatSupport.cpp
  Helpers.cpp
  JsonSupport.cpp
  Particle.cpp
  ParticleData.cpp
  Predictions.cpp
  PredictionsData.cpp
  Process.cpp
  ReferenceModels.cpp
  UncertainMass.cpp)
target_include_directories(
  HiggsPredictions PUBLIC ${PROJECT_SOURCE_DIR}/include
                          ${PROJECT_BINARY_DIR}/include)
target_compile_features(HiggsPredictions PUBLIC cxx_std_17)
set_target_properties(HiggsPredictions PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(HiggsPredictions PRIVATE HiggsUtilities)
target_compile_options(
  HiggsPredictions
  PRIVATE -Wall
          -Wextra
          -pedantic
          -pedantic-errors
          -Wmissing-declarations
          ${HiggsTools_EXTRA_COMPILE_ARGS}
          -Werror)
