add_executable(
  testPredictions
  Basics.cpp
  BsmParticle.cpp
  Channels.cpp
  Clustering.cpp
  Constraints.cpp
  EffectiveCouplings.cpp
  Helpers.cpp
  JsonSupport.cpp
  Particle.cpp
  Predictions.cpp
  Process.cpp
  ReferenceModels.cpp
  UncertainMass.cpp)
target_link_libraries(testPredictions Catch2::Catch2WithMain HiggsPredictions
                      HiggsUtilities prettyprint)
include(CTest)
include(Catch)
catch_discover_tests(testPredictions TEST_PREFIX "Predictions::" PROPERTIES
                     LABELS "Predictions")
