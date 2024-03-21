add_library(
  HiggsUtilities OBJECT
  Algorithm.cpp
  ArithmeticArray.cpp
  Format.cpp
  Json.cpp
  JsonFwd.cpp
  LinearInterpolator.cpp
  Logging.cpp
  RootFinding.cpp)
target_include_directories(HiggsUtilities PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_compile_features(HiggsUtilities PUBLIC cxx_std_17)
set_target_properties(HiggsUtilities PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(
  HiggsUtilities PUBLIC range-v3::range-v3 magic_enum::magic_enum
                        spdlog::spdlog fmt::fmt nlohmann_json::nlohmann_json)
target_compile_options(
  HiggsUtilities
  PRIVATE -Wall
          -Wextra
          -pedantic
          -pedantic-errors
          -Wsign-conversion
          ${HiggsTools_EXTRA_COMPILE_ARGS}
          -Werror)
