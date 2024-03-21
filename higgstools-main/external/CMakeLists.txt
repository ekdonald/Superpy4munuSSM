include(FetchContent)

set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(BUILD_TESTING OFF)
set(BUILD_SHARED_LIBS OFF)

FetchContent_Declare(rangev3
                     GIT_REPOSITORY https://github.com/ericniebler/range-v3.git)
FetchContent_Declare(
  magic_enum
  GIT_REPOSITORY https://github.com/Neargye/magic_enum.git
  GIT_TAG v0.7.3)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 8.0.1)

FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.9.2)
set(SPDLOG_FMT_EXTERNAL ON)

FetchContent_Declare(
  json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.10.5)
set(JSON_BuildTests OFF)
set(JSON_Install OFF)
set(JSON_ImplicitConversions OFF)
set(JSON_MultipleHeaders ON)

FetchContent_Declare(
  eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)

if(${PROJECT_NAME}_BUILD_PYTHON_MODULE)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG master)
  set(FETCH_PYBIND11 pybind11)
endif()

FetchContent_MakeAvailable(
  rangev3
  magic_enum
  fmt
  spdlog
  json
  eigen3
  ${FETCH_PYBIND11})
