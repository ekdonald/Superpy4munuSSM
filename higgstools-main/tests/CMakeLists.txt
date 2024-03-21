add_library(prettyprint INTERFACE)
target_include_directories(prettyprint INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(utilities)
add_subdirectory(predictions)
add_subdirectory(bounds)
add_subdirectory(signals)
