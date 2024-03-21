pybind11_add_module(
  _Higgs
  MODULE
  THIN_LTO
  HiggsPredictions.cpp
  HiggsSignals.cpp
  HiggsBounds.cpp
  Higgs.cpp)
target_link_libraries(
  _Higgs
  PRIVATE HiggsPredictions
          HiggsBounds
          HiggsSignals
          HiggsUtilities
          magic_enum::magic_enum
          range-v3::range-v3
          fmt::fmt)
if(HiggsTools_STATIC_CXX_STDLIB)
  target_link_options(_Higgs PRIVATE -static-libstdc++)
endif()
if(SKBUILD)
  install(TARGETS _Higgs DESTINATION .)
endif()
