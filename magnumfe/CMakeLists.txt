set(SRC
  Mesher.cpp
  SubMeshInterpolator.cpp
  DofForm.cpp
  DofAssembler.cpp
  CGDofForm.cpp
)

append_magnumfe_src(${CMAKE_CURRENT_SOURCE_DIR} "${SRC}")
