set(SRC
  assemble_scalar_product.cpp
  Mesher.cpp
  SubMeshInterpolator.cpp
  DofAssembler.cpp
  DofForm.cpp
)

append_magnumfe_src(${CMAKE_CURRENT_SOURCE_DIR} "${SRC}")
