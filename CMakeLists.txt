enable_language(Fortran)
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug")

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release")

ADD_CUSTOM_TARGET(native
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Native")
file(GLOB_RECURSE SMALL_SFEMaNS_SRC ./GAUSS_POINTS_OBJET/*.f90
  ./FEMSUB/*.f90
  ./FEMSUB_OBJET/*.f90)
include_directories($ENV{MKLROOT}/include)                                                                               
add_library(SMALL_SFEMaNS SHARED ${SMALL_SFEMaNS_SRC})
