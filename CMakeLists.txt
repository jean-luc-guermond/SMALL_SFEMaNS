enable_language(Fortran)
file(GLOB_RECURSE SMALL_SFEMaNS_SRC $ENV{SMALL_SFEMaNS_DIR}/GAUSS_POINTS_OBJET/*.f90
  $ENV{SMALL_SFEMaNS_DIR}/FEMSUB/*.f90
  $ENV{SMALL_SFEMaNS_DIR}/FEMSUB_OBJET/*.f90)
include_directories($ENV{MKLROOT}/include)                                                                               
add_library(SMALL_SFEMaNS SHARED ${SMALL_SFEMaNS_SRC})
