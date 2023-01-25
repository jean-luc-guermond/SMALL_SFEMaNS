enable_language(Fortran)
file(GLOB_RECURSE SMALL_SFEMaNS_SRC ./GAUSS_POINTS_OBJET/*.f90
  ./FEMSUB/*.f90
  ./FEMSUB_OBJET/*.f90)
include_directories($ENV{MKLROOT}/include)                                                                               
add_library(SMALL_SFEMaNS SHARED ${SMALL_SFEMaNS_SRC})
