cmake_minimum_required(VERSION 3.19)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(ESMF 8.6.0 MODULE REQUIRED)

set(CMAKE_Fortran_COMPILER "${ESMF_F90COMPILER}")

project(GenIO
        VERSION 1.0
        LANGUAGES Fortran)

find_package(OpenMP REQUIRED)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-ffree-line-length-none -Wall -Wextra -Wconversion -Wno-unused -Wno-unused-dummy-argument")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -fbacktrace -Ofast")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -traceback -Ofast")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -O0 -check all -fpe0 -ftrapuv -init=snan,arrays")
else()
  message(WARNING "${CMAKE_Fortran_COMPILER_ID} Fortran compiler will be used with default options")
endif()

add_library(GenIO genio.F90
                  genio_defaults.F90
                  genio_field.F90
                  genio_geom.F90
                  genio_util.F90
                  genio_struct.F90
                  genio_params.F90)
target_link_libraries(GenIO PUBLIC ESMF)
target_link_libraries(GenIO PUBLIC OpenMP::OpenMP_Fortran)

install(TARGETS GenIO
        EXPORT GenIO
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION mod
)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
        DESTINATION mod
        FILES_MATCHING
        PATTERN "*.mod"
        PATTERN "CMakeFiles" EXCLUDE
)
install(EXPORT GenIO DESTINATION cmake)
