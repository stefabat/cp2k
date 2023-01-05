#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(FFTW3 "Fftw")
# Check if we can use PkgConfig
find_package(PkgConfig)

# First try with pkg
if(PKG_CONFIG_FOUND)
  pkg_search_module(CP2K_FFTW3 fftw3)
  pkg_search_module(CP2K_FFTW3F fftw3f)
  pkg_search_module(CP2K_FFTW3L fftw3l)
  pkg_search_module(CP2K_FFTW3Q fftw3q)
endif()

foreach(_lib fftw3 fftw3f fftw3l fftw3q)
  if(NOT CP2K_${__lib_up}_FOUND)
    set(CP2K_${__lib_up}_PREFIX "${CP2k_FFTW3_PREFIX}")
    cp2k_find_libraries("${__lib_up}" "${_lib}")
    unset(CP2K_${__lib_up}_PREFIX CACHE)
  endif()

  # OMP variant
  foreach(_subtype "mpi" "omp" "threads")
    string(TOUPPER "${_lib}_${_subtype}" _sub_lib)

    # yeah not super deep
    set(CP2K_${_sub_lib}_PREFIX "${CP2K_FFTW3_PREFIX}")
    cp2k_find_libraries("${_sub_lib}" "${_lib}_${_subtype}")
    unset(CP2K_${_sub_lib}_PREFIX CACHE)
  endforeach()
endforeach()

if(NOT CP2K_FFTW3_INCLUDE_DIRS)
  cp2k_include_dirs(FFTW3 "fftw3.h;fftw3/fftw3.h")
endif()

if(CP2K_FFTW3_INCLUDE_DIRS)
  find_package_handle_standard_args(Fftw DEFAULT_MSG CP2K_FFTW3_INCLUDE_DIRS
                                    CP2K_FFTW3_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(Fftw DEFAULT_MSG CP2K_FFTW3_LINK_LIBRARIES)
endif()

foreach(lib_name "fftw3" "fftw3l" "fftw3q" "fftw3f")
  string(TOUPPER "${lib_name}" __lib_name_up)

  if(CP2K_${__lib_name_up}_FOUND AND NOT TARGET CP2K_FFTW3::${lib_name})
    add_library(CP2K_FFTW3::${lib_name} INTERFACE IMPORTED)
    # we do not recheck if the libraries are found when pkg_config is
    # successful.
    set_target_properties(
      CP2K_FFTW3::${lib_name}
      PROPERTIES INTERFACE_LINK_LIBRARIES
                 "${CP2K_${__lib_name_up}_LINK_LIBRARIES}")

    if(CP2K_FFTW3_INCLUDE_DIRS)
      set_target_properties(
        CP2K_FFTW3::${lib_name} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                           "${CP2K_FFTW3_INCLUDE_DIRS}")
    endif()

    foreach(sub_type "threads" "mpi" "omp")
      string(TOUPPER "${lib_name}_${sub_type}" __libs)
      if(CP2K_${__libs}_FOUND AND NOT TARGET
                                  CP2K_FFTW3::${lib_name}_${sub_type})
        add_library(CP2K_FFTW3::${lib_name}_${sub_type} INTERFACE IMPORTED)
        set_target_properties(
          CP2K_FFTW3::${lib_name}_${sub_type}
          PROPERTIES INTERFACE_LINK_LIBRARIES
                     "${CP2K_${__libs}_LINK_LIBRARIES}")
      endif()
    endforeach()
  endif()
endforeach()

set(CP2K_FFTW3_FOUND ON)
mark_as_advanced(
  CP2K_FFTW3_FOUND
  CP2K_FFTW3_PREFIX
  CP2K_FFTW3_INCLUDE_DIRS
  CP2K_FFTW3_MPI
  FFTW3_OMP
  CP2K_FFTW3_THREADS
  CP2K_FFTW3Q_OMP
  CP2K_FFTW3Q_THREADS
  CP2K_FFTW3F_MPI
  CP2K_FFTW3_OMP
  CP2K_FFTW3F_THREADS
  CP2K_FFTW3L_MPI
  CP2K_FFTW3L_OMP
  CP2K_FFTW3L_THREADS)
