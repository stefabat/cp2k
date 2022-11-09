#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(LIBXC "LibXC")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(
    CP2K_LIBXC
    REQUIRED
    IMPORTED_TARGET
    GLOBAL
    libxc>=${LibXC_FIND_VERSION}
    libxcf03
    libxcf90)
endif()

if(NOT CP2K_LIBXC_FOUND)
  foreach(_var xc xcf03 xcf90)
    string(TOUPPER LIB${_var} _var_up)
    cp2k_find_libraries(${_var_up} ${_var})
  endforeach()
endif()

if(NOT CP2K_LIBXC_INCLUDE_DIRS)
  cp2k_include_dirs(LIBXC "xc.h;libxc/xc.h")
endif()

if(CP2K_LIBXC_INCLUDE_DIRS)
  find_package_handle_standard_args(
    LibXC DEFAULT_MSG CP2K_LIBXC_FOUND CP2K_LIBXC_LINK_LIBRARIES
    CP2K_LIBXC_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(LibXC DEFAULT_MSG CP2K_LIBXC_FOUND
                                    CP2K_LIBXC_LINK_LIBRARIES)
endif()
if(CP2K_LIBXC_FOUND AND NOT TARGET CP2K_CP2K_Libxc::xc)
  add_library(CP2K_Libxc::xc INTERFACE IMPORTED)
  if(LIBXC_INCLUDE_DIRS)
    set_target_properties(
      CP2K_Libxc::xc
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBXC_INCLUDE_DIRS}"
                 INTERFACE_LINK_LIBRARIES "${CP2K_LIBXC_LINK_LIBRARIES}")
  else()
    set_target_properties(CP2K_Libxc::xc PROPERTIES INTERFACE_LINK_LIBRARIES
                                                    "${CP2K_LIBXC_LIBRARIES}")
  endif()
  add_library(CP2K_Libxc::xcf90 INTERFACE IMPORTED)
  set_target_properties(
    CP2K_Libxc::xcf90 PROPERTIES INTERFACE_LINK_LIBRARIES
                                 "${CP2K_LIBXCF90_LINK_LIBRARIES}")
  add_library(CP2K_Libxc::xcf03 INTERFACE IMPORTED)
  set_target_properties(
    CP2K_Libxc::xcf03 PROPERTIES INTERFACE_LINK_LIBRARIES
                                 "${CP2K_LIBXCF03_LINK_LIBRARIES}")
endif()

mark_as_advanced(
  CP2K_LIBXC_FOUND CP2K_LIBXC_LINK_LIBRARIES CP2K_LIBXC_INCLUDE_DIRS
  CP2K_LIBXCF03_LINK_LIBRARIES CP2K_LIBXCF90_LINK_LIBRARIES)
