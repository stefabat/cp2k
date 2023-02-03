#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ ${MPI_MODE} != "intelmpi" ] && exit 0
rm -f "${BUILDDIR}/setup_intelmpi"

INTELMPI_CFLAGS=""
INTELMPI_LDFLAGS=""
INTELMPI_LIBS=""
mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_intelmpi}" in
  __INSTALL__)
    echo "==================== Installing Intel MPI ===================="
    echo '__INSTALL__ is not supported; please manually install Intel MPI'
    exit 1
    ;;
  __SYSTEM__)
    echo "==================== Finding Intel MPI from system paths ===================="
    check_command mpiexec "intelmpi" && MPIRUN="$(command -v mpiexec)" || exit 1
    if [ "${with_intel}" != "__DONTUSE__" ]; then
      check_command mpiicc "intelmpi" && MPICC="$(command -v mpiicc)" || exit 1
      check_command mpiicpc "intelmpi" && MPICXX="$(command -v mpiicpc)" || exit 1
      check_command mpiifort "intelmpi" && MPIFC="$(command -v mpiifort)" || exit 1
    else
      check_command mpicc "intelmpi" && MPICC="$(command -v mpicc)" || exit 1
      check_command mpicxx "intelmpi" && MPICXX="$(command -v mpicxx)" || exit 1
      check_command mpifort "intelmpi" && MPIFC="$(command -v mpifort)" || exit 1
    fi
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    # include path is already handled by compiler wrapper scripts (can cause wrong mpi.mod with GNU Fortran)
    # add_include_from_paths INTELMPI_CFLAGS "mpi.h" $INCLUDE_PATHS
    add_lib_from_paths INTELMPI_LDFLAGS "libmpi.*" $LIB_PATHS
    check_lib -lmpi "intelmpi"
    check_lib -lmpicxx "intelmpi"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking INTELMPI to user paths ===================="
    pkg_install_dir="${with_intelmpi}"
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    check_command ${pkg_install_dir}/bin/mpiexec "intel" && MPIRUN="${pkg_install_dir}/bin/mpiexec" || exit 1
    if [ "${with_intel}" != "__DONTUSE__" ]; then
      check_command ${pkg_install_dir}/bin/mpiicc "intel" && MPICC="${pkg_install_dir}/bin/mpiicc" || exit 1
      check_command ${pkg_install_dir}/bin/mpiicpc "intel" && MPICXX="${pkg_install_dir}/bin/mpiicpc" || exit 1
      check_command ${pkg_install_dir}/bin/mpiifort "intel" && MPIFC="${pkg_install_dir}/bin/mpiifort" || exit 1
    else
      check_command ${pkg_install_dir}/bin/mpicc "intel" && MPICC="${pkg_install_dir}/bin/mpicc" || exit 1
      check_command ${pkg_install_dir}/bin/mpicxx "intel" && MPICXX="${pkg_install_dir}/bin/mpicxx" || exit 1
      check_command ${pkg_install_dir}/bin/mpifort "intel" && MPIFC="${pkg_install_dir}/bin/mpifort" || exit 1
    fi
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    # include path is already handled by compiler wrapper scripts (can cause wrong mpi.mod with GNU Fortran)
    #INTELMPI_CFLAGS="-I'${pkg_install_dir}/include'"
    INTELMPI_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "${with_intelmpi}" != "__DONTUSE__" ]; then
  INTELMPI_LIBS="-lmpi -lmpicxx"
  cat << EOF > "${BUILDDIR}/setup_intelmpi"
export MPI_MODE="${MPI_MODE}"
export MPIRUN="${MPIRUN}"
export MPICC="${MPICC}"
export MPICXX="${MPICXX}"
export MPIFC="${MPIFC}"
export MPIFORT="${MPIFORT}"
export MPIF77="${MPIF77}"
export INTELMPI_CFLAGS="${INTELMPI_CFLAGS}"
export INTELMPI_LDFLAGS="${INTELMPI_LDFLAGS}"
export INTELMPI_LIBS="${INTELMPI_LIBS}"
export MPI_CFLAGS="${INTELMPI_CFLAGS}"
export MPI_LDFLAGS="${INTELMPI_LDFLAGS}"
export MPI_LIBS="${INTELMPI_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__parallel -D__MPI_F08|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${INTELMPI_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${INTELMPI_LDFLAGS}|)"
export CP_LIBS="\${CP_LIBS} IF_MPI(${INTELMPI_LIBS}|)"
EOF
  if [ "${with_intelmpi}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_intelmpi"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include"
EOF
  fi
  cat "${BUILDDIR}/setup_intelmpi" >> ${SETUPFILE}
fi

load "${BUILDDIR}/setup_intelmpi"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "intelmpi"
