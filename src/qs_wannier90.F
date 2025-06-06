!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to Wannier90 code
!> \par History
!>      06.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_wannier90
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_element,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                              kpoint_init_cell_index,&
                                              kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos,&
                                              rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_env_type,&
                                              kpoint_release,&
                                              kpoint_type
   USE machine,                         ONLY: m_datum
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_berry_kpoint_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_scf_diagonalization,          ONLY: do_general_diag_kp
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
   USE wannier90,                       ONLY: wannier_setup
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_wannier90'

   TYPE berry_matrix_type
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER      :: sinmat => NULL(), cosmat => NULL()
   END TYPE berry_matrix_type

   PUBLIC :: wannier90_interface

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param input ...
!> \param logger ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE wannier90_interface(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'wannier90_interface'

      INTEGER                                            :: handle, iw
      LOGICAL                                            :: explicit
      TYPE(section_vals_type), POINTER                   :: w_input

      !--------------------------------------------------------------------------------------------!

      CALL timeset(routineN, handle)
      w_input => section_vals_get_subs_vals(section_vals=input, &
                                            subsection_name="DFT%PRINT%WANNIER90")
      CALL section_vals_get(w_input, explicit=explicit)
      IF (explicit) THEN

         iw = cp_logger_get_default_io_unit(logger)

         IF (iw > 0) THEN
            WRITE (iw, '(/,T2,A)') &
               '!-----------------------------------------------------------------------------!'
            WRITE (iw, '(T32,A)') "Interface to Wannier90"
            WRITE (iw, '(T2,A)') &
               '!-----------------------------------------------------------------------------!'
         END IF

         CALL wannier90_files(qs_env, w_input, iw)

         IF (iw > 0) THEN
            WRITE (iw, '(/,T2,A)') &
               '!--------------------------------End of Wannier90-----------------------------!'
         END IF
      END IF
      CALL timestop(handle)

   END SUBROUTINE wannier90_interface

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE wannier90_files(qs_env, input, iw)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input
      INTEGER, INTENT(IN)                                :: iw

      INTEGER, PARAMETER                                 :: num_nnmax = 12

      CHARACTER(len=2)                                   :: asym
      CHARACTER(len=20), ALLOCATABLE, DIMENSION(:)       :: atom_symbols
      CHARACTER(LEN=256)                                 :: datx
      CHARACTER(len=default_string_length)               :: filename, seed_name
      INTEGER :: i, i_rep, ib, ib1, ib2, ibs, ik, ik2, ikk, ikpgr, ispin, iunit, ix, iy, iz, k, &
         n_rep, nadd, nao, nbs, nexcl, nkp, nmo, nntot, nspins, num_atoms, num_bands, &
         num_bands_tot, num_kpts, num_wann
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: exclude_bands
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nblist, nnlist
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: nncell
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(3)                              :: mp_grid
      INTEGER, DIMENSION(:), POINTER                     :: invals
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: diis_step, do_kpoints, gamma_only, &
                                                            my_kpgrp, mygrp, spinors
      REAL(KIND=dp)                                      :: cmmn, ksign, rmmn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atoms_cart, b_latt, kpt_latt
      REAL(KIND=dp), DIMENSION(3)                        :: bvec
      REAL(KIND=dp), DIMENSION(3, 3)                     :: real_lattice, recip_lattice
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(berry_matrix_type), DIMENSION(:), POINTER     :: berry_matrix
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_mmn, matrix_struct_work
      TYPE(cp_fm_type)                                   :: fm_tmp, mmn_imag, mmn_real
      TYPE(cp_fm_type), DIMENSION(2)                     :: fmk1, fmk2
      TYPE(cp_fm_type), POINTER                          :: fmdummy, fmi, fmr
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_environment_type), POINTER                 :: qs_env_kp
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      !--------------------------------------------------------------------------------------------!

      ! add code for exclude_bands and projectors

      ! generate all arrays needed for the setup call
      CALL section_vals_val_get(input, "SEED_NAME", c_val=seed_name)
      CALL section_vals_val_get(input, "MP_GRID", i_vals=invals)
      CALL section_vals_val_get(input, "WANNIER_FUNCTIONS", i_val=num_wann)
      CALL section_vals_val_get(input, "ADDED_MOS", i_val=nadd)
      mp_grid(1:3) = invals(1:3)
      num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)
      ! excluded bands
      CALL section_vals_val_get(input, "EXCLUDE_BANDS", n_rep_val=n_rep)
      nexcl = 0
      DO i_rep = 1, n_rep
         CALL section_vals_val_get(input, "EXCLUDE_BANDS", i_rep_val=i_rep, i_vals=invals)
         nexcl = nexcl + SIZE(invals)
      END DO
      IF (nexcl > 0) THEN
         ALLOCATE (exclude_bands(nexcl))
         nexcl = 0
         DO i_rep = 1, n_rep
            CALL section_vals_val_get(input, "EXCLUDE_BANDS", i_rep_val=i_rep, i_vals=invals)
            exclude_bands(nexcl + 1:nexcl + SIZE(invals)) = invals(:)
            nexcl = nexcl + SIZE(invals)
         END DO
      END IF
      !
      ! lattice -> Angstrom
      CALL get_qs_env(qs_env, cell=cell)
      CALL get_cell(cell, h=real_lattice, h_inv=recip_lattice)
      real_lattice(1:3, 1:3) = angstrom*real_lattice(1:3, 1:3)
      recip_lattice(1:3, 1:3) = (twopi/angstrom)*TRANSPOSE(recip_lattice(1:3, 1:3))
      ! k-points
      ALLOCATE (kpt_latt(3, num_kpts))
      CALL get_qs_env(qs_env, particle_set=particle_set)
      NULLIFY (kpoint)
      CALL kpoint_create(kpoint)
      kpoint%kp_scheme = "MONKHORST-PACK"
      kpoint%symmetry = .FALSE.
      kpoint%nkp_grid(1:3) = mp_grid(1:3)
      kpoint%verbose = .FALSE.
      kpoint%full_grid = .TRUE.
      kpoint%eps_geo = 1.0e-6_dp
      kpoint%use_real_wfn = .FALSE.
      kpoint%parallel_group_size = 0
      i = 0
      DO ix = 0, mp_grid(1) - 1
         DO iy = 0, mp_grid(2) - 1
            DO iz = 0, mp_grid(3) - 1
               i = i + 1
               kpt_latt(1, i) = REAL(ix, KIND=dp)/REAL(mp_grid(1), KIND=dp)
               kpt_latt(2, i) = REAL(iy, KIND=dp)/REAL(mp_grid(2), KIND=dp)
               kpt_latt(3, i) = REAL(iz, KIND=dp)/REAL(mp_grid(3), KIND=dp)
            END DO
         END DO
      END DO
      kpoint%nkp = num_kpts
      ALLOCATE (kpoint%xkp(3, num_kpts), kpoint%wkp(num_kpts))
      kpoint%wkp(:) = 1._dp/REAL(num_kpts, KIND=dp)
      DO i = 1, num_kpts
         kpoint%xkp(1:3, i) = (angstrom/twopi)*MATMUL(recip_lattice, kpt_latt(:, i))
      END DO
      ! number of bands in calculation
      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao, nmo=num_bands_tot)
      num_bands_tot = MIN(nao, num_bands_tot + nadd)
      num_bands = num_wann
      num_atoms = SIZE(particle_set)
      ALLOCATE (atoms_cart(3, num_atoms))
      ALLOCATE (atom_symbols(num_atoms))
      DO i = 1, num_atoms
         atoms_cart(1:3, i) = particle_set(i)%r(1:3)
         CALL get_atomic_kind(particle_set(i)%atomic_kind, element_symbol=asym)
         atom_symbols(i) = asym
      END DO
      gamma_only = .FALSE.
      spinors = .FALSE.
      ! output
      ALLOCATE (nnlist(num_kpts, num_nnmax))
      ALLOCATE (nncell(3, num_kpts, num_nnmax))
      nnlist = 0
      nncell = 0
      nntot = 0

      IF (iw > 0) THEN
         ! setup
         CALL wannier_setup(mp_grid, num_kpts, real_lattice, recip_lattice, &
                            kpt_latt, nntot, nnlist, nncell, iw)
      END IF

      CALL get_qs_env(qs_env, para_env=para_env)
      CALL para_env%sum(nntot)
      CALL para_env%sum(nnlist)
      CALL para_env%sum(nncell)

      IF (para_env%is_source()) THEN
         ! Write the Wannier90 input file "seed_name.win"
         WRITE (filename, '(A,A)') TRIM(seed_name), ".win"
         CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
         !
         CALL m_datum(datx)
         WRITE (iunit, "(A)") "! Wannier90 input file generated by CP2K "
         WRITE (iunit, "(A,/)") "! Creation date "//TRIM(datx)
         !
         WRITE (iunit, "(A,I5)") "num_wann     = ", num_wann
         IF (num_bands /= num_wann) THEN
            WRITE (iunit, "(A,I5)") "num_bands    = ", num_bands
         END IF
         WRITE (iunit, "(/,A,/)") "length_unit  = bohr "
         WRITE (iunit, "(/,A,/)") "! System"
         WRITE (iunit, "(/,A)") "begin unit_cell_cart"
         WRITE (iunit, "(A)") "bohr"
         DO i = 1, 3
            WRITE (iunit, "(3F12.6)") cell%hmat(i, 1:3)
         END DO
         WRITE (iunit, "(A,/)") "end unit_cell_cart"
         WRITE (iunit, "(/,A)") "begin atoms_cart"
         DO i = 1, num_atoms
            WRITE (iunit, "(A,3F15.10)") atom_symbols(i), atoms_cart(1:3, i)
         END DO
         WRITE (iunit, "(A,/)") "end atoms_cart"
         WRITE (iunit, "(/,A,/)") "! Kpoints"
         WRITE (iunit, "(/,A,3I6/)") "mp_grid      = ", mp_grid(1:3)
         WRITE (iunit, "(A)") "begin kpoints"
         DO i = 1, num_kpts
            WRITE (iunit, "(3F12.6)") kpt_latt(1:3, i)
         END DO
         WRITE (iunit, "(A)") "end kpoints"
         CALL close_file(iunit)
      ELSE
         iunit = -1
      END IF

      ! calculate bands
      NULLIFY (qs_env_kp)
      CALL get_qs_env(qs_env, do_kpoints=do_kpoints)
      IF (do_kpoints) THEN
         ! we already do kpoints
         qs_env_kp => qs_env
      ELSE
         ! we start from gamma point only
         ALLOCATE (qs_env_kp)
         CALL create_kp_from_gamma(qs_env, qs_env_kp)
      END IF
      IF (iw > 0) THEN
         WRITE (unit=iw, FMT="(/,T2,A)") "Start K-Point Calculation ..."
      END IF
      CALL get_qs_env(qs_env=qs_env_kp, para_env=para_env, blacs_env=blacs_env)
      CALL kpoint_env_initialize(kpoint, para_env, blacs_env)
      CALL kpoint_initialize_mos(kpoint, mos, nadd)
      CALL kpoint_initialize_mo_set(kpoint)
      !
      CALL get_qs_env(qs_env=qs_env_kp, sab_orb=sab_nl, dft_control=dft_control)
      CALL kpoint_init_cell_index(kpoint, sab_nl, para_env, dft_control)
      !
      CALL get_qs_env(qs_env=qs_env_kp, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s, &
                      scf_env=scf_env, scf_control=scf_control)
      CALL do_general_diag_kp(matrix_ks, matrix_s, kpoint, scf_env, scf_control, .FALSE., diis_step)
      !
      IF (iw > 0) THEN
         WRITE (iw, '(T69,A)') "... Finished"
      END IF
      !
      ! Calculate and print Overlaps
      !
      IF (para_env%is_source()) THEN
         WRITE (filename, '(A,A)') TRIM(seed_name), ".mmn"
         CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
         CALL m_datum(datx)
         WRITE (iunit, "(A)") "! Wannier90 file generated by CP2K "//TRIM(datx)
         WRITE (iunit, "(3I8)") num_bands, num_kpts, nntot
      ELSE
         iunit = -1
      END IF
      ! create a list of unique b vectors and a table of pointers
      ! nblist(ik,i) -> +/- b_latt(1:3,x)
      ALLOCATE (nblist(num_kpts, nntot))
      ALLOCATE (b_latt(3, num_kpts*nntot))
      nblist = 0
      nbs = 0
      DO ik = 1, num_kpts
         DO i = 1, nntot
            bvec(1:3) = kpt_latt(1:3, nnlist(ik, i)) - kpt_latt(1:3, ik) + nncell(1:3, ik, i)
            ibs = 0
            DO k = 1, nbs
               IF (SUM(ABS(bvec(1:3) - b_latt(1:3, k))) < 1.e-6_dp) THEN
                  ibs = k
                  EXIT
               END IF
               IF (SUM(ABS(bvec(1:3) + b_latt(1:3, k))) < 1.e-6_dp) THEN
                  ibs = -k
                  EXIT
               END IF
            END DO
            IF (ibs /= 0) THEN
               ! old lattice vector
               nblist(ik, i) = ibs
            ELSE
               ! new lattice vector
               nbs = nbs + 1
               b_latt(1:3, nbs) = bvec(1:3)
               nblist(ik, i) = nbs
            END IF
         END DO
      END DO
      ! calculate all the operator matrices (a|bvec|b)
      ALLOCATE (berry_matrix(nbs))
      DO i = 1, nbs
         NULLIFY (berry_matrix(i)%cosmat)
         NULLIFY (berry_matrix(i)%sinmat)
         bvec(1:3) = twopi*MATMUL(TRANSPOSE(cell%h_inv(1:3, 1:3)), b_latt(1:3, i))
         CALL build_berry_kpoint_matrix(qs_env_kp, berry_matrix(i)%cosmat, &
                                        berry_matrix(i)%sinmat, bvec)
      END DO
      ! work matrices for MOs (all group)
      kp => kpoint%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1, 1), nmo=nmo)
      NULLIFY (matrix_struct_work)
      CALL cp_fm_struct_create(matrix_struct_work, nrow_global=nao, &
                               ncol_global=nmo, &
                               para_env=para_env, &
                               context=blacs_env)
      CALL cp_fm_create(fm_tmp, matrix_struct_work)
      DO i = 1, 2
         CALL cp_fm_create(fmk1(i), matrix_struct_work)
         CALL cp_fm_create(fmk2(i), matrix_struct_work)
      END DO
      ! work matrices for Mmn(k,b) integrals
      NULLIFY (matrix_struct_mmn)
      CALL cp_fm_struct_create(matrix_struct_mmn, nrow_global=nmo, &
                               ncol_global=nmo, &
                               para_env=para_env, &
                               context=blacs_env)
      CALL cp_fm_create(mmn_real, matrix_struct_mmn)
      CALL cp_fm_create(mmn_imag, matrix_struct_mmn)
      ! allocate some work matrices
      ALLOCATE (rmatrix, cmatrix)
      CALL dbcsr_create(rmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)
      !
      CALL get_kpoint_info(kpoint=kpoint, cell_to_index=cell_to_index)
      NULLIFY (fmdummy)
      nspins = dft_control%nspins
      DO ispin = 1, nspins
         ! loop over all k-points
         DO ik = 1, num_kpts
            ! get the MO coefficients for this k-point
            my_kpgrp = (ik >= kpoint%kp_range(1) .AND. ik <= kpoint%kp_range(2))
            IF (my_kpgrp) THEN
               ikk = ik - kpoint%kp_range(1) + 1
               kp => kpoint%kp_env(ikk)%kpoint_env
               CPASSERT(SIZE(kp%mos, 1) == 2)
               fmr => kp%mos(1, ispin)%mo_coeff
               fmi => kp%mos(2, ispin)%mo_coeff
               CALL cp_fm_copy_general(fmr, fmk1(1), para_env)
               CALL cp_fm_copy_general(fmi, fmk1(2), para_env)
            ELSE
               NULLIFY (fmr, fmi, kp)
               CALL cp_fm_copy_general(fmdummy, fmk1(1), para_env)
               CALL cp_fm_copy_general(fmdummy, fmk1(2), para_env)
            END IF
            ! loop over all connected neighbors
            DO i = 1, nntot
               ! get the MO coefficients for the connected k-point
               ik2 = nnlist(ik, i)
               mygrp = (ik2 >= kpoint%kp_range(1) .AND. ik2 <= kpoint%kp_range(2))
               IF (mygrp) THEN
                  ikk = ik2 - kpoint%kp_range(1) + 1
                  kp => kpoint%kp_env(ikk)%kpoint_env
                  CPASSERT(SIZE(kp%mos, 1) == 2)
                  fmr => kp%mos(1, ispin)%mo_coeff
                  fmi => kp%mos(2, ispin)%mo_coeff
                  CALL cp_fm_copy_general(fmr, fmk2(1), para_env)
                  CALL cp_fm_copy_general(fmi, fmk2(2), para_env)
               ELSE
                  NULLIFY (fmr, fmi, kp)
                  CALL cp_fm_copy_general(fmdummy, fmk2(1), para_env)
                  CALL cp_fm_copy_general(fmdummy, fmk2(2), para_env)
               END IF
               !
               ! transfer realspace overlaps to connected k-point
               ibs = nblist(ik, i)
               ksign = SIGN(1.0_dp, REAL(ibs, KIND=dp))
               ibs = ABS(ibs)
               CALL dbcsr_set(rmatrix, 0.0_dp)
               CALL dbcsr_set(cmatrix, 0.0_dp)
               CALL rskp_transform(rmatrix, cmatrix, rsmat=berry_matrix(ibs)%cosmat, ispin=1, &
                                   xkp=kpoint%xkp(1:3, ik2), cell_to_index=cell_to_index, sab_nl=sab_nl, &
                                   is_complex=.FALSE., rs_sign=ksign)
               CALL rskp_transform(cmatrix, rmatrix, rsmat=berry_matrix(ibs)%sinmat, ispin=1, &
                                   xkp=kpoint%xkp(1:3, ik2), cell_to_index=cell_to_index, sab_nl=sab_nl, &
                                   is_complex=.TRUE., rs_sign=ksign)
               !
               ! calculate M_(mn)^(k,b)
               CALL cp_dbcsr_sm_fm_multiply(rmatrix, fmk2(1), fm_tmp, nmo)
               CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, fmk1(1), fm_tmp, 0.0_dp, mmn_real)
               CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, fmk1(2), fm_tmp, 0.0_dp, mmn_imag)
               CALL cp_dbcsr_sm_fm_multiply(rmatrix, fmk2(2), fm_tmp, nmo)
               CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, fmk1(1), fm_tmp, 1.0_dp, mmn_imag)
               CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, fmk1(2), fm_tmp, -1.0_dp, mmn_real)
               CALL cp_dbcsr_sm_fm_multiply(cmatrix, fmk2(1), fm_tmp, nmo)
               CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, fmk1(1), fm_tmp, 1.0_dp, mmn_imag)
               CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, fmk1(2), fm_tmp, -1.0_dp, mmn_real)
               CALL cp_dbcsr_sm_fm_multiply(cmatrix, fmk2(2), fm_tmp, nmo)
               CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, fmk1(1), fm_tmp, -1.0_dp, mmn_real)
               CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, fmk1(2), fm_tmp, -1.0_dp, mmn_imag)
               !
               ! write to output file
               IF (para_env%is_source()) THEN
                  WRITE (iunit, "(2I8,3I5)") ik, ik2, nncell(1:3, ik, i)
               END IF
               DO ib2 = 1, nmo
                  DO ib1 = 1, nmo
                     CALL cp_fm_get_element(mmn_real, ib1, ib2, rmmn)
                     CALL cp_fm_get_element(mmn_imag, ib1, ib2, cmmn)
                     IF (para_env%is_source()) THEN
                        WRITE (iunit, "(2E30.14)") rmmn, cmmn
                     END IF
                  END DO
               END DO
               !
            END DO
         END DO
      END DO
      DO i = 1, nbs
         CALL dbcsr_deallocate_matrix_set(berry_matrix(i)%cosmat)
         CALL dbcsr_deallocate_matrix_set(berry_matrix(i)%sinmat)
      END DO
      DEALLOCATE (berry_matrix)
      CALL cp_fm_struct_release(matrix_struct_work)
      DO i = 1, 2
         CALL cp_fm_release(fmk1(i))
         CALL cp_fm_release(fmk2(i))
      END DO
      CALL cp_fm_release(fm_tmp)
      CALL cp_fm_struct_release(matrix_struct_mmn)
      CALL cp_fm_release(mmn_real)
      CALL cp_fm_release(mmn_imag)
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      !
      IF (para_env%is_source()) THEN
         CALL close_file(iunit)
      END IF
      !
      ! Calculate and print Projections
      !
      ! Print eigenvalues
      nspins = dft_control%nspins
      kp => kpoint%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1, 1), nmo=nmo)
      ALLOCATE (eigval(nmo))
      CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range, xkp=xkp)
      IF (para_env%is_source()) THEN
         WRITE (filename, '(A,A)') TRIM(seed_name), ".eig"
         CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
      ELSE
         iunit = -1
      END IF
      !
      DO ik = 1, nkp
         my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
         DO ispin = 1, nspins
            IF (my_kpgrp) THEN
               ikpgr = ik - kp_range(1) + 1
               kp => kpoint%kp_env(ikpgr)%kpoint_env
               CALL get_mo_set(kp%mos(1, ispin), eigenvalues=eigenvalues)
               eigval(1:nmo) = eigenvalues(1:nmo)
            ELSE
               eigval(1:nmo) = 0.0_dp
            END IF
            CALL kpoint%para_env_inter_kp%sum(eigval)
            eigval(1:nmo) = eigval(1:nmo)*evolt
            ! output
            IF (iunit > 0) THEN
               DO ib = 1, nmo
                  WRITE (iunit, "(2I8,F24.14)") ib, ik, eigval(ib)
               END DO
            END IF
         END DO
      END DO
      IF (para_env%is_source()) THEN
         CALL close_file(iunit)
      END IF
      !
      ! clean up
      DEALLOCATE (kpt_latt, atoms_cart, atom_symbols, eigval)
      DEALLOCATE (nnlist, nncell)
      DEALLOCATE (nblist, b_latt)
      IF (nexcl > 0) THEN
         DEALLOCATE (exclude_bands)
      END IF
      IF (do_kpoints) THEN
         NULLIFY (qs_env_kp)
      ELSE
         CALL qs_env_release(qs_env_kp)
         DEALLOCATE (qs_env_kp)
         NULLIFY (qs_env_kp)
      END IF

      CALL kpoint_release(kpoint)

   END SUBROUTINE wannier90_files

END MODULE qs_wannier90
