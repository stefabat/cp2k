!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Parametrization based on GTH pseudo potentials
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_gth
   USE arnoldi_api,                     ONLY: arnoldi_extremal
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks,&
                                              dbcsr_reserve_diag_blocks
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_comm_type
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE pao_param_fock,                  ONLY: pao_calc_U_block_fock
   USE pao_param_methods,               ONLY: pao_calc_AB_from_U,&
                                              pao_calc_grad_lnv_wrt_U
   USE pao_potentials,                  ONLY: pao_calc_gaussian
   USE pao_types,                       ONLY: pao_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              pao_potential_type,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pao_param_init_gth, pao_param_finalize_gth, pao_calc_AB_gth
   PUBLIC :: pao_param_count_gth, pao_param_initguess_gth

CONTAINS

! **************************************************************************************************
!> \brief Initialize the linear potential parametrization
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_init_gth(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init_gth'

      INTEGER                                            :: acol, arow, handle, iatom, idx, ikind, &
                                                            iterm, jatom, maxl, n, natoms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pri, col_blk_size, nterms, &
                                                            row_blk_size
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_V_term, vec_V_terms
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: pao_potentials
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      natom=natoms, &
                      matrix_s=matrix_s, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

      maxl = 0
      ALLOCATE (row_blk_size(natoms), col_blk_size(natoms), nterms(natoms))
      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL pao_param_count_gth(qs_env, ikind, nterms(iatom))
         CALL get_qs_kind(qs_kind_set(ikind), pao_potentials=pao_potentials)
         CPASSERT(SIZE(pao_potentials) == 1)
         maxl = MAX(maxl, pao_potentials(1)%maxl)
      END DO
      CALL init_orbital_pointers(maxl) ! needs to be called before gth_calc_term()

      ! allocate matrix_V_terms
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=blk_sizes_pri)
      col_blk_size = SUM(nterms)
      row_blk_size = blk_sizes_pri**2
      CALL dbcsr_create(pao%matrix_V_terms, &
                        name="PAO matrix_V_terms", &
                        dist=pao%diag_distribution, &
                        matrix_type="N", &
                        row_blk_size=row_blk_size, &
                        col_blk_size=col_blk_size)
      CALL dbcsr_reserve_diag_blocks(pao%matrix_V_terms)
      CALL dbcsr_set(pao%matrix_V_terms, 0.0_dp)

      ! calculate and store poential terms
!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,qs_env,blk_sizes_pri,natoms,nterms) &
!$OMP PRIVATE(iter,arow,acol,iatom,jatom,N,idx,vec_V_terms,block_V_term)
      CALL dbcsr_iterator_start(iter, pao%matrix_V_terms)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, vec_V_terms)
         iatom = arow; CPASSERT(arow == acol)
         n = blk_sizes_pri(iatom)
         DO jatom = 1, natoms
            IF (jatom == iatom) CYCLE ! waste some storage to simplify things later
            DO iterm = 1, nterms(jatom)
               idx = SUM(nterms(1:jatom - 1)) + iterm
               block_V_term(1:n, 1:n) => vec_V_terms(:, idx) ! map column into matrix
               CALL gth_calc_term(qs_env, block_V_term, iatom, jatom, iterm)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      IF (pao%precondition) &
         CALL pao_param_gth_preconditioner(pao, qs_env, nterms)

      DEALLOCATE (row_blk_size, col_blk_size, nterms)
      CALL timestop(handle)
   END SUBROUTINE pao_param_init_gth

! **************************************************************************************************
!> \brief Finalize the GTH potential parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_finalize_gth(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CALL dbcsr_release(pao%matrix_V_terms)
      IF (pao%precondition) THEN
         CALL dbcsr_release(pao%matrix_precon)
         CALL dbcsr_release(pao%matrix_precon_inv)
      END IF

   END SUBROUTINE pao_param_finalize_gth

! **************************************************************************************************
!> \brief Builds the preconditioner matrix_precon and matrix_precon_inv
!> \param pao ...
!> \param qs_env ...
!> \param nterms ...
! **************************************************************************************************
   SUBROUTINE pao_param_gth_preconditioner(pao, qs_env, nterms)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, DIMENSION(:), POINTER                     :: nterms

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_gth_preconditioner'

      INTEGER                                            :: acol, arow, handle, i, iatom, ioffset, &
                                                            j, jatom, joffset, m, n, natoms
      LOGICAL                                            :: arnoldi_converged, converged, found
      REAL(dp)                                           :: eval_max, eval_min
      REAL(dp), DIMENSION(:, :), POINTER                 :: block, block_overlap, block_V_term
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: matrix_gth_overlap
      TYPE(ls_scf_env_type), POINTER                     :: ls_scf_env
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, ls_scf_env=ls_scf_env)
      CALL dbcsr_get_info(pao%matrix_V_terms, group=group)
      natoms = SIZE(nterms)

      CALL dbcsr_create(matrix_gth_overlap, &
                        template=pao%matrix_V_terms, &
                        matrix_type="N", &
                        row_blk_size=nterms, &
                        col_blk_size=nterms)
      CALL dbcsr_reserve_all_blocks(matrix_gth_overlap)
      CALL dbcsr_set(matrix_gth_overlap, 0.0_dp)

      DO iatom = 1, natoms
      DO jatom = 1, natoms
         ioffset = SUM(nterms(1:iatom - 1))
         joffset = SUM(nterms(1:jatom - 1))
         n = nterms(iatom)
         m = nterms(jatom)

         ALLOCATE (block(n, m))
         block = 0.0_dp

         ! can't use OpenMP here block is a pointer and hence REDUCTION(+:block) does work
         CALL dbcsr_iterator_start(iter, pao%matrix_V_terms)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, arow, acol, block_V_term)
            CPASSERT(arow == acol)
            DO i = 1, n
            DO j = 1, m
               block(i, j) = block(i, j) + SUM(block_V_term(:, ioffset + i)*block_V_term(:, joffset + j))
            END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)

         CALL group%sum(block)

         CALL dbcsr_get_block_p(matrix=matrix_gth_overlap, row=iatom, col=jatom, block=block_overlap, found=found)
         IF (ASSOCIATED(block_overlap)) &
            block_overlap = block

         DEALLOCATE (block)
      END DO
      END DO

      !TODO: good setting for arnoldi?
      CALL arnoldi_extremal(matrix_gth_overlap, eval_max, eval_min, max_iter=100, &
                            threshold=1e-2_dp, converged=arnoldi_converged)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| GTH-preconditioner converged, min, max, max/min:", &
         arnoldi_converged, eval_min, eval_max, eval_max/eval_min

      CALL dbcsr_create(pao%matrix_precon, template=matrix_gth_overlap)
      CALL dbcsr_create(pao%matrix_precon_inv, template=matrix_gth_overlap)

      CALL matrix_sqrt_Newton_Schulz(pao%matrix_precon_inv, pao%matrix_precon, matrix_gth_overlap, &
                                     threshold=ls_scf_env%eps_filter, &
                                     order=ls_scf_env%s_sqrt_order, &
                                     max_iter_lanczos=ls_scf_env%max_iter_lanczos, &
                                     eps_lanczos=ls_scf_env%eps_lanczos, &
                                     converged=converged)
      CALL dbcsr_release(matrix_gth_overlap)

      IF (.NOT. converged) &
         CPABORT("PAO: Sqrt of GTH-preconditioner did not converge.")

      CALL timestop(handle)
   END SUBROUTINE pao_param_gth_preconditioner

! **************************************************************************************************
!> \brief Takes current matrix_X and calculates the matrices A and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param gradient ...
!> \param penalty ...
! **************************************************************************************************
   SUBROUTINE pao_calc_AB_gth(pao, qs_env, ls_scf_env, gradient, penalty)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      LOGICAL, INTENT(IN)                                :: gradient
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: penalty

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_AB_gth'

      INTEGER                                            :: handle
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_M, matrix_U

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_create(matrix_U, matrix_type="N", dist=pao%diag_distribution, template=matrix_s(1)%matrix)
      CALL dbcsr_reserve_diag_blocks(matrix_U)

      !TODO: move this condition into pao_calc_U, use matrix_N as template
      IF (gradient) THEN
         CALL pao_calc_grad_lnv_wrt_U(qs_env, ls_scf_env, matrix_M)
         CALL pao_calc_U_gth(pao, matrix_U, matrix_M, pao%matrix_G, penalty)
         CALL dbcsr_release(matrix_M)
      ELSE
         CALL pao_calc_U_gth(pao, matrix_U, penalty=penalty)
      END IF

      CALL pao_calc_AB_from_U(pao, qs_env, ls_scf_env, matrix_U)
      CALL dbcsr_release(matrix_U)
      CALL timestop(handle)
   END SUBROUTINE pao_calc_AB_gth

! **************************************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param matrix_U ...
!> \param matrix_M1 ...
!> \param matrix_G ...
!> \param penalty ...
! **************************************************************************************************
   SUBROUTINE pao_calc_U_gth(pao, matrix_U, matrix_M1, matrix_G, penalty)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(dbcsr_type)                                   :: matrix_U
      TYPE(dbcsr_type), OPTIONAL                         :: matrix_M1, matrix_G
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: penalty

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_U_gth'

      INTEGER                                            :: acol, arow, handle, iatom, idx, iterm, &
                                                            n, natoms
      INTEGER, DIMENSION(:), POINTER                     :: nterms
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: gaps
      REAL(dp), DIMENSION(:), POINTER                    :: world_G, world_X
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_G, block_M1, block_M2, block_U, &
                                                            block_V, block_V_term, block_X, &
                                                            vec_V_terms
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(pao%matrix_X, row_blk_size=nterms, group=group)
      natoms = SIZE(nterms)
      ALLOCATE (gaps(natoms))
      gaps(:) = HUGE(dp)

      ! allocate arrays for world-view
      ALLOCATE (world_X(SUM(nterms)), world_G(SUM(nterms)))
      world_X = 0.0_dp; world_G = 0.0_dp

      ! collect world_X from atomic blocks
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         idx = SUM(nterms(1:iatom - 1))
         world_X(idx + 1:idx + nterms(iatom)) = block_X(:, 1)
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL group%sum(world_X) ! sync world view across MPI ranks

      ! loop over atoms
      CALL dbcsr_iterator_start(iter, matrix_U)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_U)
         iatom = arow; CPASSERT(arow == acol)
         n = SIZE(block_U, 1)
         CALL dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=iatom, col=iatom, block=vec_V_terms, found=found)
         CPASSERT(ASSOCIATED(vec_V_terms))

         ! calculate potential V of i'th atom
         ALLOCATE (block_V(n, n))
         block_V = 0.0_dp
         DO iterm = 1, SIZE(world_X)
            block_V_term(1:n, 1:n) => vec_V_terms(:, iterm) ! map column into matrix
            block_V = block_V + world_X(iterm)*block_V_term
         END DO

         ! calculate gradient block of i'th atom
         IF (.NOT. PRESENT(matrix_G)) THEN
            CALL pao_calc_U_block_fock(pao, iatom=iatom, penalty=penalty, V=block_V, U=block_U, gap=gaps(iatom))

         ELSE ! TURNING POINT (if calc grad) ------------------------------------
            CPASSERT(PRESENT(matrix_M1))
            CALL dbcsr_get_block_p(matrix=matrix_M1, row=iatom, col=iatom, block=block_M1, found=found)
            ALLOCATE (block_M2(n, n))
            CALL pao_calc_U_block_fock(pao, iatom=iatom, penalty=penalty, V=block_V, U=block_U, &
                                       M1=block_M1, G=block_M2, gap=gaps(iatom))
            DO iterm = 1, SIZE(world_G)
               block_V_term(1:n, 1:n) => vec_V_terms(:, iterm) ! map column into matrix
               world_G(iterm) = world_G(iterm) + SUM(block_V_term*block_M2)
            END DO
            DEALLOCATE (block_M2)
         END IF
         DEALLOCATE (block_V)
      END DO
      CALL dbcsr_iterator_stop(iter)

      ! distribute world_G across atomic blocks
      IF (PRESENT(matrix_G)) THEN
         CALL group%sum(world_G) ! sync world view across MPI ranks
         CALL dbcsr_iterator_start(iter, matrix_G)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, arow, acol, block_G)
            iatom = arow; CPASSERT(arow == acol)
            idx = SUM(nterms(1:iatom - 1))
            block_G(:, 1) = world_G(idx + 1:idx + nterms(iatom))
         END DO
         CALL dbcsr_iterator_stop(iter)
      END IF

      DEALLOCATE (world_X, world_G)

      ! sum penalty energies across ranks
      IF (PRESENT(penalty)) &
         CALL group%sum(penalty)

      ! print homo-lumo gap encountered by fock-layer
      CALL group%min(gaps)
      IF (pao%iw_gap > 0) THEN
         DO iatom = 1, natoms
            WRITE (pao%iw_gap, *) "PAO| atom:", iatom, " fock gap:", gaps(iatom)
         END DO
         CALL m_flush(pao%iw_gap)
      END IF

      ! one-line summary
      IF (pao%iw > 0) THEN
         WRITE (pao%iw, "(A,E20.10,A,T71,I10)") " PAO| min_gap:", MINVAL(gaps), " for atom:", MINLOC(gaps)
      END IF

      DEALLOCATE (gaps)
      CALL timestop(handle)

   END SUBROUTINE pao_calc_U_gth

! **************************************************************************************************
!> \brief Returns the number of parameters for given atomic kind
!> \param qs_env ...
!> \param ikind ...
!> \param nparams ...
! **************************************************************************************************
   SUBROUTINE pao_param_count_gth(qs_env, ikind, nparams)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nparams

      INTEGER                                            :: max_projector, maxl, ncombis
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: pao_potentials
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), pao_potentials=pao_potentials)

      IF (SIZE(pao_potentials) /= 1) &
         CPABORT("GTH parametrization requires exactly one PAO_POTENTIAL section per KIND")

      max_projector = pao_potentials(1)%max_projector
      maxl = pao_potentials(1)%maxl

      IF (maxl < 0) &
         CPABORT("GTH parametrization requires non-negative PAO_POTENTIAL%MAXL")

      IF (max_projector < 0) &
         CPABORT("GTH parametrization requires non-negative PAO_POTENTIAL%MAX_PROJECTOR")

      IF (MOD(maxl, 2) /= 0) &
         CPABORT("GTH parametrization requires even-numbered PAO_POTENTIAL%MAXL")

      ncombis = (max_projector + 1)*(max_projector + 2)/2
      nparams = ncombis*(maxl/2 + 1)

   END SUBROUTINE pao_param_count_gth

! **************************************************************************************************
!> \brief Fills the given block_V with the requested potential term
!> \param qs_env ...
!> \param block_V ...
!> \param iatom ...
!> \param jatom ...
!> \param kterm ...
! **************************************************************************************************
   SUBROUTINE gth_calc_term(qs_env, block_V, iatom, jatom, kterm)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: block_V
      INTEGER, INTENT(IN)                                :: iatom, jatom, kterm

      INTEGER                                            :: c, ikind, jkind, lpot, max_l, min_l, &
                                                            pot_max_projector, pot_maxl
      REAL(dp), DIMENSION(3)                             :: Ra, Rab, Rb
      REAL(KIND=dp)                                      :: pot_beta
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: pao_potentials
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      ! get GTH-settings from remote atom
      CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
      CALL get_qs_kind(qs_kind_set(jkind), pao_potentials=pao_potentials)
      CPASSERT(SIZE(pao_potentials) == 1)
      pot_max_projector = pao_potentials(1)%max_projector
      pot_maxl = pao_potentials(1)%maxl
      pot_beta = pao_potentials(1)%beta

      c = 0
      outer: DO lpot = 0, pot_maxl, 2
         DO max_l = 0, pot_max_projector
         DO min_l = 0, max_l
            c = c + 1
            IF (c == kterm) EXIT outer
         END DO
         END DO
      END DO outer

      ! get basis-set of central atom
      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)

      Ra = particle_set(iatom)%r
      Rb = particle_set(jatom)%r
      Rab = pbc(ra, rb, cell)

      block_V = 0.0_dp
      CALL pao_calc_gaussian(basis_set, block_V, Rab=Rab, lpot=lpot, &
                             min_l=min_l, max_l=max_l, beta=pot_beta, weight=1.0_dp)

   END SUBROUTINE gth_calc_term

! **************************************************************************************************
!> \brief Calculate initial guess for matrix_X
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_initguess_gth(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      INTEGER                                            :: acol, arow
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X
      TYPE(dbcsr_iterator_type)                          :: iter

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao) &
!$OMP PRIVATE(iter,arow,acol,block_X)
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         CPASSERT(arow == acol)
         CPASSERT(SIZE(block_X, 2) == 1)

         ! a simplistic guess, which at least makes the atom visible to others
         block_X = 0.0_dp
         block_X(1, 1) = 0.01_dp
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

   END SUBROUTINE pao_param_initguess_gth

END MODULE pao_param_gth
