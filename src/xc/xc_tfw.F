!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the Thomas-Fermi kinetic energy functional
!>      plus the von Weizsaecker term
!> \par History
!>      JGH (26.02.2003) : OpenMP enabled
!>      fawzi (04.2004)  : adapted to the new xc interface
!> \author JGH (18.02.2002)
! **************************************************************************************************
MODULE xc_tfw
   USE cp_array_utils,                  ONLY: cp_3d_r_cp_type
   USE kinds,                           ONLY: dp
   USE xc_derivative_desc,              ONLY: deriv_norm_drho,&
                                              deriv_norm_drhoa,&
                                              deriv_norm_drhob,&
                                              deriv_rho,&
                                              deriv_rhoa,&
                                              deriv_rhob
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_get_derivative
   USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                              xc_derivative_type
   USE xc_functionals_utilities,        ONLY: set_util
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                              xc_rho_set_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
   REAL(KIND=dp), PARAMETER :: f13 = 1.0_dp/3.0_dp, &
                               f23 = 2.0_dp*f13, &
                               f43 = 4.0_dp*f13, &
                               f53 = 5.0_dp*f13

   PUBLIC :: tfw_lda_info, tfw_lda_eval, tfw_lsd_info, tfw_lsd_eval

   REAL(KIND=dp) :: cf, flda, flsd, fvw
   REAL(KIND=dp) :: eps_rho
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_tfw'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cutoff ...
! **************************************************************************************************
   SUBROUTINE tfw_init(cutoff)

      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      eps_rho = cutoff
      CALL set_util(cutoff)

      cf = 0.3_dp*(3.0_dp*pi*pi)**f23
      flda = cf
      flsd = flda*2.0_dp**f23
      fvw = 1.0_dp/72.0_dp

   END SUBROUTINE tfw_init

! **************************************************************************************************
!> \brief ...
!> \param reference ...
!> \param shortform ...
!> \param needs ...
!> \param max_deriv ...
! **************************************************************************************************
   SUBROUTINE tfw_lda_info(reference, shortform, needs, max_deriv)
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      IF (PRESENT(reference)) THEN
         reference = "Thomas-Fermi-Weizsaecker kinetic energy functional {LDA version}"
      END IF
      IF (PRESENT(shortform)) THEN
         shortform = "TF+vW kinetic energy functional {LDA}"
      END IF
      IF (PRESENT(needs)) THEN
         needs%rho = .TRUE.
         needs%rho_1_3 = .TRUE.
         needs%norm_drho = .TRUE.
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 3

   END SUBROUTINE tfw_lda_info

! **************************************************************************************************
!> \brief ...
!> \param reference ...
!> \param shortform ...
!> \param needs ...
!> \param max_deriv ...
! **************************************************************************************************
   SUBROUTINE tfw_lsd_info(reference, shortform, needs, max_deriv)
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      IF (PRESENT(reference)) THEN
         reference = "Thomas-Fermi-Weizsaecker kinetic energy functional"
      END IF
      IF (PRESENT(shortform)) THEN
         shortform = "TF+vW kinetic energy functional"
      END IF
      IF (PRESENT(needs)) THEN
         needs%rho_spin = .TRUE.
         needs%rho_spin_1_3 = .TRUE.
         needs%norm_drho = .TRUE.
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 3

   END SUBROUTINE tfw_lsd_info

! **************************************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE tfw_lda_eval(rho_set, deriv_set, order)
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, INTENT(in)                                :: order

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tfw_lda_eval'

      INTEGER                                            :: handle, npoints
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: epsilon_rho
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: s
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER :: e_0, e_ndrho, e_ndrho_ndrho, &
         e_rho, e_rho_ndrho, e_rho_ndrho_ndrho, e_rho_rho, e_rho_rho_ndrho, e_rho_rho_rho, grho, &
         r13, rho
      TYPE(xc_derivative_type), POINTER                  :: deriv

      CALL timeset(routineN, handle)

      CALL xc_rho_set_get(rho_set, rho_1_3=r13, rho=rho, &
                          norm_drho=grho, local_bounds=bo, rho_cutoff=epsilon_rho)
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)
      CALL tfw_init(epsilon_rho)

      ALLOCATE (s(npoints))
      CALL calc_s(rho, grho, s, npoints)

      IF (order >= 0) THEN
         deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_0)

         CALL tfw_u_0(rho, r13, s, e_0, npoints)
      END IF
      IF (order >= 1 .OR. order == -1) THEN
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho)
         deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_ndrho)

         CALL tfw_u_1(rho, grho, r13, s, e_rho, e_ndrho, npoints)
      END IF
      IF (order >= 2 .OR. order == -2) THEN
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho_rho)
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_norm_drho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho_ndrho)
         deriv => xc_dset_get_derivative(deriv_set, &
                                         [deriv_norm_drho, deriv_norm_drho], allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrho)

         CALL tfw_u_2(rho, grho, r13, s, e_rho_rho, e_rho_ndrho, &
                      e_ndrho_ndrho, npoints)
      END IF
      IF (order >= 3 .OR. order == -3) THEN
         deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho, deriv_rho], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_rho)
         deriv => xc_dset_get_derivative(deriv_set, &
                                         [deriv_rho, deriv_rho, deriv_norm_drho], allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_ndrho)
         deriv => xc_dset_get_derivative(deriv_set, &
                                         [deriv_rho, deriv_norm_drho, deriv_norm_drho], allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho_ndrho_ndrho)

         CALL tfw_u_3(rho, grho, r13, s, e_rho_rho_rho, e_rho_rho_ndrho, &
                      e_rho_ndrho_ndrho, npoints)
      END IF
      IF (order > 3 .OR. order < -3) THEN
         CPABORT("derivatives bigger than 3 not implemented")
      END IF

      DEALLOCATE (s)
      CALL timestop(handle)
   END SUBROUTINE tfw_lda_eval

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param s ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE calc_s(rho, grho, s, npoints)
      REAL(KIND=dp), DIMENSION(*), INTENT(in)            :: rho, grho
      REAL(KIND=dp), DIMENSION(*), INTENT(out)           :: s
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip

!$OMP     PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP     SHARED(npoints,rho,eps_rho,s,grho)
      DO ip = 1, npoints
         IF (rho(ip) < eps_rho) THEN
            s(ip) = 0.0_dp
         ELSE
            s(ip) = grho(ip)*grho(ip)/rho(ip)
         END IF
      END DO
   END SUBROUTINE calc_s

! **************************************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE tfw_lsd_eval(rho_set, deriv_set, order)
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, INTENT(in)                                :: order

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tfw_lsd_eval'
      INTEGER, DIMENSION(2), PARAMETER :: &
         norm_drho_spin_name = [deriv_norm_drhoa, deriv_norm_drhob], &
         rho_spin_name = [deriv_rhoa, deriv_rhob]

      INTEGER                                            :: handle, i, ispin, npoints
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: epsilon_rho
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: s
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: e_0, e_ndrho, e_ndrho_ndrho, e_rho, &
                                                            e_rho_ndrho, e_rho_ndrho_ndrho, &
                                                            e_rho_rho, e_rho_rho_ndrho, &
                                                            e_rho_rho_rho
      TYPE(cp_3d_r_cp_type), DIMENSION(2)                :: norm_drho, rho, rho_1_3
      TYPE(xc_derivative_type), POINTER                  :: deriv

      CALL timeset(routineN, handle)
      NULLIFY (deriv)
      DO i = 1, 2
         NULLIFY (norm_drho(i)%array, rho(i)%array, rho_1_3(i)%array)
      END DO

      CALL xc_rho_set_get(rho_set, rhoa_1_3=rho_1_3(1)%array, &
                          rhob_1_3=rho_1_3(2)%array, rhoa=rho(1)%array, &
                          rhob=rho(2)%array, norm_drhoa=norm_drho(1)%array, &
                          norm_drhob=norm_drho(2)%array, rho_cutoff=epsilon_rho, &
                          local_bounds=bo)
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)
      CALL tfw_init(epsilon_rho)

      ALLOCATE (s(npoints))

      DO ispin = 1, 2
         CALL calc_s(rho(ispin)%array, norm_drho(ispin)%array, s, npoints)

         IF (order >= 0) THEN
            deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_0)

            CALL tfw_p_0(rho(ispin)%array, &
                         rho_1_3(ispin)%array, s, e_0, npoints)
         END IF
         IF (order >= 1 .OR. order == -1) THEN
            deriv => xc_dset_get_derivative(deriv_set, [rho_spin_name(ispin)], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho)
            deriv => xc_dset_get_derivative(deriv_set, [norm_drho_spin_name(ispin)], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho)

            CALL tfw_p_1(rho(ispin)%array, norm_drho(ispin)%array, &
                         rho_1_3(ispin)%array, s, e_rho, e_ndrho, npoints)
         END IF
         IF (order >= 2 .OR. order == -2) THEN
            deriv => xc_dset_get_derivative(deriv_set, [rho_spin_name(ispin), &
                                                        rho_spin_name(ispin)], allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho)
            deriv => xc_dset_get_derivative(deriv_set, [rho_spin_name(ispin), &
                                                        norm_drho_spin_name(ispin)], allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_ndrho)
            deriv => xc_dset_get_derivative(deriv_set, [norm_drho_spin_name(ispin), &
                                                        norm_drho_spin_name(ispin)], allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrho)

            CALL tfw_p_2(rho(ispin)%array, norm_drho(ispin)%array, &
                         rho_1_3(ispin)%array, s, e_rho_rho, e_rho_ndrho, &
                         e_ndrho_ndrho, npoints)
         END IF
         IF (order >= 3 .OR. order == -3) THEN
            deriv => xc_dset_get_derivative(deriv_set, [rho_spin_name(ispin), &
                                                        rho_spin_name(ispin), rho_spin_name(ispin)], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_rho)
            deriv => xc_dset_get_derivative(deriv_set, [rho_spin_name(ispin), &
                                                        rho_spin_name(ispin), norm_drho_spin_name(ispin)], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_ndrho)
            deriv => xc_dset_get_derivative(deriv_set, [rho_spin_name(ispin), &
                                                        norm_drho_spin_name(ispin), norm_drho_spin_name(ispin)], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_ndrho_ndrho)

            CALL tfw_p_3(rho(ispin)%array, norm_drho(ispin)%array, &
                         rho_1_3(ispin)%array, s, e_rho_rho_rho, e_rho_rho_ndrho, &
                         e_rho_ndrho_ndrho, npoints)
         END IF
         IF (order > 3 .OR. order < -3) THEN
            CPABORT("derivatives bigger than 3 not implemented")
         END IF
      END DO

      DEALLOCATE (s)
      CALL timestop(handle)
   END SUBROUTINE tfw_lsd_eval

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param s ...
!> \param e_0 ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE tfw_u_0(rho, r13, s, e_0, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, r13, s
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_0
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rho,eps_rho,e_0,flda,r13,s,fvw)
      DO ip = 1, npoints

         IF (rho(ip) > eps_rho) THEN

            e_0(ip) = e_0(ip) + flda*r13(ip)*r13(ip)*rho(ip) + fvw*s(ip)

         END IF

      END DO

   END SUBROUTINE tfw_u_0

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param r13 ...
!> \param s ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE tfw_u_1(rho, grho, r13, s, e_rho, e_ndrho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, grho, r13, s
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho, e_ndrho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = f53*flda

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rho,eps_rho,e_rho,e_ndrho,grho,s,r13,f,fvw)
      DO ip = 1, npoints

         IF (rho(ip) > eps_rho) THEN

            e_rho(ip) = e_rho(ip) + f*r13(ip)*r13(ip) - fvw*s(ip)/rho(ip)
            e_ndrho(ip) = e_ndrho(ip) + 2.0_dp*fvw*grho(ip)/rho(ip)

         END IF

      END DO

   END SUBROUTINE tfw_u_1

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param r13 ...
!> \param s ...
!> \param e_rho_rho ...
!> \param e_rho_ndrho ...
!> \param e_ndrho_ndrho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE tfw_u_2(rho, grho, r13, s, e_rho_rho, e_rho_ndrho, e_ndrho_ndrho, &
                      npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, grho, r13, s
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho_rho, e_rho_ndrho, e_ndrho_ndrho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = f23*f53*flda

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rho,eps_rho,e_rho_rho,e_rho_ndrho,e_ndrho_ndrho,grho,f,fvw)
      DO ip = 1, npoints

         IF (rho(ip) > eps_rho) THEN

            e_rho_rho(ip) = e_rho_rho(ip) + f/r13(ip) + 2.0_dp*fvw*s(ip)/(rho(ip)*rho(ip))
            e_rho_ndrho(ip) = e_rho_ndrho(ip) - 2.0_dp*fvw*grho(ip)/(rho(ip)*rho(ip))
            e_ndrho_ndrho(ip) = e_ndrho_ndrho(ip) + 2.0_dp*fvw/rho(ip)

         END IF

      END DO

   END SUBROUTINE tfw_u_2

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param r13 ...
!> \param s ...
!> \param e_rho_rho_rho ...
!> \param e_rho_rho_ndrho ...
!> \param e_rho_ndrho_ndrho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE tfw_u_3(rho, grho, r13, s, e_rho_rho_rho, e_rho_rho_ndrho, &
                      e_rho_ndrho_ndrho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, grho, r13, s
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho_rho_rho, e_rho_rho_ndrho, &
                                                            e_rho_ndrho_ndrho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = -f13*f23*f53*flda

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rho,eps_rho,e_rho_rho_rho,r13,s,e_rho_rho_ndrho,e_rho_ndrho_ndrho,f,fvw)
      DO ip = 1, npoints

         IF (rho(ip) > eps_rho) THEN

            e_rho_rho_rho(ip) = e_rho_rho_rho(ip) + f/(r13(ip)*rho(ip)) &
                                - 6.0_dp*fvw*s(ip)/(rho(ip)*rho(ip)*rho(ip))
            e_rho_rho_ndrho(ip) = e_rho_rho_ndrho(ip) &
                                  + 4.0_dp*fvw*grho(ip)/(rho(ip)*rho(ip)*rho(ip))
            e_rho_ndrho_ndrho(ip) = e_rho_ndrho_ndrho(ip) &
                                    - 2.0_dp*fvw/(rho(ip)*rho(ip))
         END IF

      END DO

   END SUBROUTINE tfw_u_3

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param r13a ...
!> \param sa ...
!> \param e_0 ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE tfw_p_0(rhoa, r13a, sa, e_0, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, r13a, sa
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_0
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints, rhoa,eps_rho,e_0,r13a,sa,flsd,fvw)
      DO ip = 1, npoints

         IF (rhoa(ip) > eps_rho) THEN
            e_0(ip) = e_0(ip) + flsd*r13a(ip)*r13a(ip)*rhoa(ip) + fvw*sa(ip)
         END IF

      END DO

   END SUBROUTINE tfw_p_0

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param grhoa ...
!> \param r13a ...
!> \param sa ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE tfw_p_1(rhoa, grhoa, r13a, sa, e_rho, e_ndrho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, grhoa, r13a, sa
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho, e_ndrho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = f53*flsd

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rhoa,eps_rho,r13a,sa,fvw,grhoa,e_rho,e_ndrho,f)
      DO ip = 1, npoints

         IF (rhoa(ip) > eps_rho) THEN
            e_rho(ip) = e_rho(ip) + f*r13a(ip)*r13a(ip) - fvw*sa(ip)/rhoa(ip)
            e_ndrho(ip) = e_ndrho(ip) + 2.0_dp*fvw*grhoa(ip)/rhoa(ip)
         END IF

      END DO

   END SUBROUTINE tfw_p_1

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param grhoa ...
!> \param r13a ...
!> \param sa ...
!> \param e_rho_rho ...
!> \param e_rho_ndrho ...
!> \param e_ndrho_ndrho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE tfw_p_2(rhoa, grhoa, r13a, sa, e_rho_rho, e_rho_ndrho, &
                      e_ndrho_ndrho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, grhoa, r13a, sa
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho_rho, e_rho_ndrho, e_ndrho_ndrho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = f23*f53*flsd

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rhoa,eps_rho,e_rho_rho,f,fvw,r13a,sa,e_rho_ndrho,e_ndrho_ndrho)
      DO ip = 1, npoints

         IF (rhoa(ip) > eps_rho) THEN
            e_rho_rho(ip) = e_rho_rho(ip) &
                            + f/r13a(ip) + 2.0_dp*fvw*sa(ip)/(rhoa(ip)*rhoa(ip))
            e_rho_ndrho(ip) = e_rho_ndrho(ip) &
                              - 2.0_dp*fvw*grhoa(ip)/(rhoa(ip)*rhoa(ip))
            e_ndrho_ndrho(ip) = e_ndrho_ndrho(ip) + 2.0_dp*fvw/rhoa(ip)
         END IF

      END DO

   END SUBROUTINE tfw_p_2

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param grhoa ...
!> \param r13a ...
!> \param sa ...
!> \param e_rho_rho_rho ...
!> \param e_rho_rho_ndrho ...
!> \param e_rho_ndrho_ndrho ...
!> \param npoints ...
! **************************************************************************************************
   SUBROUTINE tfw_p_3(rhoa, grhoa, r13a, sa, e_rho_rho_rho, e_rho_rho_ndrho, &
                      e_rho_ndrho_ndrho, npoints)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, grhoa, r13a, sa
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_rho_rho_rho, e_rho_rho_ndrho, &
                                                            e_rho_ndrho_ndrho
      INTEGER, INTENT(in)                                :: npoints

      INTEGER                                            :: ip
      REAL(KIND=dp)                                      :: f

      f = -f13*f23*f53*flsd

!$OMP PARALLEL DO PRIVATE(ip) DEFAULT(NONE)&
!$OMP SHARED(npoints,rhoa,eps_rho,e_rho_rho_rho,e_rho_rho_ndrho,e_rho_ndrho_ndrho,f,fvw,sa,grhoa)
      DO ip = 1, npoints

         IF (rhoa(ip) > eps_rho) THEN
            e_rho_rho_rho(ip) = e_rho_rho_rho(ip) &
                                + f/(r13a(ip)*rhoa(ip)) &
                                - 6.0_dp*fvw*sa(ip)/(rhoa(ip)*rhoa(ip)*rhoa(ip))
            e_rho_rho_ndrho(ip) = e_rho_rho_ndrho(ip) &
                                  + 4.0_dp*fvw*grhoa(ip)/(rhoa(ip)*rhoa(ip)*rhoa(ip))
            e_rho_ndrho_ndrho(ip) = e_rho_ndrho_ndrho(ip) &
                                    - 2.0_dp*fvw/(rhoa(ip)*rhoa(ip))
         END IF

      END DO

   END SUBROUTINE tfw_p_3

END MODULE xc_tfw

