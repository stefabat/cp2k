!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief methods used in the flexible partitioning scheme
!> \par History
!>      04.2006 [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE fp_methods

   USE beta_gamma_psi,                  ONLY: psi
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_iter_string,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE fp_types,                        ONLY: fp_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fac,&
                                              maxfac,&
                                              oorootpi
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fp_eval

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fp_methods'

CONTAINS

! **************************************************************************************************
!> \brief computest the forces and the energy due to the flexible potential & bias,
!>     and writes the weights file
!> \param fp_env ...
!> \param subsys ...
!> \param cell ...
!> \par History
!>      04.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE fp_eval(fp_env, subsys, cell)
      TYPE(fp_type), POINTER                             :: fp_env
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fp_eval'

      CHARACTER(LEN=15)                                  :: tmpstr
      INTEGER                                            :: handle, i, icenter, iparticle, &
                                                            output_unit
      LOGICAL                                            :: zero_weight
      REAL(KIND=dp)                                      :: c, dcdr, kT, r, rab(3), sf, strength
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_list_type), POINTER                  :: particles_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(fp_env))
      CPASSERT(fp_env%use_fp)
      CPASSERT(ASSOCIATED(subsys))
      CALL cp_subsys_get(subsys, particles=particles_list)
      particles => particles_list%els

      ! compute the force due to the reflecting walls
      ! and count the distribution in discrete and contiguous ways
      zero_weight = .FALSE.
      fp_env%restraint_energy = 0.0_dp
      icenter = fp_env%central_atom
      strength = fp_env%strength
      fp_env%i1 = 0; fp_env%i2 = 0; fp_env%o1 = 0; fp_env%o2 = 0
      fp_env%ri1 = 0.0_dp; fp_env%ri2 = 0.0_dp; fp_env%ro1 = 0.0_dp; fp_env%ro2 = 0.0_dp
      fp_env%energy = 0.0_dp

      ! inner particles
      DO i = 1, SIZE(fp_env%inner_atoms)
         iparticle = fp_env%inner_atoms(i)
         rab = particles(iparticle)%r - particles(icenter)%r
         rab = pbc(rab, cell)
         r = SQRT(SUM(rab**2))
         ! constraint wall  (they feel to outer wall)
         IF (r > fp_env%outer_radius) THEN
            zero_weight = .TRUE.
            fp_env%restraint_energy = fp_env%restraint_energy + 0.5_dp*strength*(r - fp_env%outer_radius)**2
            sf = strength*(r - fp_env%outer_radius)/r
            particles(iparticle)%f = particles(iparticle)%f - sf*rab
            particles(icenter)%f = particles(icenter)%f + sf*rab
         END IF
         ! count the distribution
         IF (r > fp_env%inner_radius) THEN
            fp_env%i2 = fp_env%i2 + 1
         ELSE
            fp_env%i1 = fp_env%i1 + 1
         END IF
         ! smooth count the distribution
         CALL smooth_count(r, fp_env%inner_radius, fp_env%smooth_width, c, dcdr)
         fp_env%ri1 = fp_env%ri1 + c
         fp_env%ri2 = fp_env%ri2 + (1.0_dp - c)
      END DO

      ! outer particles
      DO i = 1, SIZE(fp_env%outer_atoms)
         iparticle = fp_env%outer_atoms(i)
         rab = particles(iparticle)%r - particles(icenter)%r
         rab = pbc(rab, cell)
         r = SQRT(SUM(rab**2))
         ! constraint wall (they feel the inner wall)
         IF (r < fp_env%inner_radius) THEN
            zero_weight = .TRUE.
            fp_env%restraint_energy = fp_env%restraint_energy + &
                                      0.5_dp*strength*(r - fp_env%inner_radius)**2
            sf = strength*(r - fp_env%inner_radius)/r
            particles(iparticle)%f = particles(iparticle)%f - sf*rab
            particles(icenter)%f = particles(icenter)%f + sf*rab
         END IF
         ! count the distribution
         IF (r > fp_env%outer_radius) THEN
            fp_env%o2 = fp_env%o2 + 1
         ELSE
            fp_env%o1 = fp_env%o1 + 1
         END IF
         ! smooth count the distribution
         CALL smooth_count(r, fp_env%outer_radius, fp_env%smooth_width, c, dcdr)
         fp_env%ro1 = fp_env%ro1 + c
         fp_env%ro2 = fp_env%ro2 + (1.0_dp - c)
      END DO
      fp_env%energy = fp_env%energy + fp_env%restraint_energy

      ! the combinatorial weight
      i = fp_env%i2 + fp_env%o1
      CPASSERT(i <= maxfac)
      fp_env%comb_weight = (fac(fp_env%i2)*fac(fp_env%o1))/fac(i)

      ! we can add the bias potential now.
      ! this bias has the form
      ! kT * { ln[(o1+i2)!] - ln[o1!] - ln[i2!] }
      ! where the smooth counts are used for o1 and i2
      fp_env%bias_energy = 0.0_dp
      IF (fp_env%bias) THEN
         kT = fp_env%temperature
         fp_env%bias_energy = kT*(LOG_GAMMA(fp_env%ro1 + fp_env%ri2 + 1) - &
                                  LOG_GAMMA(fp_env%ro1 + 1) - LOG_GAMMA(fp_env%ri2 + 1))

         ! and add the corresponding forces
         ! inner particles
         DO i = 1, SIZE(fp_env%inner_atoms)
            iparticle = fp_env%inner_atoms(i)
            rab = particles(iparticle)%r - particles(icenter)%r
            rab = pbc(rab, cell)
            r = SQRT(SUM(rab**2))
            CALL smooth_count(r, fp_env%inner_radius, fp_env%smooth_width, c, dcdr)
            sf = kT*(psi(fp_env%ro1 + fp_env%ri2 + 1) - psi(fp_env%ri2 + 1))*(-dcdr)/r
            particles(iparticle)%f = particles(iparticle)%f - sf*rab
            particles(icenter)%f = particles(icenter)%f + sf*rab
         END DO
         ! outer particles
         DO i = 1, SIZE(fp_env%outer_atoms)
            iparticle = fp_env%outer_atoms(i)
            rab = particles(iparticle)%r - particles(icenter)%r
            rab = pbc(rab, cell)
            r = SQRT(SUM(rab**2))
            CALL smooth_count(r, fp_env%outer_radius, fp_env%smooth_width, c, dcdr)
            sf = kT*(psi(fp_env%ro1 + fp_env%ri2 + 1) - psi(fp_env%ro1 + 1))*(dcdr)/r
            particles(iparticle)%f = particles(iparticle)%f - sf*rab
            particles(icenter)%f = particles(icenter)%f + sf*rab
         END DO
      END IF
      fp_env%energy = fp_env%energy + fp_env%bias_energy
      fp_env%bias_weight = EXP(fp_env%bias_energy/kT)

      ! if this configuration is a valid one, compute its weight
      IF (zero_weight) THEN
         fp_env%weight = 0.0_dp
      ELSE
         fp_env%weight = fp_env%comb_weight*fp_env%bias_weight
      END IF

      ! put weights and other info on file
      logger => cp_get_default_logger()
      output_unit = cp_print_key_unit_nr(logger, fp_env%print_section, "", &
                                         extension=".weights")
      IF (output_unit > 0) THEN
         tmpstr = cp_iter_string(logger%iter_info, fp_env%print_section)
         WRITE (output_unit, '(T2,A15,6(1X,F16.10),4(1X,I4),4(1X,F16.10))') &
            tmpstr, &
            fp_env%weight, fp_env%comb_weight, fp_env%bias_weight, &
            fp_env%energy, fp_env%restraint_energy, fp_env%bias_energy, &
            fp_env%i1, fp_env%i2, fp_env%o1, fp_env%o2, &
            fp_env%ri1, fp_env%ri2, fp_env%ro1, fp_env%ro2
      END IF

      CALL cp_print_key_finished_output(output_unit, logger, fp_env%print_section, &
                                        "")

      CALL timestop(handle)

   END SUBROUTINE fp_eval

! **************************************************************************************************
!> \brief counts in a smooth way (error function with width=width)
!>      if r is closer than r1. Returns 1.0 for the count=c if r<<r1
!>      and the derivative wrt r dcdr
!> \param r ...
!> \param r1 ...
!> \param width ...
!> \param c ...
!> \param dcdr ...
!> \par History
!>      04.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE smooth_count(r, r1, width, c, dcdr)
      REAL(KIND=dp), INTENT(IN)                          :: r, r1, width
      REAL(KIND=dp), INTENT(OUT)                         :: c, dcdr

      REAL(KIND=dp)                                      :: arg

      arg = (r1 - r)/width

      c = (1.0_dp + ERF(arg))/2.0_dp
      dcdr = (-oorootpi/width)*EXP(-arg**2)

   END SUBROUTINE

END MODULE fp_methods
