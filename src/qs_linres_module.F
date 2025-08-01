!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains the setup for  the calculation of properties by linear response
!>      by the application of second order density functional perturbation theory.
!>      The knowledge of the ground state energy, density and wavefunctions is assumed.
!>      Uses the self consistent approach.
!>      Properties that can be calculated : none
!> \par History
!>       created 06-2005 [MI]
!> \author MI
! **************************************************************************************************
MODULE qs_linres_module
   USE bibliography,                    ONLY: Ditler2021,&
                                              Ditler2022,&
                                              Weber2009,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type,&
                                              use_qmmm,&
                                              use_qs_force
   USE input_constants,                 ONLY: lr_current,&
                                              lr_none,&
                                              ot_precond_full_all,&
                                              ot_precond_full_kinetic,&
                                              ot_precond_full_single,&
                                              ot_precond_full_single_inverse,&
                                              ot_precond_none,&
                                              ot_precond_s_inverse
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_dcdr,                         ONLY: apt_dR,&
                                              apt_dR_localization,&
                                              dcdr_build_op_dR,&
                                              dcdr_response_dR,&
                                              prepare_per_atom
   USE qs_dcdr_utils,                   ONLY: dcdr_env_cleanup,&
                                              dcdr_env_init,&
                                              dcdr_print
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_linres_current,               ONLY: current_build_chi,&
                                              current_build_current
   USE qs_linres_current_utils,         ONLY: current_env_cleanup,&
                                              current_env_init,&
                                              current_response
   USE qs_linres_epr_nablavks,          ONLY: epr_nablavks
   USE qs_linres_epr_ownutils,          ONLY: epr_g_print,&
                                              epr_g_so,&
                                              epr_g_soo,&
                                              epr_g_zke,&
                                              epr_ind_magnetic_field
   USE qs_linres_epr_utils,             ONLY: epr_env_cleanup,&
                                              epr_env_init
   USE qs_linres_issc_utils,            ONLY: issc_env_cleanup,&
                                              issc_env_init,&
                                              issc_issc,&
                                              issc_print,&
                                              issc_response
   USE qs_linres_methods,               ONLY: linres_localize
   USE qs_linres_nmr_shift,             ONLY: nmr_shift,&
                                              nmr_shift_print
   USE qs_linres_nmr_utils,             ONLY: nmr_env_cleanup,&
                                              nmr_env_init
   USE qs_linres_op,                    ONLY: current_operators,&
                                              issc_operators,&
                                              polar_operators,&
                                              polar_operators_local,&
                                              polar_operators_local_wannier
   USE qs_linres_polar_utils,           ONLY: polar_env_init,&
                                              polar_polar,&
                                              polar_print,&
                                              polar_response
   USE qs_linres_types,                 ONLY: &
        current_env_type, dcdr_env_type, epr_env_type, get_polar_env, issc_env_type, &
        linres_control_type, nmr_env_type, polar_env_type, vcd_env_type
   USE qs_mfp,                          ONLY: mfp_aat,&
                                              mfp_build_operator_gauge_dependent,&
                                              mfp_build_operator_gauge_independent,&
                                              mfp_response
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_p_env_methods,                ONLY: p_env_create,&
                                              p_env_psi0_changed
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vcd,                          ONLY: aat_dV,&
                                              apt_dV,&
                                              prepare_per_atom_vcd,&
                                              vcd_build_op_dV,&
                                              vcd_response_dV
   USE qs_vcd_utils,                    ONLY: vcd_env_cleanup,&
                                              vcd_env_init,&
                                              vcd_print
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: linres_calculation, linres_calculation_low

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_module'

CONTAINS
! *****************************************************************************
!> \brief Calculates the derivatives of the MO coefficients dC/dV^lambda_beta
!>         wrt to nuclear velocities. The derivative is indexed by `beta`, the
!>         electric dipole operator by `alpha`.
!>        Calculates the APT and AAT in velocity form
!>               P^lambda_alpha,beta = d< mu_alpha >/dV^lambda_beta
!>               M^lambda_alpha,beta = d< m_alpha >/dV^lambda_beta
!> \param qs_env ...
!> \param p_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE vcd_linres(qs_env, p_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env

      INTEGER                                            :: beta, i, latom
      LOGICAL                                            :: mfp_is_done, mfp_repeat
      TYPE(vcd_env_type)                                 :: vcd_env

      CALL cite_reference(Ditler2022)

      ! We need the position perturbation for the velocity perturbation operator
      CALL vcd_env_init(vcd_env, qs_env)

      mfp_repeat = vcd_env%distributed_origin
      mfp_is_done = .FALSE.

      qs_env%linres_control%linres_restart = .TRUE.

      ! Iterate over the list of atoms for which we want to calculate the APTs/AATs
      !  default is all atoms.
      DO latom = 1, SIZE(vcd_env%dcdr_env%list_of_atoms)
         vcd_env%dcdr_env%lambda = vcd_env%dcdr_env%list_of_atoms(latom)

         CALL prepare_per_atom(vcd_env%dcdr_env, qs_env)
         CALL prepare_per_atom_vcd(vcd_env, qs_env)

         DO beta = 1, 3                   ! in every direction

            vcd_env%dcdr_env%beta = beta
            vcd_env%dcdr_env%deltaR(vcd_env%dcdr_env%beta, vcd_env%dcdr_env%lambda) = 1._dp

            ! Since we do the heavy lifting anyways, we might also calculate the length form APTs here
            CALL dcdr_build_op_dR(vcd_env%dcdr_env, qs_env)
            CALL dcdr_response_dR(vcd_env%dcdr_env, p_env, qs_env)
            CALL apt_dR(qs_env, vcd_env%dcdr_env)

            ! And with the position perturbation ready, we can calculate the NVP
            CALL vcd_build_op_dV(vcd_env, qs_env)
            CALL vcd_response_dV(vcd_env, p_env, qs_env)

            CALL apt_dV(vcd_env, qs_env)
            CALL aat_dV(vcd_env, qs_env)

            IF (vcd_env%do_mfp) THEN
               ! Since we came so far, we might as well calculate the MFP AATs
               ! If we use a distributed origin we need to compute the MFP response again for each
               !   atom, because the reference point changes.
               IF (.NOT. mfp_is_done .OR. mfp_repeat) THEN
                  DO i = 1, 3
                     IF (vcd_env%origin_dependent_op_mfp) THEN
                        CPWARN("Using the origin dependent MFP operator")
                        CALL mfp_build_operator_gauge_dependent(vcd_env, qs_env, i)
                     ELSE
                        CALL mfp_build_operator_gauge_independent(vcd_env, qs_env, i)
                     END IF
                     CALL mfp_response(vcd_env, p_env, qs_env, i)
                  END DO
                  mfp_is_done = .TRUE.
               END IF

               CALL mfp_aat(vcd_env, qs_env)
            END IF
         END DO ! beta

         vcd_env%dcdr_env%apt_total_dcdr(:, :, vcd_env%dcdr_env%lambda) = &
            vcd_env%dcdr_env%apt_el_dcdr(:, :, vcd_env%dcdr_env%lambda) &
            + vcd_env%dcdr_env%apt_nuc_dcdr(:, :, vcd_env%dcdr_env%lambda)

         vcd_env%apt_total_nvpt(:, :, vcd_env%dcdr_env%lambda) = &
            vcd_env%apt_el_nvpt(:, :, vcd_env%dcdr_env%lambda) + vcd_env%apt_nuc_nvpt(:, :, vcd_env%dcdr_env%lambda)

         IF (vcd_env%do_mfp) &
            vcd_env%aat_atom_mfp(:, :, vcd_env%dcdr_env%lambda) = vcd_env%aat_atom_mfp(:, :, vcd_env%dcdr_env%lambda)*4._dp

      END DO !lambda

      CALL vcd_print(vcd_env, qs_env)
      CALL vcd_env_cleanup(qs_env, vcd_env)

   END SUBROUTINE vcd_linres

! **************************************************************************************************
!> \brief Calculates the derivatives of the MO coefficients dC/dR^lambda_beta
!>         wrt to nuclear coordinates. The derivative is index by `beta`, the
!>         electric dipole operator by `alpha`.
!>        Also calculates the APT
!>               P^lambda_alpha,beta = d< mu_alpha >/dR^lambda_beta
!>        and calculates the sum rules for the APT elements.
!> \param qs_env ...
!> \param p_env ...
! **************************************************************************************************
   SUBROUTINE dcdr_linres(qs_env, p_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env

      INTEGER                                            :: beta, latom
      TYPE(dcdr_env_type)                                :: dcdr_env
      TYPE(polar_env_type), POINTER                      :: polar_env

      CALL cite_reference(Ditler2021)
      CALL dcdr_env_init(dcdr_env, qs_env)

      IF (.NOT. dcdr_env%z_matrix_method) THEN

         DO latom = 1, SIZE(dcdr_env%list_of_atoms)
            dcdr_env%lambda = dcdr_env%list_of_atoms(latom)
            CALL prepare_per_atom(dcdr_env, qs_env)

            DO beta = 1, 3                   ! in every direction
               dcdr_env%beta = beta
               dcdr_env%deltaR(dcdr_env%beta, dcdr_env%lambda) = 1._dp

               CALL dcdr_build_op_dR(dcdr_env, qs_env)
               CALL dcdr_response_dR(dcdr_env, p_env, qs_env)

               IF (.NOT. dcdr_env%localized_psi0) THEN
                  CALL apt_dR(qs_env, dcdr_env)
               ELSE IF (dcdr_env%localized_psi0) THEN
                  CALL apt_dR_localization(qs_env, dcdr_env)
               END IF

            END DO !beta

            dcdr_env%apt_total_dcdr(:, :, dcdr_env%lambda) = &
               dcdr_env%apt_el_dcdr(:, :, dcdr_env%lambda) + dcdr_env%apt_nuc_dcdr(:, :, dcdr_env%lambda)
         END DO !lambda

      ELSE

         CALL polar_env_init(qs_env)
         CALL get_qs_env(qs_env=qs_env, polar_env=polar_env)
         CALL get_polar_env(polar_env=polar_env)

         IF (.NOT. dcdr_env%localized_psi0) THEN
            CALL polar_operators_local(qs_env)
         ELSE
            CALL polar_operators_local_wannier(qs_env, dcdr_env)
         END IF

         polar_env%do_periodic = .FALSE.
         CALL polar_response(p_env, qs_env)

         DO latom = 1, SIZE(dcdr_env%list_of_atoms)
            dcdr_env%lambda = dcdr_env%list_of_atoms(latom)
            CALL prepare_per_atom(dcdr_env, qs_env)

            DO beta = 1, 3                   ! in every direction
               dcdr_env%beta = beta
               dcdr_env%deltaR(dcdr_env%beta, dcdr_env%lambda) = 1._dp

               CALL dcdr_build_op_dR(dcdr_env, qs_env)
               IF (.NOT. dcdr_env%localized_psi0) THEN
                  CALL apt_dR(qs_env, dcdr_env)
               ELSE
                  CALL apt_dR_localization(qs_env, dcdr_env)
               END IF
            END DO !beta

            dcdr_env%apt_total_dcdr(:, :, dcdr_env%lambda) = &
               dcdr_env%apt_el_dcdr(:, :, dcdr_env%lambda) + dcdr_env%apt_nuc_dcdr(:, :, dcdr_env%lambda)
         END DO !lambda

      END IF

      CALL dcdr_print(dcdr_env, qs_env)
      CALL dcdr_env_cleanup(qs_env, dcdr_env)
   END SUBROUTINE dcdr_linres

! **************************************************************************************************
!> \brief Driver for the linear response calculatios
!> \param force_env ...
!> \par History
!>      06.2005 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE linres_calculation(force_env)

      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_calculation'

      INTEGER                                            :: handle
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routineN, handle)

      NULLIFY (qs_env)

      CPASSERT(ASSOCIATED(force_env))
      CPASSERT(force_env%ref_count > 0)

      SELECT CASE (force_env%in_use)
      CASE (use_qs_force)
         CALL force_env_get(force_env, qs_env=qs_env)
      CASE (use_qmmm)
         qs_env => force_env%qmmm_env%qs_env
      CASE DEFAULT
         CPABORT("Does not recognize this force_env")
      END SELECT

      qs_env%linres_run = .TRUE.

      CALL linres_calculation_low(qs_env)

      CALL timestop(handle)

   END SUBROUTINE linres_calculation

! **************************************************************************************************
!> \brief Linear response can be called as run type or as post scf calculation
!>      Initialize the perturbation environment
!>      Define which properties is to be calculated
!>      Start up the optimization of the response density and wfn
!> \param qs_env ...
!> \par History
!>      06.2005 created [MI]
!>      02.2013 added polarizability section [SL]
!> \author MI
! **************************************************************************************************
   SUBROUTINE linres_calculation_low(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_calculation_low'

      INTEGER                                            :: every_n_step, handle, iounit
      LOGICAL                                            :: dcdr_present, do_apt_fd, epr_present, &
                                                            issc_present, lr_calculation, &
                                                            nmr_present, polar_present, vcd_present
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(section_vals_type), POINTER                   :: lr_section, prop_section

      CALL timeset(routineN, handle)

      lr_calculation = .FALSE.
      nmr_present = .FALSE.
      epr_present = .FALSE.
      issc_present = .FALSE.
      polar_present = .FALSE.
      dcdr_present = .FALSE.
      do_apt_fd = .FALSE.

      NULLIFY (dft_control, linres_control, logger, prop_section, lr_section)

      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")
      CALL section_vals_get(lr_section, explicit=lr_calculation)

      CALL section_vals_val_get(lr_section, "DCDR%APT_FD", explicit=do_apt_fd)
      IF (do_apt_fd) THEN
         CALL timestop(handle)
         RETURN
      END IF

      logger => cp_get_default_logger()

      CALL section_vals_val_get(lr_section, "EVERY_N_STEP", i_val=every_n_step)

      IF (lr_calculation .AND. MODULO(qs_env%sim_step, every_n_step) == 0) THEN
         CALL linres_init(lr_section, p_env, qs_env)
         iounit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                       extension=".linresLog")
         CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, &
                         linres_control=linres_control)

         ! The type of perturbation has not been defined yet
         linres_control%property = lr_none

         ! We do NMR or EPR, then compute the current response
         prop_section => section_vals_get_subs_vals(lr_section, "NMR")
         CALL section_vals_get(prop_section, explicit=nmr_present)
         prop_section => section_vals_get_subs_vals(lr_section, "EPR")
         CALL section_vals_get(prop_section, explicit=epr_present)

         IF (nmr_present .OR. epr_present) THEN
            CALL nmr_epr_linres(linres_control, qs_env, p_env, dft_control, &
                                nmr_present, epr_present, iounit)
         END IF

         ! We do the indirect spin-spin coupling calculation
         prop_section => section_vals_get_subs_vals(lr_section, "SPINSPIN")
         CALL section_vals_get(prop_section, explicit=issc_present)

         IF (issc_present) THEN
            CALL issc_linres(linres_control, qs_env, p_env, dft_control)
         END IF

         ! We do the polarizability calculation
         prop_section => section_vals_get_subs_vals(lr_section, "POLAR")
         CALL section_vals_get(prop_section, explicit=polar_present)
         IF (polar_present) THEN
            CALL polar_linres(qs_env, p_env)
         END IF

         ! Nuclear Position Perturbation
         prop_section => section_vals_get_subs_vals(lr_section, "dcdr")
         CALL section_vals_get(prop_section, explicit=dcdr_present)

         IF (dcdr_present) THEN
            CALL dcdr_linres(qs_env, p_env)
         END IF

         ! VCD
         prop_section => section_vals_get_subs_vals(lr_section, "VCD")
         CALL section_vals_get(prop_section, explicit=vcd_present)

         IF (vcd_present) THEN
            CALL vcd_linres(qs_env, p_env)
         END IF

         ! Other possible LR calculations can be introduced here

         CALL p_env_release(p_env)

         IF (iounit > 0) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A,/,T25,A,/,T2,A,/)") &
               REPEAT("=", 79), &
               "ENDED LINRES CALCULATION", &
               REPEAT("=", 79)
         END IF
         CALL cp_print_key_finished_output(iounit, logger, lr_section, &
                                           "PRINT%PROGRAM_RUN_INFO")
      END IF

      CALL timestop(handle)

   END SUBROUTINE linres_calculation_low

! **************************************************************************************************
!> \brief Initialize some general settings like the p_env
!>      Localize the psi0 if required
!> \param lr_section ...
!> \param p_env ...
!> \param qs_env ...
!> \par History
!>      06.2005 created [MI]
!> \author MI
!> \note
!>      - The localization should probably be always for all the occupied states
! **************************************************************************************************
   SUBROUTINE linres_init(lr_section, p_env, qs_env)

      TYPE(section_vals_type), POINTER                   :: lr_section
      TYPE(qs_p_env_type), INTENT(OUT)                   :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: iounit, ispin
      LOGICAL                                            :: do_it
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: loc_section

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                    extension=".linresLog")
      NULLIFY (dft_control, linres_control, loc_section, rho, mos, matrix_ks, rho_ao)

      ALLOCATE (linres_control)
      CALL set_qs_env(qs_env=qs_env, linres_control=linres_control)
      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, matrix_ks=matrix_ks, mos=mos, rho=rho)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ! Localized Psi0 are required when the position operator has to be defined (nmr)
      loc_section => section_vals_get_subs_vals(lr_section, "LOCALIZE")
      CALL section_vals_val_get(loc_section, "_SECTION_PARAMETERS_", &
                                l_val=linres_control%localized_psi0)
      IF (linres_control%localized_psi0) THEN
         IF (iounit > 0) THEN
            WRITE (UNIT=iounit, FMT="(/,T3,A,A)") &
               "Localization of ground state orbitals", &
               " before starting linear response calculation"
         END IF

         CALL linres_localize(qs_env, linres_control, dft_control%nspins)

         DO ispin = 1, dft_control%nspins
            CALL calculate_density_matrix(mos(ispin), rho_ao(ispin)%matrix)
         END DO
         ! ** update qs_env%rho
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
      END IF

      CALL section_vals_val_get(lr_section, "RESTART", l_val=linres_control%linres_restart)
      CALL section_vals_val_get(lr_section, "MAX_ITER", i_val=linres_control%max_iter)
      CALL section_vals_val_get(lr_section, "EPS", r_val=linres_control%eps)
      CALL section_vals_val_get(lr_section, "EPS_FILTER", r_val=linres_control%eps_filter)
      CALL section_vals_val_get(lr_section, "RESTART_EVERY", i_val=linres_control%restart_every)
      CALL section_vals_val_get(lr_section, "PRECONDITIONER", i_val=linres_control%preconditioner_type)
      CALL section_vals_val_get(lr_section, "ENERGY_GAP", r_val=linres_control%energy_gap)

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A,/,T25,A,/,T2,A,/)") &
            REPEAT("=", 79), &
            "START LINRES CALCULATION", &
            REPEAT("=", 79)

         WRITE (UNIT=iounit, FMT="(T2,A)") &
            "LINRES| Properties to be calculated:"
         CALL section_vals_val_get(lr_section, "NMR%_SECTION_PARAMETERS_", l_val=do_it)
         IF (do_it) WRITE (UNIT=iounit, FMT="(T62,A)") "NMR Chemical Shift"
         CALL section_vals_val_get(lr_section, "EPR%_SECTION_PARAMETERS_", l_val=do_it)
         IF (do_it) WRITE (UNIT=iounit, FMT="(T68,A)") "EPR g Tensor"
         CALL section_vals_val_get(lr_section, "SPINSPIN%_SECTION_PARAMETERS_", l_val=do_it)
         IF (do_it) WRITE (UNIT=iounit, FMT="(T43,A)") "Indirect spin-spin coupling constants"
         CALL section_vals_val_get(lr_section, "POLAR%_SECTION_PARAMETERS_", l_val=do_it)
         IF (do_it) WRITE (UNIT=iounit, FMT="(T57,A)") "Electric Polarizability"

         IF (linres_control%localized_psi0) WRITE (UNIT=iounit, FMT="(T2,A,T65,A)") &
            "LINRES|", " LOCALIZED PSI0"

         WRITE (UNIT=iounit, FMT="(T2,A,T60,A)") &
            "LINRES| Optimization algorithm", " Conjugate Gradients"

         SELECT CASE (linres_control%preconditioner_type)
         CASE (ot_precond_none)
            WRITE (UNIT=iounit, FMT="(T2,A,T60,A)") &
               "LINRES| Preconditioner", "                NONE"
         CASE (ot_precond_full_single)
            WRITE (UNIT=iounit, FMT="(T2,A,T60,A)") &
               "LINRES| Preconditioner", "         FULL_SINGLE"
         CASE (ot_precond_full_kinetic)
            WRITE (UNIT=iounit, FMT="(T2,A,T60,A)") &
               "LINRES| Preconditioner", "        FULL_KINETIC"
         CASE (ot_precond_s_inverse)
            WRITE (UNIT=iounit, FMT="(T2,A,T60,A)") &
               "LINRES| Preconditioner", "      FULL_S_INVERSE"
         CASE (ot_precond_full_single_inverse)
            WRITE (UNIT=iounit, FMT="(T2,A,T60,A)") &
               "LINRES| Preconditioner", " FULL_SINGLE_INVERSE"
         CASE (ot_precond_full_all)
            WRITE (UNIT=iounit, FMT="(T2,A,T60,A)") &
               "LINRES| Preconditioner", "            FULL_ALL"
         CASE DEFAULT
            CPABORT("Preconditioner NYI")
         END SELECT

         WRITE (UNIT=iounit, FMT="(T2,A,T72,ES8.1)") &
            "LINRES| EPS", linres_control%eps
         WRITE (UNIT=iounit, FMT="(T2,A,T72,I8)") &
            "LINRES| MAX_ITER", linres_control%max_iter
      END IF

      !------------------!
      ! create the p_env !
      !------------------!
      CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE., linres_control=linres_control)

      ! update the m_epsilon matrix
      CALL p_env_psi0_changed(p_env, qs_env)

      p_env%new_preconditioner = .TRUE.
      CALL cp_print_key_finished_output(iounit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

   END SUBROUTINE linres_init

! **************************************************************************************************
!> \brief ...
!> \param linres_control ...
!> \param qs_env ...
!> \param p_env ...
!> \param dft_control ...
!> \param nmr_present ...
!> \param epr_present ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE nmr_epr_linres(linres_control, qs_env, p_env, dft_control, nmr_present, epr_present, iounit)

      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      LOGICAL                                            :: nmr_present, epr_present
      INTEGER                                            :: iounit

      INTEGER                                            :: iB
      LOGICAL                                            :: do_qmmm
      TYPE(current_env_type)                             :: current_env
      TYPE(epr_env_type)                                 :: epr_env
      TYPE(nmr_env_type)                                 :: nmr_env

      linres_control%property = lr_current

      CALL cite_reference(Weber2009)

      IF (.NOT. linres_control%localized_psi0) THEN
         CALL cp_abort(__LOCATION__, &
                       "Are you sure that you want to calculate the chemical "// &
                       "shift without localized psi0?")
         CALL linres_localize(qs_env, linres_control, &
                              dft_control%nspins, centers_only=.TRUE.)
      END IF
      IF (dft_control%nspins /= 2 .AND. epr_present) THEN
         CPABORT("LSD is needed to perform a g tensor calculation!")
      END IF
      !
      !Initialize the current environment
      do_qmmm = .FALSE.
      IF (qs_env%qmmm) do_qmmm = .TRUE.
      current_env%do_qmmm = do_qmmm
      !current_env%prop='nmr'
      CALL current_env_init(current_env, qs_env)
      CALL current_operators(current_env, qs_env)
      CALL current_response(current_env, p_env, qs_env)
      !
      IF (current_env%all_pert_op_done) THEN
         !Initialize the nmr environment
         IF (nmr_present) THEN
            CALL nmr_env_init(nmr_env, qs_env)
         END IF
         !
         !Initialize the epr environment
         IF (epr_present) THEN
            CALL epr_env_init(epr_env, qs_env)
            CALL epr_g_zke(epr_env, qs_env)
            CALL epr_nablavks(epr_env, qs_env)
         END IF
         !
         ! Build the rs_gauge if needed
         !CALL current_set_gauge(current_env,qs_env)
         !
         ! Loop over field direction
         DO iB = 1, 3
            !
            ! Build current response and succeptibility
            CALL current_build_current(current_env, qs_env, iB)
            CALL current_build_chi(current_env, qs_env, iB)
            !
            ! Compute NMR shift
            IF (nmr_present) THEN
               CALL nmr_shift(nmr_env, current_env, qs_env, iB)
            END IF
            !
            ! Compute EPR
            IF (epr_present) THEN
               CALL epr_ind_magnetic_field(epr_env, current_env, qs_env, iB)
               CALL epr_g_so(epr_env, current_env, qs_env, iB)
               CALL epr_g_soo(epr_env, current_env, qs_env, iB)
            END IF
         END DO
         !
         ! Finalized the nmr environment
         IF (nmr_present) THEN
            CALL nmr_shift_print(nmr_env, current_env, qs_env)
            CALL nmr_env_cleanup(nmr_env)
         END IF
         !
         ! Finalized the epr environment
         IF (epr_present) THEN
            CALL epr_g_print(epr_env, qs_env)
            CALL epr_env_cleanup(epr_env)
         END IF
         !
      ELSE
         IF (iounit > 0) THEN
            WRITE (iounit, "(T10,A,/T20,A,/)") &
               "CURRENT: Not all responses to perturbation operators could be calculated.", &
               " Hence: NO nmr and NO epr possible."
         END IF
      END IF
      ! Finalized the current environment
      CALL current_env_cleanup(current_env)

   END SUBROUTINE nmr_epr_linres

! **************************************************************************************************
!> \brief ...
!> \param linres_control ...
!> \param qs_env ...
!> \param p_env ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE issc_linres(linres_control, qs_env, p_env, dft_control)

      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      INTEGER                                            :: iatom
      LOGICAL                                            :: do_qmmm
      TYPE(current_env_type)                             :: current_env
      TYPE(issc_env_type)                                :: issc_env

      linres_control%property = lr_current
      IF (.NOT. linres_control%localized_psi0) THEN
         CALL cp_abort(__LOCATION__, &
                       "Are you sure that you want to calculate the chemical "// &
                       "shift without localized psi0?")
         CALL linres_localize(qs_env, linres_control, &
                              dft_control%nspins, centers_only=.TRUE.)
      END IF
      !
      !Initialize the current environment
      do_qmmm = .FALSE.
      IF (qs_env%qmmm) do_qmmm = .TRUE.
      current_env%do_qmmm = do_qmmm
      !current_env%prop='issc'
      !CALL current_env_init(current_env,qs_env)
      !CALL current_response(current_env,p_env,qs_env)
      !
      !Initialize the issc environment
      CALL issc_env_init(issc_env, qs_env)
      !
      ! Loop over atoms
      DO iatom = 1, issc_env%issc_natms
         CALL issc_operators(issc_env, qs_env, iatom)
         CALL issc_response(issc_env, p_env, qs_env)
         CALL issc_issc(issc_env, qs_env, iatom)
      END DO
      !
      ! Finalized the issc environment
      CALL issc_print(issc_env, qs_env)
      CALL issc_env_cleanup(issc_env)

   END SUBROUTINE issc_linres

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \par History
!>      06.2018 polar_env integrated into qs_env (MK)
! **************************************************************************************************
   SUBROUTINE polar_linres(qs_env, p_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env

      CALL polar_env_init(qs_env)
      CALL polar_operators(qs_env)
      CALL polar_response(p_env, qs_env)
      CALL polar_polar(qs_env)
      CALL polar_print(qs_env)

   END SUBROUTINE polar_linres

END MODULE qs_linres_module
