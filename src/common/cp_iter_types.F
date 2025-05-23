!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Collection of routines to handle the iteration info
! **************************************************************************************************
MODULE cp_iter_types
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! iteration_info
   PUBLIC :: cp_iteration_info_type, &
             cp_iteration_info_create, &
             cp_iteration_info_retain, &
             cp_iteration_info_release, &
             cp_iteration_info_copy_iter

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_iter_types'
   LOGICAL, PRIVATE, PARAMETER          :: debug_this_module = .FALSE.

   ! When adding a new iteration level PLEASE update the following list with the proper name!
   CHARACTER(LEN=default_path_length), PARAMETER, PUBLIC, DIMENSION(18) :: each_possible_labels = (/ &
                                                                           "__ROOT__          ", &
                                                                           "JUST_ENERGY       ", &
                                                                           "POWELL_OPT        ", &
                                                                           "QS_SCF            ", &
                                                                           "XAS_SCF           ", &
                                                                           "MD                ", &
                                                                           "PINT              ", &
                                                                           "METADYNAMICS      ", &
                                                                           "GEO_OPT           ", &
                                                                           "ROT_OPT           ", &
                                                                           "CELL_OPT          ", &
                                                                           "BAND              ", &
                                                                           "EP_LIN_SOLVER     ", &
                                                                           "SPLINE_FIND_COEFFS", &
                                                                           "REPLICA_EVAL      ", &
                                                                           "BSSE              ", &
                                                                           "SHELL_OPT         ", &
                                                                           "TDDFT_SCF         "/)

   CHARACTER(LEN=default_path_length), PARAMETER, PUBLIC, DIMENSION(18) ::  each_desc_labels = (/ &
                                                "Iteration level for __ROOT__ (fictitious iteration level)                      ", &
                                                "Iteration level for an ENERGY/ENERGY_FORCE calculation.                        ", &
                                                "Iteration level for POWELL based optimization steps.                           ", &
                                                "Iteration level for the SCF steps.                                             ", &
                                                "Iteration level for the X-Ray Absorption Spectroscopy (XAS) SCF steps.         ", &
                                                "Iteration level for the MD steps.                                              ", &
                                                "Iteration level for the Path integral md steps.                                ", &
                                                "Iteration level for the METADYNAMICS steps (number of hills added).            ", &
                                                "Iteration level for the Geometry optimization steps.                           ", &
                                                "Iteration level for the Rotational optimization steps in the Dimer calculation.", &
                                                "Iteration level for the Cell optimization steps.                               ", &
                                                "Iteration level for the BAND calculation steps                                 ", &
                                                "Iteration level for the Energy Perturbation (EP) linear solver                 ", &
                                                "Iteration level for the solution of the coefficients of the splines            ", &
                                                "Iteration level for the evaluation of the Replica Environment                  ", &
                                                "Iteration level for the Basis Set Superposition Error (BSSE) calculation       ", &
                                                "Iteration level for the Shell-Core distances optimization steps                ", &
                                                 "Iteration level for the Time-Dependent Density Functional Theory SCF steps.    "/)

! **************************************************************************************************
!> \brief contains the information about the current state of the program
!>      to be able to decide if output is necessary
!> \author fawzi
! **************************************************************************************************
   TYPE cp_iteration_info_type
      INTEGER                              :: ref_count = -1
      INTEGER                              :: print_level = -1, n_rlevel = -1
      INTEGER, DIMENSION(:), POINTER       :: iteration => NULL()
      LOGICAL, DIMENSION(:), POINTER       :: last_iter => NULL()
      CHARACTER(len=default_string_length) :: project_name = ""
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER           :: level_name => NULL()
   END TYPE cp_iteration_info_type

CONTAINS

! **************************************************************************************************
!> \brief creates an output info object
!> \param iteration_info the object to create
!> \param project_name name of the project, used to create the filenames
!> \author fawzi
! **************************************************************************************************
   PURE SUBROUTINE cp_iteration_info_create(iteration_info, project_name)
      TYPE(cp_iteration_info_type), POINTER              :: iteration_info
      CHARACTER(len=*), INTENT(in)                       :: project_name

      ALLOCATE (iteration_info)

      iteration_info%ref_count = 1
      iteration_info%print_level = 2
      iteration_info%n_rlevel = 1
      iteration_info%project_name = project_name
      ALLOCATE (iteration_info%iteration(iteration_info%n_rlevel))
      ALLOCATE (iteration_info%level_name(iteration_info%n_rlevel))
      ALLOCATE (iteration_info%last_iter(iteration_info%n_rlevel))
      iteration_info%iteration(iteration_info%n_rlevel) = 1
      iteration_info%level_name(iteration_info%n_rlevel) = "__ROOT__"
      iteration_info%last_iter(iteration_info%n_rlevel) = .FALSE.

   END SUBROUTINE cp_iteration_info_create

! **************************************************************************************************
!> \brief retains the iteration_info (see doc/ReferenceCounting.html)
!> \param iteration_info the iteration_info to retain
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE cp_iteration_info_retain(iteration_info)
      TYPE(cp_iteration_info_type), INTENT(INOUT)        :: iteration_info

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_iteration_info_retain', &
         routineP = moduleN//':'//routineN

      IF (iteration_info%ref_count <= 0) THEN
         CPABORT(routineP//" iteration_info%ref_counf<=0")
      END IF
      iteration_info%ref_count = iteration_info%ref_count + 1
   END SUBROUTINE cp_iteration_info_retain

! **************************************************************************************************
!> \brief releases the iteration_info (see doc/ReferenceCounting.html)
!> \param iteration_info the iteration_info to release
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE cp_iteration_info_release(iteration_info)
      TYPE(cp_iteration_info_type), POINTER              :: iteration_info

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_iteration_info_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(iteration_info)) THEN
         IF (iteration_info%ref_count <= 0) THEN
            CPABORT(routineP//" iteration_info%ref_counf<=0")
         END IF
         iteration_info%ref_count = iteration_info%ref_count - 1
         IF (iteration_info%ref_count == 0) THEN
            IF (ASSOCIATED(iteration_info%iteration)) THEN
               DEALLOCATE (iteration_info%iteration)
            END IF
            IF (ASSOCIATED(iteration_info%last_iter)) THEN
               DEALLOCATE (iteration_info%last_iter)
            END IF
            IF (ASSOCIATED(iteration_info%level_name)) THEN
               DEALLOCATE (iteration_info%level_name)
            END IF
            DEALLOCATE (iteration_info)
         END IF
      END IF
   END SUBROUTINE cp_iteration_info_release

! **************************************************************************************************
!> \brief Copies iterations info of an iteration info into another iteration info
!> \param iteration_info_in the iteration_info to be copied
!> \param iteration_info_out the iteration_info results of the copy
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE cp_iteration_info_copy_iter(iteration_info_in, iteration_info_out)
      TYPE(cp_iteration_info_type), INTENT(INOUT)        :: iteration_info_in, iteration_info_out

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_iteration_info_copy_iter', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (iteration_info_in%ref_count <= 0) THEN
         CPABORT(routineP//" iteration_info_in%ref_counf<=0")
      END IF

      iteration_info_out%n_rlevel = iteration_info_in%n_rlevel

      DEALLOCATE (iteration_info_out%iteration)
      i = SIZE(iteration_info_in%iteration)
      ALLOCATE (iteration_info_out%iteration(i))
      iteration_info_out%iteration = iteration_info_in%iteration

      DEALLOCATE (iteration_info_out%last_iter)
      i = SIZE(iteration_info_in%last_iter)
      ALLOCATE (iteration_info_out%last_iter(i))
      iteration_info_out%last_iter = iteration_info_in%last_iter

      DEALLOCATE (iteration_info_out%level_name)
      i = SIZE(iteration_info_in%level_name)
      ALLOCATE (iteration_info_out%level_name(i))
      iteration_info_out%level_name = iteration_info_in%level_name

   END SUBROUTINE cp_iteration_info_copy_iter

END MODULE cp_iter_types

