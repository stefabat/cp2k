!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Handling of the Wiener process currently employed in turn of the
!>      Langevin dynamics.
!> \par History
!>      none
!> \author Matthias Krack (05.07.2005)
! **************************************************************************************************
MODULE wiener_process

   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE md_environment_types,            ONLY: get_md_env,&
                                              md_environment_type,&
                                              need_per_atom_wiener_process
   USE message_passing,                 ONLY: mp_para_env_type
   USE metadynamics_types,              ONLY: meta_env_type
   USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                              next_rng_seed,&
                                              rng_record_length,&
                                              rng_stream_type,&
                                              rng_stream_type_from_record
   USE particle_list_types,             ONLY: particle_list_type
   USE simpar_types,                    ONLY: simpar_type
   USE string_utilities,                ONLY: compress
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters in this module
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wiener_process'

   ! Public subroutines
   PUBLIC :: create_wiener_process, create_wiener_process_cv

CONTAINS

! **************************************************************************************************
!> \brief Create a Wiener process for Langevin dynamics and initialize an
!>      independent random number generator for each atom in all force
!>      environment and all the subsystems/fragments therein.
!> \param md_env ...
!> \par History
!>      Creation (06.07.2005,MK)
! **************************************************************************************************
   SUBROUTINE create_wiener_process(md_env)

      TYPE(md_environment_type), POINTER                 :: md_env

      CHARACTER(LEN=40)                                  :: name
      INTEGER                                            :: iparticle, iparticle_kind, &
                                                            iparticle_local, nparticle, &
                                                            nparticle_kind, nparticle_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: seed
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(section_vals_type), POINTER                   :: force_env_section, subsys_section, &
                                                            work_section
      TYPE(simpar_type), POINTER                         :: simpar

      NULLIFY (work_section, force_env)
      CPASSERT(ASSOCIATED(md_env))

      CALL get_md_env(md_env=md_env, force_env=force_env, para_env=para_env, &
                      simpar=simpar)

      ![NB] shouldn't the calling process know if it's needed
      IF (need_per_atom_wiener_process(md_env)) THEN

         CALL force_env_get(force_env, force_env_section=force_env_section, &
                            subsys=subsys)

         subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")

         CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds, local_particles=local_particles, &
                            particles=particles)

         nparticle_kind = atomic_kinds%n_els
         nparticle = particles%n_els

         ! Allocate the (local) data structures for the Wiener process
         ALLOCATE (local_particles%local_particle_set(nparticle_kind))

         DO iparticle_kind = 1, nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            ALLOCATE (local_particles%local_particle_set(iparticle_kind)%rng(nparticle_local))
            DO iparticle_local = 1, nparticle_local
               ALLOCATE (local_particles%local_particle_set(iparticle_kind)%rng(iparticle_local)%stream)
            END DO
         END DO

         ! Each process generates all seeds. The seed generation should be
         ! quite fast and in this way a broadcast is avoided.
         ALLOCATE (seed(3, 2, nparticle))

         ! Load initial seed (not needed for a restart)
         seed(:, :, 1) = subsys%seed(:, :)

         DO iparticle = 2, nparticle
            seed(:, :, iparticle) = next_rng_seed(seed(:, :, iparticle - 1))
         END DO

         ! Update initial seed
         subsys%seed(:, :) = next_rng_seed(seed(:, :, nparticle))

         ! Create a random number stream (Wiener process) for each particle
         DO iparticle_kind = 1, nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            DO iparticle_local = 1, nparticle_local
               iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
               WRITE (UNIT=name, FMT="(A,I8)") "Wiener process for particle", iparticle
               CALL compress(name)
               local_particles%local_particle_set(iparticle_kind)%rng(iparticle_local)% &
                  stream = rng_stream_type(name=name, distribution_type=GAUSSIAN, &
                                           extended_precision=.TRUE., seed=seed(:, :, iparticle))
            END DO
         END DO

         DEALLOCATE (seed)

         ! Possibly restart Wiener process
         NULLIFY (work_section)
         work_section => section_vals_get_subs_vals(section_vals=subsys_section, &
                                                    subsection_name="RNG_INIT")
         CALL init_local_particle_set(distribution_1d=local_particles, &
                                      nparticle_kind=nparticle_kind, &
                                      work_section=work_section)
      END IF

   END SUBROUTINE create_wiener_process

! **************************************************************************************************
!> \brief Helper routine for create_wiener_process.
!> \param distribution_1d ...
!> \param nparticle_kind ...
!> \param work_section ...
!> \par History
!>      01.2014 moved from distribution_1d_types (Ole Schuett)
! **************************************************************************************************
   SUBROUTINE init_local_particle_set(distribution_1d, nparticle_kind, &
                                      work_section)

      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      INTEGER, INTENT(in)                                :: nparticle_kind
      TYPE(section_vals_type), POINTER                   :: work_section

      CHARACTER(LEN=rng_record_length)                   :: rng_record
      INTEGER                                            :: iparticle, iparticle_kind, &
                                                            iparticle_local, nparticle_local
      LOGICAL                                            :: explicit

! -------------------------------------------------------------------------

      CPASSERT(ASSOCIATED(distribution_1d))

      IF (ASSOCIATED(work_section)) THEN
         CALL section_vals_get(work_section, explicit=explicit)
         IF (explicit) THEN
            DO iparticle_kind = 1, nparticle_kind
               nparticle_local = distribution_1d%n_el(iparticle_kind)
               DO iparticle_local = 1, nparticle_local
                  iparticle = distribution_1d%list(iparticle_kind)%array(iparticle_local)
                  IF (iparticle == distribution_1d%list(iparticle_kind)%array(iparticle_local)) THEN
                     CALL section_vals_val_get(section_vals=work_section, &
                                               keyword_name="_DEFAULT_KEYWORD_", &
                                               i_rep_val=iparticle, &
                                               c_val=rng_record)
                     distribution_1d%local_particle_set(iparticle_kind)%rng(iparticle_local)% &
                        stream = rng_stream_type_from_record(rng_record)
                  END IF
               END DO
            END DO
         END IF
      END IF

   END SUBROUTINE init_local_particle_set

! **************************************************************************************************
!> \brief Create a Wiener process for Langevin dynamics used for
!>        metadynamics and initialize an
!>        independent random number generator for each COLVAR.
!> \param meta_env ...
!> \date   01.2009
!> \author Fabio Sterpone
!>
! **************************************************************************************************
   SUBROUTINE create_wiener_process_cv(meta_env)

      TYPE(meta_env_type), POINTER                       :: meta_env

      CHARACTER(LEN=40)                                  :: name
      INTEGER                                            :: i_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: seed
      REAL(KIND=dp), DIMENSION(3, 2)                     :: initial_seed

      IF (.NOT. ASSOCIATED(meta_env)) RETURN

      initial_seed = next_rng_seed()

      ! Each process generates all seeds. The seed generation should be
      ! quite fast and in this way a broadcast is avoided.

      ALLOCATE (seed(3, 2, meta_env%n_colvar))

      seed(:, :, 1) = initial_seed
      DO i_c = 2, meta_env%n_colvar
         seed(:, :, i_c) = next_rng_seed(seed(:, :, i_c - 1))
      END DO

      ! Update initial seed
      initial_seed = next_rng_seed(seed(:, :, meta_env%n_colvar))

      ! Create a random number stream (Wiener process) for each particle
      DO i_c = 1, meta_env%n_colvar
         WRITE (UNIT=name, FMT="(A,I8)") "Wiener process for COLVAR", i_c
         CALL compress(name)
         meta_env%rng(i_c) = rng_stream_type(name=name, distribution_type=GAUSSIAN, &
                                             extended_precision=.TRUE., seed=seed(:, :, i_c))
      END DO
      DEALLOCATE (seed)

   END SUBROUTINE create_wiener_process_cv

END MODULE wiener_process
