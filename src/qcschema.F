!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The module to read/write QCSchema JSON files for interfacing CP2K with other programs
!> \par History
!>      10.2022 created [SB]
!> \author SB
! **************************************************************************************************
MODULE qcschema

   ! USE json_module,                    ONLY: json_core, json_value
   USE cp2k_info,                      ONLY: cp2k_version
   USE cp_control_types,               ONLY: dft_control_type
   USE cp_log_handling,                ONLY: cp_get_default_logger, &
                                             cp_logger_get_default_io_unit, &
                                             cp_logger_type
   USE kinds,                          ONLY: dp, default_string_length, &
                                             default_path_length
   USE mp2_types,                      ONLY: mp2_type
   USE particle_types,                 ONLY: particle_type
   USE qs_force_types,                 ONLY: qs_force_type
   USE qs_mo_types,                    ONLY: mo_set_p_type
   USE qs_scf_types,                   ONLY: qs_scf_env_type
   USE qs_energy_types,                ONLY: qs_energy_type
   USE qs_environment_types,           ONLY: get_qs_env, &
                                             qs_environment_type
   USE hdf5,                           ONLY: h5open_f, h5close_f, &
                                             h5fcreate_f, h5fopen_f, h5fclose_f, &
                                             h5screate_f, h5sclose_f, H5S_SCALAR_F, &
                                             h5gcreate_f, &
                                             h5acreate_f, h5awrite_f, h5aclose_f, &
                                             h5tset_cset_f, h5tcopy_f, h5tset_strpad_f, &
                                             h5tcreate_f, h5tenum_insert_f, h5tclose_f, &
                                             H5T_STRING, H5T_CSET_UTF8_F, H5T_ENUM_F, &
                                             H5T_STR_NULLPAD_F, H5T_NATIVE_INTEGER, &
                                             HID_T, H5F_ACC_TRUNC_F, SIZE_T
   USE iso_c_binding,                  ONLY: c_ptr, C_LOC, c_int
   USE h5lt,                           ONLY: h5ltset_attribute_string_f, &
                                             h5ltset_attribute_double_f, &
                                             h5ltset_attribute_int_f

#include "./base/base_uses.f90"


   IMPLICIT NONE

   PRIVATE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'qcschema'

   PUBLIC :: qcschema_type
   PUBLIC :: qcschema_env_create, qcschema_to_hdf5

! **************************************************************************************************
!> \brief A derived type to store the program information that generated the QCSchema file.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#provenance
! **************************************************************************************************
   TYPE qcschema_provenance
      CHARACTER(LEN=default_string_length) :: creator ! The name of the creator of this object
      CHARACTER(LEN=default_string_length) :: version ! The version of the creator of this object
      CHARACTER(LEN=default_string_length) :: routine ! The routine that was used to create this object
   END TYPE qcschema_provenance


! **************************************************************************************************
!> \brief A derived type to store the topological information of the physical system.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#topology
! **************************************************************************************************
   TYPE qcschema_topology
      ! The list of atom symbols in this topology
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE  :: symbols
      ! The XYZ coordinates (in Bohr) of the atoms as a flat array (row major order)
      REAL(KIND=dp),  DIMENSION(:), ALLOCATABLE    :: geometry
      CHARACTER(LEN=default_string_length)         :: schema_name    ! The name of this schema
      INTEGER                                      :: schema_version ! The version of this specific schema
      INTEGER                                      :: molecular_charge
      INTEGER                                      :: molecular_multiplicity
      ! Whether translation of the geometry is allowed (.False.) or not (.True.)
      LOGICAL                                      :: fix_com
      ! A list indicating whether each atom is real (.True.) or a ghost (.False.).
      ! Its length must match that of the `symbols` list
      LOGICAL, DIMENSION(:), ALLOCATABLE           :: real
      ! A list indicating the bonds between the atoms in the molecule. Each item of this list must be
      ! a tuple of three integers, indicating the first atom index in the bond, the second atom index,
      ! and finally the order of the bond
      INTEGER, DIMENSION(:,:), ALLOCATABLE         :: connectivity
      ! Whether rotation of the geometry is allowed (`False`) or not (`True`)
      LOGICAL                                      :: fix_orientation
      ! A list of user-provided information for each atom. Its length must match that of the `symbols` list
      CHARACTER(LEN=8), DIMENSION(:), ALLOCATABLE  :: atom_labels
      ! The list of multiplicities associated with each fragment
      INTEGER, DIMENSION(:), ALLOCATABLE           :: fragment_multiplicities
      ! The maximal point group symmetry at which the `geometry` should be treated
      CHARACTER(LEN=default_string_length)         :: fix_symmetry
      ! The list of charges associated with each fragment
      INTEGER, DIMENSION(:), ALLOCATABLE           :: fragment_charges
      ! The mass numbers of all atoms. If it is an unknown isotope, the value should be -1.
      ! Its length must match that of the `symbols` list
      INTEGER, DIMENSION(:), ALLOCATABLE           :: mass_numbers
      ! The (user-given) name of the molecule
      CHARACTER(LEN=default_string_length)         :: name
      ! The masses (in au) of all atoms. Canonical weights are assumed if this is not given explicitly
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: masses
      ! Any additional (user-provided) comment
      CHARACTER(LEN=default_string_length)         :: comment
      ! An instance of :class:`QCProvenance`
      TYPE(qcschema_provenance)                    :: provenance
      ! The list of fragments. Each item of this list must be a tuple of integers with variable
      ! length (greater than 1). The first number indicates the fragment index, all following numbers
      ! refer to the (0-indexed) atom indices that constitute this fragment
      INTEGER, DIMENSION(:), ALLOCATABLE           :: fragments
      ! The atomic numbers of all atoms, indicating their nuclear charge. Its length must match that
      ! of the `symbols` list
      INTEGER, DIMENSION(:), ALLOCATABLE           :: atomic_numbers
   END TYPE qcschema_topology


! **************************************************************************************************
!> \brief A derived type to store the information of a single electron shell in a basis set.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L43
! **************************************************************************************************
   TYPE qcschema_electron_shell
      ! The angular momenta of this electron shell as a list of integers
      INTEGER, DIMENSION(:), ALLOCATABLE           :: angular_momentum
      ! The type of this shell: spherical or cartesian
      CHARACTER(LEN=9)                             :: harmonic_type
      ! The exponents of this contracted shell. The official spec stores these values as strings
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: exponents
      ! The general contraction coefficients of this contracted shell
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: coefficients
   END TYPE qcschema_electron_shell


! **************************************************************************************************
!> \brief A derived type to store the information of an ECP in a basis set.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L90
! **************************************************************************************************
   TYPE qcschema_ecp
      ! The type of this potential
      CHARACTER(LEN=default_string_length)         :: ecp_type
      ! The angular momenta of this potential as a list of integers
      INTEGER, DIMENSION(:), ALLOCATABLE           :: angular_momentum
      ! The exponents of the r terms
      INTEGER, DIMENSION(:), ALLOCATABLE           :: r_exponents
      ! The exponents of the Gaussian terms
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: gaussian_exponents
      ! The general contraction coefficients of this potential
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: coefficients
   END TYPE qcschema_ecp


! **************************************************************************************************
!> \brief A derived type to store the information of a single atom/center in the basis.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L146
! **************************************************************************************************
   TYPE qcschema_center_basis
      ! The list of electronic shells for this element
      TYPE(qcschema_electron_shell), DIMENSION(:), ALLOCATABLE &
      :: electron_shells
      ! The list of effective core potentials for this element
      TYPE(qcschema_ecp), DIMENSION(:), ALLOCATABLE &
      :: ecp_potentials
      ! The number of electrons replaced by an ECP
      INTEGER                                      :: ecp_electrons
   END TYPE qcschema_center_basis


! **************************************************************************************************
!> \brief A derived type to store the information of the basis set used in the calculation.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_basis.html#basis-set-schema
! **************************************************************************************************
   TYPE qcschema_basis_set
      ! The name of the basis set
      CHARACTER(LEN=default_string_length)         :: name
      ! A dictionary mapping the keys provided by `atom_map` to their basis center data
      TYPE(qcschema_center_basis), DIMENSION(:), ALLOCATABLE &
      :: center_data
      ! The list of atomic kinds, indicating the keys used to store the basis in `center_data`
      ! Not clear if this will be of the length of the basis set size, or rather just one
      ! entry for atomic kind. E.g. only one entry for hydrogen even though there might be
      ! many hydrogen atoms in the molecule. If this is the case, then we really need a
      ! hash table for `center_data`
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE  :: atom_map
      ! The version of this specific schema
      INTEGER                                      :: schema_version
      ! The name of this schema. This value is expected to be `qcschema_basis`
      CHARACTER(LEN=default_string_length)         :: schema_name
      ! A description of this basis set
      CHARACTER(LEN=default_string_length)         :: description
   END TYPE qcschema_basis_set


! **************************************************************************************************
!> \brief A derived type to store any additional computed wavefunction properties.
!>        Matrix quantities are stored as flat, column-major arrays.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_wf.html#wavefunction-schema
! **************************************************************************************************
   TYPE qcschema_wavefunction

      ! The basis set used during the computation
      TYPE(qcschema_basis_set)                     :: basis_set

      ! These are supposedly the "names" of the quantities even
      ! though I don't understand why do we need them
      CHARACTER(LEN=default_string_length)         :: orbitals_a
      CHARACTER(LEN=default_string_length)         :: orbitals_b
      CHARACTER(LEN=default_string_length)         :: density_a
      CHARACTER(LEN=default_string_length)         :: density_b
      CHARACTER(LEN=default_string_length)         :: density_mo_a
      CHARACTER(LEN=default_string_length)         :: density_mo_b
      CHARACTER(LEN=default_string_length)         :: fock_a
      CHARACTER(LEN=default_string_length)         :: fock_b
      CHARACTER(LEN=default_string_length)         :: fock_mo_a
      CHARACTER(LEN=default_string_length)         :: fock_mo_b
      CHARACTER(LEN=default_string_length)         :: eigenvalues_a
      CHARACTER(LEN=default_string_length)         :: eigenvalues_b
      CHARACTER(LEN=default_string_length)         :: occupations_a
      CHARACTER(LEN=default_string_length)         :: occupations_b
      CHARACTER(LEN=default_string_length)         :: eri
      CHARACTER(LEN=default_string_length)         :: eri_mo_aa
      CHARACTER(LEN=default_string_length)         :: eri_mo_ab
      CHARACTER(LEN=default_string_length)         :: eri_mo_ba
      CHARACTER(LEN=default_string_length)         :: eri_mo_bb

      ! SCF quantities in AO or MO basis
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_orbitals_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_orbitals_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_density_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_density_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_density_mo_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_density_mo_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_fock_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_fock_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_fock_mo_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_fock_mo_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_coulomb_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_coulomb_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_exchange_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_exchange_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eigenvalues_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eigenvalues_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_occupations_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_occupations_b

      ! ERI should be arrays of doubles, but in Qiskit are simple strings?
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eri
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eri_mo_aa
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eri_mo_ab
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eri_mo_ba
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eri_mo_bb

      ! Quantities with localized orbitals. All `nmo` orbitals are included,
      ! even if only a subset were localized
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: localized_orbitals_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: localized_orbitals_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: localized_fock_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: localized_fock_b

      ! Nne-electron Hamiltonians in AO basis
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: h_core_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: h_core_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: h_effective_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: h_effective_b

      ! Whether the computation used restricted spin orbitals
      LOGICAL :: restricted

   END TYPE qcschema_wavefunction


! **************************************************************************************************
!> \brief A derived type to store the computed properties of the original calculation.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_props.html#properties-schema
! **************************************************************************************************
   TYPE qcschema_properties

      REAL(KIND=dp) :: return_energy = 0.0_dp

      INTEGER :: calcinfo_nbasis = 0
      INTEGER :: calcinfo_nmo    = 0
      INTEGER :: calcinfo_nalpha = 0 ! # of alpha electrons
      INTEGER :: calcinfo_nbeta  = 0 ! # of beta electrons
      INTEGER :: calcinfo_natom  = 0

      ! SCF results
      INTEGER :: scf_iterations = 0
      REAL(KIND=dp) :: scf_one_electron_energy          = 0.0_dp
      REAL(KIND=dp) :: scf_two_electron_energy          = 0.0_dp
      REAL(KIND=dp) :: nuclear_repulsion_energy         = 0.0_dp
      ! REAL(KIND=dp) :: scf_vv10_energy                  = 0.0_dp
      REAL(KIND=dp) :: scf_xc_energy                    = 0.0_dp
      REAL(KIND=dp) :: scf_dispersion_correction_energy = 0.0_dp
      REAL(KIND=dp) :: scf_total_energy                 = 0.0_dp
      ! the dipole moment is calculated on the fly and not stored,
      ! so I just don't care about it
      ! REAL(KIND=dp), DIMENSION(3) :: scf_dipole_moment  = 0.0_dp

      ! MP2 results
      REAL(KIND=dp) :: mp2_same_spin_correlation_energy     = 0.0_dp
      REAL(KIND=dp) :: mp2_opposite_spin_correlation_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_singles_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_doubles_energy = 0.0_dp
      ! these are the only two that are saved
      REAL(KIND=dp) :: mp2_correlation_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_total_energy       = 0.0_dp
      ! the dipole moment is calculated on the fly and not stored,
      ! so I just don't care about it
      ! REAL(KIND=dp), DIMENSION(3) :: mp2_dipole_moment = 0.0_dp

      ! do we really need CC stuff?
      ! REAL(KIND=dp) :: ccsd_same_spin_correlation_energy
      ! REAL(KIND=dp) :: ccsd_opposite_spin_correlation_energy
      ! REAL(KIND=dp) :: ccsd_singles_energy
      ! REAL(KIND=dp) :: ccsd_doubles_energy
      ! REAL(KIND=dp) :: ccsd_correlation_energy
      ! REAL(KIND=dp) :: ccsd_total_energy
      ! REAL(KIND=dp) :: ccsd_prt_pr_correlation_energy
      ! REAL(KIND=dp) :: ccsd_prt_pr_total_energy
      ! REAL(KIND=dp) :: ccsdt_correlation_energy
      ! REAL(KIND=dp) :: ccsdt_total_energy
      ! REAL(KIND=dp) :: ccsdtq_correlation_energy
      ! REAL(KIND=dp) :: ccsdtq_total_energy
      ! REAL(KIND=dp), DIMENSION(3) :: ccsd_dipole_moment
      ! REAL(KIND=dp), DIMENSION(3) :: ccsd_prt_pr_dipole_moment
      ! REAL(KIND=dp), DIMENSION(3) :: ccsdt_dipole_moment
      ! REAL(KIND=dp), DIMENSION(3) :: ccsdtq_dipole_moment
      ! INTEGER :: ccsd_iterations
      ! INTEGER :: ccsdt_iterations
      ! INTEGER :: ccsdtq_iterations

   END TYPE qcschema_properties


! **************************************************************************************************
!> \brief The full QCSchema output type.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#output-components
! **************************************************************************************************
   TYPE qcschema_type
      TYPE(qcschema_provenance)   :: provenance
      TYPE(qcschema_properties)   :: properties
      TYPE(qcschema_wavefunction) :: wavefunction
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: return_result
      LOGICAL :: success
      CHARACTER(LEN=default_string_length) :: driver
   END TYPE qcschema_type


   CONTAINS


! **************************************************************************************************
!> \brief Create and initialize a qcschema environment
!> \param qcschema_env the qcschema environment to populate
!> \param qs_env the qs environment with all the info of the computation
! **************************************************************************************************
   SUBROUTINE qcschema_env_create(qcschema_env, qs_env)
      TYPE(qcschema_type)                                :: qcschema_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp2_type), POINTER                            :: mp2_env

      INTEGER :: handle, output_unit, nspins, nalpha, nbeta, natom, nbasis, nmo
      REAL(KIND=dp) :: one_el_en, two_el_en, dispersion
      CHARACTER(LEN=*), PARAMETER :: routineN = 'qcschema_env_initialize'

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      ! Collect environment info
      IF (ASSOCIATED(qs_env)) THEN
         CALL get_qs_env(qs_env, dft_control=dft_control, energy=energy, &
                         force=force, &
                         mos=mos, particle_set=particle_set, scf_env=scf_env, &
                         mp2_env=mp2_env)
      ELSE
         ! print something and quit
         RETURN
      END IF

      ! driver
      ! TODO: maybe separate this into input component and
      ! import iot here from there? This would follow the official
      ! guidelines of QCSchema...
      IF (ASSOCIATED(force)) THEN
         qcschema_env%driver = 'gradient'
      ELSE
         qcschema_env%driver = 'energy'
      END IF

      ! success
      ! TODO: how to check if the calculation was succesful?
      qcschema_env%success = .true.

      ! return result
      IF (ALLOCATED(qcschema_env%return_result)) THEN
         DEALLOCATE (qcschema_env%return_result)
      END IF
      IF (qcschema_env%success) THEN
         IF (qcschema_env%driver == 'energy') THEN
            ALLOCATE (qcschema_env%return_result(1))
            qcschema_env%return_result(1) = energy%total
         ELSE
            ALLOCATE (qcschema_env%return_result(3*SIZE(particle_set)))
            ! TODO: populate with forces!!
            qcschema_env%return_result = 0.0_dp
         END IF
      ELSE
         ! TODO: error handling and abort?
      END IF

      ! QCSchema provenanve
      qcschema_env%provenance%creator = 'CP2K'
      qcschema_env%provenance%version = cp2k_version
      qcschema_env%provenance%routine = routineN

      ! basis functions info
      nbasis = mos(1)%mo_set%nao
      nmo    = mos(1)%mo_set%nmo
      natom  = SIZE(particle_set)
      ! restricted or unrestricted calculation?
      nspins = SIZE(mos)
      IF (nspins == 1) THEN
         nalpha = mos(1)%mo_set%nelectron/2
         nbeta  = nalpha
      ELSE
         nalpha = mos(1)%mo_set%nelectron
         nbeta  = mos(2)%mo_set%nelectron
      END IF

      ! *** QCSchema properties *** !

      ! general info
      qcschema_env%properties%calcinfo_nbasis = nbasis
      qcschema_env%properties%calcinfo_nmo    = nmo
      qcschema_env%properties%calcinfo_natom  = natom
      qcschema_env%properties%calcinfo_nalpha = nalpha
      qcschema_env%properties%calcinfo_nbeta  = nbeta
      qcschema_env%properties%return_energy   = energy%total

      ! SCF info
      qcschema_env%properties%scf_iterations = scf_env%iter_count
      ! one-electron energy is the sum of all core terms
      one_el_en = energy%core_overlap + energy%core_self + energy%core
      qcschema_env%properties%scf_two_electron_energy = one_el_en
      ! two-electron energy is the sum of hartree and exact exchange (if there)
      two_el_en = energy%hartree + energy%ex + energy%hartree_1c
      qcschema_env%properties%scf_one_electron_energy = two_el_en
      ! xc energy
      qcschema_env%properties%scf_xc_energy = &
      energy%exc + energy%exc_aux_fit + energy%exc1 + energy%exc1_aux_fit
      ! dispersion energy
      dispersion = energy%dispersion + energy%gcp
      qcschema_env%properties%scf_dispersion_correction_energy = dispersion
      ! qcschema_env%properties%scf_vv10_energy
      ! remove the mp2 part, in case it's non-zero
      qcschema_env%properties%scf_total_energy = energy%total

      ! Some methods of CP2K are not supported by QCSchema, let's warn the user
      IF (dft_control%smear) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' WARNING: smearing not supported in QCSchema'
      END IF
      IF (dft_control%dft_plus_u) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' WARNING: DFT+U not supported in QCSchema'
      END IF
      IF (dft_control%do_sccs) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' WARNING: SCCS not supported in QCSchema'
      END IF
      IF (qs_env%qmmm) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' WARNING: QM/MM not supported in QCSchema'
      END IF
      IF (dft_control%qs_control%mulliken_restraint) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' WARNING: Mulliken restrains not supported in QCSchema'
      END IF
      IF (dft_control%qs_control%semi_empirical) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' WARNING: semi_empirical methods not supported in QCSchema'
      END IF
      IF (dft_control%qs_control%dftb) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' WARNING: DFTB not supported in QCSchema'
      END IF
      IF (dft_control%qs_control%xtb) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' WARNING: xTB not supported in QCSchema'
      END IF

      ! MP2 info
      IF (ASSOCIATED(qs_env%mp2_env)) THEN
         ! this info is actually computed, but not stored!
         ! qcschema_env%properties%mp2_singles_energy       ! don't know what these are
         ! qcschema_env%properties%mp2_doubles_energy       ! MP2 should be only this?
         ! qcschema_env%properties%mp2_same_spin_correlation_energy
         ! qcschema_env%properties%mp2_opposite_spin_correlation_energy

         qcschema_env%properties%mp2_correlation_energy = energy%mp2
         qcschema_env%properties%mp2_total_energy = energy%total

         ! update the scf energy
         qcschema_env%properties%scf_total_energy = energy%total - energy%mp2
      END IF

      CALL timestop(handle)

   END SUBROUTINE qcschema_env_create


   SUBROUTINE qcschema_to_hdf5(qcschema_env)!, qs_env)
      TYPE(qcschema_type)                    :: qcschema_env
      ! TYPE(qs_environment_type), POINTER     :: qs_env
      TYPE(cp_logger_type), POINTER          :: logger
      INTEGER                                :: error, output_unit
      INTEGER(KIND=8)                        :: n_forces
      INTEGER(KIND=hid_t)                    :: file_id, space_id, attr_id, group_id, str_type, bool_type
      INTEGER(KIND=size_t)                   :: enum_size = 1 ! 8-bit integers in enum bool_type
      CHARACTER(LEN=default_path_length)     :: filename

      TYPE(c_ptr), TARGET                    :: vlen_str
      TYPE(c_ptr)                            :: buffer
      CHARACTER(LEN=17), TARGET              :: schema_name = 'qcschema'
      INTEGER, TARGET                        :: schema_version = 3
      INTEGER, TARGET                        :: success

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      ! initialize HDF5 Fortran API
      CALL h5open_f(error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: initialization of HDF5 interface failed'
         RETURN
      ENDIF

      filename = '/home/stefano/Downloads/test.hdf5'
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      ! create a variable-length string type
      CALL h5tcopy_f(H5T_STRING, str_type, error)
      CAll h5tset_cset_f(str_type, H5T_CSET_UTF8_F, error)
      CALL h5tset_strpad_f(str_type, H5T_STR_NULLPAD_F, error)

      ! driver
      CALL h5ltset_attribute_string_f(file_id, '/', 'driver', qcschema_env%driver, error)

      ! return result
      n_forces = SIZE(qcschema_env%return_result)
      CALL h5ltset_attribute_double_f(file_id, '/', 'return_result', &
                                      qcschema_env%return_result, n_forces, error)

      ! schema name
      ! CALL h5ltset_attribute_string_f(file_id, '/', 'schema_name', 'qc_schema_output', error)
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      vlen_str = C_LOC(schema_name)
      buffer = C_LOC(vlen_str)
      CALL h5acreate_f(file_id, 'schema_name', str_type, space_id, attr_id, error)
      CALL h5awrite_f(attr_id, str_type, buffer, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(space_id, error)

      ! schema version
      ! CALL h5ltset_attribute_int_f(file_id, '/', 'schema_version', 3, error)
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      buffer = C_LOC(schema_version)
      CALL h5acreate_f(file_id, 'schema_version', H5T_NATIVE_INTEGER, space_id, attr_id, error)
      CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, buffer, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(space_id, error)

      ! success
      ! CALL h5tenum_create_f(H5T_ENUM_F, bool_type, error)
      CALL h5tcreate_f(H5T_ENUM_F, enum_size, bool_type, error)
      CALL h5tenum_insert_f(bool_type, "FALSE", 0, error)
      CALL h5tenum_insert_f(bool_type,  "TRUE", 1, error)
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      success = 1
      buffer = C_LOC(success)
      CALL h5acreate_f(file_id, 'success', bool_type, space_id, attr_id, error)
      CALL h5awrite_f(attr_id, bool_type, buffer, error)
      CALL h5aclose_f(attr_id, error)
      CALL h5sclose_f(space_id, error)
      CALL h5tclose_f(bool_type, error)


      ! create the provenance group
      CALL h5gcreate_f(file_id, 'provenance', group_id, error)
      ! CALL h5screate_f(h5s_scalar_f, space_id, error)

      ! create the properties group
      CALL h5gcreate_f(file_id, 'properties', group_id, error)

      ! create the molecule group
      CALL h5gcreate_f(file_id, 'molecule', group_id, error)

      ! create the wavefunction group
      CALL h5gcreate_f(file_id, 'wavefunction', group_id, error)

      ! create the model group
      CALL h5gcreate_f(file_id, 'model', group_id, error)

      ! create the keywords group
      CALL h5gcreate_f(file_id, 'keywords', group_id, error)


      CALL h5fclose_f(file_id, error)

      CALL h5close_f(error)


   END SUBROUTINE qcschema_to_hdf5

   ! SUBROUTINE qcschema_env_release
   ! END SUBROUTINE qcschema_env_release

! Write to json. write only if quantities exist
   ! SUBROUTINE qcschema_to_json(qcschema_env)
   !    ! TYPE(qs_environment_type), POINTER :: qs_env
   !    TYPE(qcschema_type), POINTER :: qcschema_env
   !    TYPE(json_core)              :: json
   !    TYPE(json_value), POINTER    :: root
   !    INTEGER :: handle

   !    CHARACTER(LEN=*), PARAMETER :: routineN = 'qcschema_to_json'

   !    CALL timeset(routineN, handle)

   !    CALL json%initialize()
   !    CALL json%create_object(root,'')



   !    CALL timestop(handle)
   ! END SUBROUTINE qcschema_to_json



END MODULE qcschema
