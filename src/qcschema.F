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
   USE cp_control_types,               ONLY: dft_control_type
   USE cp2k_info,                      ONLY: cp2k_version
   USE cp_log_handling,                ONLY: cp_get_default_logger, &
                                             cp_logger_get_default_io_unit, &
                                             cp_logger_type
   USE kinds,                          ONLY: dp, default_string_length, &
                                             default_path_length, int_8
   USE mp2_types,                      ONLY: mp2_type
   USE particle_types,                 ONLY: particle_type
   USE qs_active_space_types,          ONLY: active_space_type
   USE qs_force_types,                 ONLY: qs_force_type
   USE qs_ks_types,                    ONLY: qs_ks_env_type
   USE qs_mo_types,                    ONLY: mo_set_p_type
   USE qs_scf_types,                   ONLY: qs_scf_env_type
   USE qs_energy_types,                ONLY: qs_energy_type
   USE qs_environment_types,           ONLY: get_qs_env, &
                                             qs_environment_type
   USE hdf5,                           ONLY: h5open_f, h5close_f, &
                                             h5fcreate_f, h5fopen_f, h5fclose_f, &
                                             h5screate_f, h5sclose_f, H5S_SCALAR_F, &
                                             h5gcreate_f, h5gclose_f, &
                                             h5acreate_f, h5awrite_f, h5aclose_f, &
                                             h5tset_size_f, h5tset_cset_f, h5tcopy_f, h5tset_strpad_f, &
                                             h5tcreate_f, h5tenum_insert_f, h5tclose_f, &
                                             H5T_STRING, H5T_CSET_UTF8_F, H5T_ENUM_F, &
                                             H5T_STR_NULLPAD_F, H5T_NATIVE_INTEGER, &
                                             H5T_NATIVE_DOUBLE, HID_T, H5F_ACC_TRUNC_F, &
                                             SIZE_T, H5T_C_S1, HSIZE_T
   USE iso_c_binding,                  ONLY: c_ptr, C_LOC, c_int
   USE h5lt,                           ONLY: h5ltset_attribute_string_f, &
                                             h5ltset_attribute_double_f, &
                                             h5ltmake_dataset_double_f, &
                                             h5ltset_attribute_int_f

#include "./base/base_uses.f90"


   IMPLICIT NONE

   PRIVATE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'qcschema'

   PUBLIC :: qcschema_type
   PUBLIC :: qcschema_env_create, qcschema_env_release, qcschema_to_hdf5

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
      CHARACTER(LEN=2), DIMENSION(:), POINTER      :: symbols => Null()
      ! The XYZ coordinates (in Bohr) of the atoms as a flat array (row major order)
      REAL(KIND=dp),  DIMENSION(:), POINTER        :: geometry => Null()
      CHARACTER(LEN=default_string_length)         :: schema_name    ! The name of this schema
      INTEGER                                      :: schema_version ! The version of this specific schema
      INTEGER                                      :: molecular_charge
      INTEGER                                      :: molecular_multiplicity
      ! Whether translation of the geometry is allowed (.False.) or not (.True.)
      LOGICAL                                      :: fix_com
      ! A list indicating whether each atom is real (.True.) or a ghost (.False.).
      ! Its length must match that of the `symbols` list
      LOGICAL, DIMENSION(:), POINTER               :: real => Null()
      ! A list indicating the bonds between the atoms in the molecule. Each item of this list must be
      ! a tuple of three integers, indicating the first atom index in the bond, the second atom index,
      ! and finally the order of the bond
      INTEGER, DIMENSION(:,:), POINTER             :: connectivity => Null()
      ! Whether rotation of the geometry is allowed (`False`) or not (`True`)
      LOGICAL                                      :: fix_orientation
      ! A list of user-provided information for each atom. Its length must match that of the `symbols` list
      CHARACTER(LEN=8), DIMENSION(:), POINTER      :: atom_labels => Null()
      ! The list of multiplicities associated with each fragment
      INTEGER, DIMENSION(:), POINTER               :: fragment_multiplicities => Null()
      ! The maximal point group symmetry at which the `geometry` should be treated
      CHARACTER(LEN=default_string_length)         :: fix_symmetry
      ! The list of charges associated with each fragment
      INTEGER, DIMENSION(:), POINTER               :: fragment_charges => Null()
      ! The mass numbers of all atoms. If it is an unknown isotope, the value should be -1.
      ! Its length must match that of the `symbols` list
      INTEGER, DIMENSION(:), POINTER               :: mass_numbers => Null()
      ! The (user-given) name of the molecule
      CHARACTER(LEN=default_string_length)         :: name
      ! The masses (in au) of all atoms. Canonical weights are assumed if this is not given explicitly
      REAL(KIND=dp), DIMENSION(:), POINTER         :: masses => Null()
      ! Any additional (user-provided) comment
      CHARACTER(LEN=default_string_length)         :: comment
      ! An instance of :class:`QCProvenance`
      TYPE(qcschema_provenance)                    :: provenance
      ! The list of fragments. Each item of this list must be a tuple of integers with variable
      ! length (greater than 1). The first number indicates the fragment index, all following numbers
      ! refer to the (0-indexed) atom indices that constitute this fragment
      INTEGER, DIMENSION(:), POINTER               :: fragments => Null()
      ! The atomic numbers of all atoms, indicating their nuclear charge. Its length must match that
      ! of the `symbols` list
      INTEGER, DIMENSION(:), POINTER               :: atomic_numbers => Null()
   END TYPE qcschema_topology


! **************************************************************************************************
!> \brief A derived type to store the information of a single electron shell in a basis set.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L43
! **************************************************************************************************
   TYPE qcschema_electron_shell
      ! The angular momenta of this electron shell as a list of integers
      INTEGER, DIMENSION(:), POINTER               :: angular_momentum => Null()
      ! The type of this shell: spherical or cartesian
      CHARACTER(LEN=9)                             :: harmonic_type
      ! The exponents of this contracted shell. The official spec stores these values as strings
      REAL(KIND=dp), DIMENSION(:), POINTER         :: exponents => Null()
      ! The general contraction coefficients of this contracted shell
      REAL(KIND=dp), DIMENSION(:,:), POINTER       :: coefficients => Null()
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
      INTEGER, DIMENSION(:), POINTER               :: angular_momentum => Null()
      ! The exponents of the r terms
      INTEGER, DIMENSION(:), POINTER               :: r_exponents => Null()
      ! The exponents of the Gaussian terms
      REAL(KIND=dp), DIMENSION(:), POINTER         :: gaussian_exponents => Null()
      ! The general contraction coefficients of this potential
      REAL(KIND=dp), DIMENSION(:,:), POINTER       :: coefficients => Null()
   END TYPE qcschema_ecp


! **************************************************************************************************
!> \brief A derived type to store the information of a single atom/center in the basis.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L146
! **************************************************************************************************
   TYPE qcschema_center_basis
      ! The list of electronic shells for this element
      TYPE(qcschema_electron_shell), DIMENSION(:), POINTER :: electron_shells => Null()
      ! The list of effective core potentials for this element
      TYPE(qcschema_ecp), DIMENSION(:), POINTER            :: ecp_potentials => Null()
      ! The number of electrons replaced by an ECP
      INTEGER                                              :: ecp_electrons
   END TYPE qcschema_center_basis


! **************************************************************************************************
!> \brief A derived type to store the information of the basis set used in the calculation.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_basis.html#basis-set-schema
! **************************************************************************************************
   TYPE qcschema_basis_set
      ! The name of the basis set
      CHARACTER(LEN=default_string_length) :: name
      ! A dictionary mapping the keys provided by `atom_map` to their basis center data
      TYPE(qcschema_center_basis), DIMENSION(:), POINTER :: center_data => Null()
      ! The list of atomic kinds, indicating the keys used to store the basis in `center_data`
      ! Not clear if this will be of the length of the basis set size, or rather just one
      ! entry for atomic kind. E.g. only one entry for hydrogen even though there might be
      ! many hydrogen atoms in the molecule. If this is the case, then we really need a
      ! hash table for `center_data`
      CHARACTER(LEN=2), DIMENSION(:), POINTER            :: atom_map => Null()
      ! The version of this specific schema
      INTEGER                                            :: schema_version
      ! The name of this schema. This value is expected to be `qcschema_basis`
      CHARACTER(LEN=default_string_length)               :: schema_name
      ! A description of this basis set
      CHARACTER(LEN=default_string_length)               :: description
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
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_orbitals_a    => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_orbitals_b    => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_density_a     => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_density_b     => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_density_mo_a  => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_density_mo_b  => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_fock_a        => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_fock_b        => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_fock_mo_a     => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_fock_mo_b     => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_coulomb_a     => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_coulomb_b     => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_exchange_a    => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_exchange_b    => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_eigenvalues_a => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_eigenvalues_b => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_occupations_a => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_occupations_b => Null()

      ! ERI should be arrays of doubles, but in Qiskit are simple strings?
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_eri       => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_eri_mo_aa => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_eri_mo_ab => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_eri_mo_ba => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: scf_eri_mo_bb => Null()

      ! Quantities with localized orbitals. All `nmo` orbitals are included,
      ! even if only a subset were localized
      REAL(KIND=dp), DIMENSION(:), POINTER         :: localized_orbitals_a => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: localized_orbitals_b => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: localized_fock_a     => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: localized_fock_b     => Null()

      ! Nne-electron Hamiltonians in AO basis
      REAL(KIND=dp), DIMENSION(:), POINTER         :: h_core_a      => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: h_core_b      => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: h_effective_a => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER         :: h_effective_b => Null()

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

      ! internal flags to know the type of calculation
      LOGICAL :: mp2 = .FALSE.

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
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(active_space_type), POINTER                   :: active_space_env

      INTEGER :: handle, output_unit, nspins, nalpha, nbeta, natom, nbasis, nmo
      INTEGER :: i, m
      REAL(KIND=dp) :: one_el_en, two_el_en, dispersion
      CHARACTER(LEN=*), PARAMETER :: routineN = 'qcschema_env_initialize'

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      ! Collect environment info
      IF (ASSOCIATED(qs_env)) THEN
         CALL get_qs_env(qs_env, ks_env=ks_env, energy=energy, &
                         dft_control=dft_control, force=force, &
                         mos=mos, particle_set=particle_set, &
                         scf_env=scf_env, mp2_env=mp2_env, &
                         active_space_env=active_space_env)
      ELSE
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: QS environment is not associated, QCSchema interface quitting'
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
      ! TODO: no way to get this one?
      qcschema_env%properties%nuclear_repulsion_energy = 0.0_dp

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

      ! was it a DFT calculation?
      ! IF (ASSOCIATED(qs_env%ks_env)) qcschema_env%properties%dft = .TRUE.
      ! with dispersion?
      ! IF (ASSOCIATED(qs_env%dispersion_env)) qcschema_env%properties%dispersion = .TRUE.

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
         qcschema_env%properties%mp2 = .TRUE.
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


      ! *** QCSchema wavefunction *** !
      qcschema_env%wavefunction%scf_eigenvalues_a => mos(1)%mo_set%eigenvalues
      IF (nspins == 2) THEN
         qcschema_env%wavefunction%scf_eigenvalues_b => mos(2)%mo_set%eigenvalues
      END IF

      m = (nspins*(nspins + 1))/2
      DO i = 1,m
         active_space_env%eri%eri(i)%csr_mat
      END DO



      CALL timestop(handle)

   END SUBROUTINE qcschema_env_create


   SUBROUTINE qcschema_env_release(qcschema_env)
      TYPE(qcschema_type), INTENT(INOUT)      :: qcschema_env

      IF (ALLOCATED(qcschema_env%return_result)) THEN
         DEALLOCATE (qcschema_env%return_result)
      END IF

   END SUBROUTINE qcschema_env_release

   SUBROUTINE qcschema_to_hdf5(qcschema_env)
      TYPE(qcschema_type), INTENT(IN)        :: qcschema_env
      TYPE(cp_logger_type), POINTER          :: logger
      INTEGER                                :: error, output_unit
      INTEGER(KIND=int_8)                    :: nresult
      INTEGER(KIND=hid_t)                    :: file_id, space_id, group_id
      CHARACTER(LEN=default_path_length)     :: filename
      INTEGER(KIND=hsize_t), DIMENSION(1:1)  :: dims

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

      ! !===========================================================================!
      ! Root
      ! !===========================================================================!
      ! driver
      CALL h5awrite_fixlen_string(file_id, 'driver', trim(qcschema_env%driver))
      ! return result
      nresult = SIZE(qcschema_env%return_result)
      IF (SIZE(qcschema_env%return_result) == 1) THEN
         CALL h5awrite_double_scalar(file_id, 'return_result', qcschema_env%return_result(1))
      ELSE
         CALL h5ltset_attribute_double_f(file_id, '/', 'return_result', &
                                         qcschema_env%return_result, nresult, error)
      END IF
      ! schema name
      CALL h5awrite_fixlen_string(file_id, 'schema_name', 'qcschema')
      ! schema version
      CALL h5awrite_integer_scalar(file_id, 'schema_version', 3)
      ! success
      CALL h5awrite_boolean(file_id, 'success', qcschema_env%success)

      ! !===========================================================================!
      ! Provenance
      ! !===========================================================================!
      ! create the provenance group
      CALL h5gcreate_f(file_id, 'provenance', group_id, error)
      ! populate provenance
      CALL h5awrite_fixlen_string(group_id, 'creator', trim(qcschema_env%provenance%creator))
      CALL h5awrite_fixlen_string(group_id, 'routine', trim(qcschema_env%provenance%routine))
      CALL h5awrite_fixlen_string(group_id, 'version', trim(qcschema_env%provenance%version))
      ! close provenance group
      CALL h5gclose_f(group_id, error)

      ! !===========================================================================!
      ! Properties
      ! !===========================================================================!
      ! create the properties group
      CALL h5gcreate_f(file_id, 'properties', group_id, error)
      ! populate properties
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_natom',  qcschema_env%properties%calcinfo_natom )
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nbasis', qcschema_env%properties%calcinfo_nbasis)
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nmo',    qcschema_env%properties%calcinfo_nmo   )
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nalpha', qcschema_env%properties%calcinfo_nalpha)
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nbeta',  qcschema_env%properties%calcinfo_nbeta )
      ! energies, scf, mp2, ...
      CALL h5awrite_double_scalar(group_id, 'return_energy', qcschema_env%properties%return_energy)
      CALL h5awrite_double_scalar(group_id, 'nuclear_repulsion_energy', &
                                  qcschema_env%properties%nuclear_repulsion_energy)
      CALL h5awrite_integer_scalar(group_id, 'scf_iterations', qcschema_env%properties%scf_iterations)
      CALL h5awrite_double_scalar(group_id, 'scf_one_electron_energy', &
                                  qcschema_env%properties%scf_one_electron_energy)
      CALL h5awrite_double_scalar(group_id, 'scf_two_electron_energy', &
                                  qcschema_env%properties%scf_two_electron_energy)
      CALL h5awrite_double_scalar(group_id, 'scf_total_energy', qcschema_env%properties%scf_total_energy)
      IF (qcschema_env%properties%scf_xc_energy /= 0.0_dp) THEN
         CALL h5awrite_double_scalar(group_id, 'scf_xc_energy', &
                                     qcschema_env%properties%scf_xc_energy)
      END IF
      IF (qcschema_env%properties%scf_dispersion_correction_energy /= 0.0_dp) THEN
         CALL h5awrite_double_scalar(group_id, 'scf_dispersion_correction_energy', &
                                     qcschema_env%properties%scf_dispersion_correction_energy)
      END IF
      IF (qcschema_env%properties%mp2) THEN
         CALL h5awrite_double_scalar(group_id, 'mp2_correlation_energy', &
                                     qcschema_env%properties%mp2_correlation_energy)
      END IF
      ! close properties group
      CALL h5gclose_f(group_id, error)

      ! !===========================================================================!
      ! Molecule
      ! !===========================================================================!
      ! create the molecule group
      CALL h5gcreate_f(file_id, 'molecule', group_id, error)
      ! close molecule group
      CALL h5gclose_f(group_id, error)

      ! create the wavefunction group
      CALL h5gcreate_f(file_id, 'wavefunction', group_id, error)
      dims(1) = size(qcschema_env%wavefunction%scf_eigenvalues_a, kind=hsize_t)
      CALL h5ltmake_dataset_double_f(group_id, 'scf_eigenvalues_a', 1, dims, &
                                     qcschema_env%wavefunction%scf_eigenvalues_a, error)
      IF (nspins = 2) THEN
         dims(1) = size(qcschema_env%wavefunction%scf_eigenvalues_b, kind=hsize_t)
         CALL h5ltmake_dataset_double_f(group_id, 'scf_eigenvalues_b', 1, dims, &
                                        qcschema_env%wavefunction%scf_eigenvalues_b, error)
      END IF
      ! close wavefunction group
      CALL h5gclose_f(group_id, error)

      ! create the model group
      CALL h5gcreate_f(file_id, 'model', group_id, error)
      CALL h5screate_f(h5s_scalar_f, space_id, error)
      CALL h5sclose_f(space_id, error)
      ! close model group
      CALL h5gclose_f(group_id, error)

      ! create the keywords group
      CALL h5gcreate_f(file_id, 'keywords', group_id, error)
      ! close keywords group
      CALL h5gclose_f(group_id, error)


      CALL h5fclose_f(file_id, error)

      CALL h5close_f(error)


   END SUBROUTINE qcschema_to_hdf5


   SUBROUTINE h5awrite_varlen_string(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      CHARACTER(LEN=*), INTENT(IN), TARGET   :: attr_data
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      TYPE(c_ptr), TARGET                    :: in_between_ptr
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error, output_unit

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! create a variable-length string type
      CALL h5tcopy_f(H5T_STRING, type_id, error)
      CAll h5tset_cset_f(type_id, H5T_CSET_UTF8_F, error)
      CALL h5tset_strpad_f(type_id, H5T_STR_NULLPAD_F, error)

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 attribute failed'
         RETURN
      ENDIF

      ! weird in-between pointer needed for variable-length
      ! string to a scalar dataspace
      in_between_ptr = C_LOC(attr_data)
      ! the actual pointer to be passed
      buffer = C_LOC(in_between_ptr)

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: writing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 datatype failed'
         RETURN
      ENDIF

   END SUBROUTINE h5awrite_varlen_string


   SUBROUTINE h5awrite_fixlen_string(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      CHARACTER(LEN=*), INTENT(IN), TARGET   :: attr_data
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error, output_unit

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! create a fixed-length string datatype
      CALL h5tcopy_f(H5T_C_S1, type_id, error)
      CAll h5tset_cset_f(type_id, H5T_CSET_UTF8_F, error)
      CALL h5tset_size_f(type_id, LEN(attr_data, int_8), error)

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 attribute failed'
         RETURN
      ENDIF

      ! the actual pointer to be passed
      buffer = C_LOC(attr_data)

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: writing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 datatype failed'
         RETURN
      ENDIF

   END SUBROUTINE h5awrite_fixlen_string


   SUBROUTINE h5awrite_boolean(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      LOGICAL, INTENT(IN)                    :: attr_data
      INTEGER, TARGET                        :: attr_data_to_int
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      INTEGER(KIND=size_t)                   :: enum_size = 1 ! 8-bit integers in enum bool_type
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error, output_unit

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! create the datatype
      CALL h5tcreate_f(H5T_ENUM_F, enum_size, type_id, error)
      CALL h5tenum_insert_f(type_id, "FALSE", 0, error)
      CALL h5tenum_insert_f(type_id,  "TRUE", 1, error)

      IF (attr_data) THEN
         attr_data_to_int = 1
      ELSE
         attr_data_to_int = 0
      END IF
      ! the C pointer to the actual data
      buffer = C_LOC(attr_data_to_int)

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 attribute failed'
         RETURN
      ENDIF

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: writing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 datatype failed'
         RETURN
      ENDIF

   END SUBROUTINE h5awrite_boolean


   SUBROUTINE h5awrite_integer_scalar(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      INTEGER, INTENT(IN), TARGET            :: attr_data
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error, output_unit

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! the C pointer to the actual data
      buffer = C_LOC(attr_data)

      ! set the type of data
      type_id = H5T_NATIVE_INTEGER

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 attribute failed'
         RETURN
      ENDIF

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: writing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 dataspace failed'
         RETURN
      ENDIF

   END SUBROUTINE h5awrite_integer_scalar


   SUBROUTINE h5awrite_double_scalar(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      REAL(KIND=dp), INTENT(IN), TARGET      :: attr_data
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error, output_unit

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! the C pointer to the actual data
      buffer = C_LOC(attr_data)

      ! set the type of data
      type_id = H5T_NATIVE_DOUBLE

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 attribute failed'
         RETURN
      ENDIF

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: writing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 dataspace failed'
         RETURN
      ENDIF

   END SUBROUTINE h5awrite_double_scalar


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
