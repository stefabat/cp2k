!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The types defining the QCSchema
!> \par History
!>      10.2022 created [SB]
!> \author SB
! **************************************************************************************************
MODULE qcschema_types

   USE kinds,                           ONLY: dp

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qcschema_types'


! **************************************************************************************************
!> \brief A derived type to store the program information that generated the QCSchema file.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#provenance
! **************************************************************************************************
   TYPE qcschema_provenance
      CHARACTER(len=*)     :: creator ! The name of the creator of this object
      CHARACTER(len=*)     :: version ! The version of the creator of this object
      CHARACTER(len=*)     :: routine ! The routine that was used to create this object
   END TYPE qcschema_provenance


! **************************************************************************************************
!> \brief A derived type to store the topological information of the physical system.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#topology
! **************************************************************************************************
   TYPE qcschema_topology
      ! The list of atom symbols in this topology
      CHARACTER(len=*), DIMENSION(:), ALLOCATABLE  :: symbols
      ! The XYZ coordinates (in Bohr units) of the atoms
      REAL(KIND=dp),  DIMENSION(:,:), ALLOCATABLE  :: geometry
      CHARACTER(len=*)                             :: schema_name    ! The name of this schema
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
      INTEGER, DIMENSION(:,3), ALLOCATABLE         :: connectivity
      ! Whether rotation of the geometry is allowed (`False`) or not (`True`)
      LOGICAL                                      :: fix_orientation
      ! A list of user-provided information for each atom. Its length must match that of the `symbols` list
      CHARACTER(len=*), DIMENSION(:), ALLOCATABLE  :: atom_labels
      ! The list of multiplicities associated with each fragment
      INTEGER, DIMENSION(:), ALLOCATABLE           :: fragment_multiplicities
      ! The maximal point group symmetry at which the `geometry` should be treated
      CHARACTER(len=*)                             :: fix_symmetry
      ! The list of charges associated with each fragment
      INTEGER, DIMENSION(:), ALLOCATABLE           :: fragment_charges
      ! The mass numbers of all atoms. If it is an unknown isotope, the value should be -1.
      ! Its length must match that of the `symbols` list
      INTEGER, DIMENSION(:), ALLOCATABLE           :: mass_numbers
      ! The (user-given) name of the molecule
      CHARACTER(len=*)                             :: name
      ! The masses (in au) of all atoms. Canonical weights are assumed if this is not given explicitly
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: masses
      ! Any additional (user-provided) comment
      CHARACTER(len=*)                             :: comment
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
      CHARACTER(LEN=*)                             :: harmonic_type
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
      CHARACTER(LEN=*)                             :: ecp_type
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
      CHARACTER(LEN=*)                             :: name
      ! A dictionary mapping the keys provided by `atom_map` to their basis center data
      TYPE(qcschema_center_basis), DIMENSION(:), ALLOCATABLE &
      :: center_data
      ! The list of atomic kinds, indicating the keys used to store the basis in `center_data`
      ! Not clear if this will be of the length of the basis set size, or rather just one
      ! entry for atomic kind. E.g. only one entry for hydrogen even though there might be
      ! many hydrogen atoms in the molecule. If this is the case, then we really need a
      ! hash table for `center_data`
      CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE  :: atom_map
      ! The version of this specific schema
      INTEGER                                      :: schema_version
      ! The name of this schema. This value is expected to be `qcschema_basis`
      CHARACTER(LEN=*)                             :: schema_name
      ! A description of this basis set
      CHARACTER(LEN=*)                             :: description
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
      CHARACTER(LEN=*)                             :: orbitals_a
      CHARACTER(LEN=*)                             :: orbitals_b
      CHARACTER(LEN=*)                             :: density_a
      CHARACTER(LEN=*)                             :: density_b
      CHARACTER(LEN=*)                             :: density_mo_a
      CHARACTER(LEN=*)                             :: density_mo_b
      CHARACTER(LEN=*)                             :: fock_a
      CHARACTER(LEN=*)                             :: fock_b
      CHARACTER(LEN=*)                             :: fock_mo_a
      CHARACTER(LEN=*)                             :: fock_mo_b
      CHARACTER(LEN=*)                             :: eigenvalues_a
      CHARACTER(LEN=*)                             :: eigenvalues_b
      CHARACTER(LEN=*)                             :: occupations_a
      CHARACTER(LEN=*)                             :: occupations_b
      CHARACTER(LEN=*)                             :: eri
      CHARACTER(LEN=*)                             :: eri_mo_aa
      CHARACTER(LEN=*)                             :: eri_mo_ab
      CHARACTER(LEN=*)                             :: eri_mo_ba
      CHARACTER(LEN=*)                             :: eri_mo_bb

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

   INTEGER :: calcinfo_nbasis
   INTEGER :: calcinfo_nmo
   INTEGER :: calcinfo_nalpha ! # of alpha electrons
   INTEGER :: calcinfo_nbeta  ! # of beta electrons
   INTEGER :: calcinfo_natom
   REAL(KIND=dp) :: return_energy

   ! SCF energies
   REAL(KIND=dp) :: scf_one_electron_energy
   REAL(KIND=dp) :: scf_two_electron_energy
   REAL(KIND=dp) :: nuclear_repulsion_energy
   REAL(KIND=dp) :: scf_vv10_energy
   REAL(KIND=dp) :: scf_xc_energy
   REAL(KIND=dp) :: scf_dispersion_correction_energy
   REAL(KIND=dp) :: scf_total_energy
   REAL(KIND=dp), DIMENSION(3), ALLOCATABLE :: scf_dipole_moment
   INTEGER :: scf_iterations

   REAL(KIND=dp) :: mp2_same_spin_correlation_energy
   REAL(KIND=dp) :: mp2_opposite_spin_correlation_energy
   REAL(KIND=dp) :: mp2_singles_energy
   REAL(KIND=dp) :: mp2_doubles_energy
   REAL(KIND=dp) :: mp2_correlation_energy
   REAL(KIND=dp) :: mp2_total_energy
   REAL(KIND=dp), DIMENSION(3), ALLOCATABLE :: mp2_dipole_moment

   REAL(KIND=dp) :: ccsd_same_spin_correlation_energy
   REAL(KIND=dp) :: ccsd_opposite_spin_correlation_energy
   REAL(KIND=dp) :: ccsd_singles_energy
   REAL(KIND=dp) :: ccsd_doubles_energy
   REAL(KIND=dp) :: ccsd_correlation_energy
   REAL(KIND=dp) :: ccsd_total_energy
   REAL(KIND=dp) :: ccsd_prt_pr_correlation_energy
   REAL(KIND=dp) :: ccsd_prt_pr_total_energy
   REAL(KIND=dp) :: ccsdt_correlation_energy
   REAL(KIND=dp) :: ccsdt_total_energy
   REAL(KIND=dp) :: ccsdtq_correlation_energy
   REAL(KIND=dp) :: ccsdtq_total_energy
   REAL(KIND=dp), DIMENSION(3), ALLOCATABLE :: ccsd_dipole_moment
   REAL(KIND=dp), DIMENSION(3), ALLOCATABLE :: ccsd_prt_pr_dipole_moment
   REAL(KIND=dp), DIMENSION(3), ALLOCATABLE :: ccsdt_dipole_moment
   REAL(KIND=dp), DIMENSION(3), ALLOCATABLE :: ccsdtq_dipole_moment
   INTEGER :: ccsd_iterations: int | None = None
   INTEGER :: ccsdt_iterations: int | None = None
   INTEGER :: ccsdtq_iterations: int | None = None

   END TYPE qcschema_properties

END MODULE qcschema_types
