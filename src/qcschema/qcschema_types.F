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
      CHARACTER(len=*), DIMENSION(:)      :: symbols        ! The list of atom symbols in this topology
      REAL(KIND=dp),    DIMENSION(:,:)    :: geometry       ! The XYZ coordinates (in Bohr units) of the atoms
      CHARACTER(len=*)                    :: schema_name    ! The name of this schema
      INTEGER                             :: schema_version ! The version of this specific schema
      INTEGER                             :: molecular_charge
      INTEGER                             :: molecular_multiplicity
      ! Whether translation of the geometry is allowed (.False.) or not (.True.)
      LOGICAL                             :: fix_com
      ! A list indicating whether each atom is real (.True.) or a ghost (.False.).
      ! Its length must match that of the `symbols` list
      LOGICAL, DIMENSION(:)               :: real
      ! A list indicating the bonds between the atoms in the molecule. Each item of this list must be
      ! a tuple of three integers, indicating the first atom index in the bond, the second atom index,
      ! and finally the order of the bond
      INTEGER, DIMENSION(:,3)             :: connectivity
      ! Whether rotation of the geometry is allowed (`False`) or not (`True`)
      LOGICAL                             :: fix_orientation
      ! A list of user-provided information for each atom. Its length must match that of the `symbols` list
      CHARACTER(len=*), DIMENSION(:)      :: atom_labels
      ! The list of multiplicities associated with each fragment
      INTEGER, DIMENSION(:)               :: fragment_multiplicities
      ! The maximal point group symmetry at which the `geometry` should be treated
      CHARACTER(len=*)                    :: fix_symmetry
      ! The list of charges associated with each fragment
      INTEGER, DIMENSION(:)               :: fragment_charges
      ! The mass numbers of all atoms. If it is an unknown isotope, the value should be -1.
      ! Its length must match that of the `symbols` list
      INTEGER, DIMENSION(:)               :: mass_numbers
      ! The (user-given) name of the molecule
      CHARACTER(len=*)                    :: name
      ! The masses (in au) of all atoms. Canonical weights are assumed if this is not given explicitly
      REAL(KIND=dp), DIMENSION(:)         :: masses
      ! Any additional (user-provided) comment
      CHARACTER(len=*)                    :: comment
      ! An instance of :class:`QCProvenance`
      TYPE(qcschema_provenance)           :: provenance
      ! The list of fragments. Each item of this list must be a tuple of integers with variable
      ! length (greater than 1). The first number indicates the fragment index, all following numbers
      ! refer to the (0-indexed) atom indices that constitute this fragment
      INTEGER, DIMENSION(:)               :: fragments
      ! The atomic numbers of all atoms, indicating their nuclear charge. Its length must match that
      ! of the `symbols` list
      INTEGER, DIMENSION(:)               :: atomic_numbers
   END TYPE qcschema_topology


! **************************************************************************************************
!> \brief A derived type to store the information of a single electron shell in a basis set.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L43
! **************************************************************************************************
   TYPE qcschema_electron_shell
      ! The angular momenta of this electron shell as a list of integers
      INTEGER, DIMENSION(:)               :: angular_momentum
      ! The type of this shell: spherical or cartesian
      CHARACTER(LEN=*)                    :: harmonic_type
      ! The exponents of this contracted shell. The official spec stores these values as strings
      REAL(KIND=dp), DIMENSION(:)         :: exponents
      ! The general contraction coefficients of this contracted shell
      REAL(KIND=dp), DIMENSION(:,:)       :: coefficients
   END TYPE qcschema_electron_shell


   TYPE qcschema_wavefunction

   END TYPE qcschema_wavefunction


   TYPE qcschema_properties

   END TYPE qcschema_properties

END MODULE qcschema_types
