!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!> This module contains routines necessary to operate on plane waves on GPUs
!     >  independently of the GPU platform.
!> \par History
!>      BGL (06-Mar-2008)  : Created
!>      AG  (18-May-2012)  : Refacturing:
!>                           - added explicit interfaces to C routines
!>                           - enable double precision complex transformations
!>      AG  (11-Sept-2012) : Modifications:
!>                          - use pointers if precision mapping is not required
!>                          - use OMP for mapping
!>      MT  (Jan 2022)     : Modifications
!>                          - use a generic interface for fft calls to GPUs
!>                          - Support both Nvidia and AMD GPUs. Other GPUs manufacturers
!>                            can be added easily.
!> \author Benjamin G. Levine
! **************************************************************************************************
MODULE pw_gpu
   USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE fft_tools,                       ONLY: &
        cube_transpose_1, cube_transpose_2, fft_scratch_sizes, fft_scratch_type, get_fft_scratch, &
        release_fft_scratch, x_to_yz, xz_to_yz, yz_to_x, yz_to_xz
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_zero
   USE message_passing,                 ONLY: mp_cart_type
   USE pw_grid_types,                   ONLY: FULLSPACE
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pw_gpu_r3dc1d_3d
   PUBLIC :: pw_gpu_c1dr3d_3d
   PUBLIC :: pw_gpu_r3dc1d_3d_ps
   PUBLIC :: pw_gpu_c1dr3d_3d_ps
   PUBLIC :: pw_gpu_init, pw_gpu_finalize

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_gpu'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Allocates resources on the gpu device for gpu fft acceleration
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE pw_gpu_init()
      INTEGER                                            :: dummy
      INTERFACE
         SUBROUTINE pw_gpu_init_c() BIND(C, name="pw_gpu_init")
         END SUBROUTINE pw_gpu_init_c
      END INTERFACE

      MARK_USED(dummy) ! TODO: fix fpretty
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_init_c()
#else
      ! Nothing to do.
#endif
   END SUBROUTINE pw_gpu_init

! **************************************************************************************************
!> \brief Releases resources on the gpu device for gpu fft acceleration
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE pw_gpu_finalize()
      INTEGER                                            :: dummy
      INTERFACE
         SUBROUTINE pw_gpu_finalize_c() BIND(C, name="pw_gpu_finalize")
         END SUBROUTINE pw_gpu_finalize_c
      END INTERFACE

      MARK_USED(dummy) ! TODO: fix fpretty
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_finalize_c()
#else
      ! Nothing to do.
#endif
   END SUBROUTINE pw_gpu_finalize

! **************************************************************************************************
!> \brief perform an fft followed by a gather on the gpu
!> \param pw1 ...
!> \param pw2 ...
!> \author Benjamin G Levine
! **************************************************************************************************
   SUBROUTINE pw_gpu_r3dc1d_3d(pw1, pw2)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw1
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: pw2

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_r3dc1d_3d'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwout
      INTEGER                                            :: handle, l1, l2, l3, ngpts
      INTEGER, DIMENSION(:), POINTER                     :: npts
      INTEGER, POINTER                                   :: ptr_ghatmap
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), POINTER                             :: ptr_pwin
      INTERFACE
         SUBROUTINE pw_gpu_cfffg_c(din, zout, ghatmap, npts, ngpts, scale) BIND(C, name="pw_gpu_cfffg")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: din
            TYPE(C_PTR), VALUE                           :: zout
            TYPE(C_PTR), INTENT(IN), VALUE               :: ghatmap
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: ngpts
            REAL(KIND=C_DOUBLE), INTENT(IN), VALUE       :: scale

         END SUBROUTINE pw_gpu_cfffg_c
      END INTERFACE

      CALL timeset(routineN, handle)

      scale = 1.0_dp/REAL(pw1%pw_grid%ngpts, KIND=dp)

      ngpts = SIZE(pw2%pw_grid%gsq)
      l1 = LBOUND(pw1%array, 1)
      l2 = LBOUND(pw1%array, 2)
      l3 = LBOUND(pw1%array, 3)
      npts => pw1%pw_grid%npts

      ! pointers to data arrays
      ptr_pwin => pw1%array(l1, l2, l3)
      ptr_pwout => pw2%array(1)

      ! pointer to map array
      ptr_ghatmap => pw2%pw_grid%g_hatmap(1, 1)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_cfffg_c(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, ngpts, scale)
#else
      CPABORT("Compiled without pw offloading.")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_r3dc1d_3d

! **************************************************************************************************
!> \brief perform an scatter followed by a fft on the gpu
!> \param pw1 ...
!> \param pw2 ...
!> \author Benjamin G Levine
! **************************************************************************************************
   SUBROUTINE pw_gpu_c1dr3d_3d(pw1, pw2)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: pw1
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: pw2

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_c1dr3d_3d'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin
      INTEGER                                            :: handle, l1, l2, l3, ngpts, nmaps
      INTEGER, DIMENSION(:), POINTER                     :: npts
      INTEGER, POINTER                                   :: ptr_ghatmap
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), POINTER                             :: ptr_pwout
      INTERFACE
         SUBROUTINE pw_gpu_sfffc_c(zin, dout, ghatmap, npts, ngpts, nmaps, scale) BIND(C, name="pw_gpu_sfffc")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: dout
            TYPE(C_PTR), INTENT(IN), VALUE               :: ghatmap
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: ngpts, nmaps
            REAL(KIND=C_DOUBLE), INTENT(IN), VALUE       :: scale
         END SUBROUTINE pw_gpu_sfffc_c
      END INTERFACE

      CALL timeset(routineN, handle)

      scale = 1.0_dp

      ngpts = SIZE(pw1%pw_grid%gsq)
      l1 = LBOUND(pw2%array, 1)
      l2 = LBOUND(pw2%array, 2)
      l3 = LBOUND(pw2%array, 3)
      npts => pw1%pw_grid%npts

      ! pointers to data arrays
      ptr_pwin => pw1%array(1)
      ptr_pwout => pw2%array(l1, l2, l3)

      ! pointer to map array
      nmaps = SIZE(pw1%pw_grid%g_hatmap, 2)
      ptr_ghatmap => pw1%pw_grid%g_hatmap(1, 1)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_sfffc_c(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, ngpts, nmaps, scale)
#else
      CPABORT("Compiled without pw offloading")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_c1dr3d_3d

! **************************************************************************************************
!> \brief perform an parallel fft followed by a gather on the gpu
!> \param pw1 ...
!> \param pw2 ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_r3dc1d_3d_ps(pw1, pw2)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw1
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: pw2

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_gpu_r3dc1d_3d_ps'

      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: grays, pbuf, qbuf, rbuf, sbuf
      COMPLEX(KIND=dp), DIMENSION(:, :, :), POINTER      :: tbuf
      INTEGER                                            :: g_pos, handle, lg, lmax, mg, mmax, mx2, &
                                                            mz2, n1, n2, ngpts, nmax, numtask, rp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: p2p
      INTEGER, DIMENSION(2)                              :: r_dim, r_pos
      INTEGER, DIMENSION(:), POINTER                     :: n, nloc, nyzray
      INTEGER, DIMENSION(:, :, :, :), POINTER            :: bo
      REAL(KIND=dp)                                      :: scale
      TYPE(fft_scratch_sizes)                            :: fft_scratch_size
      TYPE(fft_scratch_type), POINTER                    :: fft_scratch
      TYPE(mp_cart_type)                                 :: rs_group

      CALL timeset(routineN, handle)

      scale = 1.0_dp/REAL(pw1%pw_grid%ngpts, KIND=dp)

      ! dimensions
      n => pw1%pw_grid%npts
      nloc => pw1%pw_grid%npts_local
      grays => pw1%pw_grid%grays
      ngpts = nloc(1)*nloc(2)*nloc(3)

      !..transform
      IF (pw1%pw_grid%para%ray_distribution) THEN
         rs_group = pw1%pw_grid%para%group
         nyzray => pw1%pw_grid%para%nyzray
         bo => pw1%pw_grid%para%bo

         g_pos = rs_group%mepos
         numtask = rs_group%num_pe
         r_dim = rs_group%num_pe_cart
         r_pos = rs_group%mepos_cart

         lg = SIZE(grays, 1)
         mg = SIZE(grays, 2)
         mmax = MAX(mg, 1)
         lmax = MAX(lg, (ngpts/mmax + 1))

         ALLOCATE (p2p(0:numtask - 1))

         CALL rs_group%rank_compare(rs_group, p2p)

         rp = p2p(g_pos)
         mx2 = bo(2, 1, rp, 2) - bo(1, 1, rp, 2) + 1
         mz2 = bo(2, 3, rp, 2) - bo(1, 3, rp, 2) + 1
         n1 = MAXVAL(bo(2, 1, :, 1) - bo(1, 1, :, 1) + 1)
         n2 = MAXVAL(bo(2, 2, :, 1) - bo(1, 2, :, 1) + 1)
         nmax = MAX((2*n2)/numtask, 2)*mx2*mz2
         nmax = MAX(nmax, n1*MAXVAL(nyzray))

         fft_scratch_size%nx = nloc(1)
         fft_scratch_size%ny = nloc(2)
         fft_scratch_size%nz = nloc(3)
         fft_scratch_size%lmax = lmax
         fft_scratch_size%mmax = mmax
         fft_scratch_size%mx1 = bo(2, 1, rp, 1) - bo(1, 1, rp, 1) + 1
         fft_scratch_size%mx2 = mx2
         fft_scratch_size%my1 = bo(2, 2, rp, 1) - bo(1, 2, rp, 1) + 1
         fft_scratch_size%mz2 = mz2
         fft_scratch_size%lg = lg
         fft_scratch_size%mg = mg
         fft_scratch_size%nbx = MAXVAL(bo(2, 1, :, 2))
         fft_scratch_size%nbz = MAXVAL(bo(2, 3, :, 2))
         fft_scratch_size%mcz1 = MAXVAL(bo(2, 3, :, 1) - bo(1, 3, :, 1) + 1)
         fft_scratch_size%mcx2 = MAXVAL(bo(2, 1, :, 2) - bo(1, 1, :, 2) + 1)
         fft_scratch_size%mcz2 = MAXVAL(bo(2, 3, :, 2) - bo(1, 3, :, 2) + 1)
         fft_scratch_size%nmax = nmax
         fft_scratch_size%nmray = MAXVAL(nyzray)
         fft_scratch_size%nyzray = nyzray(g_pos)
         fft_scratch_size%rs_group = rs_group
         fft_scratch_size%g_pos = g_pos
         fft_scratch_size%r_pos = r_pos
         fft_scratch_size%r_dim = r_dim
         fft_scratch_size%numtask = numtask

         IF (r_dim(2) > 1) THEN
            !
            ! real space is distributed over x and y coordinate
            ! we have two stages of communication
            !
            IF (r_dim(1) == 1) &
               CPABORT("This processor distribution is not supported.")

            CALL get_fft_scratch(fft_scratch, tf_type=300, n=n, fft_sizes=fft_scratch_size)

            ! assign buffers
            qbuf => fft_scratch%p2buf
            rbuf => fft_scratch%p3buf
            pbuf => fft_scratch%p4buf
            sbuf => fft_scratch%p5buf

            ! FFT along z
            CALL pw_gpu_cf(pw1, qbuf)

            ! Exchange data ( transpose of matrix )
            CALL cube_transpose_2(qbuf, bo(:, :, :, 1), bo(:, :, :, 2), rbuf, fft_scratch)

            ! FFT along y
            ! use the inbuild fft-lib
            ! CALL fft_1dm(fft_scratch%fft_plan(2), rbuf, pbuf, 1.0_dp, stat)
            ! or cufft (works faster, but is only faster if plans are stored)
            CALL pw_gpu_f(rbuf, pbuf, +1, n(2), mx2*mz2)

            ! Exchange data ( transpose of matrix ) and sort
            CALL xz_to_yz(pbuf, rs_group, r_dim, g_pos, p2p, pw1%pw_grid%para%yzp, nyzray, &
                          bo(:, :, :, 2), sbuf, fft_scratch)

            ! FFT along x
            CALL pw_gpu_fg(sbuf, pw2, scale)

            CALL release_fft_scratch(fft_scratch)

         ELSE
            !
            ! real space is only distributed over x coordinate
            ! we have one stage of communication, after the transform of
            ! direction x
            !

            CALL get_fft_scratch(fft_scratch, tf_type=200, n=n, fft_sizes=fft_scratch_size)

            ! assign buffers
            tbuf => fft_scratch%tbuf
            sbuf => fft_scratch%r1buf

            ! FFT along y and z
            CALL pw_gpu_cff(pw1, tbuf)

            ! Exchange data ( transpose of matrix ) and sort
            CALL yz_to_x(tbuf, rs_group, g_pos, p2p, pw1%pw_grid%para%yzp, nyzray, &
                         bo(:, :, :, 2), sbuf, fft_scratch)

            ! FFT along x
            CALL pw_gpu_fg(sbuf, pw2, scale)

            CALL release_fft_scratch(fft_scratch)

         END IF

         DEALLOCATE (p2p)

!--------------------------------------------------------------------------
      ELSE
         CPABORT("Not implemented (no ray_distr.) in: pw_gpu_r3dc1d_3d_ps.")
      END IF

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_r3dc1d_3d_ps

! **************************************************************************************************
!> \brief perform an parallel scatter followed by a fft on the gpu
!> \param pw1 ...
!> \param pw2 ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_c1dr3d_3d_ps(pw1, pw2)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: pw1
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: pw2

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_gpu_c1dr3d_3d_ps'

      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: grays, pbuf, qbuf, rbuf, sbuf
      COMPLEX(KIND=dp), DIMENSION(:, :, :), POINTER      :: tbuf
      INTEGER                                            :: g_pos, handle, lg, lmax, mg, mmax, mx2, &
                                                            mz2, n1, n2, ngpts, nmax, numtask, rp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: p2p
      INTEGER, DIMENSION(2)                              :: r_dim, r_pos
      INTEGER, DIMENSION(:), POINTER                     :: n, nloc, nyzray
      INTEGER, DIMENSION(:, :, :, :), POINTER            :: bo
      REAL(KIND=dp)                                      :: scale
      TYPE(fft_scratch_sizes)                            :: fft_scratch_size
      TYPE(fft_scratch_type), POINTER                    :: fft_scratch
      TYPE(mp_cart_type)                                 :: rs_group

      CALL timeset(routineN, handle)

      scale = 1.0_dp

      ! dimensions
      n => pw1%pw_grid%npts
      nloc => pw1%pw_grid%npts_local
      grays => pw1%pw_grid%grays
      ngpts = nloc(1)*nloc(2)*nloc(3)

      !..transform
      IF (pw1%pw_grid%para%ray_distribution) THEN
         rs_group = pw1%pw_grid%para%group
         nyzray => pw1%pw_grid%para%nyzray
         bo => pw1%pw_grid%para%bo

         g_pos = rs_group%mepos
         numtask = rs_group%num_pe
         r_dim = rs_group%num_pe_cart
         r_pos = rs_group%mepos_cart

         lg = SIZE(grays, 1)
         mg = SIZE(grays, 2)
         mmax = MAX(mg, 1)
         lmax = MAX(lg, (ngpts/mmax + 1))

         ALLOCATE (p2p(0:numtask - 1))

         CALL rs_group%rank_compare(rs_group, p2p)

         rp = p2p(g_pos)
         mx2 = bo(2, 1, rp, 2) - bo(1, 1, rp, 2) + 1
         mz2 = bo(2, 3, rp, 2) - bo(1, 3, rp, 2) + 1
         n1 = MAXVAL(bo(2, 1, :, 1) - bo(1, 1, :, 1) + 1)
         n2 = MAXVAL(bo(2, 2, :, 1) - bo(1, 2, :, 1) + 1)
         nmax = MAX((2*n2)/numtask, 2)*mx2*mz2
         nmax = MAX(nmax, n1*MAXVAL(nyzray))

         fft_scratch_size%nx = nloc(1)
         fft_scratch_size%ny = nloc(2)
         fft_scratch_size%nz = nloc(3)
         fft_scratch_size%lmax = lmax
         fft_scratch_size%mmax = mmax
         fft_scratch_size%mx1 = bo(2, 1, rp, 1) - bo(1, 1, rp, 1) + 1
         fft_scratch_size%mx2 = mx2
         fft_scratch_size%my1 = bo(2, 2, rp, 1) - bo(1, 2, rp, 1) + 1
         fft_scratch_size%mz2 = mz2
         fft_scratch_size%lg = lg
         fft_scratch_size%mg = mg
         fft_scratch_size%nbx = MAXVAL(bo(2, 1, :, 2))
         fft_scratch_size%nbz = MAXVAL(bo(2, 3, :, 2))
         fft_scratch_size%mcz1 = MAXVAL(bo(2, 3, :, 1) - bo(1, 3, :, 1) + 1)
         fft_scratch_size%mcx2 = MAXVAL(bo(2, 1, :, 2) - bo(1, 1, :, 2) + 1)
         fft_scratch_size%mcz2 = MAXVAL(bo(2, 3, :, 2) - bo(1, 3, :, 2) + 1)
         fft_scratch_size%nmax = nmax
         fft_scratch_size%nmray = MAXVAL(nyzray)
         fft_scratch_size%nyzray = nyzray(g_pos)
         fft_scratch_size%rs_group = rs_group
         fft_scratch_size%g_pos = g_pos
         fft_scratch_size%r_pos = r_pos
         fft_scratch_size%r_dim = r_dim
         fft_scratch_size%numtask = numtask

         IF (r_dim(2) > 1) THEN
            !
            ! real space is distributed over x and y coordinate
            ! we have two stages of communication
            !
            IF (r_dim(1) == 1) &
               CPABORT("This processor distribution is not supported.")

            CALL get_fft_scratch(fft_scratch, tf_type=300, n=n, fft_sizes=fft_scratch_size)

            ! assign buffers
            pbuf => fft_scratch%p7buf
            qbuf => fft_scratch%p4buf
            rbuf => fft_scratch%p3buf
            sbuf => fft_scratch%p2buf

            ! FFT along x
            CALL pw_gpu_sf(pw1, pbuf, scale)

            ! Exchange data ( transpose of matrix ) and sort
            IF (pw1%pw_grid%grid_span /= FULLSPACE) qbuf = z_zero
            CALL yz_to_xz(pbuf, rs_group, r_dim, g_pos, p2p, pw1%pw_grid%para%yzp, nyzray, &
                          bo(:, :, :, 2), qbuf, fft_scratch)

            ! FFT along y
            ! use the inbuild fft-lib
            ! CALL fft_1dm(fft_scratch%fft_plan(5), qbuf, rbuf, 1.0_dp, stat)
            ! or cufft (works faster, but is only faster if plans are stored)
            CALL pw_gpu_f(qbuf, rbuf, -1, n(2), mx2*mz2)

            ! Exchange data ( transpose of matrix )
            IF (pw1%pw_grid%grid_span /= FULLSPACE) sbuf = z_zero

            CALL cube_transpose_1(rbuf, bo(:, :, :, 2), bo(:, :, :, 1), sbuf, fft_scratch)

            ! FFT along z
            CALL pw_gpu_fc(sbuf, pw2)

            CALL release_fft_scratch(fft_scratch)

         ELSE
            !
            ! real space is only distributed over x coordinate
            ! we have one stage of communication, after the transform of
            ! direction x
            !

            CALL get_fft_scratch(fft_scratch, tf_type=200, n=n, fft_sizes=fft_scratch_size)

            ! assign buffers
            sbuf => fft_scratch%r1buf
            tbuf => fft_scratch%tbuf

            ! FFT along x
            CALL pw_gpu_sf(pw1, sbuf, scale)

            ! Exchange data ( transpose of matrix ) and sort
            IF (pw1%pw_grid%grid_span /= FULLSPACE) tbuf = z_zero
            CALL x_to_yz(sbuf, rs_group, g_pos, p2p, pw1%pw_grid%para%yzp, nyzray, &
                         bo(:, :, :, 2), tbuf, fft_scratch)

            ! FFT along y and z
            CALL pw_gpu_ffc(tbuf, pw2)

            CALL release_fft_scratch(fft_scratch)

         END IF

         DEALLOCATE (p2p)

!--------------------------------------------------------------------------
      ELSE
         CPABORT("Not implemented (no ray_distr.) in: pw_gpu_c1dr3d_3d_ps.")
      END IF

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_c1dr3d_3d_ps

! **************************************************************************************************
!> \brief perform a parallel real_to_complex copy followed by a 2D-FFT on the gpu
!> \param pw1 ...
!> \param pwbuf ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_cff(pw1, pwbuf)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw1
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT), TARGET                           :: pwbuf

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_cff'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwout
      INTEGER                                            :: handle, l1, l2, l3
      INTEGER, DIMENSION(:), POINTER                     :: npts
      REAL(KIND=dp), POINTER                             :: ptr_pwin
      INTERFACE
         SUBROUTINE pw_gpu_cff_c(din, zout, npts) BIND(C, name="pw_gpu_cff")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: din
            TYPE(C_PTR), VALUE                           :: zout
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
         END SUBROUTINE pw_gpu_cff_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ! dimensions
      npts => pw1%pw_grid%npts_local
      l1 = LBOUND(pw1%array, 1)
      l2 = LBOUND(pw1%array, 2)
      l3 = LBOUND(pw1%array, 3)

      ! pointers to data arrays
      ptr_pwin => pw1%array(l1, l2, l3)
      ptr_pwout => pwbuf(1, 1, 1)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_cff_c(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)
#else
      CPABORT("Compiled without pw offloading")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_cff

! **************************************************************************************************
!> \brief perform a parallel 2D-FFT followed by a complex_to_real copy on the gpu
!> \param pwbuf ...
!> \param pw2 ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_ffc(pwbuf, pw2)
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         TARGET                                          :: pwbuf
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw2

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_ffc'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin
      INTEGER                                            :: handle, l1, l2, l3
      INTEGER, DIMENSION(:), POINTER                     :: npts
      REAL(KIND=dp), POINTER                             :: ptr_pwout
      INTERFACE
         SUBROUTINE pw_gpu_ffc_c(zin, dout, npts) BIND(C, name="pw_gpu_ffc")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: dout
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
         END SUBROUTINE pw_gpu_ffc_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ! dimensions
      npts => pw2%pw_grid%npts_local
      l1 = LBOUND(pw2%array, 1)
      l2 = LBOUND(pw2%array, 2)
      l3 = LBOUND(pw2%array, 3)

      ! pointers to data arrays
      ptr_pwin => pwbuf(1, 1, 1)
      ptr_pwout => pw2%array(l1, l2, l3)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_ffc_c(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)
#else
      CPABORT("Compiled without pw offloading")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_ffc

! **************************************************************************************************
!> \brief perform a parallel real_to_complex copy followed by a 1D-FFT on the gpu
!> \param pw1 ...
!> \param pwbuf ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_cf(pw1, pwbuf)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw1
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: pwbuf

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_cf'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwout
      INTEGER                                            :: handle, l1, l2, l3
      INTEGER, DIMENSION(:), POINTER                     :: npts
      REAL(KIND=dp), POINTER                             :: ptr_pwin
      INTERFACE
         SUBROUTINE pw_gpu_cf_c(din, zout, npts) BIND(C, name="pw_gpu_cf")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: din
            TYPE(C_PTR), VALUE                           :: zout
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
         END SUBROUTINE pw_gpu_cf_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ! dimensions
      npts => pw1%pw_grid%npts_local
      l1 = LBOUND(pw1%array, 1)
      l2 = LBOUND(pw1%array, 2)
      l3 = LBOUND(pw1%array, 3)

      ! pointers to data arrays
      ptr_pwin => pw1%array(l1, l2, l3)
      ptr_pwout => pwbuf(1, 1)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_cf_c(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)
#else
      CPABORT("Compiled without pw offloading")
#endif
      CALL timestop(handle)
   END SUBROUTINE pw_gpu_cf

! **************************************************************************************************
!> \brief perform a parallel 1D-FFT followed by a complex_to_real copy on the gpu
!> \param pwbuf ...
!> \param pw2 ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_fc(pwbuf, pw2)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         TARGET                                          :: pwbuf
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw2

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_fc'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin
      INTEGER                                            :: handle, l1, l2, l3
      INTEGER, DIMENSION(:), POINTER                     :: npts
      REAL(KIND=dp), POINTER                             :: ptr_pwout
      INTERFACE
         SUBROUTINE pw_gpu_fc_c(zin, dout, npts) BIND(C, name="pw_gpu_fc")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: dout
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
         END SUBROUTINE pw_gpu_fc_c
      END INTERFACE

      CALL timeset(routineN, handle)

      npts => pw2%pw_grid%npts_local
      l1 = LBOUND(pw2%array, 1)
      l2 = LBOUND(pw2%array, 2)
      l3 = LBOUND(pw2%array, 3)

      ! pointers to data arrays
      ptr_pwin => pwbuf(1, 1)
      ptr_pwout => pw2%array(l1, l2, l3)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_fc_c(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)
#else
      CPABORT("Compiled without pw offloading")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_fc

! **************************************************************************************************
!> \brief perform a parallel 1D-FFT on the gpu
!> \param pwbuf1 ...
!> \param pwbuf2 ...
!> \param dir ...
!> \param n ...
!> \param m ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_f(pwbuf1, pwbuf2, dir, n, m)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         TARGET                                          :: pwbuf1
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: pwbuf2
      INTEGER, INTENT(IN)                                :: dir, n, m

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_f'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin, ptr_pwout
      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE pw_gpu_f_c(zin, zout, dir, n, m) BIND(C, name="pw_gpu_f")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: zout
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: dir, n, m
         END SUBROUTINE pw_gpu_f_c
      END INTERFACE

      CALL timeset(routineN, handle)

      IF (n*m /= 0) THEN
         ! pointers to data arrays
         ptr_pwin => pwbuf1(1, 1)
         ptr_pwout => pwbuf2(1, 1)

         ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
         CALL pw_gpu_f_c(c_loc(ptr_pwin), c_loc(ptr_pwout), dir, n, m)
#else
         MARK_USED(dir)
         CPABORT("Compiled without pw offloading")
#endif
      END IF

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_f
! **************************************************************************************************
!> \brief perform a parallel 1D-FFT followed by a gather on the gpu
!> \param pwbuf ...
!> \param pw2 ...
!> \param scale ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_fg(pwbuf, pw2, scale)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         TARGET                                          :: pwbuf
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: pw2
      REAL(KIND=dp), INTENT(IN)                          :: scale

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_fg'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin, ptr_pwout
      INTEGER                                            :: handle, mg, mmax, ngpts
      INTEGER, DIMENSION(:), POINTER                     :: npts
      INTEGER, POINTER                                   :: ptr_ghatmap
      INTERFACE
         SUBROUTINE pw_gpu_fg_c(zin, zout, ghatmap, npts, mmax, ngpts, scale) BIND(C, name="pw_gpu_fg")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: zout
            TYPE(C_PTR), INTENT(IN), VALUE               :: ghatmap
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: mmax, ngpts
            REAL(KIND=C_DOUBLE), INTENT(IN), VALUE       :: scale

         END SUBROUTINE pw_gpu_fg_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ngpts = SIZE(pw2%pw_grid%gsq)
      npts => pw2%pw_grid%npts

      IF ((npts(1) /= 0) .AND. (ngpts /= 0)) THEN
         mg = SIZE(pw2%pw_grid%grays, 2)
         mmax = MAX(mg, 1)

         ! pointers to data arrays
         ptr_pwin => pwbuf(1, 1)
         ptr_pwout => pw2%array(1)

         ! pointer to map array
         ptr_ghatmap => pw2%pw_grid%g_hatmap(1, 1)

         ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
         CALL pw_gpu_fg_c(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, mmax, ngpts, scale)
#else
         MARK_USED(scale)
         CPABORT("Compiled without pw offloading")
#endif
      END IF

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_fg

! **************************************************************************************************
!> \brief perform a parallel scatter followed by a 1D-FFT on the gpu
!> \param pw1 ...
!> \param pwbuf ...
!> \param scale ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_sf(pw1, pwbuf, scale)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: pw1
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: pwbuf
      REAL(KIND=dp), INTENT(IN)                          :: scale

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_sf'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin, ptr_pwout
      INTEGER                                            :: handle, mg, mmax, ngpts, nmaps
      INTEGER, DIMENSION(:), POINTER                     :: npts
      INTEGER, POINTER                                   :: ptr_ghatmap
      INTERFACE
         SUBROUTINE pw_gpu_sf_c(zin, zout, ghatmap, npts, mmax, ngpts, nmaps, scale) BIND(C, name="pw_gpu_sf")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: zout
            TYPE(C_PTR), INTENT(IN), VALUE               :: ghatmap
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: mmax, ngpts, nmaps
            REAL(KIND=C_DOUBLE), INTENT(IN), VALUE       :: scale

         END SUBROUTINE pw_gpu_sf_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ngpts = SIZE(pw1%pw_grid%gsq)
      npts => pw1%pw_grid%npts

      IF ((npts(1) /= 0) .AND. (ngpts /= 0)) THEN
         mg = SIZE(pw1%pw_grid%grays, 2)
         mmax = MAX(mg, 1)

         ! pointers to data arrays
         ptr_pwin => pw1%array(1)
         ptr_pwout => pwbuf(1, 1)

         ! pointer to map array
         nmaps = SIZE(pw1%pw_grid%g_hatmap, 2)
         ptr_ghatmap => pw1%pw_grid%g_hatmap(1, 1)

         ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
         CALL pw_gpu_sf_c(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, mmax, ngpts, nmaps, scale)
#else
         MARK_USED(scale)
         CPABORT("Compiled without pw offloading")
#endif
      END IF

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_sf

END MODULE pw_gpu

