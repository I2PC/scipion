      program diag_arpack
c
c "NMA calculation based on ARPACK library by Florence Tama & Osamu Miyashita (AICS).
c
c Lehoucq, R.B., Sorensen, D.C., Yang, C. (1998). ARPACK Users Guide: Solution of Large-Scale Eigenvalue Problems with Implicitly Restarted Arnoldi Methods. Philadelphia: SIAM. ISBN 978-0-89871-407-4."

c
c-----------------------------------------------------------------------
c
c     %------------------------------------------------------%
c     | Storage Declarations:                                |
c     |                                                      |
c     | The maximum dimensions for all arrays are            |
c     | set here to accommodate a problem size of            |
c     | N .le. MAXN                                          |
c     |                                                      |
c     | NEV is the number of eigenvalues requested.          |
c     |     See specifications for ARPACK usage below.       |
c     |                                                      |
c     | NCV is the largest number of basis vectors that will |
c     |     be used in the Implicitly Restarted Arnoldi      |
c     |     Process.  Work per major iteration is            |
c     |     proportional to N*NCV*NCV.                       |
c     |                                                      |
c     | You must set:                                        |
c     |                                                      |
c     | MAXNZ:  Maximum number of nonzero elements           |
c     | MAXN:   Maximum dimension of the A allowed.          |
c     | MAXNEV: Maximum NEV allowed.                         |
c     | MAXNCV: Maximum NCV allowed.                         |
c     %------------------------------------------------------%
c
      integer          maxnz, maxn, maxnev, maxncv, ldv
      parameter       (maxnz=90000000,maxn=50000,maxnev=200,maxncv=1000, 
     $                 ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      Double precision
     &                 v(ldv,maxncv), workl(maxncv*(maxncv+8)),
     &                 workd(3*maxn), d(maxncv,2), resid(maxn),
     &                 ax(maxn)
      logical          select(maxncv)
      integer          iparam(11), ipntr(11)

      integer          n_nz
      integer          element_i(maxnz), element_j(maxnz)
      double precision element_s(maxnz)
c     
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, ierr,
     &                 j, ishfts, maxitr, mode1, nconv
      logical          rvec
      Double precision      
     &                 tol, sigma
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                 zero
      parameter        (zero = 0.0D+0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision           
     &                 dnrm2
      external         dnrm2, daxpy
c
c     %-----------------------------%
c     | routines used |
c     %-----------------------------%
c
      Double precision matrix_tr
      integer          matrix_dim
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic        abs
      intrinsic        cpu_time
      real*8           sec

      namelist /inputs/ nev, ncv, tol
c     
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting msaupd = 1.                    |
c     %-------------------------------------------------%
c
      include 'debug.h'
      ndigit = -3
      logfil = 6
      msgets = 0
      msaitr = 0 
      msapps = 0
      msaupd = 1
      msaup2 = 0
      mseigt = 0
      mseupd = 0
c     
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
c
      call read_matrix(n_nz, element_i,element_j,element_s)
      print *, 'number of data record', n_nz
      n = matrix_dim(n_nz, element_i, element_j, element_s)
      print *, 'dimension of matrix  ', n
      print *, 'trace of the matrix  ',
     $     matrix_tr(n_nz, element_i, element_j, element_s)
      call cpu_time(sec)
      write(*,*) 'cpu time for reading the matrix', sec

      
c
c     %-----------------------------------------------%
c     |                                               | 
c     | Specifications for ARPACK usage are set       | 
c     | below:                                        |
c     |                                               |
c     |    1) NEV = 4  asks for 4 eigenvalues to be   |  
c     |       computed.                               | 
c     |                                               |
c     |    2) NCV = 20 sets the length of the Arnoldi |
c     |       factorization                           |
c     |                                               |
c     |    3) This is a standard problem              |
c     |         (indicated by bmat  = 'I')            |
c     |                                               |
c     |    4) Ask for the NEV eigenvalues of          |
c     |       largest magnitude                       |
c     |         (indicated by which = 'LM')           |
c     |       See documentation in DSAUPD for the     |
c     |       other options SM, LA, SA, LI, SI.       | 
c     |                                               |
c     | Note: NEV and NCV must satisfy the following  |
c     | conditions:                                   |
c     |              NEV <= MAXNEV                    |
c     |          NEV + 1 <= NCV <= MAXNCV             |
c     %-----------------------------------------------%
c
c     %-----------------------------------------------%
c     |               default parameters              |
c     %-----------------------------------------------%
      nev  =  56
      ncv  = 100
      tol  = 0.0

c     %-----------------------------------------------%
c     |               read parameters                 |
c     %-----------------------------------------------%
      open(unit=20,file='diag_arpack.in', status='OLD',form='FORMATTED')
      read(unit=20,nml=inputs)
      write(*,*) 'number of eigenvalues to be computed: nev = ', nev
      write(*,*) 'length of the Arnoldi factorization:  ncv = ', ncv
      write(*,*) 'the stopping criterion:               tol = ', tol


      bmat  = 'I'
      which = 'SM'
c
c     The following condition was put by Ruben Nogales @ CNB to avoid
c     the "Error with _saupd, info = -3" problem  
c      if ( ncv .gt. n) then
c          ncv = n
c      end if
c
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SSIMP: N=', n, ' is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SSIMP: NEV=', nev,
     $        ' is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SSIMP: NCV=', ncv,
     $        ' is greater than MAXNCV '
         go to 9000
      end if

c
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling DSAUPD                    |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from DSAUPD. (See usage below.)                |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to DSAUPD.                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID.)  | 
c     |                                                     |
c     | The work array WORKL is used in DSAUPD as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl = ncv*(ncv+8)

      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting PARAM(1) = 1).              |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | DSAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300 
      mode1  = 1
c
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode1
c
c     %------------------------------------------------%
c     | M A I N   L O O P (Reverse communication loop) |
c     %------------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DSAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %--------------------------------------%
c           | Perform matrix vector multiplication |
c           |              y <--- OP*x             |
c           | The user should supply his/her own   |
c           | matrix vector multiplication routine |
c           | here that takes workd(ipntr(1)) as   |
c           | the input, and return the result to  |
c           | workd(ipntr(2)).                     |
c           %--------------------------------------%
c
            call av(n, workd(ipntr(1)), workd(ipntr(2)),
     $           n_nz, element_i, element_j, element_s)
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if 
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        |                                           |
c        | The routine DSEUPD now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1.)                                   |
c        |                                           |
c        %-------------------------------------------%
c           
          rvec = .true.
c
          call dseupd ( rvec, 'All', select, d, v, ldv, sigma, 
     &         bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &         iparam, ipntr, workd, workl, lworkl, ierr )
c
c         %----------------------------------------------%
c         | Eigenvalues are returned in the first column |
c         | of the two dimensional array D and the       |
c         | corresponding eigenvectors are returned in   |
c         | the first NCONV (=IPARAM(5)) columns of the  |
c         | two dimensional array V if requested.        |
c         | Otherwise, an orthogonal basis for the       |
c         | invariant subspace corresponding to the      |
c         | eigenvalues in D is returned in V.           |
c         %----------------------------------------------%
c
          if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of DSEUPD. |
c            %------------------------------------%
c
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
c
          else
c
             nconv =  iparam(5)
             do 20 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%
c
                call av(n, v(1,j), ax,
     $               n_nz, element_i, element_j, element_s)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
c
 20          continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
             call dmout(6, nconv, 2, d, maxncv, -6,
     &            'Ritz values and relative residuals')

c            %-----------------------------%
c            | Output vector to UNIT 11    |
c            %-----------------------------%

             do j=1, nconv
                do i=1, n
                write(11,*) v(i,j)
                end do
             end do
          end if
c
c         %-------------------------------------------%
c         | Print additional convergence information. |
c         %-------------------------------------------%
c
          if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
          else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit',
     &                ' Arnoldi update, try increasing NCV.'
             print *, ' '
          end if      
c
          print *, ' '
          print *, ' _SSIMP '
          print *, ' ====== '
          print *, ' '
          print *, ' Size of the matrix is ', n
          print *, ' The number of Ritz values requested is ', nev
          print *, ' The number of Arnoldi vectors generated',
     &             ' (NCV) is ', ncv
          print *, ' What portion of the spectrum: ', which
          print *, ' The number of converged Ritz values is ', 
     &               nconv 
          print *, ' The number of Implicit Arnoldi update',
     &             ' iterations taken is ', iparam(3)
          print *, ' The number of OP*x is ', iparam(9)
          print *, ' The convergence criterion is ', tol
          print *, ' '
c
      end if
c
c     %---------------------------%
c     | Done with program dssimp. |
c     %---------------------------%
c
 9000 continue
c
      end
c 
c ------------------------------------------------------------------
c     matrix vector subroutine
c
c     Computes w <--- OP*v, where OP is the n*n by n*n block 
c     tridiagonal matrix
c
      subroutine av(n, v, w, n_nz, ei, ej, es)
      integer          n, n_nz
      integer          ei(n_nz), ej(n_nz)
      double precision es(n_nz)
      double precision v(n), w(n)
c     
      integer k, i, j
      
      do k = 1, n
         w(k)=0
      enddo
c$$$      print *,'w', w(1:n)
c$$$      print *,'v', v(1:n)
      do k = 1, n_nz
         i = ei(k)
         j = ej(k)
         if( i .eq. j ) then
            w(i) = w(i)+es(k)*v(i)
         else
            w(i) = w(i)+es(k)*v(j)
            w(j) = w(j)+es(k)*v(i)
         endif
      enddo

      return
      end
c
c-------------------------------------------------------------------
      subroutine read_matrix(n_nz, ei, ej, es)
c
      integer         n_nz, ei(*), ej(*)
      double precision            es(*)
      
      integer         m

c      open(UNIT=10,FILE='matrice.sdij',STATUS='OLD',FORM='UNFORMATTED')
       open(UNIT=10,FILE='matrice.sdijf',STATUS='OLD',FORM='FORMATTED')

      m=1
      do
         read(10,*,END=900) ei(m), ej(m), es(m)
c         read(10,END=900) ei(m), ej(m), es(m)
         m=m+1
      end do
 900  continue
c     finish to read data
      
      n_nz= m
      return
      end


      
      function matrix_dim( n_nz, ei, ej, es )
      integer          matrix_dim
      integer          n_nz, ei(*), ej(*)
      double precision es(*)
      
      integer          m, max

      max=0
      do m=1, n_nz
         if( max .lt. ei(m) ) then
            max = ei(m)
         endif
         if( max .lt. ej(m) ) then
            max = ej(m)
         endif
      enddo
      matrix_dim = max

      end

      function matrix_tr( n_nz, ei, ej, es )
      double precision matrix_tr
      integer          n_nz, ei(*), ej(*)
      double precision es(*)
      
      integer          m, max
      double precision tr

      tr=0
      do m=1, n_nz
         if( ei(m) .eq. ej(m) ) then
            tr = tr + es(m)
         endif
      enddo
      matrix_tr = tr

      end

      
