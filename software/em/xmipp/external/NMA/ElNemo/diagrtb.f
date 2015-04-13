      Program diagrtb
c======================================================================

c     DIAGonalisation of a matrix, using the RTB approximation.

c     The method rests upon the approximation that the
c     required matrix eigenvectors are well described by 
c     rigid body motions (Rotations and Translations)
c     of Blocks of particles buildt with sets of consecutive
c     ones, as found in the corresponding coordinate file.
c     (See references below for more details)

c     In the case of proteins, typically, one amino-acid is put in one 
c     block. 

c======================================================================

c     INPUT: 
c     ******
c     A parameter, or command, file, called diagrtb.dat.  

c     Note that each run of diagrtb produces a diagrtb.dat_run file,
c     where all parameter values are given (and commented).
c     diagrtb.dat_run can be used as a diagrtb.dat file, 
c     in further runs. 

c     Among the parameters: 

c    *The name of a file with the coordinates of the system, 
c     in FREE or PDB (Protein Data Bank) format.

c     PDB  format: Masses can be given in the Bfactors column.
c     Free format: x, y, z, mass, block-number.

c    *The name of a file with the matrix, assumed to be:
c     a) real, square and symmetrical.
c     b) in the following format: i, j, i-j-matrix-element

c     OUTPUT: 
c     *******
c     The eigenvalues and eigenvectors of the matrix, 
c     in x,y,z "CERFACS" format.
 
c     In the present state of the program, the matrix
c     is first re-written in a block-by-block form.

c     ATTENTION: 
c     **********
c     Temporary diagrtb_work.* files are created during each 
c     program run.

c     So, avoid running several diagrtb simultaneously, in the same 
c     directory !

c......................................................................

c     MEMORY LIMITS:
c     **************

c     The computer memory required (LIWORK and LRWORK) can be large, 
c     in particular when the size of the largest block is large.
     
c     NATMAX    = maximum number of particles in the system.
c     LLWORK    = size of the LOGICAL working array.

c     LIWORK    = size of the INTEGER working array.
c     LRWORK    = size of the DOUBLE PRECISION working array.
 
      implicit         none
      integer          liwork, llwork, lrwork, natmax
 
      parameter        (natmax =   55000)
      parameter        (llwork =3*natmax)

c     ===================================
c      parameter        (LIWORK = 5000000)
c      parameter        (LRWORK =32000000)
      parameter        (LIWORK = 100000000)
      parameter        (LRWORK =200000000)
c     ====================================
c......................................................................

c     REFERENCES: 
c     ***********

c     1) P. Durand, G. Trinquier, Y.H. Sanejouand (1994):
c    "A new approach for determining low-frequency
c     normal modes in macromolecules", 
c     Biopolymers vol.34, p759-771.
 
c     2) F. Tama, F.X. Gadea, O. Marques, Y.H. Sanejouand (2000):
c    "Building-block approach for determining 
c     low-frequency normal modes of macromolecules", 
c     Proteins: Structure, Function, and Genetics vol.41(1), p1-7.

c     3) G. Li, Q. Cui (2002):
c    "A coarse-grained normal mode approach for macromolecules: an
c     efficient implementation and application to Ca(2+)-ATPase",
c     Biophys. J. vol.83, p2457-2474.

c     In the later paper, RTB is called BNM (for Block-Normal-Modes). 
c     It is more general (flexible blocks are allowed). 
c     It is now available in the CHARMm package (since version 32). 

c......................................................................

c     1994: version 1.0, by Georges Trinquier, Philippe Durand
c           and Yves-Henri Sanejouand.
c     2001: version 2.0, by Florence Tama, Florent-Xavier Gadea, 
c           Osni Marques and Yves-Henri Sanejouand.

c     Further modifications (2002-2004): by Yves-Henri Sanejouand.

c     In case of problem, feel free to contact: 
c     Yves-Henri.Sanejouand@ens-lyon.fr

c......................................................................

      logical          LWORK(llwork)
      integer          IWORK(liwork)
      double precision RWORK(lrwork)
 
c     nbrt: Number of degrees of freedom of each block.
c          (3 translations and 3 rotations)

      integer   nbrt, nmotsmax, nresmax      
      parameter (nbrt=6,nresmax=natmax,nmotsmax=100)

      character cformat*12, clean*4, cstatus*12, eige*4, feig*128, 
     .          field*128, fmtx*128, fpdb*128, macro*4, mass*4, 
     .          namfich*128, mots(nmotsmax)*128, motinp*128, program*9, 
     .          progrer*12, progrwn*12, resname(natmax)*4, string*128, 
     .          ssu(natmax)*1, ssusel*1, sstr*4, typf*4, version*30
      logical   alive, qerror, qexist, qinter, qstop, qverbos
      integer   i, iiwkm1, iiwkm2, iiwkove, iiwkp1, iiwkp2, iiwkp3, 
     .          iiwrk1, iiwrk2, iiwrk3, irwkd1, irwkd2, irwkd3, irwkd4, 
     .          irwkm1, irwkm2, irwkm3, irwkm4, irwkm5, iiwkmx, irwkove, 
     .          irwkp1, irwkp2, irwkmx, irwrk1, irwrk2, irwrk3, irwrk4, 
     .          irwrk5, irwrk6, irwrk7, irwrk8, irwrk9, k, klist, ksep, 
     .          lmot, lnomeig, lnommtx, lnompdb, natblocs, natom, nb, 
     .          nddres, nddres2, ndim, nmots, nrbl, nunit, nvec, prtlev, 
     .          undat, uninp
 
c-----------------------------------------------------------------------
      version=' Version 2.52, November 2004.'
      program=' Diagrtb>'
 
      write(6,'(2A)') program,
     .    ' Diagonalizes a matrix, using the RTB/BNM approximation.'
      write(6,'(2A)') program, version

      progrwn='%Diagrtb-Wn>'
      progrer='%Diagrtb-Er>'

c.... Defaults .........................................................
c     Can be modified in 'diagrtb.dat'.
 
c     Number of residues per block:
      nrbl = 1
      
c     Substructuring:
c     --------------
c    (Chain breaks interrupt blocks in RESIdue case)

c     SECO: Whole secondary structures are put in blocks
c     DOMA: Whole domains are put in blocks
c     SUBU: Whole subunits are put in blocks
c     RESI: Only the number of residues per block is taken into account
c     NONE:       "                       "                   "
 
c     Structures are defined in the coordinate file, by the
c     HELIX, SHEET, LOOP or DOMAIN keywords 
c    (the two laters are not standard).
c     Non-specified portions of the sequence are split
c     into blocks, according to the NRBL value.
      
      sstr = 'NONE'

c     Atomic masses:
c     PDB : Look for masses in the B-factors column of the pdb file.
c     CONS: Set them to a constant value (mass=1).
      mass = 'CONS'
 
c     Number of eigenvectors to compute:
      nvec = 56
      
c     LOWE: The lowest-eigenvalues ones.
c     HIGH: The highest-eigenvalues ones.
      eige = 'LOWE'
 
c     ALL: Work files are deleted at the end of the job.
c     NO : They are kept.
      clean= 'ALL'
 
c     PDB default filename with the coordinates of the system.
      fpdb = 'structure.pdb'
 
c     Matrix default filename with the non-zero elements.
      fmtx = 'pdbmat.sdijf'
      typf = 'FREE'
 
c     Eigenvector output filename.
      feig = 'diagrtb.eigenfacs'

c     PRINTing level:
      prtlev=0

C.... read options from file diagrtb.dat ..............................

      nunit=10

      namfich='diagrtb.dat'
      undat=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"

      qverbos=.false.
      if (prtlev.gt.1) qverbos=.true.

      call openam(namfich,cformat,cstatus,undat,qverbos,
     .            qinter,qexist)

      if (.not.qexist) then
          write(6,'(2A)') progrwn,
     .  ' No diagrtb.dat file. Defaults assumed for all options.'
          goto 15
      else
          write(6,'(2A)') program,
     .  ' Options to be read in diagrtb.dat file.'
      endif
 
 5    continue 
      read (undat,end=10,err=20,fmt='(a)') string

      if (string(1:1).eq.'!'.or.string(1:1).eq.'*'.or.
     .    string(1:1).eq.'#') goto 5

      ksep=index(string,'=')

      motinp=' '
      if (ksep.gt.0) then 
          motinp=string(1:ksep)
          field=string(ksep+1:80)
          k=index(field,'!')
          if (k.gt.0) field=field(1:k-1)
          call stringcl(field,lmot)
          klist=-1
          if (field(1:1).eq.'?') klist=1
      else
          write(6,'(/2A/A)') progrwn,
     .  ' No separator (=) in command ligne:',
     .    string
          write(6,'(2A)') progrwn,' This ligne is skipped.'
          goto 5
      endif
      call mintomaj(motinp)

      macro='XXXX'
      if (index(motinp,'NRBL').gt.0.or.
     .    index(motinp,'BLOC').gt.0) macro='NRBL'
      if (index(motinp,'NVEC').gt.0.or.
     .    index(motinp,' VECT').gt.0) macro='NVEC'
      if (index(motinp,'FPDB').gt.0.or.
     .    index(motinp,'COOR').gt.0) macro='FPDB'
      if (index(motinp,'FMTX').gt.0.or.
     .    index(motinp,'MATR').gt.0) macro='FMTX'
      if (index(motinp,'FEIG').gt.0.or.
     .    index(motinp,' OUTP').gt.0) macro='FEIG'
      if (index(motinp,'MASS').gt.0) macro='MASS'
      if (index(motinp,'SUBS').gt.0) macro='SUBS'
      if (index(motinp,'NVAL').gt.0.or.
     .    index(motinp,' FREQ').gt.0) macro='EIGE'
      if (index(motinp,'CLEA').gt.0) macro='CLEA'
      if (index(motinp,'PRIN').gt.0) macro='PRTL'
      if (index(motinp,'FORM').gt.0.and.index(motinp,'MATR').gt.0)
     .    macro='FORM'

      if (macro.ne.'FMTX'.and.macro.ne.'FPDB'.and.
     .    macro.ne.'FEIG') call mintomaj(field)

      if      ( macro.eq.'NRBL' ) then
                read (field,'(i16)') nrbl 
      else if ( macro.eq.'NVEC' ) then
                read (field,'(i16)') nvec 
      else if ( macro.eq.'PRTL' ) then
                read (field,'(i16)') prtlev
      else if ( macro.eq.'FPDB' ) then
                read (field,'(a)') fpdb 
      else if ( macro.eq.'FMTX' ) then
                read (field,'(a)') fmtx 
      else if ( macro.eq.'FORM' ) then
                read (field,'(a)') typf 
                call mintomaj(typf)
                if (typf.ne.'FREE'.and.typf.ne.'BINA') then
                    write(6,'(3A)') program,
     .            ' Matrix FORMat: ',typf 
                    if (klist.le.0) write(6,'(2A)') progrwn,
     .                ' This is not a known keyword.'
                    write(6,'(2A)') program,
     .            ' Valid options are: FREE, BINAry.'
                    write(6,'(2A)') program,' Default assumed (FREE).'
                    typf='FREE'
                endif
      else if ( macro.eq.'FEIG' ) then
                read (field,'(a)') feig 
      else if ( macro.eq.'MASS' ) then
                read (field,'(a)') mass 
                call mintomaj(mass)
                if (mass.ne.'PDB'.and.mass.ne.'CONS'.and.
     .              mass.ne.'COOR') then
                    write(6,'(3A)') program,
     .            ' Origin of MASS values: ',mass 
                    if (klist.le.0) write(6,'(2A)') progrwn,
     .                ' This is not a known keyword.'
                    write(6,'(2A)') program,
     .            ' Valid options are: CONStant, COORdinate, PDB.'
                    write(6,'(2A)') program,' Default assumed (CONS).'
                    mass='CONS'
                endif
      else if ( macro.eq.'SUBS' ) then
                read (field,'(a)') sstr 
                call mintomaj(sstr)
                if (sstr.eq.'RESI') sstr='NONE'
                if (sstr.ne.'SECO'.and.sstr.ne.'SUBU'.and.
     .              sstr.ne.'DOMA'.and.sstr.ne.'NONE') then
                    write(6,'(3A)') program,
     .            ' Type of SUBStructures: ',sstr 
                    if (klist.le.0) write(6,'(2A)') progrwn,
     .                ' This is not a known keyword.'
                    write(6,'(2A)') program,
     .            ' Valid options are:'//
     .            ' RESIdues, SECOndary, SUBUnits, DOMAins, NONE.'
                    write(6,'(2A)') program,' Default assumed (NONE).'
                    sstr='NONE'
                endif
      else if ( macro.eq.'EIGE' ) then
                read (field,'(a)') eige 
                call mintomaj(eige)
                if (eige.ne.'LOWE'.and.eige.ne.'HIGH') then
                    write(6,'(3A)') program,
     .            ' Selected EIGEnvalues : ',eige 
                    if (klist.le.0) write(6,'(2A)') progrwn,
     .                ' This is not a known keyword.'
                    write(6,'(2A)') program,
     .            ' Valid options are: LOWEst, HIGHest.'
                    write(6,'(2A)') program,' Default assumed (LOWE).'
                    eige='LOWE'
                endif
      else if ( macro.eq.'CLEA' ) then
                read (field,'(a)') clean
                call mintomaj(clean)
                if (clean.ne.'ALL'.and.clean(1:2).ne.'NO') then
                    write(6,'(3A)') program,
     .            ' Temp. files CLEAning : ',clean
                    if (klist.le.0) write(6,'(2A)') progrwn,
     .                ' This is not a known keyword.'
                    write(6,'(2A)') program,
     .            ' Valid options are: ALL, NOne.'
                    write(6,'(2A)') program,' Default assumed (ALL).'
                    clean='ALL'
                endif
      else    
                write(6,'(/2A/A)') progrwn,
     .        ' No known keyword in ligne:',string
                write(6,'(2A)') progrwn,
     .        ' This command ligne is skipped.'
      end if
 
      goto 5
   10 close (10)
   15 continue
      close(undat)
 
      call stringcl(fpdb,lnompdb)
      call stringcl(fmtx,lnommtx)
      call stringcl(feig,lnomeig)

      if (prtlev.le.0) then
      write(6,'(/2A)')  program,' Main options:'
      else
      write(6,'(/2A)')  program,' Options taken into account:'
      endif

      write(6,'(/2A)')   ' MATRix filename          = ',fmtx(1:lnommtx)
      write(6,'(2A)')   ' COORdinates filename     = ',fpdb(1:lnompdb)
      write(6,'(2A)')   ' Eigenvector OUTPut file  = ',feig(1:lnomeig)
      write(6,'(A,I10)') ' Nb of VECTors required   = ',nvec
      write(6,'(A,6X,A)') ' EigeNVALues chosen       = ',eige
      if (sstr.ne.'NONE'.or.prtlev.gt.0)
     .write(6,'(A,6X,A)')   ' Type of SUBStructuring   = ',sstr
      write(6,'(A,I10)') ' Nb of residues per BLOck = ',nrbl
      write(6,'(A,6X,A)')   ' Origin of MASS values    = ',mass
      if (typf.ne.'FREE'.or.prtlev.gt.0)
     .write(6,'(A,6X,A)')   ' MATRix FORMat            = ',typf 
      if (prtlev.gt.0) then
      write(6,'(A,6X,A)')   ' Temporary files cleaning = ',clean
      write(6,'(A,I10)') ' Output PRINting level    = ',prtlev
      endif

c     Sauvegarde du fichier de commandes complet:

      uninp=nunit
      nunit=nunit+1
      namfich='diagrtb.dat_run'
      cformat="FORMATTED"
      cstatus="ove"
      call openam(namfich,cformat,cstatus,uninp,qverbos,
     .     qinter,qexist)

      write(uninp,'(2A)') 
     .'! This file can be modified and used as a command file',
     .' (named diagrtb.dat) for diagrtb.'

      write(uninp,'(2A)')   ' MATRix filename            = ',
     .      fmtx(1:lnommtx)
      write(uninp,'(2A)')   ' COORdinates filename       = ',
     .      fpdb(1:lnompdb)
      write(uninp,'(2A)')   ' Eigenvector OUTPut filename= ',
     .      feig(1:lnomeig)
      write(uninp,'(A,I10)') ' Nb of VECTors required     = ',nvec
      write(uninp,'(A,6X,2A)') ' EigeNVALues chosen         = ',eige,
     .            '   ! LOWEst, HIGHest.'
      write(uninp,'(A,6X,2A)')   ' Type of SUBStructuring     = ',sstr,
     .            '   ! RESIdues, SECOndary, SUBUnits, DOMAins, NONE.'
      write(uninp,'(A,I10)') ' Nb of residues per BLOck   = ',nrbl
      write(uninp,'(A,6X,2A)')   ' Origin of MASS values      = ',mass,
     .            '   ! CONStant, COORdinate, PDB.'
      write(uninp,'(A,6X,2A)')   ' Temporary files cleaning   = ',clean,
     .            '   ! ALL, NOne.'
      write(uninp,'(A,6X,2A)')   ' MATRix FORMat              = ',typf, 
     .            '   ! FREE, BINAry.'
      write(uninp,'(A,I10,A)') ' Output PRINting level      = ',prtlev,
     .            '   ! =1: More detailed; =2: Debug level.'
      close(uninp)

c     On recherche l'information/sous-unite:
 
      call string_split(fpdb,lnompdb,":",
     .                  mots,nmotsmax,nmots)
      call stringcl(mots(1),lmot)
 
      if (nmots.gt.1) then
          call stringcl(mots(2),lmot)
          ssusel=mots(2)
          write(6,'(/3A)') program,' Selected (pdb) subunit: ',ssusel
          if (nmots.gt.2.or.lmot.gt.1) then
              write(6,'(4A)') progrwn,
     .      ' The end of name: ',
     .        fpdb(1:lnompdb),', was not understood.'
          endif
      else
          ssusel='*'
      endif
      fpdb=mots(1)
      call stringcl(fpdb,lnompdb)

C     Premiers tests:

      qerror=.false.

      inquire(file=fpdb,exist=qexist)

      if (.not.qexist) then
      write(6,'(/2A)') progrer,' Coordinate file not found. Required !'
      qerror=.true.
      endif

      if (nvec.le.0) then
      write(6,'(/2A)') progrer,' No eigenvector required...'
      qerror=.true.
      endif

      if (qerror) stop '*Command error. Check diagrtb.dat file*'

C.... pointers for BLOCPDB .............................................
 
      iiwkove=-1
      irwkove=-1
      qstop=.false.

      if (prtlev.gt.1)
     .write(6,'(/2A)') program,' Memory allocation for Blocpdb.'

      iiwrk1 = 1
      iiwrk2 = iiwrk1 + nresmax
      iiwrk3 = iiwrk2 + natmax
      iiwkmx = iiwrk3 + natmax
      if ( iiwkmx.gt.liwork ) then
          write(6,'(/A,2(A,I12))') progrer,' IIWKMX= ',IIWKMX,
     .  ' Maximum allowed is LIWORK= ',liwork
          qstop=.true.
      endif
      
      irwrk1 = 1
      irwrk2 = irwrk1 + 3*natmax
      irwkmx = irwrk2 + 3*natmax
      if ( irwkmx.gt.lrwork ) then
          write(6,'(/A,2(A,I12))') progrer,' IRWKMX= ',IRWKMX,
     .  ' Maximum allowed is LRWORK= ',lrwork
          qstop=.true.
      endif

      if (qstop) then
          write(6,'(/2A/A/A)') progrer,
     .  ' Not enough memory allowed for working arrays. Sorry.',
     .  ' Recompile DIAGRTB with smaller NATMAX and/or NRESMAX.'
          stop '*Working arrays allocation error*'
      endif
 
      CALL BLOCPDB ( fpdb, ssusel, sstr, mass, natmax, nrbl, nresmax, 
     .               IWORK(iiwrk2), IWORK(iiwrk3), IWORK(iiwrk1), 
     .               RWORK(irwrk1), RWORK(irwrk2), 
     .               resname, ssu, natblocs, natom, nb, prtlev )

C.... pointers for PREPMAT .............................................
 
C     Note: The first nresmax entries of IWORK (on output from BLOCPDB)
C           are required in PREPMAT, so we don't change iiwrk2.
 
      if (prtlev.gt.1)
     .write(6,'(/2A)') program,' Memory allocation for Prepmat.'
      iiwkp1 = iiwrk1
      iiwkp2 = iiwrk2
      iiwkp3 = iiwkp2 + 9*natblocs**2
      iiwkmx = iiwkp3 + 9*natblocs**2
      if ( iiwkmx.gt.liwork ) then
          if ( iiwkmx.gt.iiwkove ) iiwkove=iiwkmx
          qstop=.true.
      endif
 
      irwkp1 = 1
      irwkp2 = irwkp1 + 9*natblocs*natom
      irwkmx = irwkp2 + 9*natblocs*natblocs
      if ( irwkmx.gt.lrwork ) then
          if ( irwkmx.gt.irwkove ) irwkove=irwkmx
          qstop=.true.
      endif
 
C.... pointers for RTB .................................................
 
      if (prtlev.gt.1)
     .write(6,'(2A)') program,' Memory allocation for RTB.'
      nddres = 3*natblocs
      nddres2 = nddres**2
 
      iiwrk1 = 1
      iiwrk2 = iiwrk1 + nddres2
      iiwrk3 = iiwrk2 + nddres2 
      iiwkmx = iiwrk3 + nb
      if ( iiwkmx.gt.liwork ) then
          if ( iiwkmx.gt.iiwkove ) iiwkove=iiwkmx
          qstop=.true.
      endif
 
      irwrk1 = 1
      irwrk2 = irwrk1 + 3*natom
      irwrk3 = irwrk2 + 3*natom
      irwrk4 = irwrk3 + nbrt*nddres
      irwrk5 = irwrk4 + nddres2
      irwrk6 = irwrk5 + nddres*nbrt 
      irwrk7 = irwrk6 + nddres2
      irwrk8 = irwrk7 + 36*nb**2 
      irwrk9 = irwrk8 + nbrt*nddres*nb
      irwkmx = irwrk9 + nbrt*nbrt
      if ( irwkmx.gt.lrwork ) then
          if ( irwkmx.gt.irwkove ) irwkove=irwkmx
          qstop=.true.
      endif
 
C.... pointers for DIAGSTD .............................................
 
      if (prtlev.gt.1)
     .write(6,'(2A)') program,' Memory allocation for Diagstd.'
      ndim = 6*nb
 
      iiwkmx = ndim
      if ( iiwkmx.gt.liwork ) then
          if ( iiwkmx.gt.iiwkove ) iiwkove=iiwkmx
          qstop=.true.
      endif

      irwkd1 = 1
      irwkd2 = irwkd1 + ndim*ndim
      irwkd3 = irwkd2 + ndim
      irwkd4 = irwkd3 + ndim
      irwkmx = irwkd4 + ndim
      if ( irwkmx.gt.lrwork ) then
          if ( irwkmx.gt.irwkove ) irwkove=irwkmx
          qstop=.true.
      endif
 
      if (ndim.lt.nvec) then
          nvec = ndim
          write(6,'(/A,I10,A)') progrwn,
     .    nvec,' eigenvectors, only, can be determined.'
      endif
 
C.... pointers for RTBTOMODES ..........................................
 
      if (prtlev.gt.1)
     .write(6,'(2A)') program,' Memory allocation for RTB_to_modes.'
      iiwkm1 = 1
      iiwkm2 = iiwkm1 + nb
      iiwkmx = iiwkm2 + nvec
      if ( iiwkmx.gt.liwork ) then
          if ( iiwkmx.gt.iiwkove ) iiwkove=iiwkmx
          qstop=.true.
      endif
 
      irwkm1 = 1
      irwkm2 = irwkm1 + 3*natom*nvec
      irwkm3 = irwkm2 + nvec
      irwkm4 = irwkm3 + 6*nb*nvec
      irwkm5 = irwkm4 + nvec
      irwkmx = irwkm5 + nbrt*nddres*nb
      if ( irwkmx.gt.lrwork ) then
          if ( irwkmx.gt.irwkove ) irwkove=irwkmx
          qstop=.true.
      endif
 
      if (qstop) then
          if ( iiwkove.gt.liwork ) then
          write(6,'(/A,2(A,I12))') progrer,' IIWKMX up to: ',IIWKove,
     .  ' Maximum allowed is LIWORK= ',liwork
          endif
          if ( irwkove.gt.lrwork ) then
          write(6,'(/A,2(A,I12))') progrer,' IRWKMX up to: ',IRWKove,
     .  ' Maximum allowed is LRWORK= ',lrwork
          endif
          write(6,'(/2A/A/A)') progrer,
     .  ' Not enough memory allowed for working arrays. Sorry.',
     .  ' Lower the number of blocks, the sizes of the largest ones...',
     .  ' Or recompile DIAGRTB with larger WORKing arrays.'
          stop '*Working arrays allocation error*'
      endif

      inquire(file=fmtx,exist=qexist)
      if (.not.qexist) then
          write(6,'(/2A)') progrer,' Matrix not found. Required !'
          stop '*Nothing to do*'
      endif

      CALL PREPMAT ( fmtx, typf, natblocs, natom, nb, 
     &               LWORK, IWORK(iiwkp2), IWORK(iiwkp3), IWORK(iiwkp1), 
     &               RWORK(irwkp1), RWORK(irwkp2), prtlev )
 
      CALL RTB     ( natom, nb, nbrt, nddres, nddres2, IWORK(iiwrk1),
     &               IWORK(iiwrk2), IWORK(iiwrk3), RWORK(irwrk1),
     &               RWORK(irwrk2), RWORK(irwrk3), RWORK(irwrk4),
     &               RWORK(irwrk5), RWORK(irwrk6), RWORK(irwrk7),
     &               RWORK(irwrk8), RWORK(irwrk9), prtlev )
 
      CALL DIAGSTD ( ndim, nvec, eige, IWORK, 
     &               RWORK(irwkd1), RWORK(irwkd2),
     &               RWORK(irwkd3), RWORK(irwkd4), prtlev )

      CALL RTBTOMODES ( natom, nb, nbrt, nddres, nvec, IWORK(iiwkm1),
     &                  IWORK(iiwkm2), RWORK(irwkm1), RWORK(irwkm2),
     &                  RWORK(irwkm3), RWORK(irwkm4), RWORK(irwkm5),
     &                  prtlev, feig )
 
C.... clean up .........................................................
 
      if (clean.eq.'ALL') then
      inquire (file='diagrtb_work.xyzm',exist=alive)
      if ( alive ) then
         open  (10,file='diagrtb_work.xyzm')
         close (10,status='delete')
      end if
      inquire (file='diagrtb_work.blocs',exist=alive)
      if ( alive ) then
         open  (10,file='diagrtb_work.blocs')
         close (10,status='delete')
      end if
      inquire (file='diagrtb_work.matblocs',exist=alive)
      if ( alive ) then
         open  (10,file='diagrtb_work.matblocs')
         close (10,status='delete')
      end if
      inquire (file='diagrtb_work.sdijb',exist=alive)
      if ( alive ) then
         open  (10,file='diagrtb_work.sdijb')
         close (10,status='delete')
      end if
      inquire (file='diagrtb_work.eigenfacs',exist=alive)
      if ( alive ) then
         open  (10,file='diagrtb_work.eigenfacs')
         close (10,status='delete')
      end if
      end if
 
C.... End of diagrtb ..................................................
 
      write(6,'(/2A)') program,' Normal end.'
      stop
   20 stop '*Command error. Check diagrtb.dat file*'
      end
c----------------------------------------------------------------------
      SUBROUTINE BLOCPDB ( fpdb, ssusel, sstr, mass, natmax, nbb, 
     &                     nblocmx, resi, numbloc, tabb, amass, corlin, 
     &                     resname, ssu, natbmax, natom, nb, prtlev )
      implicit none

c     Scalar Arguments
 
      character        fpdb*(*), mass*(*), sstr*(*), ssusel*(*)
      integer          natbmax, natmax, natom, nb, nbb, nblocmx, prtlev 
 
c     Array Arguments
 
      integer          numbloc(natmax), resi(natmax), tabb(nblocmx)
      double precision amass(3*natmax), corlin(3*natmax), unknown 
      character        resname(natmax)*4, ssu(natmax)*1
      parameter       (unknown=9999.999999)
 
c     Purpose:
c     =======
 
c     -------------------------------------------------------------
c     Atomic coordinates are read and block lengths are determined,
c     as a consequence of the number of residues per block and/or
c     of the substructuring requirements (SSTR command).
c     -------------------------------------------------------------
c    (The peptidic bond between residues is simply split)
 
c     Arguments
c     =========
 
c     fpdb      : Pdb filename (with the x,y,z coordinates)
c     sstr      : keyword for substructuring
c     mass      : keyword for choosing masses 
c     natmax    : maximum number of atoms in the system.
c     nbb       : number of residues per block.
c     nblocmx   : maximum number of blocks.
c     resi      : residue number each atom belongs to.
c     numbloc   : block number each atom belong to.
c     tabb      : block lengths.
c     amass     : coordinate masses.
c     corlin    : coordinates (x1, y1, z1, x2, y2, ...).
c     resname   : residue name each atom belong to.
c     ssu       : subunit name each atom belong to.
c     natbmax   : maximum number of atoms found in blocks.
c     natom     : number of atoms in the system.
c     nb        : number of blocks the system is split into.
  
c-----------------------------------------------------------------------
 
c     Local variables
 
      character aa1*4, aa2*4, atname*4, cformat*12, cstatus*12, 
     .          lign80*80, namfil*20, program*9, progrer*12, progrwn*12, 
     .          ssu1*1, ssu2*1, ssuch*4
      logical   qerror, qexist, qinterr, qmass, qok, qpdb, qverbos
      integer   at, i, ii, j, k, lnompdb, natbloc, natbmin, ndat, 
     .          nl, noff, nres, nresb, nssu, num1, num2, nunit, 
     .          uninp, unrtb
      double precision rave, rdev, rmax, rmin, tot
 
C-----------------------------------------------------------------------
 
      program=' Blocpdb>'
      progrwn='%Blocpdb-Wn>'
      progrer='%Blocpdb-Er>'

      if (prtlev.gt.1)
     .    write(6,'(/2A)') program,' Entering in.'
 
      if (sstr.eq.'NONE'.and.nbb.le.0) then
          write(6,'(/A,I6,A)') progrer,nbb,
     .  ' residues per bloc required. Not a positive number !'
          stop '*Command error. Check diagrtb.dat file*'
      endif
 
      nunit=10

c     Fichier de sortie
c     -----------------
      namfil='diagrtb_work.xyzm'
      cformat='UNFORMATTED'
      cstatus='ove'
      unrtb=nunit
      nunit=nunit+1

      qverbos=.false.
      if (prtlev.gt.1) qverbos=.true.

      call openam(namfil,cformat,cstatus,unrtb,qverbos,
     .            qinterr,qexist)
      if (qinterr) stop

c     ----------------------------
c     The coordinate file is read:
c     ----------------------------

      cformat='FORMATTED'
      cstatus='OLD'
      uninp=nunit
      nunit=nunit+1

      if (prtlev.gt.1) then
        call stringcl(fpdb,lnompdb)
        write(6,'(/4A)') program,
     .' Coordinate file ',fpdb(1:lnompdb),' to be opened.'
      endif

      call openam(fpdb,cformat,cstatus,uninp,qverbos,
     .            qinterr,qexist)
      if (qinterr) stop

c     Format pdb pour les coordonnees ?

      nl=0
      qpdb=.false.
 120  continue
      read(uninp,'(A)',end=130) lign80
      if (lign80(1:5).eq.'ATOM '.or.lign80(1:6).eq.'HETATM') then
          qpdb=.true.
          goto 130
      else
          nl=nl+1
      endif
      goto 120
 130  continue
      rewind(uninp)

      qmass=.true.
      if (mass.eq.'CONS') qmass=.false.
 
      do i=1,natmax
         ii=3*i-2
         corlin(ii)=unknown
         corlin(ii+1)=unknown
         corlin(ii+2)=unknown
         amass(ii)=unknown
         amass(ii+1)=unknown
         amass(ii+2)=unknown
         resi(i)=i
         ssu(i)='*'
      enddo

      if (.not.qpdb) then
          if (nl.eq.0) then
              write(6,'(/2A)') progrer,' Empty coordinate file.'
              stop '*I/O error*'
          endif
          if (prtlev.gt.1)
     .    write(6,'(2A)') program,
     .  ' Coordinate file in Free format.'

          call readrlin(uninp,corlin,amass,resi,natmax,natom,
     .         ndat,qerror,prtlev)
          if (qerror) stop '*I/O error*'

          if (qmass.and.ndat.lt.4) then
              write(6,'(/2A)') progrer,
     .      ' Masses have not been all found, as they should.'
              qmass=.false.
          else if (mass.eq.'PDB'.and..not.qpdb) then
              write(6,'(2A)') progrwn,
     .      ' PDB keyword instead of COOR, for masses ? '
              mass='COOR' 
          endif
          goto 35
      endif

      if (prtlev.gt.1)
     .write(6,'(/2A)') program,' Coordinate file in PDB format.'

      i=0
      noff=0
 25   continue
      read(uninp,'(A)',end=30) lign80
      if (lign80(1:5).ne.'ATOM '.and.lign80(1:6).ne.'HETATM') goto 25
      i=i+1
      if (i.gt.natmax) then
          write(6,'(/2A,I6,A)') progrer,' More than ',natmax,
     .        ' atoms, the maximum allowed. Sorry.'
          stop '*Too large system*'
      endif

c     Coordinates are read according to the "pdb" format.
c     Masses can be in the B-factors column.

      ii=3*i-2
      read(lign80,'(6X,I5,1X,A4,1X,A4,A1,I4,4X,F8.3,F8.3,
     .     F8.3,1X,F5.2,F6.2,6X,A4)')
     .     at,atname,resname(i),ssu(i),resi(i),corlin(ii),corlin(ii+1),
     .     corlin(ii+2),tot,amass(ii),ssuch

      if (ssu(i).ne.ssusel.and.ssusel.ne.'*') then
          noff=noff+1
          goto 25
      endif

      if (amass(ii).le.0.and.qmass) then
          write(6,'(2A,I6,A,F8.2)') progrwn,' Atom ',i,
     .  ' Mass read in B-factor column= ',amass(ii)
          qmass=.false.
      endif
      amass(ii+1)=amass(ii)
      amass(ii+2)=amass(ii)
      goto 25
 30   continue
      natom=i
 
      write(6,'(/A,I6,A)') program,natom,' atoms picked in pdb file.'
      if (noff.gt.0)
     .write(6,'(A,I6,A)') program,noff,
     .  ' skipped, because they belong to other subunit(s).'
 35   continue

c     Tests:

      if (natom.lt.3) then
          write(6,'(/2A)') progrer,
     .  ' Not enough atoms found in file. Nothing done.'
          stop '*Too small system*'
      endif

      if (qmass) then
          write(6,'(/2A)') program,
     .        ' Masses were read in coordinate file.'

          write(6,'(/2A)') program,' Mass statistics: '
          call vecstat(amass,3*natom,rmin,rmax,rave,rdev)
          write(6,'(4(A,F12.6))') 
     .  ' <m>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

          if (rmin.le.0.d0) then
              write(6,'(/2A)') progrer,
     .      ' Negative or null masses found !'
              qmass=.false.
          endif
      endif

      if (.not.qmass) then
         do i=1,3*natom
            amass(i)=1.d0
         enddo
         write(6,'(2A)') program,'   All masses set to unity.'
         if (mass.eq.'COOR'.or.mass.eq.'PDB') then
             write(6,'(/2A)') progrwn,
     .     ' If the matrix was mass-weighted,'//
     .     ' wrong results are expected.'
         endif
      endif

      if (sstr.eq.'SUBU'.and..not.qpdb) then
          write(6,'(/2A)') progrer,
     .  ' SUBUnit information can be found in PDB files only.'
          stop '*Not consistent commands*'
      endif

c     ---------------------------------------------------------
c     Creation du fichier diagrtb_work.xyzm
c     Ce fichier contient:
c          1- les coordonnes (x1, y1, z1, x2, y2, z2, ...)
c          2- les masses associees
c     ---------------------------------------------------------

      if (prtlev.gt.1) 
     .write(6,'(2A)') program,' Coordinate file is rewritten.'
      write(unrtb) natom
      write(unrtb) (corlin(ii),ii=1,3*natom)
      write(unrtb) (amass(ii),ii=1,3*natom)
      close(unrtb)
 
c     --------------------------------------------------------
c     Substructuring information is picked in coordinate file:
c     --------------------------------------------------------

      nb=0
      do i=1,natom
         numbloc(i)=0
      enddo

      if (prtlev.gt.0.or.(sstr.eq.'DOMA'.or.sstr.eq.'SECO')) 
     .write(6,'(/2A)') program,' Substructuring:'
      if (prtlev.gt.1) then
      if (sstr.eq.'DOMA') then
      write(6,'(2A)') program,
     .  ' DOMAIN keyword expected in coordinate file.' 
      else if (sstr.eq.'SECO') then
      write(6,'(2A)') program,
     .  ' HELIX or SHEET keyword expected in coordinate file.' 
      endif
      endif

      if (sstr.eq.'SECO'.or.sstr.eq.'DOMA') then
          rewind(uninp)
 40       continue
          read(uninp,'(A)',end=45) lign80
          if ((sstr.eq.'DOMA'.and.lign80(1:6).ne.'DOMAIN').or.
     .        (sstr.eq.'SECO'.and.lign80(1:5).ne.'HELIX'.and.
     .         lign80(1:5).ne.'SHEET')) goto 40

          nb=nb+1

          if (lign80(1:5).eq.'HELIX'.or.sstr.eq.'DOMA') 
     .        read(lign80,'(15X,A4,A1,I5,2X,A4,A1,I5)',err=70) 
     .        aa1,ssu1,num1,aa2,ssu2,num2
          if (lign80(1:5).eq.'SHEET') 
     .        read(lign80,'(17X,A4,A1,I4,2X,A4,A1,I4)',err=70) 
     .        aa1,ssu1,num1,aa2,ssu2,num2

          write(6,'(6A,I5,3A,I5)') program,' ',lign80(1:6),' from ',
     .          aa1,ssu1,num1,' to ',aa2,ssu2,num2

          if (ssu1.ne.ssu2.or.num2.lt.num1) then
              write(6,'(/2A)') progrer,' Not coherent.'
              goto 70
          endif

          natbloc=0
          do i=1,natom
             if ((ssu(i).eq.'*'.or.ssu1.eq.ssu(i)).and.
     .          resi(i).ge.num1.and.resi(i).le.num2) then
                if (numbloc(i).gt.0) then
                   write(6,'(2A,I6,A)') progrwn,' Atom ',i,
     .           ' is already assigned block ',numbloc(i)
                endif
                natbloc=natbloc+1
                numbloc(i)=nb
             endif
          enddo

          if (natbloc.eq.0) then
              write(6,'(2A)') progrwn,
     .      ' Substructuring information ignored.'
              nb=nb-1
          endif
          goto 40

c         Split the loops into blocks: 

 45       continue
          if (prtlev.gt.0.or.nb.le.0) then
          if (nb.le.0) then
          write(6,'(/2A)') progrwn,
     .  ' Substructuring information not found. NONE assumed.'
          else
          if (sstr.eq.'SECO') then
          write(6,'(A,I6,A)') program,nb,
     .        ' secondary structure elements, each in its own block.'
          else
          write(6,'(A,I6,A)') program,nb,
     .        ' domains, each in its own block.'
          endif
          write(6,'(A,I6,A)') program,nbb,
     .  ' residue(s) per block, otherwise.'
          endif
          endif
       
          do i=1,natom

c            Fin d'une serie de blocs:

             if (i.ne.natom.and.numbloc(i+1).gt.0) then
               if (numbloc(i).eq.0) then
                   if (natbloc.eq.0) then
                      write(6,'(2A)') progrwn,
     .              ' Substructuring information ignored.'
                      nb=nb-1
                   endif
               else
                   goto 50
               endif
             elseif (i.eq.1) then
               nresb=1
               natbloc=0
               nb=nb+1
             endif

c            Debut d'une serie de blocs:

             if (i.ne.natom.and.numbloc(i).gt.0.and.
     .           numbloc(i+1).eq.0) then
                 nresb=0
                 natbloc=0
                 nb=nb+1
                 goto 50
             endif

             if (i.gt.1.and.resi(i).ne.resi(i-1)) then
                 nresb=nresb+1
             endif

c           "Saut" de numerotation => Fin de bloc.

             qok=.true.
             if (i.gt.1) 
     .       qok=(resi(i).eq.resi(i-1).and.ssu(i).eq.ssu(i-1)).or.
     .           (resi(i)-resi(i-1).eq.1.and.ssu(i).eq.ssu(i-1))

             if (nresb.gt.nbb.or..not.qok) then
                 if (natbloc.eq.0) then
                     write(6,'(2A)') progrwn,
     .             ' Substructuring information ignored.'
                     nb=nb-1
                 endif
                 natbloc=1
                 nresb=1
                 nb=nb+1
             elseif (i.eq.natom) then
                 if (natbloc.eq.0) then
                     write(6,'(2A)') progrwn,
     .             ' Substructuring information ignored.'
                     nb=nb-1
                 endif
             else
                 natbloc=natbloc+1
             endif

             numbloc(i)=nb

 50          continue
          enddo

          goto 100
 70       continue
          write(6,'(/2A/A)') progrer,
     .  ' While reading substructuring information in :',lign80
          stop '*Wrong data*'

c     Subunits are supposed to be written one after the other:

      elseif (sstr.eq.'SUBU') then
          nssu=0
          natbloc=0

          do i=1,natom

c           "Saut-arriere" de numerotation => Fin de bloc.

             qok=.true.
             if (i.gt.1) 
     .       qok=resi(i).gt.resi(i-1).and.ssu(i).eq.ssu(i-1)

             if (.not.qok) then
                nssu=nssu+1
                if (nb.gt.0) then
                    if (natbloc.eq.0) then
                        write(6,'(2A)') progrwn,
     .                ' Substructuring information ignored.'
                        nb=nb-1
                    endif
                endif
                nb=nb+1
                natbloc=0
             endif
             natbloc=natbloc+1
             numbloc(i)=nb
          enddo

          if (prtlev.gt.0)
     .    write(6,'(A,I6,A)') program,nssu,' subunits.'
 
c     sstr='NONE'
c     ----------
c     nbb residue(s) per block:
c    (Residues are supposed to be written one after the other)

      else
          if (prtlev.gt.0)
     .    write(6,'(A,I6,A)') program,nbb,' residue(s) per block.'
          nb=1
          nres=1
          nresb=1
          natbloc=0

          do i=1,natom
             if (i.gt.1.and.resi(i).ne.resi(i-1)) then
                 nres=nres+1
                 nresb=nresb+1
             endif

c           "Saut" de numerotation => Fin de bloc.

             qok=.true.
             if (i.gt.1) 
     .       qok=(resi(i).eq.resi(i-1).and.ssu(i).eq.ssu(i-1)).or.
     .           (resi(i)-resi(i-1).eq.1.and.ssu(i).eq.ssu(i-1))

             if (nresb.gt.nbb.or..not.qok) then
                 if (natbloc.eq.0) then
                     write(6,'(2A)') progrwn,
     .             ' Substructuring information ignored.'
                     nb=nb-1
                 endif
                 natbloc=1
                 nresb=1
                 nb=nb+1
             else
                 natbloc=natbloc+1
             endif
             numbloc(i)=nb
          enddo
          if (natbloc.eq.0) then
              write(6,'(2A)') progrwn,
     .      ' Substructuring information ignored.'
              nb=nb-1
          endif

          if (prtlev.gt.0)
     .    write(6,'(A,I6,A)') program,nres,' residues.'
      endif

 100  continue
      close(uninp)

c     ---------------------
c     Length of each block:    
c     ---------------------
      k=1
      natbloc=1

      do i=2,natom
         if (i.eq.natom.or.numbloc(i).ne.numbloc(i-1)) then
             if (i.eq.natom) natbloc=natbloc+1

             if (natbloc.lt.3) then
                 write(6,'(/A,I2,A,I4,A)') progrwn,
     .           natbloc,' atoms in block ',k,
     .         ' i.e., less than what is required for a rigid body.'

                 if (i.eq.natom) then
                     write(6,'(2A,I6,2A,I6)') progrwn,
     .             ' Last atom in this block is the ',i,
     .              'th, in residue ',ssu(i),resi(i)
                     write(6,'(2A)') progrwn,
     .             ' It is merged with last block.'

                     k=k-1
                     TABB(k)=TABB(k)+3*natbloc
                     if (k.eq.1.or.TABB(k).lt.natbmin) natbmin=TABB(k)
                     if (k.eq.1.or.TABB(k).gt.natbmax) natbmax=TABB(k) 

                 elseif (ssu(i).ne.ssu(i-1)) then
                     write(6,'(2A,I6,2A,I6)') progrwn,
     .             ' Last atom in this block is the ',i-1,
     .              'th, in residue ',ssu(i-1),resi(i-1)
                     write(6,'(2A)') progrwn,
     .             ' It is merged with the previous one.'

                     k=k-1
                     TABB(k)=TABB(k)+3*natbloc

                     if (k.eq.1.or.TABB(k).lt.natbmin) natbmin=TABB(k)
                     if (k.eq.1.or.TABB(k).gt.natbmax) natbmax=TABB(k) 

                     if (prtlev.gt.1) 
     .               write(6,'(A,I6,A,I6)') program,TABB(k)/3,
     .             ' atoms in block ',k 

                     k=k+1
                     natbloc=1
                 else
                     write(6,'(2A,I6,2A,I6)') progrwn,
     .             ' Last atom in this block is the ',i-1,
     .              'th, in residue ',ssu(i-1),resi(i-1)
                     write(6,'(2A)') progrwn,
     .             ' It will be merged with next block.'
                     natbloc=natbloc+1
                 endif
             else
                 TABB(k)=3*natbloc

                 if (k.eq.1.or.TABB(k).lt.natbmin) natbmin=TABB(k)
                 if (k.eq.1.or.TABB(k).gt.natbmax) natbmax=TABB(k) 

                 if (prtlev.gt.1) then 
                 write(6,'(A,I6,A,I6)') program,natbloc,
     .         ' atoms in block ',k 
                 write(6,'(A,I6)') ' Block first atom: ',i-natbloc
                 endif

                 if (i.ne.natom) then
                     k=k+1
                     natbloc=1
                 endif
             endif
         else
             natbloc=natbloc+1
         endif
      enddo
      nb=k
      natbmax=natbmax/3
      natbmin=natbmin/3
 
      write(6,'(/A,I6,A)') program,nb,' blocks.'

      if (prtlev.gt.0) then
      write(6,'(/2A,I6,A)') program,' At most, ',
     .   natbmax,' atoms in each of them.'
      write(6,'(2A,I6,A)') program,' At least,',
     .   natbmin,' atoms in each of them.'
      endif

      if (prtlev.gt.1) 
     .write(6,'(/2A)') program,' Normal end of Blocpdb.'

      RETURN
      END

c----------------------------------------------------------------------

      SUBROUTINE DIAGSTD(ndim,nvecout,eige,evord,amat,ev,evsort,work,
     .           prtlev)

c     ------------------------------------------------
c     Diagonalisation d'une matrice reelle symetrique.
c     ------------------------------------------------
c     Format de la matrice lue: i, j, element-ij-non-nul.
c     CERFBIN -Matrice binaire   (creuse)  : diagrtb_work.sdijb
c     Valeurs et vecteurs propres -> diagrtb_work.eigenfacs
c    (format CERFACS aussi) 
c
c     Routine de diagonalisation: TQLI (EISPACK).
c    -simple, et dans le domaine public...
c
c     --------------------------------
c     Ndim: maximum size of the matrix
c     --------------------------------
 
      implicit none

      logical qcrois, qexist, qinterr, qverbos
      integer evord(*), i, ii, j, jj, k, natom, nbig, nbelow, ndim, 
     .        nord, nredond, ntrace, nvec, nvecout, nunit, nzero, 
     .        prtlev, rdunit, unmess, unmodes
      double precision amat(ndim,*), bigzero, eigvalmin, 
     .       ev(*), evsort(*), matrd, trace, work(*)
      character cformat*20, cstatus*20, eige*4, matrice*20, nomfich*40, 
     .       program*9, progrer*12, progrwn*12
 
      program=' Diagstd>'
      progrer='%Diagstd-Er>'
      progrwn='%Diagstd-Wn>'
 
      if (prtlev.gt.1)
     .write(6,'(/2A)') program,' Entering in.'

c     Sortie standard: 
      unmess=6
      nvec=ndim
 
c     Detection de la matrice d'entree:
c     --------------------------------
      qverbos=.false.
      if (prtlev.gt.1) qverbos=.true.
      
      nunit=10
      rdunit=nunit
      nunit=nunit+1
      cformat='UNFORMATTED'
      cstatus='OLD'
      nomfich='diagrtb_work.sdijb'

      call openam(nomfich,cformat,cstatus,
     .     rdunit,qverbos,qinterr,qexist)

      if (qexist) then
          matrice='CERFBIN'
      else
          matrice='NONE'
          write(unmess,'(/2A)') progrer,
     .' Projected matrix not found. Expected name: diagrtb_work.sdijb '
          stop '*Projected matrix lost*'
      endif
 
      if (prtlev.gt.1)
     .write(unmess,'(3A)') program,
     .    ' Projected matrix to be read from file: ',nomfich
 
c     Lecture matrice d'entree (CERFACS, CERFBIN):
c     --------------------------------------------
 
      if (prtlev.gt.1)
     .write(unmess,'(2A)') program,
     .     ' CERFACS matrix format.'
 
c     1) Ordre de la matrice, nombre de lignes:
c     -----------------------------------------
      k=0
      nord=0
  90  continue
      if (matrice.eq.'CERFACS') then
      read(rdunit,*,end=100) i,j
      else
      read(rdunit,end=100) i,j
      endif
      k=k+1
      if (i.le.0.or.j.le.0) then
          write(unmess,'(/2A,I9,2(A,I6))')
     .    progrer,' in ligne: ',k,' I= ',i,' J= ',j
          stop '*Wrong matrix*'
      endif
      if (i.gt.nord) nord=i
      if (j.gt.nord) nord=j
      goto 90
 100  continue
 
      if (prtlev.gt.0)
     .write(unmess,'(2A,I9)')
     .     program,' Projected matrix order =',nord
      if (prtlev.gt.1)
     .write(unmess,'(2A,I9)')
     .     program,' Nb of non-zero elements:',k
 
      if (nord.gt.ndim) then
          write(unmess,'(/2A)')
     .    progrer,' Matrix can not be read.'
          if (nord.gt.ndim) write(unmess,'(2(A,I9))')
     .   ' Nord=  ',nord,' > Ndim=  ',ndim
          stop '*Wrong matrix*'
      endif
 
c     2) Lecture de la matrice:
c     -------------------------
      rewind(rdunit)
 
      nredond=0
      ntrace=0
      trace=0.d0
      nbig=0
      do i=1,nord
        do j=1,nord
         amat(i,j)=0.d0
        enddo
      enddo

      do jj=1,k
         if (matrice.eq.'CERFACS') then
         read(rdunit,*,err=95) i,j,matrd
         else
         read(rdunit,err=95) i,j,matrd
         endif
 
         if (dabs(matrd).gt.0.d0) then
             amat(i,j)=matrd
             amat(j,i)=matrd
             if (i.eq.j) then 
                trace=trace+matrd
                ntrace=ntrace+1
             endif
             if (matrd.gt.1E+10) then
                 nbig=nbig+1
                 if (nbig.lt.10) then
                     write(unmess,'(2A,2I12,A,G12.3)') 
     .               progrwn,' Element: ',i,j,' = ',matrd
                 else 
                     if (nbig.eq.10) write(unmess,*) '...'
                 endif
             endif
         else
             nredond=nredond+1
         endif
      enddo
      goto 105
  95  continue
      write(unmess,'(/2A,I6)')
     .     progrer,' while reading ligne ',k
      write(unmess,'(2I6,F16.8)') ' i, j, matrd= ',i,j,matrd
      stop '*Wrong matrix*'
 105  continue
 
      if (nredond.gt.0)
     .write(unmess,'(2A,I9)') progrwn,
     .    ' Nb of matrix elements found twice: ',nredond
      if (nbig.gt.0)
     .write(unmess,'(2A,I9)') progrwn,
     .    ' Nb of elements    > 1E+10 :',nbig
      if (prtlev.gt.1)
     .write(6,'(2A,F12.4)') program,' Projected matrix trace = ',trace
      if (nord-ntrace.ne.0)
     .write(unmess,'(2A,I11)') progrwn,
     .    ' Nb on zero elements on the trace:',nord-ntrace
 
c     Diagonalisation:
c     ----------------
      nomfich='diagrtb_work.eigenfacs'
      cformat='FORMATTED'
      cstatus='ove'
      unmodes=nunit
      nunit=nunit+1
      call openam(nomfich,cformat,cstatus,
     .     unmodes,qverbos,qinterr,qexist)
 
      if (prtlev.gt.1)
     .write(unmess,'(2A)') program,' Diagonalization.'
 
      if (nvec.gt.nord) nvec=nord
      write(unmess,'(/A,I7,A)') program,
     .      nvec,' eigenvectors are computed. '
      if (nvecout.gt.nvec) nvecout=nvec

      if (prtlev.gt.1)
     .write(unmess,'(A,I7,A)') program,
     .    nvecout,' of them to be saved. '
 
c     Initialisations:
      do i=1,ndim
         ev(i)=0.d0
      enddo
 
c     Eigenvalues/Matrix Diagonalization

c     The following routines (from the original EISPACK library) 
c     perform a diagonalization of a real symmetric matrix based 
c     on the QL algorithm. 

      CALL TRED2(amat,nord,ndim,ev,work)
      CALL TQLI(ev,work,nord,ndim,amat)
 
      trace=0.d0
      do i=1,nvec
         trace=trace+ev(i)
      enddo
      if (prtlev.gt.0)
     .write(unmess,'(/2A,F12.4)') program,
     .     ' Sum of eigenvalues =',trace
 
c     Compter les modes a Valeur propre "nulle":
c     a) Valeur propre la plus proche de la nullite.
c     b) Valeur propres "proches" de cette valeur (a un facteur 100 pres).
 
      do i=1,nvec
         if (i.eq.1.or.dabs(ev(i)).lt.eigvalmin) 
     .   eigvalmin=dabs(ev(i))
      enddo 
      if (prtlev.gt.1)
     .write(6,'(2A,F13.6)') program,
     .' Best zero-eigenvalue found  :',eigvalmin
 
      nbelow=0
      nzero=0
      do i=1,nvec
         if (dabs(ev(i)).lt.100*eigvalmin) then
             if (dabs(ev(i)).gt.bigzero) bigzero=dabs(ev(i))
             nzero=nzero+1
         endif
         if (ev(i).lt.-eigvalmin) then
             nbelow=nbelow+1
         endif
      enddo 
 
      if (prtlev.gt.0.or.nzero.ne.6) 
     .write(6,'(/A,I7,A,F13.7)') program,nzero,
     . ' zero-eigenvalues, that is, below or equal to:',bigzero
      if (nbelow.gt.0) 
     .write(6,'(/A,I7,A,F10.4)') program,nbelow,
     . ' negative eigenvalues, that is, below: ',-eigvalmin

      if (nzero.ne.6.or.nbelow.gt.0) then
          if (nbelow.gt.0.or.nzero.lt.6) then
              write(6,'(/2A)') progrwn,
     .      ' The studied structure is not a minimum-energy one ?'
              write(6,'(A)') 
     .'(Otherwise: coordinates, masses and matrix are not consistent ?)'
          else
              write(6,'(2A)') progrwn,
     .     ' Six expected. Parts of the structure interact too little ?'
          endif
      endif

c     Par ordre croissant ou decroissant:
      qcrois=.true.
      if (eige.eq.'HIGH') qcrois=.false.

      call trier(ev,nvec,ndim,evsort,evord,qcrois)

      write(unmess,'(/2A/(5F15.7))') program,
     .    ' Selected eigenvalues: ',(ev(evord(i)),i=1,nvecout)

      if (prtlev.gt.0)
     .WRITE(unmess,'(/2A/(5F15.7))') program,
     .    ' Frequencies (cm-1, '//
     .     'if the input matrix is a hessian in CHARMM units):',
     .    (sqrt(dabs(ev(evord(i))))*108.591365,i=1,nvecout)
 
c     Ecriture des modes normaux au format 'CERFACS':
c     -----------------------------------------------
      do j=1,nvecout
         i=evord(j)
         write(unmodes,'(A,I5,7X,A,1PG12.4)') 
     .       ' VECTOR',j,'VALUE',ev(i)
         write(unmodes,'(1X,35(1H-))') 
         write(unmodes,'(3(1PG12.4))') 
     .        (amat(k,i),k=1,nord)
      enddo
      close(unmodes)
 
      if (prtlev.gt.1)
     .write(unmess,'(/2A)')
     .      program,' Normal end.'
 
      return
      end
c----------------------------------------------------------------------
      SUBROUTINE TRED2(A,N,NP,D,E)

c     Reduce the matrix to tridiagonal form.

      integer i, j, k, l, n, np
      double precision A(NP,NP), D(NP), E(NP), f, g, h, hh, 
     .  scale

      IF(N.GT.1)THEN
        DO 18 I=N,2,-1
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.
      E(1)=0.
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.
            A(J,I)=0.
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END
c----------------------------------------------------------------------
      SUBROUTINE TQLI(D,E,N,NP,Z)

c     Finds the eigenvalues and eigenvectors of a tridiagonal matrix:

      integer i, iter, j, k, l, m, n, np
      double precision b, c, D(NP), dd, E(NP), f, g, p, r, s, 
     .       Z(NP,NP)

      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END
c----------------------------------------------------------------------
      SUBROUTINE PREPMAT ( fmtx, typf, natblocs, natom, nb, qtrace,
     &                     indexi, indexj, tabb, a, hw, prtlev)
 
c     Scalar Arguments
 
      character        fmtx*(*), typf*(*)
      integer          natblocs, natom, nb 
 
c     Array Arguments
 
      logical          qtrace(3*natom)
      integer          indexi(9*natblocs*natblocs),
     &                 indexj(9*natblocs*natblocs), 
     &                 prtlev, tabb(nb)
      double precision a(3*natblocs,3*natom), hw(9*natblocs*natblocs)
 
c     Purpose:
c     =======
 
c     --------------------------------------------
c     Matrix preparation: it is split into blocks.
c     --------------------------------------------
c
c     Arguments
c     =========
c
c     natblocs  : maximum number of atoms found in a block.
c     natom     : number of atoms in the system.
c     nb        : number of blocks the system is split into.
c     qtrace    : flag to remember that a given diagonal element is known.
c     indexi    : i-index of non-zero elements in block hw.
c     indexj    : j-index of non-zero elements in block hw.
c     tabb      : block lengths.
c     a         : band-matrix corresponding to a given block.
c     hw        : block being filled.
c 
C-----------------------------------------------------------------------
c
c     Local variables
 
      logical          qinit
      integer          bloc, i, ii, iii, imax, j, jbloc, k, kk, 
     &                 nblocs, nbnnul, nempty, ni1, ni2, nj1, nj2, 
     &                 nlign, nn, nread, ntotn1, tbloc, tt
      double precision toto, trace
      character        program*9, progrer*12, progrwn*12
c
C-----------------------------------------------------------------------
      program=' Prepmat>'
      progrwn='%Prepmat-Wn>'
      progrer='%Prepmat-Er>'
 
      if (prtlev.gt.1)
     .write(6,'(/2A)') program,' Entering in.'

      nblocs=nb
 
c     Fichiers d'entree 
c     -----------------
      if (typf.eq.'FREE') then
      open(unit=50,file=fmtx,status='old',form='formatted')
      else
      open(unit=50,file=fmtx,status='old',form='unformatted')
      endif
      
c     Fichiers de sortie
c     ------------------
      open(unit=51,file='diagrtb_work.matblocs',status='unknown',
     .     form='unformatted')

c    --------------------------------------------------------------
c    Creation du fichier diagrtb_work.matblocs => matrice ecrite par bloc
c    Ce fichier contient
c       1- la longueur de chaque bloc
c       2- les elements de la matrice correspondant a chaque bloc 
c    --------------------------------------------------------------

      if (prtlev.gt.1)
     .write(6,'(2A)') program,' Rewriting of the matrix begins.'

      write(51) natom,NB,(TABB(k),k=1,NB)
      nlign=1
      nempty=0
      
c     ----------------------------      
c     Initialisation compteur
c     ----------------------------

      NI1=1
      NTOTN1=0
      nn=1
      NI2=TABB(nn)
      qinit=.true.
      TBLOC=0
 
      imax=-1
      trace=0.d0
      do i=1,natom
         qtrace(i)=.false.
      enddo

      nread=1
 50   continue

c     NI1,NI2 bornes du residue vertical
c     ----------------------------------

c        Matrice-bande A(i,j)=0.d0:
c        -------------------------
	 if (qinit) then
         do i=1,TABB(nn)
            kk=NTOTN1+1
            do j=kk,3*natom
               A(i,j)=0.d0
            enddo
         enddo
         endif

         if (typf.eq.'FREE') then
         READ(50,*,end=99,err=900) i,j,toto
         else
         READ(50,end=99,err=900) i,j,toto
         endif
         nread=nread+1

         if (i.gt.3*natom.or.j.gt.3*natom) then
             write(6,'(/2A,I6,A,I6,A/A,I6,A)') progrer,
     .     ' Matrix element ',i,' j= ',j,
     .     ' found.',' More than ',natom,' atoms in matrix file !'
             stop '*Wrong matrix or wrong coordinate file*'
         endif

         if (i.eq.j) then
         if (.not.qtrace(i)) then
             trace=trace+toto
             qtrace(i)=.true. 
         endif
         endif

         if (i.gt.imax) imax=i
         if (j.gt.imax) imax=j

         IF (i.le.NI2) then
            iii=i-NTOTN1
            A(iii,j)=toto
	    qinit=.false.
         endif
 
c        Decoupage de la matrice-bande lue,
c        et sortie des blocs, a raison d'un bloc par ligne:
 
         if (i.gt.NI2) then
c           On a depasse NI2, on remonte d'un cran:
            backspace 50
	    NJ1=NI1
	    NJ2=NI2

c    NJ1,NJ2 borne du residue en horizontal
c    --------------------------------------
c    au depart NJ1 et NJ2 = NI1 et NI2
c    TBLOC pour TABB 
c    a chaque fois il faut decaler de un pour avoir
c    la longueur du bon bloc

	 tt=0
	 TBLOC=TBLOC+1
	 BLOC=TBLOC
 
c    Boucle sur les blocs pour une bande NI1-NI2
c    NI1 et NI2 constant dans la boucle
c    -------------------------------------------
         do JBLOC=1,NB
	    nbnnul=0
            ii=0
	    tt=tt+1
c
c    Si premier bloc, seulement une partie donc on
c    cree la partie inferieure
c    ici NJ1 et NJ2 sont egaux a NI1 et NI2
c    ---------------------------------------------

         if (tt.eq.1) then
               do i=1,TABB(nn)
                  do j=kk,NI2   
                     A(j-NTOTN1,i+NTOTN1)=A(i,j)
                  enddo
               enddo
               do i=1,TABB(nn)
                  do j=NJ1,NJ2
                     if (A(i,j).ne.0.d0) then
                        ii=ii+1
                        HW(ii)=0.d0
                        nbnnul=ii
                        HW(ii)=A(i,j)
                        indexi(ii)=i+NTOTN1
                        indexj(ii)=j
                     endif
                  enddo
              enddo
         if (nbnnul.gt.9*natblocs*natblocs) then
             write(6,'(/2A/2(A,I12))') progrer,
     .     ' Too many matrix elements in a single block: ',
     .     ' Nbnnul= ',nbnnul,' Max= ',9*natblocs*natblocs
             stop '*Too large block*'
         endif
         if (nbnnul.le.0) nempty=nempty+1
         nlign=nlign+1
         write(51) nbnnul,(indexi(ii),indexj(ii),HW(ii),ii=1,nbnnul) 
         endif
c
c    Pour les autres blocs de la bande
c    NJ1 et NJ2 different de NI1 et NI2
c    on decale au bloc suivant
c    bloc suivant=> bloc+1 sachant que la valeur de depart
c    pour bloc est tbloc
c    ----------------------------------
         if (tt.gt.1) then
	       BLOC=BLOC+1
	       NJ1=NJ2+1
	       NJ2=NJ2+TABB(BLOC)
               do i=1,TABB(nn)
                  do j=NJ1,NJ2
                     if (A(i,j).ne.0.d0) then
                        ii=ii+1
                        HW(ii)=0.d0
                        nbnnul=ii
                        HW(ii)=A(i,j)
                        indexi(ii)=i+NTOTN1
                        indexj(ii)=j
                     endif
                  enddo
               enddo
         if (nbnnul.gt.9*natblocs*natblocs) then
             write(6,'(/2A/2(A,I12))') progrer,
     .     ' Too many matrix elements in a single block: ',
     .     ' Nbnnul= ',nbnnul,' Max= ',9*natblocs*natblocs
             stop '*Too large block*'
         endif
         if (nbnnul.le.0) nempty=nempty+1
         nlign=nlign+1
          write(51) nbnnul,(indexi(ii),indexj(ii),HW(ii),ii=1,nbnnul)
	 endif
        enddo
 
c       fin boucle sur les blocs => initialisation de A
c       necessaire donc qinit=vrai
c       on aura un bloc en moins => nb=nb-1
c       changement des bornes verticales
c       ------------------------------------------------------
	qinit=.true.
        NTOTN1=NTOTN1+TABB(nn)
	NB=NB-1
	nn=nn+1
	NI1=NI2+1
	NI2=NI2+TABB(nn)
       endif

      goto 50
 99   continue

      if (prtlev.gt.0) 
     .write(6,'(/A,I13,A)') program,nread,' matrix lines read.'

      if (prtlev.gt.0)
     .write(6,'(/2A,I10)') program,' Matrix order    = ',imax
      write(6,'(/2A,F15.4)') program,' Matrix trace    = ',trace

      if (imax/3.ne.natom) then
          write(6,'(/A,I6,A/A,I6,A)') progrer,imax,
     .  ' coordinates found in this matrix, ',
     .  ' instead of ',natom,'*3, as expected. Not the right one ?'
          stop '*Wrong matrix or wrong coordinate file*'
      endif

      if (prtlev.gt.1)
     .write(6,'(2A,2I7,F15.4)') program,' Last element read: ',i,j,toto

            iii=i-NTOTN1
            A(iii,j)=toto
            NJ1=NI1
            NJ2=NI2
	    ii=0
               do i=1,TABB(nn)
                  do j=kk,NJ2
                     A(j-NTOTN1,i+NTOTN1)=A(i,j)
                  enddo
               enddo
               do i=1,TABB(nn)
                  do j=NJ1,NJ2
                     if (A(i,j).ne.0.d0) then
                        ii=ii+1
                        HW(ii)=0.d0
                        nbnnul=ii
                        HW(ii)=A(i,j)
                        indexi(ii)=i+NTOTN1
                        indexj(ii)=j
                     endif
                  enddo
               enddo
         if (nbnnul.gt.9*natblocs*natblocs) then
             write(6,'(/2A/2(A,I12))') progrer,
     .     ' Too many matrix elements in a single block: ',
     .     ' Nbnnul= ',nbnnul,' Max= ',9*natblocs*natblocs
             stop '*Too large block*'
         endif
         if (nbnnul.le.0) nempty=nempty+1

      nlign=nlign+1
      write(51) nbnnul,(indexi(ii),indexj(ii),HW(ii),ii=1,nbnnul)
 
      nb=nblocs 
      if (prtlev.gt.1) then
          write(6,'(A,I13,A)') program,nlign,' lines saved.'
          write(6,'(A,I13,A)') program,nempty,' empty lines.'
      endif

      if (nlign.ne.nb*(nb+1)/2+1) then
          write(6,'(/A,I6,A)') progrer,
     .    nb*(nb+1)/2+1,' lines are expected on output !'
          stop '*Unexpected big problem*'
      else if (prtlev.gt.1) then
          write(6,'(2A)') program, 
     .  ' Number of lines on output is as expected.'
      endif
 
      close(50)
      close(51)

      if (prtlev.gt.1) 
     .    write(6,'(2A)') program,' Normal end of Prepmat.'

      RETURN
 900  write(6,'(/2A,I9)') progrer,' While reading matrix ligne ',nread
      stop '*Wrong matrix*'
      END
c----------------------------------------------------------------------
      SUBROUTINE RTB ( natom, nb, nbrt, nddres, nddres2, indexi, indexj,
     &           n, amass, corlin, d, h, hd, hh, mat, rt, s, prtlev )
c
c     Scalar Arguments
      integer          natom, nb, nbrt, nddres, nddres2
c
c     Array Arguments
      integer          indexi(nddres2), indexj(nddres2), n(nb), prtlev
      double precision amass(3*natom), corlin(3*natom), 
     &                 d(nbrt,nddres), h(nddres2), hd(nddres,nbrt),
     &                 hh(nddres,nddres), mat(6*nb,6*nb),
     &                 rt(nbrt,nddres,nb), s(nbrt,nbrt) 
c
c     Purpose:
c     =======
c
c     ------------------
c     Matrix projection.
c     ------------------
c
c     Lit les fichiers produits par "prepmat", c'est-a-dire:
c     Les coordonnes et les masses ("diagrtb_work.xyzm")
c     La matrice re-ecrite bloc par bloc ("diagrtb_work.matblocs")
c
c     Ecrit la matrice projetee ("diagrtb_work.sdijb")
c    -format i,j, element-non-nul.
c
c     Decrit les rotation-translations de
c     chaque element de chaque bloc ("diagrtb_work.blocs")
c
c     A faire: taux de remplissage de la matrice projetee.
c
c     Arguments
c     =========
c
c     natom   : number of atoms in the system.
c     nb      : number of blocks the system is split into.
c     nbrt    : number of degrees of freedom of each block.
c     nddres  : size of largest block.
c     nddres2 : nddres**2
c     indexi  : i-index of non-zero elements in block hw.
c     indexj  : j-index of non-zero elements in block hw.
c     n       : block lengths.
c     amass   : masses.
c     corlin  : coordinates.
c     d       : tranlation and rotation vectors of each block.
c     h       : block being considered.
c     hd      : working matrix for projection of a given block.
c     hh      : working matrix for projection of a given block.
c     mat     : projected matrix.
c     rt      : tranlation and rotation vectors of all blocks.
c     s       : working matrix for projection of a given block.
c
C-----------------------------------------------------------------------
c
c     Local variables
c
      logical          qok
      integer          i, ibloc, ii, indxi, indxj, j, jbloc, k, kk, kr,
     &                 mmm, natcur, nbcur, nbnnul, nempty, ni, nj,
     &                 nlign, nnskip, nnuls, nskip, ntot, ntoti, ntotj
      double precision amasstot, trace, xg, yg, zg
      character        program*5, progrer*8, progrwn*8
C
C     ------------------------- CHANGEMENT D'ORIGINE DES COORDONNEES CARTESIENNES
 
      program=' RTB>'
      progrwn='%RTB-Wn>'
      progrer='%RTB-Er>'

      if (prtlev.gt.1)
     .write(6,'(/2A)') program,' Entering in.'

c     Fichier d'entree
c     ----------------

      open(unit=74,file='diagrtb_work.xyzm',status='old',
     .     form='unformatted')
      open(unit=51,file='diagrtb_work.matblocs',status='old',
     .     form='unformatted')
 
c     Fichier de sortie
c     -----------------

      open(unit=52,file='diagrtb_work.sdijb',status='unknown',
     .     form='unformatted')
      open(unit=35,file='diagrtb_work.blocs',status='unknown',
     .     form='unformatted')

c     Lectures + Tests:
c     -----------------

      natcur=natom
      read(74) natom

      if (prtlev.gt.1)
     .write(6,'(2A,I6)') program,
     .  ' Number of atoms found in temporary coordinate file:',natom
      if (natom.ne.natcur) then
          write(6,'(/A,I6,A)') progrer,natcur,' atoms... up to now.'
          stop '*Unexpected big problem*'
      endif
       
      read(74) (corlin(ii),ii=1,3*natom)
      read(74) (amass(ii),ii=1,3*natom)
 
      qok=.true.
       AMASSTOT = 0.D0 
       DO I=1,3*natom
	 AMASSTOT = AMASSTOT + AMASS(I)
         if (amass(i).lt.0) qok=.false.
       ENDDO      
      AMASSTOT = AMASSTOT / 3.D0 

      if (.not.qok) then
          write(6,'(2A)') progrer,
     .  ' Masses with negative values read in temporary file.'
          stop '*Unexpected big problem*'
      endif
      if (prtlev.gt.1)
     .write(6,'(2A,F16.4)') program,' Total mass = ',AMASSTOT
C
      natcur=natom
      nbcur=nb
      read(51) natom,NB,(N(k),k=1,NB)

      if (prtlev.gt.1)
     .write(6,'(2A,I6)') program,
     .  ' Number of atoms found in matrix:',natom
      if (natom.ne.natcur) then
          write(6,'(2A)') progrer,
     .  ' Pdb file and matrix are not consistent any more !'
          stop '*Unexpected big problem*'
      endif
 
      if (prtlev.gt.1)
     .write(6,'(2A,I6)') program,' Number of blocks =',nbcur
      if (nb.ne.nbcur) then
          write(6,'(/A,I6,A)') progrer,nbcur,' blocks... up to now !'
          stop '*Unexpected big problem*'
      endif
      write(35) natom,NB,(N(k),k=1,NB)

      natcur=0
      do i=1,nb
         natcur=natcur+n(i)
      enddo
      if (natcur.ne.3*natom) then
          write(6,'(2A,I6,A,I6)') progrer,
     .  ' Sum of block lengths: ',natcur,
     .  ' instead of ',3*natom
          stop '*Unexpected big problem*'
      endif

      nlign=1
      nempty=0
      nskip=0

C     <------------------------------------------ DEBUT DE BOUCLE SUR LES BLOCS
C
      if (prtlev.gt.1)
     .write(6,'(2A)') program,' Projection begins.'
      kk=0
      NTOT=0

      DO IBLOC=1,NB  
C
C     REMPLISSAGE DES BLOCS
C
      nbnnul=0

      read(51,end=895,err=900) 
     .     nbnnul,(indexi(ii),indexj(ii),H(ii),ii=1,nbnnul)
      nlign=nlign+1

      if (nbnnul.gt.nddres2) then
          write(6,'(A,I6,A,I6,A,I6)') progrer,
     .    nbnnul,' elements in bloc ',ibloc,
     .  ' Maximum allowed: ',nddres2
          stop '*Too large block*'
      endif
      if (nbnnul.le.0) then
          nempty=nempty+1
      else
      do i=1,nbnnul
         if (indexi(i).le.NTOT.or.indexj(i).le.NTOT) then
             write(6,'(A,5(A,I6))') progrer,
     .     ' Ntot=',ntot,' but i= ',indexi(i),' j= ',
     .       indexj(i),' for element ',i,' of bloc ',ibloc 
             stop '*Unexpected big problem*'
         endif
      enddo
      endif
     
c     <------------------------------------------On se place au bloc suivant 

      kk=kk+1
      do kr=1,NB-kk
      read(51,end=905,err=905) nnskip
      if (nnskip.le.0) nempty=nempty+1
      nskip=nskip+1
      enddo
 
c     <------------------------------------------Initialisation du bloc

      if (n(ibloc).gt.nddres.or.n(ibloc).le.0) then
          write(6,'(A,I6,A,I6,A,I6)') progrer,' N(ibloc)= ',n(ibloc),
     .  ' for bloc ',ibloc,
     .  ' Maximum allowed: ',nddres
          stop '*Too large block*'
      endif

      do i=1,N(IBLOC)
         do j=1,N(IBLOC)
         HH(i,j)=0.d0
         enddo
      enddo

c     <------------------------------------------Remplissage du bloc

      do ii=1,nbnnul
      i=indexi(ii)-NTOT
      j=indexj(ii)-NTOT
      HH(i,j)=H(ii)
      enddo

C 
C     CHANGEMENT D'ORIGINE DANS CHAQUE BLOC (CDG DU BLOC)
c     ---------------------------------------------------------

        AMASSTOT = 0.D0
        DO I=1,N(IBLOC)
         AMASSTOT = AMASSTOT + AMASS(NTOT+I)
        ENDDO      
	AMASSTOT=AMASSTOT/3.d0

C
        XG = 0.D0
        YG = 0.D0  
        ZG = 0.D0

        DO I=0,( (N(IBLOC)/3)-1 )
         XG = XG + AMASS(NTOT+3*I+1)*corlin(NTOT+3*I+1)
         YG = YG + AMASS(NTOT+3*I+2)*corlin(NTOT+3*I+2)
         ZG = ZG + AMASS(NTOT+3*I+3)*corlin(NTOT+3*I+3)
        ENDDO

       XG = XG / AMASSTOT
       YG = YG / AMASSTOT
       ZG = ZG / AMASSTOT
C
       DO I=0,( (N(IBLOC)/3)-1)      
        corlin(NTOT+3*I+1)= corlin(NTOT+3*I+1) - XG
        corlin(NTOT+3*I+2)= corlin(NTOT+3*I+2) - YG
        corlin(NTOT+3*I+3)= corlin(NTOT+3*I+3) - ZG
       ENDDO

C
C     PROJECTION DES TRANSLATIONS ET ROTATIONS --------------------------
C

      DO I=1,6
       DO J=1,N(IBLOC)
        D(I,J) = 0.D0 
       ENDDO
      ENDDO
C
      DO I=1, N(IBLOC)/3
       II=NTOT 
       D(1,1+3*(I-1))= DSQRT(AMASS(II+1+3*(I-1)))
       D(2,2+3*(I-1))= DSQRT(AMASS(II+2+3*(I-1)))
       D(3,3+3*(I-1))= DSQRT(AMASS(II+3+3*(I-1)))
       D(4,2+3*(I-1))=-DSQRT(AMASS(II+2+3*(I-1)))*corlin(II+3+3*(I-1))
       D(4,3+3*(I-1))= DSQRT(AMASS(II+3+3*(I-1)))*corlin(II+2+3*(I-1))
       D(5,1+3*(I-1))= DSQRT(AMASS(II+1+3*(I-1)))*corlin(II+3+3*(I-1))
       D(5,3+3*(I-1))=-DSQRT(AMASS(II+3+3*(I-1)))*corlin(II+1+3*(I-1))
       D(6,1+3*(I-1))=-DSQRT(AMASS(II+1+3*(I-1)))*corlin(II+2+3*(I-1))
       D(6,2+3*(I-1))= DSQRT(AMASS(II+2+3*(I-1)))*corlin(II+1+3*(I-1))
      ENDDO
C

      MMM=6

      CALL SCHMIDT(MMM,N(IBLOC),D ) 
C
      do i=1,6
         do j=1,N(IBLOC)
            RT(i,j,ibloc)=D(i,j)
	    write(35) i,j,ibloc,rt(i,j,ibloc)
         enddo
      enddo

      NTOT=NTOT+N(IBLOC)
      ENDDO
C   <-------------------------------------------FIN DE BOUCLE SUR LES BLOCS

      close(51)

      open(unit=51,file='diagrtb_work.matblocs',status='old',
     .     form='unformatted')

      read(51) natom,NB,(N(k),k=1,NB)

      ni=0
      nj=0
      indxi=0
      indxj=0
      NTOTI=0
      NTOTJ=0

      do IBLOC=1,NB
         do JBLOC=IBLOC,NB

         read(51,end=910,err=910) 
     .   nbnnul,(indexi(ii),indexj(ii),H(ii),ii=1,nbnnul)


         do i=1,N(IBLOC)
              do j=1,N(JBLOC)
               HH(i,j)=0.d0
              enddo
         enddo

         do ii=1,nbnnul
              i=indexi(ii)-NTOTI
              j=indexj(ii)-NTOTJ
              HH(i,j)=H(ii)
         enddo


         do j=1,N(IBLOC)
              do i=1,6
                HD(j,i)=0.d0
                 do k=1,N(JBLOC)
                  HD(j,i)=HD(j,i)+HH(j,k)*RT(i,k,jbloc)   
                 enddo
              enddo
         enddo

         do i=1,6
             do j=1,6
	       s(i,j)=0.d0
               do k=1,N(IBLOC)
                  s(i,j)=s(i,j)+RT(i,k,ibloc)*HD(k,j)
               enddo
                  ni=i+indxi
                  nj=j+indxj
                  mat(ni,nj)=0.d0
		  mat(ni,nj)=s(i,j)
             enddo
         enddo


      NTOTJ=NTOTJ+N(JBLOC)

      indxj=indxj+6

      ENDDO

      indxi=indxi+6
      indxj=indxi 
      NTOTI=NTOTI+N(IBLOC)
      NTOTJ=NTOTI

      ENDDO

c     Ecriture de la matrice projetee:

      if (prtlev.gt.1) 
     .write(6,'(2A)') program,' Projected matrix is being saved.'

      trace=0.d0
      nnuls=0
      do i=1,6*nb
	 do j=i,6*nb
         if (mat(i,j).ne.0.d0) then
            nnuls=nnuls+1
            write(52) i,j,mat(i,j)
            if (i.eq.j) trace=trace+mat(i,j)
         endif
	 enddo
      enddo

      close (35)
      close (51)
      close (52)
      close (74)

      if (prtlev.gt.0) 
     .write(6,'(2A,F12.4)') program,' Projected matrix trace = ',trace

      if (prtlev.gt.1) then
      write(6,'(A,I10,A)') program,nnuls,' non-zero elements.'
      write(6,'(2A)') program,' Normal end of RTB.'
      endif

      return
 895  continue
      write(6,'(2A,I3,A,I3)') progrer,
     .    ' End-of-file while reading I-bloc: ',ibloc,
     .    ' nbnnul= ',nbnnul
      goto 990
 900  continue
      write(6,'(2A,I3,A,I3)') progrer,
     .    ' Error while reading I-bloc: ',ibloc,
     .    ' nbnnul= ',nbnnul
      goto 990
 905  continue
      write(6,'(A,3(A,I3))') progrer,
     .    ' Error after reading I-bloc: ',ibloc, 
     .    ' NB= ',NB,' kk= ',kk 
      goto 990
 910  continue 
      write(6,'(2A,I3,A,I3)') progrer,
     .    ' Error while reading J-bloc: ',jbloc,
     .    ' for band: ',ibloc
 990  continue
      write(6,'(A,I6,A)') program,nlign,' lines read.'
      write(6,'(A,I6,A)') program,nempty,' empty lines found.'
      write(6,'(A,I6,A)') program,nskip,' lines skipped.'
      stop '*Matrix i/o error*'

      END

c -----------------------------------------------
      SUBROUTINE SCHMIDT(M,N,C)
c ----------------------------------------------- 
C
C     ORTHOGONALISATION PAR GRAM-SCHMIDT
C
      integer   I,J,K,M,N
      double precision aaa, anorm, C(6,*), REC(6,6)
C
      ANORM = 0.D0
      DO I=1,N
       ANORM = ANORM + C(1,I)*C(1,I)
      ENDDO
      ANORM = 1.D0 / ( DSQRT(ANORM) )
C
      DO I=1,N
       C(1,I) = ANORM * C(1,I)
      ENDDO
C 
      DO  I=2,M      
       DO J=1,I-1
        REC(J,I) = 0.D0
        DO K=1,N
         REC(J,I) = REC(J,I) + C(J,K)*C(I,K)
        ENDDO
       ENDDO
       DO K=1,N
        AAA = 0.D0
        DO J=1,I-1
         AAA = AAA + C(J,K)*REC(J,I)
        ENDDO
        C(I,K) = C(I,K) - AAA
       ENDDO
C
       ANORM = 0.D0
       DO K=1,N
        ANORM = ANORM + C(I,K)*C(I,K)
       ENDDO 
       ANORM = 1.D0 / ( DSQRT(ANORM) )
C
       DO K=1,N
        C(I,K) = ANORM*C(I,K)
       ENDDO
C
      ENDDO
      RETURN
      END  
c----------------------------------------------------------------------
      SUBROUTINE RTBTOMODES ( natom, nb, nbrt, nddres, nvec, n, numvec1,
     &                        covec, freq1, matvec1, norm, rt, prtlev,
     &                        feig )
c
c     Scalar Arguments
c
      integer          natom, nb, nbrt, nddres, nvec
c
c     Array Arguments
c
      integer          n(*), numvec1(*), prtlev
      double precision covec(3*natom,nvec), freq1(nvec),
     &                 matvec1(6*nb,nvec), norm(nvec),
     &                 rt(nbrt,nddres,nb)
      character        feig*(*)
c
c     Purpose:
c     =======
c
c     -----------------------------------------------------------
c     From eigenvectors in block-rotation-translation coordinates
c     to eigenvectors in cartesian coordinates.
c     -----------------------------------------------------------
c
c     Input : matrix file 'diagrtb_work.eigenfacs'
c     Output: matrix file feig
c
c     Arguments
c     =========
c
c     natom     : maximum number of atoms in the system.
c     nresmax   : maximum number of residues and of blocks.
c     ntab      : working array.
c     resi      : residue number each atom belongs to.
c     tabb      : block lengths.
c     amass     : coordinate masses.
c     corlin    : coordinates.
c     massat    : atomic masses.
c     xat       : atomic x-coordinate.
c     yat       : atomic y-coordinate.
c     zat       : atomic z-coordinate.
c     natblocs  : maximum number of atoms found in blocks.
c     natom     : number of atoms in the system.
c     nb        : number of blocks the system is split into.
c
C-----------------------------------------------------------------------
c
c     Local variables
c
      character*40     cformat, cstatus, nomeig1, nomeig2 
      logical          qinter, qexist, qverbos
      integer          i, ibloc, ii, ivec, j, k, natcur, nbcur,
     .                 nbread, nddl1,
     .                 ntot, nunit, uneig1, uneig2
      double precision dvec
      character        program*14, progrer*17, progrwn*17
c
C-----------------------------------------------------------------------
 
      program=' Rtb_to_modes>'
      progrer='%Rtb_to_modes-Er>'
      progrer='%Rtb_to_modes-Wn>'

      if (prtlev.gt.1)
     .write(6,'(/2A)') program,' Entering in.'

c     Fichier en lecture
c      -------------------

      nunit=10
      open (unit=35,file='diagrtb_work.blocs',status='old',
     .      form='unformatted')
 
      natcur=natom
      nbcur=nb
      read(35) natom,nb,(n(i),i=1,nb)

      if (prtlev.gt.1)
     .write(6,'(/2A,I6)') program,
     .  ' Number of atoms in temporary block-file = ',natom
      if (natom.ne.natcur) then
          write(6,'(/A,I6,A)') progrer,natcur,' atoms... up to now.'
          stop '*Unexpected big problem*'
      endif

      if (prtlev.gt.1)
     .write(6,'(2A,I6)') program,' Number of blocs = ',nb
      if (nb.ne.nbcur) then
          write(6,'(/A,I6,A)') progrer,nbcur,' blocks... up to now.'
          stop '*Unexpected big problem*'
      endif

      qverbos=.false.
      if (prtlev.gt.1) qverbos=.true.

c     Matrice-input:

      nomeig1="diagrtb_work.eigenfacs"
      uneig1=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"
      call openam(nomeig1,cformat,cstatus,uneig1,qverbos,
     .            qinter,qexist)
      if (qinter.or..not.qexist) stop
 
c     Matrice-output:

      nomeig2=feig
      uneig2=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomeig2,cformat,cstatus,uneig2,qverbos,
     .            qinter,qexist)
      if (qinter) stop

c     Lecture des fichiers:
c     --------------------

      call rdmodfacs(uneig1,6*nb,nvec,numvec1,freq1,
     .               matvec1,nddl1,prtlev)
 
      if (prtlev.gt.1)
     .write(6,'(/A,I5,A,I6,A)') program,
     .nvec,' vectors, with ',nddl1,' coordinates in vector file.'

      if (nddl1.ne.6*nb) then
          write(6,'(/2A,I6,A)') progrer,' Vectors of length: ',6*nb,
     .  ' were expected. Block and vector files are not consistent.'
          stop '*Unexpected big problem*'
      endif
c
      if (prtlev.gt.1) then
      do ivec=1,nvec
         norm(ivec)=0.d0
         do i=1,nddl1
         dvec=matvec1(i,ivec)**2
         norm(ivec)=norm(ivec)+dvec
         enddo
         norm(ivec)=dsqrt(norm(ivec))
      enddo
      write(6,'(/A)') 
     .' Norm of eigenvectors in projected coordinates (one expected):'
       write(6,'(5F8.5)') (norm(ivec),ivec=1,nvec)
      endif
c
      if (prtlev.gt.1)
     .write(6,'(/2A)') program,' RTB block-file is being read.'
c
      nbread=0
      do i=1,nbrt
       do j=1,nddres
         do k=1,nb
           rt(i,j,k)=0.d0
         enddo
       enddo
      enddo
c
   90 continue
      read(35,end=99) i,j,ibloc,rt(i,j,ibloc)
      nbread=nbread+1
      if (i.gt.nbrt.or.i.le.0) then
          write(6,'(A,I6,A,I1,A,I4)') progrer,i,
     .  ' rigid body ddl instead of ',nbrt,' for block ',ibloc
          stop '*Unexpected big problem*'
      endif
      if (j.gt.nddres.or.j.le.0) then
          write(6,'(A,I6,A,I6,A,I6)') progrer,j,
     .  ' ddl for block ',ibloc,'. Maximum allowed is: ',nddres
          stop '*Unexpected big problem*'
      endif
      if (ibloc.gt.nb.or.ibloc.le.0) then
          write(6,'(A,I4,A,I4)') progrer,ibloc,
     .  ' blocks at least. Maximum allowed is: ',nb
          stop '*Unexpected big problem*'
      endif
      goto 90
   99 continue
 
      if (prtlev.gt.1) 
     .write(6,'(A,I9,A)') program,nbread,' lines found in RTB file.'
 
      do ivec=1,nvec
	 do i=1,3*natom
	    covec(i,ivec)=0.d0
         enddo
      enddo
 
      do ivec=1,nvec
       ntot=0
       ii=0
       do ibloc=1,nb
	 do j=1,n(ibloc)
	   covec(j+ntot,ivec)=0.d0
	   do i=1,6
	   covec(j+ntot,ivec)=covec(j+ntot,ivec)+
     .     rt(i,j,ibloc)*matvec1(i+ii,ivec)
	   enddo
          enddo
	   ii=ii+6
           ntot=ntot+n(ibloc)
       enddo
       enddo
 
      if (prtlev.gt.1) then
      write(6,'(/A)') 
     .' Norm of eigenvectors in cartesian coordinates (one expected):'
      do ivec=1,nvec
         norm(ivec)=0.d0
         do i=1,natom
         ii=3*i-2
         dvec=covec(ii,ivec)**2+covec(ii+1,ivec)**2
     .       +covec(ii+2,ivec)**2
         norm(ivec)=norm(ivec)+dvec
         enddo
         norm(ivec)=dsqrt(norm(ivec))
       enddo
       write(6,'(5F8.5)') (norm(ivec),ivec=1,nvec)
      endif
 
      if (prtlev.gt.1) then
      write(6,'(/A)') 
     .' Orthogonality of first eigenvectors (zero expected):'
      do i=2,min(10,nvec)
         do j=1,i-1
            norm(j)=0.d0
            do k=1,natom
               ii=3*k-2
               norm(j)=norm(j)+
     .         covec(ii,i)*covec(ii,j)+covec(ii+1,i)*covec(ii+1,j)+
     .         covec(ii+2,i)*covec(ii+2,j)
            enddo
         enddo
         write(6,'(A,I3,A,10F6.3)') 
     .       ' Vector ',i,':',(norm(j),j=1,i-1)
      enddo
      endif

      do i=1,nvec
         write(uneig2,'(A,I5,7X,A,1PG12.4)')
     .       ' VECTOR',i,'VALUE',freq1(i)
         write(uneig2,'(1X,35(1H-))')
         write(uneig2,'(3(1PG12.4))')
     .  (covec(ii,i),covec(ii+1,i),covec(ii+2,i),ii=1,3*natom,3)
      enddo
 
      close(35)
      close(uneig1)
      close(uneig2)

      write(6,'(/A,I6,A)') program,nvec,
     .    ' eigenvectors saved.'

      if (prtlev.gt.1) 
     .write(6,'(/2A)') program,' Normal end.'
      return
      end
c
      subroutine rdmodfacs(uneig,nddlmax,nvec,numvec,freq,
     .           matvec,nddl,prtlev)
c
c     Lecture de modes au format "CERFACS".
c
cI/O:
      integer numvec(*), nddlmax, nvec, nddl, uneig
      double precision freq(*), matvec(nddlmax,*)
cLocal:
      integer nmotsmax
      parameter(nmotsmax=100)
      integer nerr, ivec, indnm_cfacs, nmots,
     .        i, ii, k, prtlev
      logical qfound, qold, qfirst
      character*1 carnum
      character program*11, progrer*14, progrwn*14
      character*132 lign132, mots(nmotsmax)
cDefaut:
      program=' Rdmodfacs>'
      progrwn='%Rdmodfacs-Wn>'
      progrer='%Rdmodfacs-Er>'

      if (prtlev.gt.1)
     .    write(6,'(/2A)') program,' Entering in.'

      nerr=0
      qfirst=.true.
      qold=.false.
      qfound=.false.
 100  continue
      read (uneig,'(A)',end=300,err=110) lign132
      goto 120
 110  continue
      nerr=nerr+1
 120  continue
 
      qfound=qfound.or.
     .      (index(lign132,' value ').gt.0.and.
     .       index(lign132,' vector ').gt.0.and.
     .       index(lign132,' residual ').le.0)
      qold=qold.or.
     .      (index(lign132,' VALUE ').gt.0.and.
     .       index(lign132,' VECTOR ').gt.0)
 
      if (.not.qfound.and..not.qold) goto 100
c________________________________________
c
c     Lecture des frequences des modes :
c________________________________________
 
      if (qfirst) then
          if (prtlev.gt.1) then
          if (qold) then
          write(6,'(/A)') 
     .  ' Rdmodfacs> Old Blzpack file format detected.'
          else
          write(6,'(/A)') 
     .  ' Rdmodfacs> Blzpack file format detected.'
          endif
          endif
          qfirst=.false.
      endif
c
      ivec=0
 250  continue
      ivec=ivec+1
      if (ivec.gt.nvec) then
          write(6,'(/A,I5,A)') 
     .  '%Rdmodfacs-Wn> More than ',nvec,' vectors in file.'
          return
      endif
c
      read(lign132,'(7X,I5,12X,G12.4)',end=240,err=240)
     .     numvec(ivec), freq(ivec)
 
      goto 255
 240  continue
      write(6,'(/3A)')
     .    '%Rdmodfacs-W> Pb with ligne: ',lign132(1:36),'...'
 255  continue
c
      if (prtlev.gt.1) then
      write(6,'(/A,I6)')
     .' Rdmodfacs> Eigenvector number:',
     .  numvec(ivec)
      write(6,'(A,1PG12.4)')
     .' Rdmodfacs> Corresponding eigenvalue:',
     .  freq(ivec)
      endif
c
      if (numvec(ivec).le.0)
     .    write(6,'(/A/A)')
     .  '%Rdmodfacs-W> Vector number was expected in:',
     .    lign132
c
      read(uneig,'(A)',end=230,err=230) lign132
 230  continue
      read(lign132,'(1X,A1)',end=232,err=232) carnum
 232  continue
      if ((qfound.and.carnum.ne.'=').or.
     .    (qold.and.carnum.ne.'-')) then
          write(6,'(2A/A)')
     .       ' %Rdmodfacs-Warning> Unexpected character ',
     .       ' in second column of line:',
     .    lign132
      endif
c____________________________________________________
c
c     2) Lecture des coordonnees des modes:
c        Format libre.
c____________________________________________________
c
      k=0
 257  continue
      if (k.gt.nddlmax) then
          write(6,'(/A,I6,A,I5)') 
     .  '%Rdmodfacs-Err> More than ',nddlmax,
     .  ' coordinates for vector ',ivec
          return
      endif
c
      read(uneig,'(A)',end=300,err=270) lign132
c
c     Nombre de coordonnees par ligne:
      call string_split(lign132,132,' ',
     .                  mots,nmotsmax,nmots)
c
      if (lign132.eq.' ') then
          read(uneig,'(A)',end=300,err=260) lign132
      else if (.not.qold.or.index(lign132,' VALUE ').le.0) then
          read(lign132,*,end=258)
     .   (matvec(k+ii,ivec),ii=1,nmots)
          k=k+nmots
          goto 257
 258      continue
      endif
      nddl=k
c
 260  continue
      indnm_cfacs=index(lign132,'       VALUE')
      if (indnm_cfacs.le.0)
     .    indnm_cfacs=index(lign132,'       value')
      if (indnm_cfacs.gt.0) then
          goto 250
      else
          write(6,'(A,A/A/A)')
     .  ' Rdmodfacs> Modes read.',
     .  ' Item VALUE not found in ligne :',lign132
          goto 300
      endif
c
 270  continue
      write(6,'(A,I6)')
     .   ' %Rdmodfacs-Er> While reading coordinate ',
     .      i
      stop '*Wrong matrix*'
c
 220  continue
c*****Ligne suivante de la lecture du fichier des modes en cours :
c
      goto 100
c
c     Fin de la lecture du fichier des modes :
c
 300  continue
      return
      end
c
      subroutine string_split(chaine,taille,delimiteur,
     .                        souschaine,nbremax,nbre)
c
c     "Chaine" est coupee en "nbre" "souschaine" de part et d'autre du
c     "delimiteur"
c I/O:
      integer taille, nbremax, nbre
      character*(*) chaine, souschaine(*), delimiteur
c Local:
      integer icar, iprev
c
      nbre=1
      iprev=1
      souschaine(1)=chaine
      do icar=1,taille
         if (chaine(icar:icar).eq.delimiteur) then
            if (icar-1.ge.iprev) then
               souschaine(nbre)=chaine(iprev:icar-1)
               nbre=nbre+1
               if (nbre.le.nbremax) then
                  if (icar+1.le.taille.and.
     .               chaine(icar+1:taille).ne.' ') then
                     souschaine(nbre)=chaine(icar+1:taille) 
                  else
                     nbre=nbre-1
                     return
                  endif
               else
                  write(6,'(A,I6,A/A)') 
     .               ' %String_split-Err: more than ',nbremax,
     .               ' substrings in : ',chaine
                  return
               endif
            endif
            iprev=icar+1
         endif
      enddo
c
      return
      end

      subroutine mintomaj(chaine)
 
c     Les caracteres minuscules sont mis en MAJUSCULES.
c     Les autres ne sont pas touches.
 
      character*(*) chaine
c Local:
      integer icar, ilettre, taille
      character*26  carmaj, carmin
 
      carmin='qwertyuiopasdfghjklzxcvbnm'
      carmaj='QWERTYUIOPASDFGHJKLZXCVBNM'
 
      taille=len(chaine)
      if (taille.le.0) return

      do icar=1,taille
         ilettre=index(carmin,chaine(icar:icar))
         if (ilettre.gt.0) then
             chaine(icar:icar)=carmaj(ilettre:ilettre)
         endif
      enddo
 
      return
      end
c----------------------------------------------------------------------
      SUBROUTINE openam(namfil,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c
c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT,
c     a priori suite a une interrogation...
c
c     input:
c        namfil: nom du fichier a ouvrir. 
c        "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.
c     output: 
c        qexist: flag / existence du fichier 
c        qinterr: Pas de nom pour le fichier cherche.
c I/O:
      logical qinterr, qverbos, qexist
      integer unit
      character*(*) namfil, cformat, cstatus
c Local
      character*132 ordrunix
c begin:
      if (cstatus.eq.'old') cstatus='OLD'
      if (cstatus.eq.'new') cstatus='NEW'
      if (cstatus.eq.'ove') cstatus='OVE'
      if (cstatus.eq.'unknown') cstatus='UNKNOWN'
 
      qinterr=.false.
      qexist=.false.
 
      if (namfil.eq.' ') then 
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> No filename.'
          return
      endif
 
      if (namfil.eq.'stop'.or.namfil.eq.'end'                         
     &    .or.namfil.eq.'fin'.or.namfil.eq.'quit') then 
         write(6,'(A)') 'Openam> Program is stopping on user request.'
         stop                                                                   
      endif 
 
c     Checks if filename is consistent with the opening:
 
      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          if (qverbos) write(6,'(A)') '%Openam-Err> File not found.'
          return
      endif
 
      if (qexist.and.cstatus.eq.'NEW') then
         write(6,'(/A)') 
     .      '%Openam-Err> This file exists:',namfil
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
                                                                    
      if (qverbos) then
         write(6,'(/A,I6,A)')
     .           ' Openam> file on opening on unit ',unit,':'
         write(6,*) namfil
      endif
      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)                
         
      return  
      end

      subroutine trier(y,npoint,nmax,ysort,iord,qcrois)
c
c     Tri par ordre croissant (qcrois=T) ou non.

      implicit none
      logical qcrois
      integer i, icur, iord(*), j, nmax, npoint
      double precision y(*), ycur, ysort(*)
      character progrer*10

      progrer='%Trier-Er>'

      if (npoint.gt.nmax) then
          write(6,'(/A,I9,A,I9,A)') progrer,npoint,
     .  ' points to be sorted, i.e., more than ',nmax,' Sorry.'
          stop '*Unexpected big problem*'
      endif

      do i=1,npoint
         ysort(i)=y(i)
         iord(i)=i
      enddo

      do i=1,npoint
        do j=1,npoint
          if (qcrois) then
            if (ysort(i).lt.ysort(j)) then
                ycur=ysort(i)
                icur=iord(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
                iord(i)=iord(j)
                iord(j)=icur
            endif
          else
            if (ysort(i).gt.ysort(j)) then
                ycur=ysort(i)
                icur=iord(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
                iord(i)=iord(j)
                iord(j)=icur
            endif
          endif
        enddo
      enddo

      return
      end
 
      subroutine stringcl(chaine,nonblancs)
c
c     Les caracteres "blancs" de la CHAINE sont retires (a gauche et au milieu).
c     L'entier NONBLANCS donne la position du dernier caractere.
c
c I/O:
      integer nonblancs
      character*(*) chaine
c Local:
      integer icar, ncar, taille
c Begin:
      taille=len(chaine)
c
      if (index(chaine(1:taille),' ').le.0) then
          nonblancs=taille
          return
      endif
c
c*****Nettoyage des blancs a gauche.
c     Premier non-blanc:
c
      do icar=1,taille
         if (chaine(icar:icar).ne.' ') goto 150
      enddo
 150  continue
      chaine=chaine(icar:taille)
c
c*****Nettoyage des blancs au milieu.
c
          icar=1
          ncar=1
 170      continue
          icar=icar+1
          ncar=ncar+1
          if (chaine(icar:icar).eq.' ') then
              chaine=chaine(1:icar-1)//chaine(icar+1:taille) 
              icar=icar-1
          endif
          if (ncar.lt.taille-1) goto 170
c
      nonblancs=index(chaine,' ')-1
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine vecstat(vect,nmax,rmin,rmax,rave,rdev)

c     Statistics for vector Vect(NMAX):
c     minimum, maximum, average (rave) and standard deviation (rdev).

cI/O:
      integer nmax
      double precision rave, rdev, rmax, rmin, vect(*)
cLocal:
      integer i
      character program*9, progrer*11, progrwn*11

cBegin:
      program=' Vecstat>'
      progrer='%Vecstat-Er>'
      progrwn='%Vecstat-Wn>'

      rave=0.d0
      rdev=0.d0
      rmin=-9999.d0
      rmax=9999.d0

      if (nmax.le.0) then
          write(6,'(/2A)') progrer,' Zero-length vector.'
          return
      endif

      do i=1,nmax
         if (vect(i).gt.rmax.or.i.eq.1) rmax=vect(i)
         if (vect(i).lt.rmin.or.i.eq.1) rmin=vect(i)
         rave=rave+vect(i)
         rdev=rdev+vect(i)**2.0
      enddo

      rave=rave/dfloat(nmax) 
      rdev=rdev/dfloat(nmax)-rave*rave
      if (rdev.gt.0.d0) rdev=dsqrt(rdev)

      return
      end
c-----------------------------------------------------------------------------
      subroutine readrlin(uninp,x,w,ic,nmax,ncoor,ndat,qerror,prtlev)

c     Reads at most NMAX coordinates in free format. 
c     They are stored in a linear array.

c     Either:
c     x, y, z
c     or:
c     x, y, z, w
c     or:
c     x, y, z, w, ic

c     If first word in ligne is not a number, the whole ligne is
c     assumed to be a title or a commentary.

cI/O:
      logical qerror
      integer ic(*), ncoor, ndat, nmax, prtlev, uninp
      double precision w(*), x(*), xc, yc, zc
cLocal:
      integer nmotsmax
      parameter(nmotsmax=255)

      integer i, lmot, nchi, nlmax, nl, nlu, nmots, stats(0:nmotsmax)
      double precision rlu
      character chiffres*15, lignlg*(nmotsmax), 
     .        mots(nmotsmax)*(nmotsmax), program*9, progrer*11, 
     .        progrwn*11

cBegin:
      program=' Readxyz>'
      progrer='%Readxyz-Er>'
      progrwn='%Readxyz-Wn>'

      chiffres='1234567890.eE+-'

      do i=1,nmotsmax
         stats(i)=0
      enddo

c     Lecture ligne a ligne:
c     ----------------------

      if (prtlev.gt.1) write(6,'(/2A)') program,
     .  ' Comments, or lines with less than three numbers: '

      qerror=.false.
      ncoor=0
      nl=0
 100  continue
      read(uninp,'(A)',end=200) lignlg 
      nl=nl+1
      call string_split(lignlg,nmotsmax," ",mots,nmotsmax,nmots)

      nfound=0
      do i=1,nmots
         call stringcl(mots(i),lmot)
         if (lmot.le.0) goto 150

c        Commentaire ?
         if (mots(i)(1:1).eq.'!') goto 150

c        Chiffre ?

         do k=1,lmot
            if (index(chiffres,mots(i)(k:k)).le.0) goto 110
         enddo

         nfound=nfound+1

         if (nfound.le.4) then
             read(mots(i)(1:lmot),*,err=110) rlu
             if (nfound.eq.1) xc=rlu
             if (nfound.eq.2) yc=rlu
             if (nfound.eq.3) zc=rlu
             if (nfound.eq.4) wc=rlu
         else if (nfound.eq.5) then 
             read(mots(i)(1:lmot),*,err=110) nlu
         endif

c        Mot suivant:
 110     continue

c        Le premier mot n'est pas un chiffre => ligne de commentaires

         if (nfound.eq.0) goto 150
      enddo
 150  continue

c     Stockage des coordonnees:
c     -------------------------

      stats(nfound)=stats(nfound)+1

      if (nfound.ge.3) then
          ncoor=ncoor+1
          if (ncoor.le.nmax) then
              n3=3*ncoor-2
              x(n3)  =xc
              x(n3+1)=yc
              x(n3+2)=zc
              if (nfound.eq.4) then
                  w(n3)=wc
                  w(n3+1)=wc
                  w(n3+2)=wc
              else
                  w(n3)=wc
                  w(n3+1)=wc
                  w(n3+2)=wc
                  ic(ncoor)=nlu 
              endif
          else
              write(6,'(/2A,I9,A)') progrer,' More than ',
     .        nmax,' particles in file.'
              write(6,'(2A)') progrer,
     .      ' Please increase program memory limits (Sorry for that).'
              stop '*Too large system*'
          endif
      else
          if (prtlev.gt.1) then
              write(6,'(2A)') lignlg(1:72),'...'
          endif
      endif
      
c     Ligne suivante:
      goto 100
 200  continue
      if (prtlev.gt.1.and.nl.eq.ncoor) write(6,'(2A)') program,' None.'

      write(6,'(/2A,I7)') program,
     .' Number of particles in file (with x,y,z coordinates): ',ncoor
 
      if (ncoor.eq.0) then
          write(6,'(/2A)') progrer,' No coordinate found in file.'
          qerror=.true.
      endif

      nchi=0
      ndat=0
      nlmax=0
      do i=1,nmotsmax
         if (stats(i).gt.0) nchi=nchi+1
         if (stats(i).gt.nlmax) then
             nlmax=stats(i)
             ndat=i
         endif
      enddo

      do i=0,nmotsmax
         if (stats(i).gt.0.and.(prtlev.gt.1.or.nchi.gt.1)) then
            write(6,'(A,I6,A,I7,A)') program,i,
     .    ' numbers found in ',stats(i),' lines.'  
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------------
