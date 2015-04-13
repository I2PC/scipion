      program pdbmat
      implicit none
      integer natmax, nresmx, nvoismx
 
   
      parameter( NATMAX=55000 )
      parameter( NRESMX=55000 )
      parameter( NVOISMX=1000000 )
    

c     YHS-Nov-1996: Version 1.00.
c     Versions released at http://ecole-modelisation.free.fr/modes.html
c     YHS-Mar-2001: Version 3.31. 
c     YHS-Feb-2004: Version 3.50. 
c     YHS-Feb-2008: Version 3.73. 
c     Version used by the ELNEMO Web site (http://www.elnemo.org):
c     YHS-Feb-2004: Version 3.46.
c
      integer nmotsmax, ntopmax
      parameter(ntopmax=10*natmax,nmotsmax=100)
 
      integer fatres(nresmx+1), 
     .        i, idmax, idres(nresmx), ii, imax, imin, ires, 
     .        iresat(natmax), iseed, ivois(nvoismx),
     .        j, jangle(ntopmax), jat, jbond(ntopmax), jdihe(ntopmax), 
     .        jj, jsomang(ntopmax), jvois(nvoismx), 
     .        k, kcom, kk, klist, 
     .        l, ll, lmot, lnom, lnomlst, lnommls, lnommtx, lnompdb, 
     .        lnomvmd,
     .        namax, namin, nangle(ntopmax), nangles, natom, 
     .        nbig, nbmax, nbmin, nbond(ntopmax), nbonds, ndat, 
     .        ndihe(ntopmax), ndihs, nl, nmax, nmin, nmots, nntr, 
     .        nnzero, nres, nunit, nunknown, nvois, nvoisat(natmax), 
     .        prtlev, uninp, unlst, unmol, unout, unpdb, unrsd, unvmd
      double precision cutbnd, cutoff, ddf, der2(3,3*natmax), 
     .        dist, dist2, dmax, dmin, dmoy, drms, 
     .        elemnt, elmax, fvois(natmax), 
     .        kangle, kbond, kdihe, kfce, kij, knonb, kvois(natmax),
     .        levelshft, massat(natmax), nmoy, nrms,
     .        random, rave, rbig, rdev, rinput, rkh, rmax, rmin, rsmall, 
     .        rx, ry, rz, 
     .        trace, unknown, xat(natmax), yat(natmax), zat(natmax) 
      logical qbinary, qerror, qexist, qfread, qinter, qlist, 
     .        qmasse, qmtx, qok, qpdb, qvois(natmax)
      character atonam(natmax)*4, cformat*32, csep*1, cstatus*32, 
     .        lign80*80, motinp*80, mots(nmotsmax)*132, 
     .        nomfich*64, nomlst*64, nommls*64, nommtx*64, nompdb*64, 
     .        nomvmd*64,
     .        program*8, progrer*11, progrwn*11, 
     .        residus_standards*132, residus_stshort*21, 
     .        resnam(natmax)*4, segid(natmax)*4, ssunam(natmax)*1, 
     .        ssusel*1, typbond*80, typmas*80, typout*80, version*32
      parameter(rbig=1e10,rsmall=1e-10,unknown=9999.d9)
c.......................................................................
      version=' Version 3.73, February 2008.'
c.......................................................................
      idmax=21
      residus_standards='   ILE PHE TRP LEU CYS VAL MET TYR ALA HIS '//
     .                     'GLY THR SER PRO ARG GLN ASN ASP GLU LYS '
      residus_stshort='IFWLCVMYAHGTSPRQNDEKX'
 
      program=' Pdbmat>'
      progrer='%Pdbmat-Er>'
      progrwn='%Pdbmat-Wn>'

      write(6,'(2A)') program,
     .' Computes the Hessian matrix, using an Elastic Network Model.'
      write(6,'(2A)') program,version

c     ===============================================
c     Ouverture et Lecture du fichier d'instructions:
c     ===============================================
 
c     Default values:
      nompdb='pdbmat.ent'
      cutoff=10.d0
c     Neighbor list:
      nomlst='NONE'
      qlist=.false.
c     Typical distance for Hinsen's weigth:
c    (negative value means 1/rkh=0)
      rkh=-1.d0
      qfread=.false.
      nommls='NONE'
      nomvmd='NONE'
      typmas='CONS'
      qmasse=.false.
      knonb=1.0d0
c     Topological terms (like in ref.6):
      typbond='NONE'
c     Distance-cutoff value for defining covalent (chemical) bonds:
      cutbnd=4.0d0
      kangle=0.0d0
      kbond=1000.0d0
      kdihe=0.0d0
c     Others:
      typout='  FREE'
      qmtx=.false.
      qbinary=.false.
      prtlev=0
      levelshft=1e-8
c     Used only with a levelshift, to add some noise (not useful ?).
      iseed=27041961

      nunit=10
      uninp=nunit
      nunit=nunit+1
      nomfich='pdbmat.dat'
      cformat="FORMATTED"
      cstatus="old"
      call openam(nomfich,cformat,cstatus,uninp,.false.,
     .     qinter,qexist)
      if (qinter.or..not.qexist) then 
          write(6,'(/2A/(2A))') progrwn,
     .  ' No pdbmat.dat command file found.',
     .    progrwn,' Defaults assumed for all options. ',
     .    progrwn,' See the pdbmat.dat_run file if you need an example.'
          goto 110
      else
          write(6,'(/2A)') program,
     .  ' Options to be read in pdbmat.dat file.'
      endif
 
 50   continue
      read(uninp,'(A)',end=100) lign80
 
      kcom=index(lign80,'!')
      k=index(lign80,'=') 
      motinp=' '
      if (k.gt.0.and.(kcom.le.0.or.kcom.gt.k)) then 
          motinp=lign80(1:k)
      else
          if (k.le.0.and.kcom.gt.1) then
          write(6,'(/2A/A)') progrwn,
     .  ' No separator (=) in command ligne:',
     .    lign80
          write(6,'(2A)') progrwn,' This ligne is skipped.'
          endif
          goto 50
      endif
      call mintomaj(motinp)

      kcom=index(lign80(k+1:80),'!')
      if (kcom.gt.0) lign80(k+kcom:80)=' '
      klist=index(lign80(k+1:80),'?')

      if (index(motinp,' FILENAME').gt.0.or.
     .    index(motinp,' SCRIPT').gt.0) then 
          if (index(motinp,'MATRI').gt.0) then
              nommtx=lign80(k+1:80)
              qmtx=.true.
          elseif (index(motinp,'LIST ').gt.0.or.
     .            index(motinp,'NEIGH').gt.0) then
              nomlst=lign80(k+1:80)
              qlist=.true.
          elseif (index(motinp,'MOLS').gt.0) then
              nommls=lign80(k+1:80)
          elseif (index(motinp,'VMD').gt.0) then
              nomvmd=lign80(k+1:80)
          else
              nompdb=lign80(k+1:80)
          endif
      else if (index(motinp,' DEFINITION').gt.0) then
          typbond=lign80(k+1:80)
          call mintomaj(typbond)
          call stringcl(typbond,lnom)
          if (typbond(1:3).eq.'ALL') then 
              typbond=' ALL'
          else if (typbond(1:3).eq.'NON') then 
              typbond='NONE'
          else if (typbond(1:3).eq.'CON') then
              typbond='CONSECUTIF'
          else
              write(6,'(/3A)') progrwn,' Bond definition :',
     .        typbond(1:4)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: NONe, ALL, CONsecutive.'
              write(6,'(A)') ' Default assumed.'
              typbond='NONE'
          endif
      else if (index(motinp,'MASS').gt.0) then
          typmas=lign80(k+1:80)
          call mintomaj(typmas)
          call stringcl(typmas,lnom)
          if (typmas(1:3).eq.'PDB'.or.typmas(1:3).eq.'COO') then
              qmasse=.true.
              typmas='COOR'
          else if (typmas(1:3).ne.'CON') then
              write(6,'(/3A)') progrwn,' Origin of mass values :',
     .        typmas(1:3)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: CONstant, COOr, PDB.'
              write(6,'(A)') ' Default assumed.'
              qmasse=.false.
              typmas='CONS'
          endif
      else if (index(motinp,'FORMAT').gt.0) then
          typout=lign80(k+1:80)
          call mintomaj(typout)
          call stringcl(typout,lnom)
          if (typout(1:1).eq.'B'.or.typout(1:1).eq.'U') then
              qbinary=.true.
              typout='BINARY'
          else if (typout(1:1).ne.'F') then
              write(6,'(/3A)') progrwn,' Kind of matrix format :',
     .        typout(1:1)
              if (klist.le.0)
     .        write(6,'(2A)') progrwn,' This is not a known keyword.'
              write(6,'(2A)') progrwn,
     .      ' Valid options are: Free, Binary, Formatted, Unformatted.'
              write(6,'(A)') ' Default assumed.'
              qbinary=.false.
              typout='  FREE'
          else
              qbinary=.false.
              typout='  FREE'
          endif
      else 
          qok=.false.
          read(lign80(k+1:80),*,end=90,err=90) rinput
          if (index(motinp,'SHIFT ').gt.0) then 
               qok=.true.
               levelshft=rinput
          else if (index(motinp,'CUTOF').gt.0.or.
     .             index(motinp,'DISTANCE').gt.0) then
               qok=.true.
               cutoff=rinput
          else if (index(motinp,'HINSEN').gt.0) then
               qok=.true.
               rkh=rinput
          else if (index(motinp,'INTERAC').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   knonb=rinput
               endif
               if (index(motinp,' CUTOF').gt.0.or.
     .             index(motinp,' DIST').gt.0) then
                   qok=.true.
                   cutoff=rinput
               endif 
          else if (index(motinp,'BOND').gt.0.and.
     .        (index(motinp,' FORCE ').gt.0.or.
     .         index(motinp,' CONST').gt.0)) then
               qok=.true.
               kbond=rinput
          else if (index(motinp,' LENGTH').gt.0) then
               qok=.true.
               cutbnd=rinput
          else if (index(motinp,'PRINT').gt.0) then
               qok=.true.
               prtlev=int(rinput)
          else if (index(motinp,'ANGLE').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   kangle=rinput
               endif
          else if (index(motinp,'DIHE').gt.0) then
               if (index(motinp,' FORCE ').gt.0.or.
     .             index(motinp,' CONST').gt.0) then
                   qok=.true.
                   kdihe=rinput
               endif
          endif
  90      continue
          if (.not.qok) then
               write(6,'(/2A/A)') progrwn,
     .       ' No known or incomplete set of keywords in ligne:',
     .         motinp
               write(6,'(2A)') progrwn,
     .       ' This command ligne is skipped.'
          endif
      endif
      goto 50
 
 100  continue
      close(uninp)
 110  continue

      call stringcl(nompdb,lnompdb)
      call stringcl(nomlst,lnomlst)
      call stringcl(nommls,lnommls)
      call stringcl(nomvmd,lnomvmd)
      if (nomlst.eq.'none'.or.nomlst.eq.'NONE') qlist=.false.
      if (nommls.eq.'none') nommls='NONE'
      if (nomvmd.eq.'none') nomvmd='NONE'
      
      if (.not.qmtx) then
      if (qbinary) then
        nommtx="pdbmat.sdijb"
      else
        nommtx="pdbmat.sdijf"
      endif
      endif
      call stringcl(nommtx,lnommtx)

c     Resume des commandes:
c     ---------------------

      write(6,'(/3A)') program,' Coordinate filename     = ',
     .      nompdb(1:lnompdb)
      if (qlist) then
      write(6,'(3A)') program,' Neighbor-list filename  = ',
     .      nomlst(1:lnomlst)
      else
      write(6,'(/2A,F10.2)') program,
     .        ' Distance cutoff         = ',cutoff
      endif
      if (rkh.gt.0.d0)
     .write(6,'(8X,A,F10.2)') " Hinsen's typical range  = ",rkh
      write(6,'(A,F10.2)') 
     .'         Force constant          = ',knonb
      if (typbond.ne.'NONE') then
      write(6,'(A,6X,A)') 
     .'         Kind of bond definition = ',typbond(1:4)
      write(6,'(A,F10.2)') 
     .'         Maximum bond length     = ',cutbnd,
     .'         Bond force constant     = ',kbond,
     .'         Angle force constant    = ',kangle,
     .'         Dihedral force constant = ',kdihe
      endif
      write(6,'(A,6X,A)') 
     .'         Origin of mass values   = ',typmas(1:4)
      write(6,'(3A)') program,' Matrix filename         = ',
     .      nommtx(1:lnommtx)
      if (prtlev.gt.0) then
      write(6,'(2A,1PG10.1)') program,
     .        ' Levelshift              = ',levelshft
      write(6,'(A,3X,I7)') 
     .'         PRINTing level          = ',prtlev
      endif
      if (nommls.ne.'NONE')
     .write(6,'(3A)') program,' Molscript filename      = ',
     .      nommls(1:lnommls)
      if (nomvmd.ne.'NONE')
     .write(6,'(3A)') program,' VMD script filename     = ',
     .      nomvmd(1:lnomvmd)

c     Sauvegarde du fichier de commandes complet:
c     -------------------------------------------

      uninp=nunit
      nunit=nunit+1
      nomfich='pdbmat.dat_run'
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,cformat,cstatus,uninp,.false.,
     .     qinter,qexist)

c     Pour etre plus clair:
      if (typbond.eq.'NONE') then
          cutbnd=0.d0
          kbond=0.d0
          kangle=0.d0
          kdihe=0.d0
      endif
      if (qlist) cutoff=-1

      write(uninp,'(2A)') 
     .'! This file can be modified and used as a command file',
     .' (named pdbmat.dat) for pdbmat.'
      write(uninp,'(2A)') ' Coordinate FILENAME        = ',
     .      nompdb(1:lnompdb)
      write(uninp,'(2A)') ' MATRIx FILENAME            = ',
     .      nommtx(1:lnommtx)
      write(uninp,'(A,F10.3,A)') ' INTERACtion DISTance CUTOF = ',
     .      cutoff,' ! For defining the list of interacting atoms.'
      write(uninp,'(A,F10.3,A)') ' INTERACtion FORCE CONStant = ',
     .      knonb,' ! For specifying frequency units.'
      if (prtlev.gt.0.or.rkh.gt.0.d0)
     .write(uninp,'(A,F10.3,A)') " HINSEN's typical range     = ",
     .      rkh,' ! Force constant weighting (if negative: none).'
      if (prtlev.gt.0.or.nomlst(1:lnomlst).ne.'NONE')
     .write(uninp,'(3A)') ' NEIGHbor-list FILENAME     = ',
     .      nomlst(1:lnomlst),
     .  ' ! For defining this list yourself.'
      write(uninp,'(A,6X,2A)') ' Origin of MASS values      = ',
     .      typmas(1:4),' ! CONstant, or from COOrdinate file.'
      write(uninp,'(A,8X,I2,A)') ' Output PRINTing level      = ',
     .      prtlev,' ! =1: more detailled. =2: debug level.'
c     Not often used:
      if (prtlev.gt.0.or.nommls(1:lnommls).ne.'NONE') 
     .write(uninp,'(3A)') ' MOLScript command FILEname = ',
     .      nommls(1:lnommls),
     .  ' ! To draw the network with Molscript.'
      if (prtlev.gt.0.or.nomvmd(1:lnomvmd).ne.'NONE') 
     .write(uninp,'(3A)') ' VMD command FILEname       = ',
     .      nomvmd(1:lnomvmd),
     .  ' ! vmd -e this-file (to visualize the network with VMD).'
c     Rarely used:
      if (prtlev.gt.0.or.typbond.ne.'NONE') then
      write(uninp,'(A,6X,2A)') ' Bond DEFINITION            = ',
     .      typbond(1:4),' ! NONe, ALL, or between CONsecutive atoms.'
      write(uninp,'(A,F10.3)') ' Maximum bond LENGTH        = ',cutbnd
      write(uninp,'(A,F10.3)') ' BOND FORCE CONStant        = ',kbond
      write(uninp,'(A,F10.3)') ' ANGLE FORCE CONStant       = ',kangle
      write(uninp,'(A,F10.3)') ' DIHEdral FORCE CONStant    = ',kdihe
      write(uninp,'(A,1PG10.1,A)') ' LevelSHIFT                 = ',
     .      levelshft,
     .  ' ! Non-zero value often required (numerical reasons).'
      write(uninp,'(A,4X,2A)') ' Matrix FORMAT              = ',
     .      typout(1:6),' ! Free, or Binary, matrix saved.'
      endif
      close(uninp)

c     Tests:
      if ((.not.qlist.and.cutoff.lt.0.d0.and.rkh.lt.0.d0).or.
     .    knonb.lt.0.d0.or.
     .   (typbond(1:4).ne.'NONE'.and.(cutbnd.lt.0.d0.or.kbond.lt.0.d0
     .   .or.kangle.lt.0.d0.or.kdihe.lt.0.d0))) then
          write(6,'(/2A)') progrer,
     .  ' Distances and force constants can not have negative values !' 
          stop '*Commands are not consistent*'
      endif

c     On recherche l'information/sous-unite:
 
      call string_split(nompdb,lnompdb,":",
     .                  mots,nmotsmax,nmots)
      call stringcl(mots(1),lnom)
 
      if (nmots.gt.1) then
          call stringcl(mots(2),lnom)
          ssusel=mots(nmots)
          write(6,'(3A)') program,' Subunit to be selected: ',ssusel
          if (nmots.gt.2) then
              write(6,'(4A)') progrwn,' The end of filename, ',
     .        nompdb(1:lnompdb),', was not understood.'
          endif
      else
          ssusel=' '
      endif
      nompdb=mots(1)
      call stringcl(nompdb,lnompdb)
c                                          
c     Lecture du fichier de coordonnees:
c     ==================================

      if (prtlev.gt.0)
     .  write(6,'(/(4A))') program,
     .' Coordinate file ',nompdb(1:lnompdb),' to be opened.'

      unpdb=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"
      call openam(nompdb,cformat,cstatus,unpdb,.true.,
     .     qinter,qexist)
      if (qinter) stop '*No readable coordinate file found*'

c     Format pdb ?

      nl=0
      qpdb=.false.
 120  continue
      read(unpdb,'(A)',end=130) lign80
      if (lign80(1:5).eq.'ATOM '.or.lign80(1:6).eq.'HETATM') then
          qpdb=.true.
          goto 130
      else
          nl=nl+1
      endif
      goto 120
 130  continue
      rewind(unpdb)

      do i=1,natmax
         xat(i)=unknown
         yat(i)=unknown
         zat(i)=unknown
         massat(i)=unknown
         iresat(i)=i
      enddo

      if (qpdb) then
          write(6,'(/2A)') program,
     .  ' Coordinate file in PDB format.'
          call rdatompdb(unpdb,ssusel,xat,yat,zat,massat,
     .         atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .         fatres,nresmx,nres,qerror,prtlev)
          if (natom.eq.nres) then
              write(6,'(/2A)') program,' Study of a standard ENM model.'
          else
              write(6,'(/2A)') progrwn,
     .      ' Study of a several-atom-per-residue model '//
     .      '(this is not that standard).'
          endif
      else
          if (nl.eq.0) then
              write(6,'(/2A)') progrer,' Empty coordinate file.'
              stop
          endif
          write(6,'(/2A)') program,
     .  ' Coordinate file in Free format.'

          call readxyz(unpdb,xat,yat,zat,massat,iresat,natmax,natom,
     .         ndat,qerror,prtlev)

          if (qmasse.and.ndat.lt.4) then
              write(6,'(/2A)') progrer,
     .      ' Masses were not all found, as expected.'
              qmasse=.false.
          endif
      endif

c     Tests:
c     ======

      if (qerror) stop
      if (natom.le.1) then
          write(6,'(2A)') progrer,
     .  ' Not enough atoms found in file. Nothing done.'
          stop
      endif

      write(6,'(/2A)') program,' Coordinate statistics: '

      call vecstat(xat,natom,rmin,rmax,rave,rdev)
      write(6,'(4(A,F12.6))') 
     .' <x>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

      call vecstat(yat,natom,rmin,rmax,rave,rdev)
      write(6,'(4(A,F12.6))') 
     .' <y>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

      call vecstat(zat,natom,rmin,rmax,rave,rdev)
      write(6,'(4(A,F12.6))') 
     .' <z>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

      if (qmasse) then
          write(6,'(/2A)') program,' Mass statistics: '
          call vecstat(massat,natom,rmin,rmax,rave,rdev)
          write(6,'(4(A,F12.6))') 
     .  ' <m>= ',rave,' +/- ',rdev,' From: ',rmin,' To: ',rmax

          if (rmin.le.0.d0) then
              write(6,'(2A)') progrer,
     .      ' Negative or null masses found !'
              qmasse=.false.
          endif
      endif

      if (.not.qmasse) then
          write(6,'(2A)') program,' Masses are all set to one.'
          do i=1,natom
             massat(i)=1.d0
          enddo
      endif
 
c     Test/identification des residus.
 
      if (qpdb) then
      nunknown=0
      do i=1,nres
         ires=fatres(i)
         idres(i)=index(residus_standards,resnam(ires))/4
         if (idres(i).le.0) then
             nunknown=nunknown+1
             if (nunknown.lt.10) then
                 write(6,'(4A)') progrwn," residue:'",
     .           resnam(ires),"' is not a well known amino-acid."
                 idres(i)=idmax
             else if (nunknown.eq.10) then
                 write(6,'(2A)') progrwn,' ........'
                 idres(i)=idmax
             endif
         endif
      enddo
      if (nunknown.gt.0) 
     .write(6,'(/A,I6,A)') progrwn,nunknown,' residue(s) not known.'
      endif
 
c     Bonds i-j and j-i are stored, because below the matrix is
c     calculated and saved three lines at a time.
 
      nbonds=0
      if (typbond.ne.' ALL'.and.typbond.ne.'CONSECUTIF') goto 200
 
      if (typbond.eq.' ALL') then
         nbmax=0
         nbmin=999
         imax=-1
         imin=-1
         k=1
         nbond(1)=1
         do i=1,natom
            do j=1,natom
               if (i.ne.j) then
                   rx=xat(i)-xat(j)
                   ry=yat(i)-yat(j)
                   rz=zat(i)-zat(j)
                   dist=dsqrt(rx*rx + ry*ry + rz*rz)
                   if (dist.le.cutbnd) then
                       jbond(k)=j
                       k=k+1
                   endif
               endif
            enddo
            nbond(i+1)=k
            if (nbond(i+1)-nbond(i).gt.nbmax) then
                nbmax=nbond(i+1)-nbond(i)
                imax=i
            endif
            if (nbond(i+1)-nbond(i).lt.nbmin) then
                nbmin=nbond(i+1)-nbond(i)
                imin=i
            endif
            if (k-1.gt.ntopmax) then
                write(6,'(/2A,I12)') progrer,
     .        ' Too many bonds. Maximum is: ',ntopmax
                stop
            endif
         enddo
         nbonds=k-1
      else if (typbond.eq.'CONSECUTIF') then
 
c        On fait attention aux distances...
c        Il peut y avoir plusieurs molecules,
c        plusieurs chaines, dans le systeme.
 
         k=1
         do i=1,natom
            nbond(i)=k
            if (i.gt.1) then
            j=i-1
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist=dsqrt(rx*rx + ry*ry + rz*rz)
            if (dist.le.cutbnd) then
                jbond(k)=j
                k=k+1
            endif
            endif
            if (i.lt.natom) then
            j=i+1
            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist=dsqrt(rx*rx + ry*ry + rz*rz)
            if (dist.le.cutbnd) then
                jbond(k)=j
                k=k+1
            endif
            endif
            if (k.gt.ntopmax) then
                write(6,'(/2A,I12)') progrer,
     .        ' Too many bonds. Maximum is: ',ntopmax
                stop
            endif
         enddo
         nbond(natom+1)=k
         imax=2
         imin=1
         nbmin=1
         nbmax=2
         nbonds=k
      endif
 
      if (nbonds.eq.0) then
          write(6,'(/2A/)') progrwn,' No bond found.'
          goto 200
      endif
 
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  nbonds/2,' covalent bonds, i.e.,',
     .  float(nbonds)/float(2*natom),' per atom.'

      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)') 
     .'         Maximum number found =',nbmax,' for atom ',imax,
     .'         Minimum number found =',nbmin,' for atom ',imin
 
 
      nangles=0
      if (kangle.le.0.d0) then
          write(6,'(/2A)') program,' BOND but no ANGLe energy term.'
          goto 200
      endif

      namax=0
      namin=9999
      imax=-1
      imin=-1
      ii=1
      nangle(1)=1
      do i=1,natom
         if (nbond(i+1).gt.nbond(i)) then
         do jj=nbond(i),nbond(i+1)-1
            j=jbond(jj)
            if (nbond(j+1).gt.nbond(j)) then
            do kk=nbond(j),nbond(j+1)-1
               k=jbond(kk)
               if (k.ne.i) then
                   jangle(ii)=k
                   jsomang(ii)=j
                   ii=ii+1
               endif
            enddo
            endif
         enddo
         endif
         nangle(i+1)=ii
         if (nangle(i+1)-nangle(i).gt.namax) then
             namax=nangle(i+1)-nangle(i)
             imax=i
         endif
         if (nangle(i+1)-nangle(i).lt.namin) then
             namin=nangle(i+1)-nangle(i)
             imin=i
         endif
         if (ii.gt.ntopmax) then
             write(6,'(/2A,I12)') progrer,
     .     ' Too many angles. Maximum is: ',ntopmax
             stop
         endif
      enddo
      nangles=ii-1
 
      if (nangles.eq.0) then
          write(6,'(/2A/)') progrwn,' No bond-angle found.'
          goto 200
      endif
 
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  nangles/2,' valence angles, i.e.,',
     .  float(nangles)/float(2*natom),' per atom.'

      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)') 
     .'         Maximum number found =',namax,' for atom ',imax,
     .'         Minimum number found =',namin,' for atom ',imin

 
      ndihs=0
      if (kdihe.le.0.d0) then
         write(6,'(/2A)') program,' ANGLes but no DIHEdral energy term.'
         goto 200
      endif

      namax=0
      namin=9999
      imax=-1
      imin=-1
      ii=1
      ndihe(1)=1
c     For each atom:
      do i=1,natom
c        All angles where it is the first atom:
         if (nangle(i+1).gt.nangle(i)) then
         do jj=nangle(i),nangle(i+1)-1
c           For the "top-atom" of this angle:
            j=jsomang(jj)
c           All angles where it is the first atom:
            if (nangle(j+1).gt.nangle(j)) then
            do kk=nangle(j),nangle(j+1)-1
               k=jangle(kk)
               if (k.ne.i.and.j.ne.jsomang(kk).and.i.ne.jsomang(kk))then
                   jdihe(ii)=k
                   ii=ii+1
               endif
            enddo
            endif
         enddo
         endif
         ndihe(i+1)=ii
         if (ndihe(i+1)-ndihe(i).gt.namax) then
             namax=ndihe(i+1)-ndihe(i)
             imax=i
         endif
         if (ndihe(i+1)-ndihe(i).lt.namin) then
             namin=ndihe(i+1)-ndihe(i)
             imin=i
         endif
         if (ii.gt.ntopmax) then
             write(6,'(/2A,I12)') progrer,
     .     ' Too many dihes. Maximum is: ',ntopmax
             stop
         endif
      enddo
      ndihs=ii-1
 
      if (ndihs.eq.0) then
          write(6,'(/2A/)') progrwn,' No bond-dihedral found.'
          goto 200
      endif
 
      if (prtlev.gt.0)
     .write(6,'(A,I6,A,F5.2,A)') program,
     .  ndihs/2,' dihedral angles, i.e.,',
     .  float(ndihs)/float(2*natom),' per atom.'

      if (prtlev.gt.1)
     .write(6,'(A,I3,A,I6)') 
     .'         Maximum number found =',namax,' for atom ',imax,
     .'         Minimum number found =',namin,' for atom ',imin

c     ==============================
c     Matrice des derivees secondes:
c     ==============================                                
 200  continue

c     Lecture eventuelle de la liste des voisins:
c     -------------------------------------------
c     Formats possibles: 
c     atom-number atom-number 
c     atom-number atom-number force-constant

      if (qlist) then
      write(6,'(/2A)') program,' Neighbor list to be read.'
      unlst=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"
      call openam(nomlst,cformat,cstatus,unlst,.true.,
     .     qinter,qexist)

      nvois=0
      ndat=0
 210  continue
      read(unlst,'(A)',end=230,err=220) lign80

c     Les commentaires ne sont pas pris en compte:
      kcom=index(lign80,'!') 
      if (kcom.le.0) kcom=index(lign80,'#')
      if (kcom.gt.0) then
          if (kcom.eq.1) then
              write(6,'(2A)') lign80(1:50),' ...'
              goto 210
          else
              lign80=lign80(1:kcom-1)
          endif
      endif

c     Plusieurs separateurs sont possibles: , ; ou blanc
      csep=","
      k=index(lign80,csep)
      if (k.le.0) then
          csep=";"
          k=index(lign80,csep)
      endif
      if (k.le.0) csep=" "

      call string_split(lign80,80,csep,
     .     mots,nmotsmax,nmots)
      if (nmots.lt.2) goto 225 
      if (ndat.eq.0) then
          ndat=min(nmots,3)
          if (ndat.eq.3) qfread=.true.
      endif
      if (nmots.lt.ndat) then
          write(6,'(/A,I3,A/A)') progrer,ndat,
     .  ' data per ligne until ligne: ',lign80
          goto 220
      endif

c     Lecture de: i, j (kij le cas echeant).

      read(mots(1),*) ii
      read(mots(2),*) jj
      if (qfread) read(mots(3),*) kfce

      if (ii.le.0.or.jj.le.0) then
          write(6,'(/2A)') progrer,
     .        ' Null or negative atom number found.'
          goto 225
      endif
      if (ii.gt.natom.or.jj.gt.natom) then
          write(6,'(/2A,I6,A,I6,A)') progrer,
     .  ' Atom number: ',ii,' or ',jj,
     .  ' larger than the number of atoms.'
          stop '*Wrong file*'
      endif

      nvois=nvois+1
      if (nvois.gt.nvoismx) then
          write(6,'(/2A,I6,A)') progrer,' More than ',nvoismx,
     .  ' pairs of neighbors, the maximum allowed. Sorry.'
          stop '*Recompile with larger array*'
      endif

      ivois(nvois)=ii
      jvois(nvois)=jj
      if (qfread) fvois(nvois)=kfce

c     Ligne suivante:
      goto 210

c     Probleme de lecture:
 220  continue
      write(6,'(2A,I6)') progrer,
     .    ' While reading neigbors pair: ',nvois+1
      stop '*Wrong or corrupted file*'

 225  continue
      write(6,'(2A/A)') progrer,
     .' No (or wrong) pair of atom numbers found in ligne: ',lign80
      stop '*Wrong or corrupted file*'

c     Fin du fichier:
 230  continue     
      write(6,'(/A,I6,A)') program,nvois,' pairs of neighbors.' 
      endif

c     Coordonnees et masses utilisees, sauvegardees:
c     ----------------------------------------------

      unout=nunit
      nunit=nunit+1
      nomfich="pdbmat.xyzm"
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,cformat,cstatus,unout,.true.,
     .     qinter,qexist)

      do i=1,natom
         write(unout,'(4(1PG20.12),I9)')  
     .   xat(i), yat(i), zat(i), massat(i), iresat(i)
      enddo
      close(unout) 
      if (prtlev.gt.0)
     .write(6,'(2A)') program,
     .    ' Coordinates and masses considered are saved.'

c     Fichier de commandes pour VMD:
c     ------------------------------

      unvmd=-1
      if (nomvmd.ne.'NONE') then
      unvmd=nunit
      nunit=nunit+1
      nomfich=nomvmd
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,cformat,cstatus,unvmd,.true.,
     .     qinter,qexist)

      write(unvmd,'(A)') '#!/usr/local/bin/vmd'
      write(unvmd,'(A)') '# script for VMD (Visual Molecular Dynamics)'
      write(unvmd,'(A)') '# Goal: visualizing the elastic network'
      write(unvmd,'(A)') '# Type: vmd -e this-file'
      write(unvmd,'(A)') 'color Display {Background} white'
      write(unvmd,'(A)') 'mol new'
      write(unvmd,'(A)') 'draw color black'
      endif

c     Fichier de commandes pour Molscript:
c     ------------------------------------

      unmol=-1
      if (nommls.ne.'NONE') then
      unmol=nunit
      nunit=nunit+1
      nomfich=nommls
      cformat="FORMATTED"
      cstatus="ove"
      call openam(nomfich,cformat,cstatus,unmol,.true.,
     .     qinter,qexist)

      write(unmol,'(A)') '! Script for Molscript (Kraulis, 1993)'
      write(unmol,'(A)') '! Goal: visualizing the elastic network'
      write(unmol,'(A)')  ' set bonddistance 99.0 ;'
      endif

c     Matrice:
c     --------

      if (qbinary) then
        if (.not.qmtx) nommtx="pdbmat.sdijb"
        cformat="UNFORMATTED"
      else
        if (.not.qmtx) nommtx="pdbmat.sdijf"
        cformat="FORMATTED"
      endif
      unout=nunit
      nunit=nunit+1
      cstatus="ove"
      call openam(nommtx,cformat,cstatus,unout,.true.,
     .     qinter,qexist)

c     ========================================
c     Les atomes sont tous lies deux a deux,
c     par un potentiel "universel" (M.Tirion).
c     ========================================
 
      elmax=0.d0
      trace=0.d0
      dmin=0.d0
      dmax=0.d0
      dmoy=0.d0
      drms=0.d0
      nnzero=0
      nntr=0
      nbig=0
      ll=0

      do i=1,natom
         ii=3*i-2
         nvoisat(i)=0

c        Liste eventuelle des voisins de i:
c        ----------------------------------

         if (qlist) then
             do j=1,natom
                qvois(j)=.false.
                kvois(j)=0.d0
             enddo
             do j=1,nvois
                if (ivois(j).eq.i) then
                    qvois(jvois(j))=.true.
                    if (qfread) kvois(jvois(j))=fvois(j)
                endif
                if (jvois(j).eq.i) then
                    qvois(ivois(j))=.true.
                    if (qfread) kvois(ivois(j))=fvois(j)
                endif
             enddo
         endif
 
c        On calcule trois lignes de la matrice a la fois:
c        -----------------------------------------------
         do j=1,3*natom
            der2(1,j)=0.d0
            der2(2,j)=0.d0
            der2(3,j)=0.d0
         enddo
 
         do j=1,natom
            if (.not.qlist.or.(qlist.and.qvois(j))) then
            if (i.ne.j) then
            jj=3*j-2
            kij=knonb
            if (qfread) kij=kvois(j) 

            rx=xat(i)-xat(j)
            ry=yat(i)-yat(j)
            rz=zat(i)-zat(j)
            dist2=rx*rx + ry*ry + rz*rz
            dist=dsqrt(dist2)
 
            if (dist.lt.rsmall) then
                write(6,'(/2A,1PG10.4,A/2(I6,2A,I6,A,1X,2A))') 
     .          progrer,' Too small distance = ',dist,
     .        ' between following atoms.',
     .          i,': ',resnam(i),iresat(i),ssunam(i),atonam(i),' and ',
     .          j,': ',resnam(j),iresat(j),ssunam(j),atonam(j)
                stop '*Wrong coordinates*'
            endif

            if (rkh.gt.0.d0) kij=kij*exp(-(dist/rkh)**2.d0)

c           Constantes de force topologiques:
 
            if (nbonds.gt.0) then
            if (nbond(i+1).gt.nbond(i)) then
                do k=nbond(i),nbond(i+1)-1
                   if (jbond(k).eq.j) then
                       kij=kbond
                       goto 300
                   endif
                enddo
            endif
            else
            goto 300
            endif
 
            if (nangles.gt.0) then
            if (nangle(i+1).gt.nangle(i)) then
                do k=nangle(i),nangle(i+1)-1
                   if (jangle(k).eq.j) then
                       kij=kangle
                       goto 300
                   endif
                enddo
            endif
            else
            goto 300
            endif

            if (ndihs.gt.0) then
            if (ndihe(i+1).gt.ndihe(i)) then
                do k=ndihe(i),ndihe(i+1)-1
                   if (jdihe(k).eq.j) then
                       kij=kdihe
                       goto 300
                   endif
                enddo
            endif
            else
            goto 300
            endif

 300        continue
 
c           Calcul des elements: (potentiel harmonique)
c           -------------------------------------------
            if (dist.le.cutoff.or.qlist.or.
     .         (cutoff.le.0.d0.and.rkh.gt.0.d0)) then

                ll=ll+1
                nvoisat(i)=nvoisat(i)+1
                if (j.gt.i) then
                   if (unvmd.gt.0) then
                   write(unvmd,'(A,3F12.4,A,3F12.4,A)') 'draw line {',
     .             xat(i),yat(i),zat(i),'} {',xat(j),yat(j),zat(j),'}'
                   endif
                   if (unmol.gt.0) then
                   write(unmol,'(A,I6,3A)')  
     .           ' bonds require in residue ',iresat(i),', atom ',
     .             atonam(i),' and in molecule mol1 ',
     .           ' require in residue ',iresat(j),', atom ',
     .             atonam(j),' and in molecule mol1 ; '
                   endif
c                  Potentiel harmonique: 1/eval*knonb*(d - rval)**eval
                endif

                if (ll.eq.1.or.dist.lt.dmin) dmin=dist
                if (ll.eq.1.or.dist.gt.dmax) dmax=dist
                dmoy=dmoy+dist
                drms=drms+dist2

c               Elements diagonaux des blocs i et j:
c               -----------------------------------
                ddf=kij/dist2
                elemnt=rx*rx*ddf
                der2(1,ii)=der2(1,ii)+elemnt
                der2(1,jj)=der2(1,jj)-elemnt
                elemnt=ry*ry*ddf
                der2(2,ii+1)=der2(2,ii+1)+elemnt
                der2(2,jj+1)=der2(2,jj+1)-elemnt
                elemnt=rz*rz*ddf
                der2(3,ii+2)=der2(3,ii+2)+elemnt
                der2(3,jj+2)=der2(3,jj+2)-elemnt
 
c               Elements extra-diagonaux des deux blocs:
c               ---------------------------------------
                elemnt=rx*ry*ddf
                der2(1,ii+1)=der2(1,ii+1)+elemnt
                der2(2,ii)=der2(2,ii)+elemnt
                der2(1,jj+1)=der2(1,jj+1)-elemnt
                der2(2,jj)=der2(2,jj)-elemnt
                elemnt=rx*rz*ddf
                der2(1,ii+2)=der2(1,ii+2)+elemnt
                der2(3,ii)=der2(3,ii)+elemnt
                der2(1,jj+2)=der2(1,jj+2)-elemnt
                der2(3,jj)=der2(3,jj)-elemnt
                elemnt=ry*rz*ddf
                der2(2,ii+2)=der2(2,ii+2)+elemnt
                der2(3,ii+1)=der2(3,ii+1)+elemnt
                der2(2,jj+2)=der2(2,jj+2)-elemnt
                der2(3,jj+1)=der2(3,jj+1)-elemnt
            endif
            endif
            endif
         enddo
 
c        Sortie de la matrice-bande calculee:
c        -----------------------------------
c       (Uniquement la demi-matrice superieure)
 
c        Level-shift, pour eviter les zeros numeriques,
c        lors de la diagonalisation a venir 
c       (la minimisation est parfaite, par definition).
c        Le hasard est la pour lever la degenerescence
c        des six valeurs propres nulles, et differencier
c        rotations et translations.
 
         der2(1,ii)  =der2(1,ii)   + levelshft*random(iseed)
         der2(2,ii+1)=der2(2,ii+1) + levelshft*random(iseed)
         der2(3,ii+2)=der2(3,ii+2) + levelshft*random(iseed)
 
         do j=ii,3*natom
            jat=(j-1)/3+1
            if (der2(1,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout) 
     .          ii,j,der2(1,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I10,1PG20.12)') 
     .          ii,j,der2(1,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(1,j)).gt.rbig)  nbig=nbig+1
                if (dabs(der2(1,j)).gt.elmax) elmax=dabs(der2(1,j))
            endif
         enddo
         do j=ii+1,3*natom
            jat=(j-1)/3+1
            if (der2(2,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout) 
     .          ii+1,j,der2(2,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I10,1PG20.12)') 
     .          ii+1,j,der2(2,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(2,j)).gt.rbig) nbig=nbig+1
                if (dabs(der2(2,j)).gt.elmax) elmax=dabs(der2(2,j))
            endif
         enddo
         do j=ii+2,3*natom
            jat=(j-1)/3+1
            if (der2(3,j).ne.0.d0) then
                nnzero=nnzero+1
                if (qbinary) then
                write(unout) 
     .          ii+2,j,der2(3,j)/dsqrt(massat(i)*massat(jat))
                else
                write(unout,'(2I10,1PG20.12)') 
     .          ii+2,j,der2(3,j)/dsqrt(massat(i)*massat(jat))
                endif
                if (dabs(der2(3,j)).gt.rbig) nbig=nbig+1
                if (dabs(der2(3,j)).gt.elmax) elmax=dabs(der2(3,j))
            endif
         enddo
         elemnt=(der2(1,ii)+der2(2,ii+1)+der2(3,ii+2))/massat(i)
         if (elemnt.eq.0.d0) then
             write(6,'(2A,I6,A)') progrwn,
     .     ' Atom ',i,' has a null second derivatives...'
         else
             nntr=nntr+1
         endif
         trace=trace+elemnt
      enddo
      close(unout)

      if (unvmd.gt.0) then
          write(unvmd,'(2A)') 'mol load pdb ',nompdb(1:lnompdb)
          close(unvmd)
      endif
      if (unmol.gt.0) close(unmol)
      if (unrsd.gt.0) close(unrsd)
 
      nmoy=0.d0
      nrms=0.d0
      nmin=natom
      nmax=0
      do i=1,natom
         if (nvoisat(i).gt.nmax) nmax=nvoisat(i)
         if (nvoisat(i).lt.nmin) nmin=nvoisat(i)
         nmoy=nmoy+nvoisat(i)
         nrms=nrms+nvoisat(i)**2.d0
      enddo
      nmoy=nmoy/float(natom)
      nrms=nrms/float(natom)-nmoy**2.d0
      if (nrms.gt.0.d0) nrms=dsqrt(nrms)
 
      if (ll.eq.0) then
          write(6,'(/2A,I12,A)') progrer,
     .  ' No atom-atom interaction found. Too short cutoff ?'
          stop '*Empty matrix*'
      else
          dmoy=dmoy/float(ll)
          drms=drms/float(ll)-dmoy**2.d0
          if (drms.gt.0.d0) drms=dsqrt(drms)
      endif

      if (prtlev.gt.0)
     .write(6,'(/2A)') program,' Matrix statistics:'

      write(6,'(/2A,F8.4,A)') program,' The matrix is ',
     .  100.d0*dfloat(nnzero)/dfloat(3*natom*(3*natom+1)/2),' % Filled.'
      write(6,'(A,I12,A)') program,nnzero,'  non-zero elements.'

      if (prtlev.gt.0) then
      write(6,'(A,I12,A)') program,ll/2,' atom-atom interactions.'
      write(6,'(/2A,F9.2,A,F9.2/(A,I6))') program,
     .        ' Number per atom= ',nmoy,' +/- ',nrms,
     .'         Maximum number = ',nmax,
     .'         Minimum number = ',nmin
      write(6,'(/2A,F9.2,A,F9.2/(A,F9.2))') program,
     .        ' Average dist.  = ',dmoy,' +/- ',drms,
     .'         Maximum dist.  = ',dmax,
     .'         Minimum dist.  = ',dmin
      endif
 
      write(6,'(/2A,1PG12.6)') program,' Matrix trace   = ',trace

      if (prtlev.gt.0) then
      write(6,'(2A,1PG12.6)') program,' Larger element = ',elmax
      write(6,'(A,I6,A,1PG8.1)') program,
     .      nbig,' elements larger than +/- ',rbig
      endif
 
      write(6,'(/2A)') program,' Hessian matrix ready.'
      write(6,'(2A)') program,
     .' To diagonalize it (and get the modes),'//
     .' you may use diagstd, blzpack, diagrtb...'

      write(6,'(/2A)') program,' Normal end.'
      stop
      end
c-----------------------------------------------------------------------
      subroutine stringcl(chaine,nonblancs)
c
c     Les caracteres "blancs" de la CHAINE sont retires 
c    (a gauche et au milieu).
c     L'entier NONBLANCS donne la position du dernier caractere.
c
c     YHS-Jan-95, Toulouse.
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
c-----------------------------------------------------------------------
      subroutine string_split(chaine,taille,delimiteur,
     .                        souschaine,nbremax,nbre)
c
c     "Chaine" est coupee en "nbre" "souschaine" de part et d'autre du
c     "delimiteur"
c      YHS-Sep-93, Uppsala
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
c-----------------------------------------------------------------------
      subroutine openam(namfil,cformat,cstatus,unit,qverbos,
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
      integer lnom
      character*132 ordrunix
c begin:
      if (cstatus.eq.'old') cstatus='OLD'
      if (cstatus.eq.'new') cstatus='NEW'
      if (cstatus.eq.'ove') cstatus='OVE'
      if (cstatus.eq.'unknown') cstatus='UNKNOWN'
c
      qinterr=.false.
      qexist=.false.
c
      if (namfil.eq.' ') then 
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> No filename.'
          return
      endif
c
      if (namfil.eq.'stop'.or.namfil.eq.'end'                         
     &    .or.namfil.eq.'fin'.or.namfil.eq.'quit') then 
         write(6,'(2A)') 'Openam> Program is stopping on user request.'
         stop                                                                   
      endif 
 
c     Checks if filename is consistent with the opening:

      call stringcl(namfil,lnom)

      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          if (qverbos) then
              write(6,'(/2A)') '%Openam-Err> File: ',namfil(1:lnom)
              write(6,'(A)') 
     .      ' Expected in the current directory, but not found.'
          endif
          return
      endif

      if (qexist.and.cstatus.eq.'NEW') then
         write(6,'(/2A)') 
     .      '%Openam-Err> This file exists:',namfil(1:lnom)
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'

      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)                

      if (qverbos) then
         write(6,'(/2A)') ' Openam> File opened: ',namfil(1:lnom)
      endif

      return                                                                       
      end
c-----------------------------------------------------------------------
      subroutine rdatompdb(unpdb,ssusel,xat,yat,zat,binfo,
     .           atonam,iresat,resnam,ssunam,segid,natmax,natom,
     .           fatres,nresmx,nres,qerror,prtlev)
 
c     Lecture ligne a ligne d'un fichier pdb.
c     Uniquement les lignes commencant par 'ATOM'.
c     Uniquement ceux de la sous-unite selectionnee.
 
c     fatres(i): numero du premier atome du residu i.
 
c     YHS-nov-1996: version 1.0 (Toulouse).
c     YHS-sep-2004: version 1.1 (Lyon).
 
      implicit none
cI/O:
      integer unpdb, natmax, iresat(*), natom,  lnom,
     .        nresmx, nres, fatres(*), prtlev
      double precision xat(*), yat(*), zat(*), binfo(*)
      logical qerror
      character*4 atonam(*), resnam(*), segid(*)
      character*1 ssusel, ssunam(*)
cLocal:
      integer iatom, irs, irsprev, nerr, ntit,
     .        i, j, k, ii
      double precision bfact, x, y, z
      character*1  ssu
      character*4  ren, segat
      character*5  atncur
      character*80 lign80
cBegin:
      if (prtlev.gt.0)
     .write(6,'(/A)') ' Rdatompdb> Reading pdb file.'

      qerror=.false.
      nerr=0
 
      irsprev=-1
      ntit=0
      nres=0
      iatom=0
 105  continue   
      read(unpdb,'(A)',end=200,err=110) lign80 
  
      goto 120                                
 110  continue
      nerr=nerr+1                            
 
 120  continue                              
c     if (lign80(1:4).eq.'ATOM') then
      if (lign80(1:4).eq.'ATOM'.or.lign80(1:6).eq.'HETATM') then
      read(lign80,'(12X,A4,1X,A4,A1,I4,4X,3F8.3,6X,F6.2,6X,A4)',
     .            end=130,err=130) 
     .            atncur, ren, ssu, irs, x, y, z, 
     .            bfact, segat 
 130  continue
      if (iatom.lt.natmax) then
          if (ssu.eq.ssusel.or.ssusel.eq.' ') then
          iatom=iatom+1
          xat(iatom)=x
          yat(iatom)=y
          zat(iatom)=z
          binfo(iatom)=bfact
 
          call stringcl(atncur,lnom)
          atonam(iatom)=atncur
          call stringcl(ren,lnom)
          resnam(iatom)=ren
          iresat(iatom)=irs
          ssunam(iatom)=ssu
          segid(iatom)=segat
 
          if (irs.ne.irsprev) then
              nres=nres+1
              if (nres.gt.nresmx) then
                  write(6,'(A/A,I6)') 
     .          '%Rdatompdb-Er> Too many residues in this file.',
     .          ' Maximum allowed is = ',nresmx
                  stop
              endif
              irsprev=irs
              fatres(nres)=iatom
          endif
          endif
      else
          write(6,'(A/A,I6)') 
     .      '%Rdatompdb-Er> Too many atoms in this file.',
     .      ' Maximum allowed is = ',natmax
          stop
      endif
      else if (lign80(1:6).eq.'REMARK'.and.prtlev.gt.0) then
          ntit=ntit+1
          if (ntit.le.10) then
              write(6,'(A)') lign80
          else if (ntit.eq.11) then
              write(6,'(A)') ' .... '
          endif
      endif
 
c     2) Ligne suivante du fichier pdb :
 
      goto 105
 
c     3) Fin de la lecture du fichier pdb :
 
 200  continue 
      if (prtlev.gt.1) then
      write(6,*) 'Rdatompdb> End of file reached.'
      write(6,*) 'Rdatompdb> Number of I/O errors: ',
     .            nerr
      endif
 
      natom=iatom
      fatres(nres+1)=natom+1
      irs=0
      if (natom.gt.0) irs=iresat(natom)
 
      write(6,'(/(A,I6))') 
     .' Rdatompdb> Number of residues found = ',nres,
     .'            First residue number     = ',iresat(1),
     .'            Last  residue number     = ',irs,
     .'            Number of atoms found    = ',natom
      if (prtlev.gt.0)
     .write(6,'(A,F8.1)') 
     .'            Mean number per residue  = ',float(natom)/float(nres)
 
      if (natom.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Er> No atom found in file.'
          qerror=.true.
      endif
      if (nres.eq.0) then
          write(6,'(A)')
     .  '%Rdatompdb-Er> No residue found in file.'
          qerror=.true.
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine mintomaj(chaine)
 
c     Les caracteres minuscules sont mis en MAJUSCULES.
c     Les autres ne sont pas touches.
 
c     YHS-Oct-98: Premiere version (Toulouse).
c     YHS-Sep-03: Dernieres modifications (Lyon).
 
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
c-----------------------------------------------------------------------
      subroutine readxyz(uninp,x,y,z,w,ic,nmax,ncoor,ndat,qerror,prtlev)

c     Reads at most NMAX coordinates in free format. 

c     Either:
c     x, y, z
c     or:
c     x, y, z, w
c     or:
c     x, y, z, w, ic

c     If first word in ligne is not a number, the whole ligne is
c     assumed to be a title or a commentary.

c     YHS-Sep-03: First version (Lyon).

cI/O:
      logical qerror
      integer ic(*), ncoor, ndat, nmax, prtlev, uninp
      double precision w(*), x(*), xc, y(*), yc, z(*), zc
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
     .  ' Comments, or lignes with less than three numbers: '

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
              x(ncoor)=xc
              y(ncoor)=yc
              z(ncoor)=zc
              if (nfound.eq.4) then
                  w(ncoor)=wc
              else
                  w(ncoor)=wc
                  ic(ncoor)=nlu 
              endif
          else
              write(6,'(/2A,I9,A)') progrer,' More than ',
     .        nmax,' particles in file.'
              write(6,'(2A)') progrer,
     .      ' Please increase program memory limits (Sorry for that).'
              stop
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
     .    ' numbers found in ',stats(i),' lignes.'  
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine vecstat(vect,nmax,rmin,rmax,rave,rdev)

c     Statistics for vector Vect(NMAX):
c     minimum, maximum, average (rave) and standard deviation (rdev).

c     YHS-Sep-03: First version (Lyon).

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
          write(6,'(2A)') progrer,' Zero-length vector.'
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
C-----------------------------------------------------------------------
      REAL*8 FUNCTION RANDOM(ISEED)
C-----------------------------------------------------------------------
C     RANDOM NUMBER GENERATOR: UNIFORM DISTRIBUTION (0,1)
C     ISEED: SEED FOR GENERATOR. ON THE FIRST CALL THIS HAS TO
C     HAVE A VALUE IN THE EXCLUSIVE RANGE (1, 2147483647)
C     AND WILL BE REPLACED BY A NEW VALUE TO BE USED IN
C     FOLLOWING CALL.
C
C     REF: Lewis, P.A.W., Goodman, A.S. & Miller, J.M. (1969)
C     "Pseudo-random number generator for the System/360", IBM
C     Systems Journal 8, 136.
C
C     This is a "high-quality" machine independent generator.
C     INTEGERS are supposed to be 32 bits or more.
C     The same algorithm is used as the basic IMSL generator.
C
C     Author: Lennart Nilsson
C
      implicit none
      INTEGER ISEED
      REAL*8 DSEED,DIVIS,DENOM,MULTIP
      DATA  DIVIS/2147483647.D0/
      DATA  DENOM /2147483711.D0/
      DATA  MULTIP/16807.D0/
C
      IF(ISEED.LE.1) ISEED=314159
      DSEED=MULTIP*ISEED
      DSEED=MOD(DSEED,DIVIS)
      RANDOM=DSEED/DENOM
      ISEED=DSEED
C
      RETURN
      END
