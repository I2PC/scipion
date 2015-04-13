      program Chkmod
c
c     Pour un ensemble de vecteurs CERFACS.
c     Collectivite des modes=f(frequence)
c
c     YHS-Mai-2001: Premiere version, from rms_mode, v1.16.
c
      implicit none
      integer  natmax, nvecmx
      parameter(natmax=55000,nvecmx=200)
      integer i, ii, ivec, j, jj, k,
     .        lnomeig, natom, nddl,
     .        numvec(natmax), nunit, nvec,
     .        uneig, unout 
      double precision coll, freq(3*natmax), 
     .       matvec(3*natmax,nvecmx),
     .       norme, normtot, rsmall
      logical qinter, qexist, qok
      character*8  program
      character*10 cformat, cstatus
      character*11 progrwn, progrer
      character*40 version
      character*64 namfil, nomeig
c
      version=' Version 1.00, Bordeaux.'
      program=' Chkmod>'
      progrwn='%Chkmod-Wn>'
      progrer='%Chkmod-Er>'
c
      rsmall=1e-4
c
      write(6,'(2A)') program,version
c
      call getnam('Eigenvector filename ?',nomeig,lnomeig,qok)
      if (.not.qok) stop
c
      nunit=10
      uneig=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="old"
      call openam(nomeig,cformat,cstatus,uneig,.true.,
     .            qinter,qexist)
      if (qinter.or..not.qexist) stop
c
c     Lecture:
c     -------
      call rdmodfacs(uneig,3*natmax,nvecmx,numvec,freq,
     .               matvec,nddl,nvec)
      natom=nddl/3
c
      write(6,'(/A,I5,A,I6,A)') program,
     .      nvec,' vectors, ',nddl,' coordinates in file.'
      write(6,'(A,5X,A,I6,A)') program,
     .    ' That is: ',natom,' cartesian points.'
c
      if (nvec.le.0.or.natom.le.0) then
          write(6,'(2A)') progrer,' Wrong vector file.'
          stop
      endif
c
c     Sortie:
c     ------
      namfil='Chkmod.res'
      unout=nunit
      nunit=nunit+1
      cformat="FORMATTED"
      cstatus="ove"
      call openam(namfil,cformat,cstatus,unout,.true.,
     .            qinter,qexist)
c
      write(6,'(2A)') program,
     .' Collectivity=f(frequency) to be written in this file.'
c
      do ivec=1,nvec
c
      normtot=0.d0
      do i=1,natom
         ii=3*i-2
         norme=matvec(ii,ivec)**2.d0+
     .         matvec(ii+1,ivec)**2.d0+
     .         matvec(ii+2,ivec)**2.d0
         normtot=normtot+norme
      enddo
c
      if (normtot.le.rsmall) then
         write(6,'(2A,I5,A)') progrer,
     . ' Vector ',ivec,' has a null norm.'
         stop
      endif 
      if (dabs(normtot-1.d0).gt.rsmall) then
         write(6,'(2A,I5,A,F8.4,A)') progrwn,
     . ' Norm of vector ',ivec,' is: ',normtot,
     .' (instead of 1.0000). '
      endif 
c     
      coll=0.d0
      do i=1,natom
         ii=3*i-2
         norme=matvec(ii,ivec)**2.d0+
     .         matvec(ii+1,ivec)**2.d0+
     .         matvec(ii+2,ivec)**2.d0
         norme=norme/normtot
         coll=coll+norme*log(norme)
      enddo
      coll=exp(-coll)/dfloat(natom)
c
      write(unout,'(2F12.4)') freq(ivec), coll
c
c     Vecteur suivant: 
      enddo
c
      write(6,'(/2A)') program,' Normal end.'
c
      stop
      end
c
c     Questions/reponses.
c     Serie de routines utilitaires:
c     getnam
c     Version 1.3, Bordeaux.
c
      subroutine getnam(message,nomlu,lnomlu,qok)
c
c     NOMLU obtenu en reponse au MESSAGE.
c     NTRYMX essais en cas de probleme.
c     YHS-oct-96
c
      implicit none
cI/O:
      integer lnomlu
      logical qok
      character*(*) message, nomlu
cLocal:
      integer ntry, ntrymx
cBegin:
      ntrymx=5
c
      qok=.false.
      ntry=0
c
 100  continue
      ntry=ntry+1
      if (ntry.ge.ntrymx) return
c
      write(6,'(A,A)') ' Getnam> ',message
      read(5,'(A)',end=200,err=100) nomlu
c
      call stringcl(nomlu,lnomlu)
      write(6,'(A,A)') ' Getnam> ',nomlu(1:lnomlu)
c
      qok=.true.
      return
 200  continue
      return
      end 
c
      subroutine rdmodfacs(uneig,nddlmax,nvecmx,numvec,freq,
     .           matvec,nddl,nvec)
c
c     Lecture de modes CERFACS.
c     Devra remplacer rdcerfacs.
c     Difference: comptage de l'ordre de la matrice.
c    (et pas des atomes)
c
c     Premieres versions (rdcerfacs):
c     YHS-Nov-1996.
c     Dernieres modifications:
c     YHS-Jan-2001.
cI/O:
      integer numvec(*), nvecmx, nddlmax, nvec, nddl, uneig
      double precision freq(*), matvec(nddlmax,*)
cLocal:
      integer nmotsmax
      parameter(nmotsmax=100)
      integer nerr, ivec, indnm_cfacs, nmots,
     .        i, ii, j, jj, k, kk
      double precision wtofreq
      logical qfound, qold, qfirst
      character*1 carnum
      character*132 lign132, mots(nmotsmax)
cDefaut:
c     Facteur de conversion (2*pi*f)**2 -> f (cm-1):
      wtofreq=108.586698
c
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
c
      qfound=qfound.or.
     .      (index(lign132,' value ').gt.0.and.
     .       index(lign132,' vector ').gt.0.and.
     .       index(lign132,' residual ').le.0)
      qold=qold.or.
     .      (index(lign132,' VALUE ').gt.0.and.
     .       index(lign132,' VECTOR ').gt.0)
c
      if (.not.qfound.and..not.qold) goto 100
c________________________________________
c
c     Lecture des frequences des modes :
c________________________________________
c
      if (qfirst) then
          if (qold) then
          write(6,'(/A)') 
     .  ' Rdmodfacs> Old Blzpack file format detected.'
          else
          write(6,'(/A)') 
     .  ' Rdmodfacs> Blzpack file format detected.'
          endif
          qfirst=.false.
      endif
c
      ivec=0
      nvec=0
 250  continue
      ivec=ivec+1
      if (ivec.gt.nvecmx) then
          write(6,'(/A,I5,A)') 
     .  '%Rdmodfacs-Err> More than ',nvecmx,' vectors in file.'
          return
      endif
c
      read(lign132,'(7X,I5,12X,G12.4)',end=240,err=240)
     .     numvec(ivec), freq(ivec)
      freq(ivec)=wtofreq*dsqrt(abs(freq(ivec)))
c
      goto 255
 240  continue
      write(6,'(/3A)')
     .    '%Rdmodfacs-W> Pb with ligne: ',lign132(1:36),'...'
 255  continue
c
      nvec=ivec
      write(6,'(/A,I6)')
     .    ' Rdmodfacs> Numero du vecteur CERFACS en lecture:',
     .      numvec(ivec)
      write(6,'(A,1PG12.4)')
     .    ' Rdmodfacs> Frequence du vecteur en lecture:',
     .      freq(ivec)
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
c     2) Lecture des coordonnees des modes CERFACS :
c        Format libre.
c____________________________________________________
c
      k=0
 257  continue
      if (k.ge.nddlmax) then
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
     .  ' Rdmodfacs: Lecture des modes CERFACS terminee.',
     .  ' Item VALUE non trouve dans la ligne:',lign132
          goto 300
      endif
c
 270  continue
      write(6,'(A,I6,A)')
     .   ' %Rdmodfacs-Error: durant la lecture de la coordonnee ',
     .      i,' du mode.'
      stop
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
c     Lecture de noms de fichier + Ouverture.
c     Serie de routines utilitaires:
c     openam
c     stringcl
c     string_split
c     Version 1.11, Bordeaux.
c
      SUBROUTINE openam(namfil,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c
c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT.
c
c     input:
c        namfil: nom du fichier a ouvrir. 
c        "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.
c     output: 
c        qexist: flag / existence du fichier 
c        qinterr: Pas de nom pour le fichier cherche.
c
c     YHS-oct-1993: Premiere version.
c     YHS-jan-2000: Derniere modification.
c I/O:
      logical qinterr, qverbos, qexist
      integer unit
      character*10 cformat, cstatus
      character*64 namfil
c Local
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
         write(*,*) 'Openam> Program is stopping on user request.'
         stop                                                                   
      endif 
c
c     Checks if filename is consistent with the opening:
c
      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> File not found.'
          return
      endif
c
      if (qexist.and.cstatus.eq.'NEW') then
         write(*,'(/A)') 
     .      '%Openam-Err> This file exists:',namfil
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
c                                                                   
      if (qverbos) then
         write(*,'(/A,I6,A)')
     .           ' Openam> file on opening on unit ',unit,':'
         write(*,*) namfil
      endif
      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)                
c        
      return                                                                       
      end
c
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
c
      subroutine stringcl(chaine,nonblancs)
c
c     Les caracteres "blancs" de la CHAINE sont retires (a gauche et au milieu).
c     L'entier NONBLANCS donne la position du dernier caractere.
c
c     YHS-Jan-95, Toulouse.
c     YHS-Oct-00, Bordeaux.
c I/O:
      integer nonblancs
      character*(*) chaine
c Local:
      integer icar, ncar, taille
c Begin:
      nonblancs=0
      taille=len(chaine)
      if (taille.le.0) return
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
      icar=taille
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
