  program mkinfo

    use rddatamod

    implicit none

    integer*8                            :: ncoo,nzero,nqc,nsta
    integer*8, parameter                 :: maxdim=60,maxset=7000
    integer*8, dimension(maxdim)         :: iimag
    integer*8, dimension(maxset,maxdim)  :: ignore
    real*8                               :: e0
    real*8, dimension(maxdim)            :: xcoo0,mass,peig,freq
    real*8, dimension(maxdim,maxdim)     :: pvec,invpvec
    real*8, dimension(maxset,maxdim)     :: qcoo,xcoo,energy
    character(len=80)                    :: ain,adat,aset,acalc,&
                                            aprog,ainfo,asymfile
    character(len=80), dimension(maxset) :: aqc
    character(len=16), dimension(maxdim) :: asym
 
!------------------------------------------------------------------
! Read input file name
!------------------------------------------------------------------
    ain=''
    call getarg(1,ain)
    if (ain.eq.'') then
       write(6,'(/,a,/)') 'Filename not specified!'
       STOP
    endif

!------------------------------------------------------------------
! Read input file
!------------------------------------------------------------------
    call rdinp(ain,adat,aset,e0,acalc,aprog,nsta,asymfile)

!------------------------------------------------------------------
! Read the data file
!------------------------------------------------------------------
    call rddat(adat,ncoo,pvec,peig,iimag,xcoo0,maxdim,mass)

!------------------------------------------------------------------
! Discard any unwanted modes and frequency scale the remaining 
! modes
!------------------------------------------------------------------
    call scaleq(peig,pvec,invpvec,nzero,ncoo,maxdim,freq,iimag,&
         mass)

!------------------------------------------------------------------
! Read the set file
!------------------------------------------------------------------
    call rdset(aset,aqc,maxset,nqc,ignore,maxdim)

!------------------------------------------------------------------
! Read the mode symmetries
!------------------------------------------------------------------
    call rdsym(asymfile,ncoo,nzero,asym,maxdim)

!------------------------------------------------------------------
! Read the QC files
!------------------------------------------------------------------
    call rdqc(aqc,maxset,nqc,acalc,aprog,xcoo,maxdim,ncoo,nsta,&
         energy,e0,qcoo,invpvec,nzero,xcoo0)

!------------------------------------------------------------------
! Write the info file
!------------------------------------------------------------------
    call wrinfo(ain,ainfo,ncoo,nzero,nsta,e0,xcoo0,maxdim,mass,&
         asym,freq,pvec,nqc,energy,maxset,qcoo,xcoo,ignore)

    STOP

  contains
      
!##################################################################

    subroutine rdinp(ain,adat,aset,e0,acalc,aprog,nsta,asymfile)

      use iomod
      
      implicit none
      
      integer*8         :: i,nsta
      real*8            :: e0
      character(len=80) :: ain,adat,aset,acalc,aprog,asymfile

      integer*8                                   :: iin,ic,iz,&
                                                     ierr
      integer*8, parameter                        :: maxkey=48,&
                                                     maxkeylen=200
      integer*8, dimension(maxkey)                :: lc      
      character(len=maxkeylen), dimension(maxkey) :: keyword,&
           keyorig
      character(len=240)                          :: inptit

!------------------------------------------------------------------
! Set defaults
!------------------------------------------------------------------
      adat=''
      aset=''
      acalc=''
      aprog=''
      asymfile=''
      nsta=0
      e0=0.0d0

!------------------------------------------------------------------
! Open the input file
!------------------------------------------------------------------
      iin=20
      open(iin,file=ain,form='formatted',status='old')

!------------------------------------------------------------------
! Read keywords
!------------------------------------------------------------------
10    continue
      call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,maxkey)

      i=1
15    continue
      if (keyword(i) .eq. 'end-input') then
         goto 20

      else if (keyword(i).eq.'datafile') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') adat
         else
            write(6,'(/,a,/)') &
                 'No argument given with the datafile keyword'
            STOP
         endif

      else if (keyword(i).eq.'files') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') aset
         else
            write(6,'(/,a,/)') &
                 'No argument given with the files keyword'
            STOP
         endif

      else if (keyword(i).eq.'e0') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),*) e0
         else
            write(6,'(/,a,/)') &
                 'No argument given with the e0 keyword'
            STOP
         endif
         
      else if (keyword(i).eq.'calc_type') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') acalc
         else
            write(6,'(/,a,/)') &
                 'No argument given with the calc_type keyword'
            STOP
         endif

      else if (keyword(i).eq.'qc_prog') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') aprog
         else
            write(6,'(/,a,/)') &
                 'No argument given with the qc_prog keyword'
            STOP
         endif

      else if (keyword(i).eq.'nstates') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),*) nsta
         else
            write(6,'(/,a,/)') &
                 'No argument given with the nstates keyword'
            STOP
         endif

      else if (keyword(i).eq.'symmetries') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') asymfile
         else
            write(6,'(/,a,/)') &
                 'No argument given with the symmetries keyword'
            STOP
         endif         

      else
         write(6,'(/,2(a,2x),/)') 'Unknown keyword:',keyword(i)
         STOP
      endif

!------------------------------------------------------------------
! Check if more keywords are to be read, else read new line
!------------------------------------------------------------------
      if (i.lt.ic) then
         i=i+1
         go to 15
      else
         go to 10
      endif
      
20    continue
      close(iin)

!------------------------------------------------------------------
! Check that all the required information has been given
!------------------------------------------------------------------
      if (adat.eq.'') then
         write(6,'(/,a,/)') 'Datafile not specified!'
         STOP
      endif
      
      if (aset.eq.'') then
         write(6,'(/,a,/)') 'Set file not specified!'
         STOP
      endif

      if (e0.eq.0) then
         write(6,'(/,a,/)') 'e0 not specified!'
         STOP
      endif

      if (acalc.eq.'') then
         write(6,'(/,a,/)') 'Calculation type not specified!'
         STOP
      endif

      if (aprog.eq.'') then
         write(6,'(/,a,/)') 'QC program type not specified!'
         STOP
      endif

      if (nsta.eq.0) then
         write(6,'(/,a,/)') &
              'Number of states program type not specified!'
         STOP
      endif

      if (asymfile.eq.'') then
         write(6,'(/,a,/)') 'Symmetry file not specified!'
         STOP
      endif

      return
    end subroutine rdinp

!##################################################################

    subroutine scaleq(peig,pvec,invpvec,nzero,ncoo,maxdim,freq,&
         iimag,mass)

      implicit none

      integer*8                        :: nzero,maxdim,ncoo,i,j,k
      integer*8, dimension(ncoo)       :: itmp1
      integer*8, dimension(maxdim)     :: iimag
      real*8, dimension(maxdim,maxdim) :: pvec,invpvec
      real*8, dimension(maxdim)        :: peig,freq,mass
      real*8, parameter                :: tol=1.0d0
      real*8                           :: ftmp
      real*8, dimension(ncoo)          :: ftmp1
      real*8, dimension(ncoo,ncoo)     :: ftmp2

!------------------------------------------------------------------
! Determine the no. of zero frequencies
!------------------------------------------------------------------
      nzero=0
      do i=1,ncoo
         ftmp=sqrt(peig(i))*5140.66d0
         if (ftmp.lt.tol) nzero=nzero+1
      enddo

!------------------------------------------------------------------
! Discard unwanted modes
!------------------------------------------------------------------
      ftmp1=0.0d0
      ftmp2=0.0d0
      itmp1=0

      k=0
      do i=nzero+1,ncoo
         k=k+1
         ftmp1(k)=peig(i)
         itmp1(k)=iimag(i)
         do j=1,ncoo
            ftmp2(j,k)=pvec(j,i)
         enddo
      enddo

      peig=0.0d0
      pvec=0.0d0
      iimag=0
      do i=1,ncoo-nzero
         peig(i)=ftmp1(i)
         freq(i)=sqrt(peig(i))
         iimag(i)=itmp1(i)
         do j=1,ncoo
            pvec(j,i)=ftmp2(j,i)
         enddo
      enddo

!------------------------------------------------------------------
! Frequency scale (in eV) and mass scale (in amu) the remaining 
! modes and calculate the inverse transformation matrix
!------------------------------------------------------------------
      invpvec=0.0d0
      do i=1,ncoo-nzero
         ftmp=freq(i)*0.6373641d0
         do j=1,ncoo
            invpvec(i,j)=pvec(j,i)*(15.4644*sqrt(ftmp*mass(j)))
            pvec(j,i)=pvec(j,i)/(15.4644*sqrt(ftmp*mass(j)))
         enddo
      enddo

      return
    end subroutine scaleq

!##################################################################

    subroutine rdset(aset,aqc,maxset,nqc,ignore,maxdim)

      use iomod

      implicit none

      integer*8                            :: maxset,nqc,i,maxdim,s
      integer*8, dimension(maxset,maxdim)  :: ignore
      character(len=80)                    :: aset,string
      character(len=80), dimension(maxset) :: aqc

      integer*8                                   :: iin,ic,iz,&
                                                     ierr
      integer*8, parameter                        :: maxkey=48,&
                                                     maxkeylen=200
      integer*8, dimension(maxkey)                :: lc      
      character(len=maxkeylen), dimension(maxkey) :: keyword,&
           keyorig
      character(len=240)                          :: inptit

!------------------------------------------------------------------
! Set defaults
!------------------------------------------------------------------
      ignore=0

!------------------------------------------------------------------
! Open the set file
!------------------------------------------------------------------
      iin=20
      open(iin,file=aset,form='formatted',status='old')

!------------------------------------------------------------------
! Read the set file
!------------------------------------------------------------------
      aqc=''
      nqc=0
10    continue

      call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,maxkey)

!      read(unit,'(a)') string

      if (keyword(1).ne.'end-files') then
         nqc=nqc+1
         i=1
         ! Read file name         
         if (i.eq.1) read(keyword(i),'(a)') aqc(nqc)
         ! Read any additional keywordds
20       continue
         if (i.lt.ic) then
            i=i+1
            if (keyword(i).eq.'ignore') then
               i=i+2
               read(keyword(i),*) s
               ignore(nqc,s)=1
25             continue
               if (keyword(i+1).eq.',') then
                  i=i+2
                  read(keyword(i),*) s
                  ignore(nqc,s)=1
                  goto 25
               endif
            endif
            goto 20
         endif
         goto 10
      endif

!------------------------------------------------------------------
! Close the set file
!------------------------------------------------------------------
      close(iin)

      return
    end subroutine rdset

!##################################################################

    subroutine rdqc(aqc,maxset,nqc,acalc,aprog,xcoo,maxdim,ncoo,&
         nsta,energy,e0,qcoo,invpvec,nzero,xcoo0)

      implicit none

      integer*8                            :: maxset,nqc,n,unit,&
                                              maxdim,ncoo,i,j,&
                                              nsta,nzero
      real*8, dimension(maxset,maxdim)     :: xcoo,qcoo,energy
      real*8, dimension(maxdim,maxdim)     :: invpvec
      real*8, dimension(maxdim)            :: xcoo0
      real*8                               :: e0
      character(len=80), dimension(maxset) :: aqc
      character(len=80)                    :: acalc,aprog

      real*8               :: r0,rcurr,tmass
      real*8, dimension(3) :: com1,com2

!------------------------------------------------------------------
! Loop over the QC files and read the energies and Cartesian
! coordinates
!------------------------------------------------------------------
      xcoo=0.0d0
      energy=0.0d0
      unit=20
      do n=1,nqc

         ! Open file
         open(unit,file=aqc(n),form='formatted',status='old')

         ! Read Cartesian coordinates
         call rdcoord(n,unit,acalc,aprog,xcoo,maxdim,maxset,ncoo,xcoo0)

         ! Calculate the coordinates in terms of the projected
         ! normal modes
         call x2prq(xcoo(n,1:ncoo),xcoo0(1:ncoo),qcoo(n,1:ncoo),&
              ncoo,nzero,invpvec(1:ncoo-nzero,1:ncoo))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! NMP
!
!         do i=1,ncoo
!            qcoo(n,i)=0.0d0
!         enddo
!
!         r0=1.4509340364524157d0
!
!         rcurr=0.0d0
!         do i=1,3
!            rcurr=rcurr+(xcoo(n,3-3+i)-xcoo(n,30-3+i))**2
!         enddo
!         rcurr=sqrt(rcurr)
!
!         qcoo(n,1)=(rcurr-r0)/0.529177249d0
!        
!
! Pyrrole dimer, T-shaped, R
         do i=1,ncoo
            qcoo(n,i)=0.0d0
         enddo

         r0=1.012488002483254

         rcurr=0.0d0
         do i=1,3
            rcurr=rcurr+(xcoo(n,27-3+i)-xcoo(n,30-3+i))**2
         enddo
         rcurr=sqrt(rcurr)

         qcoo(n,1)=(rcurr-r0)/0.529177249d0

!! Pyrrole dimer, stacked, R
!         do i=1,ncoo
!            qcoo(n,i)=0.0d0
!         enddo
!
!         r0=1.012488002483254
!
!         rcurr=0.0d0
!         do i=1,3
!            rcurr=rcurr+(xcoo(n,15-3+i)-xcoo(n,18-3+i))**2
!         enddo
!         rcurr=sqrt(rcurr)
!
!         qcoo(n,1)=(rcurr-r0)/0.529177249d0
!
!
! Pyrrole dimer, rcom
!
!
!         ! Reference
!         tmass=0.0d0
!         ! CoM 1
!         com1=0.0d0
!         do i=1,10
!            do j=1,3
!               com1(j)=com1(j)+mass(i*3)*xcoo0(i*3-3+j)
!            enddo
!            tmass=tmass+mass(i*3)
!         enddo         
!         ! CoM 2
!         com2=0.0d0
!         do i=11,20
!            do j=1,3
!               com2(j)=com2(j)+mass(i*3)*xcoo0(i*3-3+j)
!            enddo
!            tmass=tmass+mass(i*3)
!         enddo
!         com1=com1/tmass
!         com2=com2/tmass
!
!         ! r0
!         r0=0.0d0
!         do i=1,3
!            r0=r0+(com1(i)-com2(i))**2
!         enddo
!         r0=sqrt(r0)
!
!         ! Current
!         tmass=0.0d0
!         do i=1,ncoo
!            qcoo(n,i)=0.0d0
!         enddo
!
!         ! CoM 1
!         com1=0.0d0
!         do i=1,10
!            do j=1,3
!               com1(j)=com1(j)+mass(i*3)*xcoo(n,i*3-3+j)
!            enddo
!            tmass=tmass+mass(i*3)
!         enddo
!
!         
!         ! CoM 2
!         com2=0.0d0
!         do i=11,20
!            do j=1,3
!               com2(j)=com2(j)+mass(i*3)*xcoo(n,i*3-3+j)
!            enddo
!            tmass=tmass+mass(i*3)
!         enddo
!
!         com1=com1/tmass
!         com2=com2/tmass
!
!         ! rcom
!         rcurr=0.0d0
!         do i=1,3
!            rcurr=rcurr+(com1(i)-com2(i))**2
!         enddo
!         rcurr=sqrt(rcurr)
!
!
!         qcoo(n,1)=(rcurr-r0)/0.529177249d0
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! Read the state energies
         call rdenergy(n,unit,nsta,energy,maxdim,maxset,acalc,&
              aprog,e0)

         ! Close file
         close(unit)

      enddo

    return
    end subroutine rdqc

!##################################################################

    subroutine rdcoord(n,unit,acalc,aprog,xcoo,maxdim,maxset,ncoo,&
         xcoo0)

      use eckartmod

      implicit none
      
      integer*8                        :: unit,maxdim,maxset,ncoo,&
                                          n,i,j,k,l,sx,sy,sz
      real*8, dimension(maxset,maxdim) :: xcoo
      real*8, dimension(maxdim)        :: xcoo0
      real*8, dimension(ncoo,ncoo)     :: tmat
      real*8, dimension(ncoo)          :: ftmp,tmp,rotxcoo
      real*8                           :: high,sum,minus
      character(len=80)                :: string,acalc,aprog

!------------------------------------------------------------------
! Read the coordinates
!------------------------------------------------------------------
      ! Molpro
      if (aprog.eq.'molpro') then
         call coord_mol(n,unit,xcoo,maxdim,maxset,ncoo)

      else if (aprog.eq.'turbomole') then
         call coord_turbomole(n,unit,xcoo,maxdim,maxset,ncoo)
      else
         write(6,'(/,2(2x,a),/)') "File type not supported:",&
              aprog(1:len_trim(aprog))
         STOP
      endif

!------------------------------------------------------------------
! Rotate to the Eckart frame
!------------------------------------------------------------------
      ! Calculate the transformation matrix
      call ecktransmat(xcoo(n,1:ncoo),ncoo,mass(1:ncoo),tmat)

      ! Rotate the coordinates
      ftmp=0.0d0
      do i=1,ncoo
         do j=1,ncoo
            ftmp(i)=ftmp(i)+tmat(i,j)*xcoo(n,j)
         enddo
      enddo
      do i=1,ncoo
         xcoo(n,i)=ftmp(i)
      enddo

      ! Take the combination of +/-x, +/-y and +/-z-components that 
      ! yields the maximum coincidence with the reference geometry
      high=999999999.9
      do i=1,2
         do j=1,2
            do k=1,2
               sx=(-1)**i
               sy=(-1)**j
               sz=(-1)**k
               do l=1,ncoo/3
                  tmp(3*l-2)=sx*xcoo(n,3*l-2)
                  tmp(3*l-1)=sy*xcoo(n,3*l-1)
                  tmp(3*l)=sz*xcoo(n,3*l)
               enddo
               sum=0.0d0
               do l=1,ncoo
                  minus=tmp(l)-xcoo0(l)
                  minus=abs(minus)
                  sum=sum+minus
               enddo
               if (sum.lt.high) then
                  high=sum
                  do l=1,ncoo
                     rotxcoo(l)=tmp(l)
                  enddo
!                  bsx=sx
!                  bsy=sy
!                  bsz=sz
               endif
            enddo
         enddo
      enddo
      do i=1,ncoo
         xcoo(n,i)=rotxcoo(i)
      enddo

      return
    end subroutine rdcoord

!##################################################################

    subroutine coord_mol(n,unit,xcoo,maxdim,maxset,ncoo)

      implicit none
      
      integer*8                        :: n,unit,maxdim,maxset,&
                                          ncoo,i,ngeom
      real*8, dimension(maxset,maxdim) :: xcoo
      character(len=80)                :: string

!------------------------------------------------------------------
! Determine the number of geometries
!------------------------------------------------------------------
      rewind(unit)
      ngeom=0
5     continue
      read(unit,'(a)') string
      if (string(2:25).ne.'Variable memory released') then
         if (string(2:19).eq.'ATOMIC COORDINATES') ngeom=ngeom+1
         goto 5
      endif

!------------------------------------------------------------------
! Read to coordinates
!------------------------------------------------------------------
      rewind(unit)
      i=0
10    continue
      read(unit,'(a)') string
      if (string(2:19).eq.'ATOMIC COORDINATES') i=i+1
      if (i.lt.ngeom) goto 10

      do i=1,3
         read(unit,*)
      enddo

!------------------------------------------------------------------
! Read the coordinates and convert to Angstrom
!------------------------------------------------------------------
      do i=1,ncoo/3
         read(unit,'(18x,3(F15.9))') &
              xcoo(n,i*3-2),xcoo(n,i*3-1),xcoo(n,i*3)
      enddo

      do i=1,ncoo
         xcoo(n,i)=xcoo(n,i)*0.529177249d0
      enddo

      return
    end subroutine coord_mol

!##################################################################

    subroutine rdenergy(n,unit,nsta,energy,maxdim,maxset,acalc,&
         aprog,e0)

      implicit none

      integer*8                        :: n,unit,nsta,maxdim,maxset
      real*8, dimension(maxset,maxdim) :: energy
      real*8                           :: e0
      character(len=80)                :: acalc,aprog

!------------------------------------------------------------------
! Read the energies
!------------------------------------------------------------------
      ! Molpro
      if (aprog.eq.'molpro') then
         ! EOM
         if (acalc.eq.'eom') then
            call eom_mol(n,unit,nsta,e0,energy,maxdim,maxset)
         else if (acalc.eq.'cas') then
            call cas_mol(n,unit,nsta,e0,energy,maxdim,maxset)
         else
            write(6,'(/,2(2x,a),/)') &
                 "Calculation type not supported:",&
                 acalc(1:len_trim(acalc))
         endif

      ! Turbomole
      else if  (aprog.eq.'turbomole') then
         ! ADC
         if (acalc.eq.'adc') then
            call adc_turbomole(n,unit,nsta,e0,energy,maxdim,maxset)
         else
            write(6,'(/,2(2x,a),/)') &
                 "Calculation type not supported:",&
                 acalc(1:len_trim(acalc))
         endif

      else
         write(6,'(/,2(2x,a),/)') "File type not supported:",&
              aprog(1:len_trim(aprog))
         STOP
      endif

!------------------------------------------------------------------
! Sort the energies
!------------------------------------------------------------------
      call sorten(energy(n,1:nsta),nsta)

      return
    end subroutine rdenergy

!##################################################################

    subroutine coord_turbomole(n,unit,xcoo,maxdim,maxset,ncoo)

      implicit none
      
      integer*8                        :: n,unit,maxdim,maxset,&
                                          ncoo,i
      real*8, dimension(maxset,maxdim) :: xcoo
      character(len=80)                :: string

!------------------------------------------------------------------
! Read to coordinates
!------------------------------------------------------------------
10    continue
      read(unit,'(a)') string
      if (string(15:32).ne.'atomic coordinates') goto 10

!------------------------------------------------------------------
! Read the coordinates and convert to Angstrom
!------------------------------------------------------------------
      do i=1,ncoo/3
         read(unit,'(1x,3(F14.8))') &
              xcoo(n,i*3-2),xcoo(n,i*3-1),xcoo(n,i*3)
      enddo

      do i=1,ncoo
         xcoo(n,i)=xcoo(n,i)*0.529177249d0
      enddo

      return
    end subroutine coord_turbomole

!##################################################################

    subroutine eom_mol(n,unit,nsta,e0,energy,maxdim,maxset)

      implicit none

      integer*8                        :: n,unit,nsta,maxdim,&
                                          maxset,i,k,itmp
      real*8                           :: e0
      real*8, dimension(maxset,maxdim) :: energy
      character(len=80)                :: string

!------------------------------------------------------------------
! Read the ground state energy
!------------------------------------------------------------------
10    continue
      read(unit,'(a)') string
      if (string(3:19).ne.'CCSD total energy') goto 10
      
      read(string,'(36x,F18.12)') energy(n,1)

!------------------------------------------------------------------
! Read the excited state energies
!------------------------------------------------------------------
15    continue
      read(unit,'(a)') string
      if (string(3:15).ne.'Final Results') goto 15

      do i=1,4
         read(unit,*)
      enddo

      do i=2,nsta
         read(unit,'(24x,F14.8)') energy(n,i)
      enddo

!------------------------------------------------------------------
! Take energies relative to e0 and convert to eV
!------------------------------------------------------------------
      do i=1,nsta
         energy(n,i)=(energy(n,i)-e0)*27.2116d0
      enddo

      return
    end subroutine eom_mol

!##################################################################

    subroutine cas_mol(n,unit,nsta,e0,energy,maxdim,maxset)

      implicit none

      integer*8                        :: n,unit,nsta,maxdim,&
                                          maxset,ngeom,k
      real*8                           :: e0
      real*8, dimension(maxset,maxdim) :: energy
      real*8, dimension(maxdim)        :: entmp
      character(len=80)                :: string

!------------------------------------------------------------------
! Determine the number of geometries
!------------------------------------------------------------------
      rewind(unit)
      ngeom=0
10    continue
      read(unit,'(a)') string
      if (string(2:25).ne.'Variable memory released') then
         if (string(2:19).eq.'ATOMIC COORDINATES') ngeom=ngeom+1
         goto 10
      endif

!------------------------------------------------------------------
! Read to the last MULTI calculation
!------------------------------------------------------------------
      rewind(unit)
      k=0
15    continue
      read(unit,'(a)') string
      if (string(1:16).eq.'1PROGRAM * MULTI') k=k+1
      if (k.lt.ngeom) goto 15

!------------------------------------------------------------------
! Read the state energies
!------------------------------------------------------------------
      entmp=99999.0d0
      k=0
20    continue
      read(unit,'(a)') string
      if (string(2:15).ne.'State-averaged') then
         if (string(19:24).eq.'Energy') then
            k=k+1
            read(string,'(36x,F18.12)') entmp(k)
         endif
         goto 20
      endif

!------------------------------------------------------------------
! Sort the state energies
!------------------------------------------------------------------
      call sorten(entmp,nsta)

!------------------------------------------------------------------
! Take energies relative to e0 and convert to eV
!------------------------------------------------------------------
      do k=1,nsta
         energy(n,k)=(entmp(k)-e0)*27.2116d0
      enddo

      return
    end subroutine cas_mol

!##################################################################

    subroutine adc_turbomole(n,unit,nsta,e0,energy,maxdim,maxset)

      implicit none

      integer*8                        :: n,unit,nsta,maxdim,&
                                          maxset,i,j,k,ibreak,itmp
      real*8                           :: e0,ftmp,swap
      real*8, dimension(maxset,maxdim) :: energy
      character(len=80)                :: string
      logical(kind=4)                  :: lsym

!------------------------------------------------------------------
! Read the ground state energy
!------------------------------------------------------------------
10    continue
      read(unit,'(a)') string
      if (string(10:25).ne.'Final MP2 energy') goto 10
      
      read(string,'(53x,F15.10)') energy(n,1)
      
!------------------------------------------------------------------
! Read the excited state energies
!------------------------------------------------------------------
15    continue
      read(unit,'(a)') string
      if (string(67:72).ne.'||T1||') goto 15
      
      do i=1,3
         read(unit,*)
      enddo

      ! Check to see if symmetry has been used
      k=0
      lsym=.false.
20    continue
      read(unit,'(a)') string
      if (string(2:4).ne.'+==') then
         if (string(2:4).eq.'+--') then
            lsym=.true.
            k=k+1
            ibreak=k
         else
            k=k+1
         endif
         goto 20
      endif
      do i=1,k+1
         backspace(unit)
      enddo

      ! Read the excitation energies
      if (lsym) then
         itmp=nsta+1
      else
         itmp=nsta
      endif

      k=1
      do i=2,itmp
         if (i.ne.ibreak+1) then
            k=k+1            
            read(unit,'(28x,F9.7)') ftmp
            energy(n,k)=energy(n,1)+ftmp
         else
            read(unit,*)
         endif
      enddo

!------------------------------------------------------------------
! Take energies relative to e0 and convert to eV
!------------------------------------------------------------------
      do i=1,nsta
         energy(n,i)=(energy(n,i)-e0)*27.2116d0
      enddo

!------------------------------------------------------------------
! Sort the energies
!------------------------------------------------------------------
      do i=1,nsta-1
         do j=i+1,nsta
            if (energy(n,i).gt.energy(n,j)) then
               swap=energy(n,i)
               energy(n,i)=energy(n,j)
               energy(n,j)=swap
            endif
         enddo
      enddo

      return
    end subroutine adc_turbomole

!##################################################################

    subroutine x2prq(xcoo,xcoo0,qcoo,ncoo,nzero,invpvec)

      implicit none
      
      integer*8                          :: ncoo,nzero,i,j
      real*8, dimension(ncoo)            :: xcoo,xcoo0,qcoo
      real*8, dimension(ncoo-nzero,ncoo) :: invpvec

      qcoo=0.0d0
      do i=1,ncoo-nzero
         do j=1,ncoo
            qcoo(i)=qcoo(i)+invpvec(i,j)*(xcoo(j)-xcoo0(j))
         enddo
      enddo

      return
    end subroutine x2prq

!##################################################################

    subroutine sorten(energy,nsta)

      implicit none
      
      integer*8               :: nsta,i,j
      real*8, dimension(nsta) :: energy
      real*8                  :: ftmp

      do i=1,nsta-1
         do j=i+1,nsta
            if (energy(i).gt.energy(j)) then
               ftmp=energy(i)
               energy(i)=energy(j)
               energy(j)=ftmp
            endif
         enddo
      enddo

      return
    end subroutine sorten

!##################################################################

    subroutine wrinfo(ain,ainfo,ncoo,nzero,nsta,e0,xcoo0,maxdim,&
         mass,asym,freq,pvec,nqc,energy,maxset,qcoo,xcoo,ignore)

      implicit none

      integer*8                            :: i,j,k,m,unit,ncoo,&
                                              nzero,nsta,nqc,j1,&
                                              j2,maxdim,nblk,n,&
                                              maxset,s
      integer*8, dimension(maxset,maxdim)  :: ignore
      real*8                               :: e0
      real*8, dimension(maxdim)            :: xcoo0,mass,freq
      real*8, dimension(maxdim,maxdim)     :: pvec
      real*8, dimension(maxset,maxdim)     :: energy,qcoo,xcoo
      character(len=80)                    :: ain,ainfo
      character(len=16), dimension(maxdim) :: asym
      character(len=8)                     :: label

!------------------------------------------------------------------
! Open the .info file
!------------------------------------------------------------------
      ainfo=''
      k=index(ain,'.')
      write(ainfo(1:k),'(a)') ain(1:k)
      write(ainfo(k+1:k+4),'(a)') 'info'

      unit=20
      open(unit,file=ainfo,form='formatted',status='unknown')

!------------------------------------------------------------------
! Version number (dummy number)
!------------------------------------------------------------------
      write(unit,'(a,f8.5)') '# VCHam Version :',10.10050
      write(unit,'(a)') '# All energies and frequencies in eV'
      write(unit,'(a)') ' '  

!------------------------------------------------------------------
! System information
!------------------------------------------------------------------
      write(unit,'(a)') '#system_dimensions'
      write(unit,'(a)') '#internal_coordinates'
      write(unit,'(i3)') ncoo-nzero
      write(unit,'(a)')'#trivial_modes'
      write(unit,'(i3)') nzero
      write(unit,'(a)') '#num_states'
      write(unit,'(i3)') nsta
      write(unit,'(a)') '#states'
      write(unit,'(20i3)') (i, i=1,nsta)
      write(unit,'(a)') '#end_system_dimensions'

      write(unit,*)
      write(unit,'(a)') &
           '#general_information from the ground state'
      write(unit,'(a)') '#zero_of_energy'
      write(unit,'(F20.10)') e0
      write(unit,'(a)') '#equilibrium_state_geometry'
      do i=1,ncoo,3
         write(unit,'(3F15.10)') (xcoo0(j),j=i,i+2)
      enddo
      write(unit,'(a)') '#masses'
      do i=1,ncoo,3
         write(unit,'(F15.10)') mass(i)
      enddo

!------------------------------------------------------------------
! Normal mode information
!------------------------------------------------------------------
      write(unit,'(a)') '#vibrational_frequencies'
      do m=1,ncoo-nzero
         write(unit,'(i5,2x,a,F20.10)') m,asym(m),&
              freq(m)*0.6373641d0
      enddo

      write(unit,'(a)') '#normal_modes'
      nblk=((ncoo-1)/5)+1
      do m=1,ncoo-nzero
         do i=1,nblk
            j1=(i-1)*5+1
            j2=min(j1+4,ncoo)
            write(unit,'(5F15.10)') (pvec(j,m),j=j1,j2)
         enddo
         write(unit,*) 
      enddo

      write(unit,'(a)') '#end_general_information'
      write(unit,*) 

!------------------------------------------------------------------
! Coordinates and energies
!------------------------------------------------------------------
      write(unit,'(a)') '#dataset_section'

      nblk=((ncoo-nzero-1)/5)+1
      k=0
      do n=1,nqc

         do s=1,nsta

            if (ignore(n,s).eq.1) cycle

            k=k+1

            call getlabel(k,label)
            write(unit,'(a)') '#dataset_'//label//'    '
            write(unit,'(a)') '#energy'
            write(unit,'(f20.10)') energy(n,s)

            write(unit,'(a)') '#state'
            write(unit,'(i3)') s

            write(unit,'(a)') '#point_in_q'
            do i=1,nblk
               j1=(i-1)*5+1
               j2=min(j1+4,ncoo-nzero)
               write(unit,'(5F15.10)') (qcoo(n,j),j=j1,j2) 
            enddo

            write(unit,'(a)') '#point_in_cartesian_coordinates'
            do i=1,ncoo,3
               write(unit,'(3F16.10)') (xcoo(n,j),j=i,i+2)
            enddo

            write(unit,'(a)') '#end_dataset_'//label//'    '

         enddo
      enddo

!------------------------------------------------------------------
! End
!------------------------------------------------------------------
      write(unit,'(a)') '#end_dataset_section'

      return
    end subroutine wrinfo

!##################################################################

    subroutine rdsym(asymfile,ncoo,nzero,asym,maxdim)

      implicit none

      integer*8                            :: ncoo,nzero,unit,&
                                              maxdim,k
      character(len=80)                    :: asymfile,string
      character(len=16), dimension(maxdim) :: asym

!------------------------------------------------------------------
! Open file
!------------------------------------------------------------------
      unit=20
      open(unit,file=asymfile,form='formatted',status='old')

!------------------------------------------------------------------
! Read mode symmetries
!------------------------------------------------------------------
      asym=''
      k=0
10    continue
      read(unit,'(a)') string
      if (string.ne.'end-file') then
         k=k+1
         read(string,'(a)') asym(k)
         goto 10
      endif

!------------------------------------------------------------------
! Close file
!------------------------------------------------------------------
      close(unit)

!------------------------------------------------------------------
! Check that the correct no. symmetry labels have been given
!------------------------------------------------------------------
      if (k.ne.ncoo-nzero) then
         write(6,'(/,a,/,2(a,2x,i2,/))') &
              'No. symmetry labels given is not correct:',&
              'No. modes:',ncoo-nzero,&
              'No. labels given',k
         STOP
      endif

      return
    end subroutine rdsym

!##################################################################

    subroutine getlabel(n,string)
        
      implicit none
      
      integer*8, intent(in)         :: n
      character(len=8), intent(out) :: string
      
      string='     '
      if (n.le.9) then
         write(string(1:1),'(i1)') n
      else if (n.le.99) then
         write(string(1:2),'(i2)') n
      else if (n.le.999) then
         write(string(1:3),'(i3)') n
      else if (n.le.9999) then
         write(string(1:4),'(i4)') n
      else if (n.le.99999) then
         write(string(1:5),'(i5)') n
      else if (n.le.999999) then
         write(string(1:6),'(i6)') n
      else if (n.le.9999999) then
         write(string(1:7),'(i7)') n
      else if (n.le.99999999) then
         write(string(1:8),'(i8)') n
      endif
        
      return
    end subroutine getlabel

!##################################################################

  end program mkinfo
