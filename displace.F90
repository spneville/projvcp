  program displace

    use rddatamod

    implicit none

    integer*8                            :: iord,ncoo,nqc,icut,&
                                            npnts,nzero
    integer*8, parameter                 :: maxdim=60
    integer*8, parameter                 :: maxlen=200
    integer*8, dimension(maxlen)         :: geomflag
    integer*8, dimension(maxdim)         :: iimag
    real*8, dimension(maxdim*3)          :: xcoo0,qi,qf
    real*8                               :: diff
    real*8, dimension(maxdim)            :: peig,mass,freq
    real*8, dimension(maxdim,maxdim)     :: pvec,invpvec
    character(len=80)                    :: ain,aqcfile,astem,adat,&
                                            aprog
    character(len=2), dimension(maxdim)  :: aatm
    character(len=90), dimension(maxlen) :: aqc
    logical(kind=4)                      :: lcutall
    
!------------------------------------------------------------------
! Read input file name
!------------------------------------------------------------------
    ain=''
    call getarg(1,ain)
    
    if (ain.eq.'') then
       write(6,'(/,a,/)') 'Input file name not specified!'
       STOP
    endif

!------------------------------------------------------------------
! Read input file
!------------------------------------------------------------------
    call rdinp(ain,iord,aqcfile,ncoo,aatm,maxdim,xcoo0,diff,astem,&
         icut,npnts,adat,lcutall,aprog)

!------------------------------------------------------------------
! Read QC file
!------------------------------------------------------------------
    call rdqcfile(aqcfile,nqc,aqc,maxlen,geomflag)

!------------------------------------------------------------------
! Read the projected normal modes and frequencies
!------------------------------------------------------------------
    call rddat(adat,ncoo,pvec,peig,iimag,xcoo0,maxdim,mass)

!------------------------------------------------------------------
! Discard any unwanted modes and frequency scale the remaining 
! modes
!------------------------------------------------------------------
    call scaleq(peig,pvec,invpvec,nzero,ncoo,maxdim,freq,iimag,&
         mass)

!------------------------------------------------------------------
! Create the QC input files
!------------------------------------------------------------------
    call wrqcfiles(iord,xcoo0,maxdim,maxlen,ncoo,diff,geomflag,&
         aatm,aqc,nqc,astem,icut,npnts,nzero,pvec,lcutall,aprog)

    STOP
    
  contains

!##################################################################

    subroutine rdinp(ain,iord,aqcfile,ncoo,aatm,maxdim,xcoo0,&
         diff,astem,icut,npnts,adat,lcutall,aprog)

      use iomod
      
      implicit none
      
      integer*8                           :: iord,i,k,ncoo,maxdim,&
                                             icut,npnts
      real*8, dimension(maxdim*3)         :: xcoo0,qi,qf
      real*8                              :: diff
      character(len=80)                   :: ain,aqcfile,string,&
                                             astem,adat,aprog
      character(len=2), dimension(maxdim) :: aatm
      logical(kind=4)                     :: lcutall

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
      iord=0
      icut=0
      ncoo=0
      npnts=0
      aqcfile=''
      aatm=''
      astem=''
      adat=''
      aprog=''
      xcoo0=0.0d0
      diff=0.0d0
      qi=0.0d0
      qf=0.0d0
      lcutall=.false.

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

      else if (keyword(i).eq.'type') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            if (keyword(i).eq.'quadratic') then
               iord=2
            else if (keyword(i).eq.'quartic') then
               iord=4
            else if (keyword(i).eq.'cut') then
               iord=999
            else
               write(6,'(/,2(a,2x),/)') &
                    'Unknown job type:',keyword(i)
               STOP
            endif
         else
            write(6,'(/,a,/)') &
                 'No argument given with the type keyword'
            STOP
         endif

      else if (keyword(i).eq.'qc_file') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') aqcfile
         else
            write(6,'(/,a,/)') &
                 'No argument given with the qc_file keyword'
            STOP
         endif

      else if (keyword(i).eq.'reference_geometry') then
         k=0
16       continue
         call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,&
              maxkey)
         if (keyword(1).ne.'end-reference_geometry') then
            k=k+1
            read(keyword(1),'(a)') aatm(k)
            read(keyword(2),*) xcoo0(k*3-3+1)
            read(keyword(3),*) xcoo0(k*3-3+2)
            read(keyword(4),*) xcoo0(k*3-3+3)
            goto 16
         endif
         ! Set ncoo
         ncoo=k*3
         ! Convert atom labels to uppercase
         do k=1,ncoo/3
            call uppercase(aatm(k),2)
         enddo

      else if (keyword(i).eq.'step_size') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),*) diff
         else
            write(6,'(/,a,/)') &
                 'No argument given with the diff keyword'
            STOP
         endif

      else if (keyword(i).eq.'file_stem') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') astem
         else
            write(6,'(/,a,/)') &
                 'No argument given with the file_stem keyword'
            STOP
         endif

      else if (keyword(i).eq.'file_type') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') aprog
         else
            write(6,'(/,a,/)') &
                 'No argument given with the file_stem keyword'
            STOP
         endif

      else if (keyword(i).eq.'mode') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            if (keyword(i).eq.'all') then
               lcutall=.true.
            else
               read(keyword(i),*) icut
            endif
         else
            write(6,'(/,a,/)') &
                 'No argument given with the mode keyword'
            STOP
         endif

      else if (keyword(i).eq.'npoints') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),*) npnts
         else
            write(6,'(/,a,/)') &
                 'No argument given with the npoints keyword'
            STOP
         endif

      else if (keyword(i).eq.'datafile') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') adat
         else
            write(6,'(/,a,/)') &
                 'No argument given with the datafile keyword'
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
! Check that all required information has been specified
!------------------------------------------------------------------
      if (iord.eq.0) then
         write(6,'(/,a,/)') 'Job type has not been specified!'
         STOP
      endif

      if (aqcfile.eq.'') then
         write(6,'(/,a,/)') 'QC file has not been specified!'
         STOP
      endif

      if (ncoo.eq.0) then
         write(6,'(/,a,/)') &
              'Reference geometry has not been specified!'
         STOP
      endif

      if (diff.eq.0) then
         write(6,'(/,a,/)') &
              'Step size has not been specified!'
         STOP
      endif

      if (astem.eq.'') then
         write(6,'(/,a,/)') 'File stem has not been specified!'
         STOP
      endif

      if (iord.eq.999.and.icut.eq.0) then
         if (.not.lcutall) then
            write(6,'(/,a,/)') 'Cut mode has not been specified!'
            STOP
         endif
      endif

      if (iord.eq.999.and.npnts.eq.0) then
         write(6,'(/,a,/)') 'No. points has not been specified!'
         STOP
      endif

      if (iord.eq.999.and.adat.eq.'') then
         write(6,'(/,a,/)') 'Data file has not been specified!'
         STOP
      endif

      return
    end subroutine rdinp

!##################################################################

    subroutine uppercase(string,dim)

        implicit none

        integer*4          :: dim,j
        character(len=dim) :: string
        
        do j = 1,len(string)
           if(string(j:j).ge."a".and.string(j:j).le."z")&
                then
              string(j:j)=achar(iachar(string(j:j))-32)
           end if
        enddo

      return
    end subroutine uppercase

!##################################################################

    subroutine rdqcfile(aqcfile,nqc,aqc,maxlen,geomflag)

      implicit none

      integer*8                            :: maxlen,nqc,unit,i,k
      integer*8, dimension(maxlen)         :: geomflag
      character(len=80)                    :: aqcfile,string
      character(len=90), dimension(maxlen) :: aqc

!------------------------------------------------------------------
! Initialise arrays
!------------------------------------------------------------------
      geomflag=0

!------------------------------------------------------------------
! Open file
!------------------------------------------------------------------
      unit=20
      open(unit,file=aqcfile,form='formatted',status='old')

!------------------------------------------------------------------
! Read file and save to the aqc array
!------------------------------------------------------------------
      nqc=0
10    continue
      read(unit,'(a)') string
      
      if (string(1:8).ne.'end-file') then
         nqc=nqc+1
         read(string,'(a)') aqc(nqc)       
         k=len_trim(string)
         do i=1,k
            if (string(i:i).eq.'$') then
               if (string(i+1:i+4).eq.'geom') geomflag(nqc)=1
            endif
         enddo
         goto 10
      endif

!------------------------------------------------------------------
! Close file
!------------------------------------------------------------------
      close(unit)

      return
    end subroutine rdqcfile

!##################################################################

    subroutine wrqcfiles(iord,xcoo0,maxdim,maxlen,ncoo,diff,&
         geomflag,aatm,aqc,nqc,astem,icut,npnts,nzero,pvec,&
         lcutall,aprog)

      implicit none

      integer*8                            :: iord,maxdim,maxlen,&
                                              ncoo,nqc,icut,npnts,&
                                              nzero
      integer*8, dimension(maxlen)         :: geomflag
      real*8, dimension(maxdim*3)          :: xcoo0
      real*8                               :: diff
      real*8, dimension(maxdim,maxdim)     :: pvec
      character(len=2), dimension(maxdim)  :: aatm
      character(len=90), dimension(maxlen) :: aqc
      character(len=80)                    :: astem,aprog
      logical(kind=4)                      :: lcutall

!------------------------------------------------------------------
! Input files for the calculation of Hessians
!------------------------------------------------------------------
      if (iord.eq.2) call makeinp2(xcoo0,maxdim,maxlen,ncoo,&
           geomflag,diff,aatm,aqc,nqc,astem,aprog)

!------------------------------------------------------------------
! Input files for the calculation of quartic terms
!------------------------------------------------------------------
      if (iord.eq.4) call makeinp4(xcoo0,maxdim,maxlen,ncoo,&
           geomflag,diff,aatm,aqc,nqc,astem,aprog)

!------------------------------------------------------------------
! Cuts along the projected normal modes
!------------------------------------------------------------------
      if (iord.eq.999) call makecut(pvec,maxdim,icut,npnts,&
           nzero,diff,ncoo,xcoo0,nqc,geomflag,aqc,maxlen,astem,&
           lcutall,aprog)

      return
    end subroutine wrqcfiles

!##################################################################

    subroutine makeinp2(xcoo0,maxdim,maxlen,ncoo,geomflag,diff,&
         aatm,aqc,nqc,astem,aprog)
      
      implicit none
      
      integer*8                            :: ncoo,maxdim,maxlen,&
                                              s,m,n,dirm,k,ilbl,&
                                              jlbl,unit,i,j,l,nqc,&
                                              dirn
      integer*8, dimension(maxlen)         :: geomflag
      real*8, dimension(maxdim*3)          :: xcoo0
      real*8, dimension(ncoo)              :: xcoo
      real*8                               :: diff
      character(len=80)                    :: atmp,astem,aprog
      character(len=90), dimension(maxlen) :: aqc
      character(len=2), dimension(maxlen)  :: aatm

      unit=20

! Loop over DOFs
      do m=1,ncoo
         do n=m,ncoo

! Loop over negative and positive directions
            do k=1,2
               dirm=(-1)**k
               
               do l=1,2
                  dirn=(-1)**l

                  ! Skip the creation of unecessary input
                  if (n.eq.m.and.k.ne.l) goto 100

! Write the input for the current state/DOF/displacement
                 
                  ! Write filename
                  atmp=''
                  ilbl=len_trim(astem)
                  if (m.lt.10) then
                     write(atmp(1:ilbl+4),'(2a,i1)') astem(1:ilbl),'_x0',m
                  else
                     ilbl=len_trim(astem)
                     write(atmp(1:ilbl+4),'(2a,i2)') astem(1:ilbl),'_x',m
                  endif
                    
                  ilbl=len_trim(atmp)+1
                  if (dirm.lt.0) then
                     write(atmp(ilbl:ilbl),'(a1)') 'l'
                  else
                     write(atmp(ilbl:ilbl),'(a1)') 'r'
                  endif

                  ilbl=len_trim(atmp)+1
                  if (n.lt.10) then
                     write(atmp(ilbl:ilbl+3),'(a3,i1)') '_x0',n
                  else
                     write(atmp(ilbl:ilbl+3),'(a2,i2)') '_x',n
                  endif
                 
                  ilbl=len_trim(atmp)+1
                  if (dirn.lt.0) then
                     write(atmp(ilbl:ilbl+1),'(a5)') 'l'
                  else
                     write(atmp(ilbl:ilbl+1),'(a5)') 'r'
                  endif

                  call wrending(atmp,aprog)

                  ! Open file
                  open(unit,file=atmp,form='formatted',status='unknown')

                  ! Calculate the current coordinates
                  do i=1,ncoo
                     xcoo(i)=xcoo0(i)
                     if (i.eq.m) xcoo(i)=xcoo(i)+dirm*diff
                     if (i.eq.n) xcoo(i)=xcoo(i)+dirn*diff
                     if (i.eq.m.and.m.eq.n) xcoo(i)=xcoo(i)-dirn*diff
                  enddo

                  ! Write to file
                  do i=1,nqc
                     if (geomflag(i).eq.1) then
                        do j=1,ncoo/3
                           write(unit,'(a2,3(F11.8,2x))') aatm(j),&
                                xcoo(j*3-2),xcoo(j*3-1),xcoo(j*3)
                        enddo
                     else
                        ilbl=len_trim(aqc(i))
                        write(unit,'(a)') aqc(i)(1:ilbl)
                     endif
                  enddo
                  
                  ! Close file
                  close(unit)
                  
100               continue
                  
               enddo
            enddo
            
         enddo
      enddo

      return
    end subroutine makeinp2

!##################################################################

    subroutine makeinp4(xcoo0,maxdim,maxlen,ncoo,geomflag,diff,&
         aatm,aqc,nqc,astem,aprog)

      implicit none

      integer*8                            :: ncoo,maxdim,maxlen,&
                                              nqc,n,dir,ilbl,k,&
                                              unit1,unit2,i,j
      integer*8, dimension(maxlen)         :: geomflag
      real*8, dimension(maxdim*3)          :: xcoo0
      real*8, dimension(ncoo)              :: xcoo1,xcoo2
      real*8                               :: diff
      character(len=80)                    :: atmp1,atmp2,astem,aprog
      character(len=90), dimension(maxlen) :: aqc
      character(len=2), dimension(maxlen)  :: aatm

      unit1=20
      unit2=30

! Loop over DOFs
      do n=1,ncoo

! Loop over negative and positive directions
          do k=1,2
             dir=(-1)**k
             
! Write the input for the current state/DOF/displacement

             ! Write filenames
             atmp1=''
             atmp2=''
             
             ilbl=len_trim(astem)

             if (n.lt.10) then
                write(atmp1(1:ilbl+7),'(2a,i1)') &
                     astem(1:ilbl),'_d1_x0',n
                write(atmp2(1:ilbl+7),'(2a,i1)') &
                     astem(1:ilbl),'_d2_x0',n
             else
                write(atmp1(1:ilbl+7),'(2a,i2)') &
                     astem(1:ilbl),'_d1_x',n
                write(atmp2(1:ilbl+7),'(2a,i2)') &
                     astem(1:ilbl),'_d2_x',n
             endif

             ilbl=len_trim(atmp1)+1
             if (dir.lt.0) then
                write(atmp1(ilbl:ilbl+1),'(a)') 'l'
                write(atmp2(ilbl:ilbl+1),'(a)') 'l'
             else
                write(atmp1(ilbl:ilbl+1),'(a)') 'r'
                write(atmp2(ilbl:ilbl+1),'(a)') 'r'
             endif

             call wrending(atmp1,aprog)
             call wrending(atmp2,aprog)

             ! Open file
             open(unit1,file=atmp1,form='formatted',status='unknown')
             open(unit2,file=atmp2,form='formatted',status='unknown')


             ! Calculate the current coordinates
             do i=1,ncoo
                xcoo1(i)=xcoo0(i)
                xcoo2(i)=xcoo0(i)
                if (i.eq.n) then
                   xcoo1(i)=xcoo1(i)+dir*diff
                   xcoo2(i)=xcoo2(i)+2*dir*diff
                endif
             enddo
             
             ! Write to files
             do i=1,nqc
                if (geomflag(i).eq.1) then
                   do j=1,ncoo/3
                      write(unit1,'(a2,3(F11.8,2x))') aatm(j),&
                           xcoo1(j*3-2),xcoo1(j*3-1),xcoo1(j*3)
                      write(unit2,'(a2,3(F11.8,2x))') aatm(j),&
                           xcoo2(j*3-2),xcoo2(j*3-1),xcoo2(j*3)
                   enddo
                else
                   ilbl=len_trim(aqc(i))
                   write(unit1,'(a)') aqc(i)(1:ilbl)
                   write(unit2,'(a)') aqc(i)(1:ilbl)
                endif
             enddo

             ! Close files
             close(unit1)
             close(unit2)

          enddo

      enddo

      return
    end subroutine makeinp4

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
! Discard the unwanted zero-frequency modes
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

    subroutine makecut(pvec,maxdim,icut,npnts,nzero,diff,ncoo,&
         xcoo0,nqc,geomflag,aqc,maxlen,astem,lcutall,aprog)

      implicit none

      integer*8                            :: maxdim,icut,npnts,&
                                              nzero,ncoo,i,j,k,&
                                              ilbl,nqc,maxlen,&
                                              unit,mi,mf,n
      integer*8, dimension(maxlen)         :: geomflag
      real*8, dimension(maxdim,maxdim)     :: pvec
      real*8, dimension(maxdim*3)          :: xcoo0
      real*8                               :: diff
      real*8, dimension(ncoo-nzero)        :: q
      real*8, dimension(ncoo)              :: x
      character(len=90), dimension(maxlen) :: aqc
      character(len=80)                    :: atmp,astem,aprog
      logical(kind=4)                      :: lcutall

      x=0.0d0
      

      unit=20

      if (lcutall) then
         mi=1
         mf=ncoo-nzero
      else
         mi=icut
         mf=icut
      endif

      ! Loop over modes
      do n=mi,mf

         ! Loop over the points
         do i=-npnts,npnts

            ! Calculate the current Cartesian coordinates
            q=0.0d0
            q(n)=diff*i
            do j=1,ncoo
               x(j)=xcoo0(j)
               do k=1,ncoo-nzero
                  x(j)=x(j)+pvec(j,k)*q(k)
               enddo
            enddo
            
            ! Open the output file
            atmp=''
            ilbl=len_trim(astem)
            write(atmp(1:ilbl),'(a)') astem(1:ilbl)
            ilbl=ilbl+1
            write(atmp(ilbl:ilbl+1),'(a)') '_q'
            ilbl=ilbl+2
            
            if (n.lt.10) then
               write(atmp(ilbl:ilbl+2),'(a1,i1,a1)') '0',n,'_'
            else
               write(atmp(ilbl:ilbl+2),'(i2,a1)') n,'_'
            endif
            ilbl=ilbl+3
            
            if (abs(i).lt.10.and.i.ne.0) then
               write(atmp(ilbl:ilbl+1),'(a,i1)') '0',abs(i)
            else if (abs(i).ge.10) then
               write(atmp(ilbl:ilbl+1),'(i2)') abs(i)
            else
               write(atmp(ilbl:ilbl+1),'(a)') '00'
            endif
            ilbl=ilbl+2
            
            if (i.lt.0) then
               write(atmp(ilbl:ilbl+1),'(a)') 'l'
            else if (i.gt.0) then
               write(atmp(ilbl:ilbl+1),'(a)') 'r'
            else
               write(atmp(ilbl:ilbl+1),'(a)') '0'
            endif
            
            call wrending(atmp,aprog)

            open(unit,file=atmp,form='formatted',status='unknown')
            
            ! Write to the output file
            do j=1,nqc
               if (geomflag(j).eq.1) then
                  do k=1,ncoo/3
                     write(unit,'(a2,3(F11.8,2x))') aatm(k),&
                          x(k*3-2),x(k*3-1),x(k*3)
                  enddo
               else
                  ilbl=len_trim(aqc(j))
                  write(unit,'(a)') aqc(j)(1:ilbl)
               endif
            enddo
            
            ! Close the output file
            close(unit)
            
         enddo

      enddo

      return
    end subroutine makecut

!##################################################################

    subroutine wrending(atmp,aprog)

      implicit none
      
      integer*8         :: k
      character(len=80) :: atmp,aprog

      k=len_trim(atmp)+1

      if (aprog.eq.'xyz') then
         write(atmp(k:k+4),'(a4)') '.xyz'
      else if (aprog.eq.'gaussian') then
         write(atmp(k:k+4),'(a4)') '.com'
      else if (aprog.eq.'molpro') then
         write(atmp(k:k+4),'(a4)') '.inp'
      else if (aprog.eq.'') then
         write(atmp(k:k+4),'(a4)') '.inp'
      else
         write(6,'(/,2(2x,a),/)') 'Unknown file type:',aprog
         STOP
      endif

      return
    end subroutine wrending

!##################################################################

  end program displace
