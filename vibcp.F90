  program vibcp

    use rddatamod

    implicit none

    integer*8                                     :: ncoo,ngrad,&
                                                     nnact,nhess,&
                                                     nquart,nsta,&
                                                     nzero
    integer*8, parameter                          :: maxdim=60
    integer*8, dimension(maxdim)                  :: iimag,&
                                                     zeromodes
    real*8, dimension(maxdim)                     :: peig,xcoo0,&
                                                     qcxcoo,e0,&
                                                     mass,freq,&
                                                     enact
    real*8, dimension(maxdim,maxdim)              :: pvec,grad,&
                                                     quart
    real*8, dimension(maxdim,maxdim,maxdim)       :: nact,hess
    real*8                                        :: qdiff,qdiff2
    character(len=80)                             :: ain,adat,&
                                                     agrad,anact,&
                                                     ahess,aoper,&
                                                     aquart
    character(len=80), dimension(maxdim*2)        :: gradfile
    character(len=80), dimension(maxdim)          :: nactfile
    character(len=80), dimension(maxdim*maxdim*4) :: hessfile
    character(len=80), dimension(maxdim*4)        :: quartfile
    logical(kind=4)                               :: lgrad,lnact,&
                                                     lhess,lquart

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
    call rdinp(ain,adat,agrad,anact,ahess,nsta,qcxcoo,maxdim,e0,&
         zeromodes,nzero,aoper,enact,aquart,qdiff,qdiff2,lgrad,&
         lnact,lhess,lquart)
    
!------------------------------------------------------------------
! Read the data file
!------------------------------------------------------------------
    call rddat(adat,ncoo,pvec,peig,iimag,xcoo0,maxdim,mass)

!------------------------------------------------------------------
! Discard any unwanted modes and frequency scale the remaining 
! modes
!------------------------------------------------------------------
    call scaleq(pvec,peig,freq,zeromodes,iimag,ncoo,maxdim,nzero,&
         mass)

!------------------------------------------------------------------
! Read the set files
!------------------------------------------------------------------
    call rdsets(agrad,anact,ahess,gradfile,nactfile,hessfile,&
         maxdim,ngrad,nnact,nhess,aquart,quartfile,nquart,lgrad,&
         lnact,lhess,lquart)

!------------------------------------------------------------------
! Read the QC output files
!------------------------------------------------------------------
    call rdqc(gradfile,nactfile,hessfile,ngrad,nnact,nhess,maxdim,&
         nsta,ncoo,grad,nact,hess,qcxcoo,e0,mass,quartfile,nquart,&
         quart,qdiff,qdiff2,lgrad,lnact,lhess,lquart)

!------------------------------------------------------------------
! Rotate the gradient vectors, NACT vectors and Hessian matrices
! to the same frame as the reference geometry
!------------------------------------------------------------------
   call rotate(grad,nact,hess,maxdim,ncoo,xcoo0,qcxcoo,mass,nsta,&
        quart,lgrad,lnact,lhess,lquart)

!------------------------------------------------------------------
! Transform the gradient vectors, NACT vectors and Hessian 
! matrices to be in terms of the eigenvectors of the ground state
! projected Hessian
!------------------------------------------------------------------
   call transform(grad,nact,hess,maxdim,ncoo,nsta,pvec,nzero,&
        quart,lgrad,lnact,lhess,lquart)

!------------------------------------------------------------------
! Output the vibronic coupling coefficients
!------------------------------------------------------------------
   call outcoeffs(grad,nact,hess,e0,maxdim,ncoo,nsta,freq,peig,&
        iimag,nzero,aoper,enact,quart,lgrad,lnact,lhess,lquart)

    STOP

    contains

!##################################################################

    subroutine rdinp(ain,adat,agrad,anact,ahess,nsta,qcxcoo,&
         maxdim,e0,zeromodes,nzero,aoper,enact,aquart,qdiff,&
         qdiff2,lgrad,lnact,lhess,lquart)

      use iomod

      implicit none

      integer*8                    :: i,j,k,nsta,maxdim,nzero
      integer*8, dimension(maxdim) :: zeromodes
      real*8, dimension(maxdim)    :: qcxcoo,e0,enact
      real*8                       :: qdiff,qdiff2,ftmp
      character(len=80)            :: ain,adat,agrad,anact,ahess,&
                                      aoper,aquart,atmp
      logical(kind=4)              :: lgrad,lnact,lhess,lquart

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
      agrad=''
      anact=''
      ahess=''
      aoper=''
      aquart=''
      nsta=0
      nzero=0
      e0=0.0d0
      enact=0.0d0
      qdiff=0.0d0
      qdiff2=0.0d0
      zeromodes=0

      lgrad=.false.
      lnact=.false.
      lhess=.false.
      lquart=.false.

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

      else if (keyword(i).eq.'modes') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') adat
         else
            write(6,'(/,a,/)') &
                 'No argument given with the modes keyword'
            STOP
         endif

      else if (keyword(i).eq.'gradient_files') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') agrad
            lgrad=.true.
         else
            write(6,'(/,a,/)') &
                 'No argument given with the gradient_files keyword'
            STOP
         endif

      else if (keyword(i).eq.'nact_files') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') anact
            lnact=.true.
         else
            write(6,'(/,a,/)') &
                 'No argument given with the nact_files keyword'
            STOP
         endif

      else if (keyword(i).eq.'hessian_files') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') ahess
            lhess=.true.
         else
            write(6,'(/,a,/)') &
                 'No argument given with the hessian_files keyword'
            STOP
         endif

      else if (keyword(i).eq.'quartic_files') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),'(a)') aquart
            lquart=.true.
         else
            write(6,'(/,a,/)') &
                 'No argument given with the quartic_files keyword'
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

      else if (keyword(i).eq.'qc_reference_geometry') then
         k=0
16       continue
         call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,&
              maxkey)
         if (keyword(1).ne.'end-qc_reference_geometry') then
            do j=1,3
               k=k+1
               read(keyword(j),*) qcxcoo(k)
            enddo
            goto 16
         endif
         i=ic

      else if (keyword(i).eq.'reference_energies') then
         k=0
17       continue
         call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,&
              maxkey)
         if (keyword(1).ne.'end-reference_energies') then
            k=k+1
            read(keyword(1),*) e0(k)            
            goto 17
         endif
         i=ic

      else if (keyword(i).eq.'zero_modes') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            read(keyword(i),*) k
            zeromodes(k)=1
            nzero=nzero+1
18          continue
            if (keyword(i+1).eq.',') then
               i=i+2
               read(keyword(i),*) k
               zeromodes(k)=1
               nzero=nzero+1
               goto 18
            endif
         else
            write(6,'(/,a,/)') &
                 'No argument given with the zero_modes keyword'
            STOP
         endif

      else if (keyword(i).eq.'oper_name') then
         if (keyword(i+1).eq.'=') then
            i=i+2
            aoper=keyword(i)
         else
            write(6,'(/,a,/)') &
                 'No argument given with the oper_name keyword'
            STOP
         endif

      else if (keyword(i).eq.'nact_energies') then
         k=0
19       continue
         call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,&
              maxkey)
         if (keyword(1).ne.'end-nact_energies') then
            k=k+1
            read(keyword(1),*) enact(k)
            goto 19
         endif
         i=ic

      else if (keyword(i).eq.'quartic_displacements') then
         if (keyword(i+1).eq.'=') then
            i=i+2            
            read(keyword(i),*) qdiff
            i=i+2
            read(keyword(i),*) qdiff2
            if (qdiff.gt.qdiff2) then
               ftmp=qdiff
               qdiff=qdiff2
               qdiff2=qdiff
            endif
            if (keyword(i+1).eq.',') then
               i=i+2
               read(keyword(i),'(a)') atmp
               if (atmp.eq.'angstrom') then
                  qdiff=qdiff/0.529177249d0
                  qdiff2=qdiff2/0.529177249d0
               else
                  write(6,'(/,2(a,2x),/)') &
                       'Unknown conversion keyword:',atmp
                  STOP
               endif
            endif
         else
            write(6,'(/,a,/)') &
                 'No argument given with the quartic_displacements keyword'
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
         write(6,'(/,a,/)') 'Data file not specified!'
         STOP
      endif

      if (lgrad.and.agrad.eq.'') then
         write(6,'(/,a,/)') 'Gradient set file not specified!'
         STOP
      endif

      if (lnact.and.anact.eq.'') then
         write(6,'(/,a,/)') 'NACT set file not specified!'
         STOP
      endif

      if (lhess.and.ahess.eq.'') then
         write(6,'(/,a,/)') 'Hessian set file not specified!'
         STOP
      endif

      if (lquart.and.aquart.eq.'') then
         write(6,'(/,a,/)') 'Quartic set file not specified!'
         STOP
      endif

      if (nsta.eq.0) then
         write(6,'(/,a,/)') 'No. states not specified!'
         STOP
      endif

      if (e0(1).eq.0) then
         write(6,'(/,a,/)') 'Reference energy not specified!'
         STOP
      endif
      
      if (enact(1).eq.0) then
         do i=1,maxdim
            enact(i)=e0(i)
         enddo
      endif

      if (aoper.eq.'') aoper='vibcp.op'

      return
    end subroutine rdinp

!##################################################################

    subroutine scaleq(pvec,peig,freq,zeromodes,iimag,ncoo,maxdim,&
         nzero,mass)

      implicit none

      integer*8                        :: maxdim,ncoo,nzero,i,j,k
      integer*8, dimension(maxdim)     :: zeromodes,iimag
      integer*8, dimension(ncoo)       :: itmp1
      real*8                           :: ftmp
      real*8, dimension(maxdim)        :: peig,freq,mass
      real*8, dimension(maxdim,maxdim) :: pvec
      real*8, dimension(ncoo)          :: ftmp1
      real*8, dimension(ncoo,ncoo)     :: ftmp2

!------------------------------------------------------------------
! Discard unwanted modes
!------------------------------------------------------------------
      ftmp1=0.0d0
      ftmp2=0.0d0
      itmp1=0

      k=0
      do i=1,ncoo
         if (zeromodes(i).ne.1) then
            k=k+1
            ftmp1(k)=peig(i)
            itmp1(k)=iimag(i)
            do j=1,ncoo
               ftmp2(j,k)=pvec(j,i)
            enddo
         endif
      enddo

      peig=0.0d0
      pvec=0.0d0
      iimag=0
      do i=1,ncoo
         peig(i)=ftmp1(i)
         freq(i)=sqrt(peig(i))
         iimag(i)=itmp1(i)
         do j=1,ncoo
            pvec(i,j)=ftmp2(i,j)
         enddo
      enddo

!------------------------------------------------------------------
! Frequency scale the remaining modes
!------------------------------------------------------------------
      do i=1,ncoo-nzero
         ftmp=freq(i)*0.6373641d0
         do j=1,ncoo
            pvec(j,i)=pvec(j,i)/(15.4644*sqrt(ftmp)*sqrt(mass(j)))
         enddo
      enddo

      return
    end subroutine scaleq

!##################################################################

    subroutine rdsets(agrad,anact,ahess,gradfile,nactfile,&
         hessfile,maxdim,ngrad,nnact,nhess,aquart,quartfile,&
         nquart,lgrad,lnact,lhess,lquart)
      
      use iomod

      implicit none

      integer*8                              :: i,maxdim,ngrad,&
                                                nnact,nhess,nquart
      character(len=80)                      :: agrad,anact,ahess,&
                                                aquart
      character(len=80), dimension(maxdim*3) :: gradfile
      character(len=80), dimension(maxdim)   :: nactfile
      character(len=80), dimension(maxdim*maxdim*4) :: hessfile
      character(len=80), dimension(maxdim*4) :: quartfile
      logical(kind=4)                        :: lgrad,lnact,lhess,&
                                                lquart

      integer*8                                   :: iin,ic,iz,&
                                                     ierr
      integer*8, parameter                        :: maxkey=48,&
                                                     maxkeylen=200
      integer*8, dimension(maxkey)                :: lc      
      character(len=maxkeylen), dimension(maxkey) :: keyword,&
                                                     keyorig
      character(len=240)                          :: inptit

!------------------------------------------------------------------
! Read gradient filenames
!------------------------------------------------------------------
      if (lgrad) then
         iin=20
         open(iin,file=agrad,form='formatted',status='old')
         
         ngrad=0
10       continue
         call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,maxkey)
         if (keyword(1).ne.'end-files') then
            ngrad=ngrad+1
            read(keyword(1),'(a)') gradfile(ngrad)
            goto 10
         endif
         
         close(iin)
      endif

!------------------------------------------------------------------
! Read NACT filenames
!------------------------------------------------------------------
      if (lnact) then
         open(iin,file=anact,form='formatted',status='old')

         nnact=0
20       continue
         call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,maxkey)
         if (keyword(1).ne.'end-files') then
            nnact=nnact+1
            read(keyword(1),'(a)') nactfile(nnact)
            goto 20
         endif
         
         close(iin)
      endif

!------------------------------------------------------------------
! Read Hessian filenames
!------------------------------------------------------------------
      if (lhess) then
         open(iin,file=ahess,form='formatted',status='old')

         nhess=0
30       continue
         call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,maxkey)
         if (keyword(1).ne.'end-files') then
            nhess=nhess+1
            read(keyword(1),'(a)') hessfile(nhess)
            goto 30
         endif
         
         close(iin)
      endif

!------------------------------------------------------------------
! Read Quartic filenames
!------------------------------------------------------------------
      if (lquart) then
         open(iin,file=aquart,form='formatted',status='old')
         
         nquart=0
40       continue
         call rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,maxkey)
         if (keyword(1).ne.'end-files') then
            nquart=nquart+1
            read(keyword(1),'(a)') quartfile(nquart)
            goto 40
         endif
         
         close(iin)
      endif

      return
    end subroutine rdsets

!##################################################################

    subroutine rdqc(gradfile,nactfile,hessfile,ngrad,nnact,nhess,&
         maxdim,nsta,ncoo,grad,nact,hess,qcxcoo,e0,mass,quartfile,&
         nquart,quart,qdiff,qdiff2,lgrad,lnact,lhess,lquart)

      implicit none

      integer*8                               :: maxdim,ngrad,&
                                                 nnact,nhess,i,j,&
                                                 k,nsta,ncoo,nquart
      real*8, dimension(maxdim)               :: qcxcoo,e0,mass
      real*8, dimension(maxdim,maxdim)        :: grad,quart
      real*8, dimension(maxdim,maxdim,maxdim) :: nact,hess
      real*8, dimension(ncoo,ncoo,nsta)       :: pospos,negneg,&
                                                 posneg,negpos
      real*8, dimension(ncoo,nsta)            :: pos,neg,pos2,neg2
      real*8                                  :: diff,qdiff,&
                                                 qdiff2
      character(len=80), dimension(maxdim*3)  :: gradfile
      character(len=80), dimension(maxdim)    :: nactfile
      character(len=80), dimension(maxdim*maxdim*4) :: hessfile
      character(len=80), dimension(maxdim*4)  :: quartfile
      logical(kind=4)                         :: lgrad,lnact,&
                                                 lhess,lquart

!------------------------------------------------------------------
! Initialise arrays
!------------------------------------------------------------------
      grad=0.0d0
      nact=0.0d0
      hess=0.0d0
      pos=0.0d0
      neg=0.0d0
      pospos=0.0d0
      negneg=0.0d0
      posneg=0.0d0
      negpos=0.0d0
      
!------------------------------------------------------------------
! Read gradient files
!------------------------------------------------------------------
      if (lgrad) then
         write(6,'(/,3x,a)') 'Reading the gradient files...'
         do i=1,ngrad
            call rdgrad(gradfile(i),nsta,ncoo,grad,maxdim)
         enddo
         
!         ! Mass-weight the gradients
!         do i=1,nsta
!            do j=1,ncoo
!               grad(i,j)=grad(i,j)/sqrt(mass(j))
!            enddo
!         enddo
      endif

!------------------------------------------------------------------
! Read the NACT files
!------------------------------------------------------------------
      if (lnact) then
         write(6,'(/,3x,a)') 'Reading the NACT files...'
         do i=1,nnact
            call rdnact(nactfile(i),nsta,ncoo,nact,maxdim,qcxcoo)
         enddo
         
!         ! Mass-weight the NACTs
!         do i=1,nsta
!            do j=1,nsta
!               do k=1,ncoo
!                  nact(i,j,k)=nact(i,j,k)/sqrt(mass(k))
!               enddo
!            enddo
!         enddo
      endif

!------------------------------------------------------------------
! Read the Quartic files
!------------------------------------------------------------------
      if (lquart) then
         write(6,'(/,3x,a)') 'Reading the quartic files...'
         do i=1,nquart
            call rdquart(quartfile(i),nsta,ncoo,maxdim,qcxcoo,qdiff,&
                 qdiff2,e0,pos,neg,pos2,neg2)
         enddo
         
         ! Calculate the 4th-order derivatives of the energies
         do i=1,nsta
            do j=1,ncoo            
               quart(i,j)=pos2(j,i)+neg2(j,i)-4*pos(j,i)-4*neg(j,i)&
                    +6*(e0(i)-e0(1))
               quart(i,j)=quart(i,j)/(qdiff**4)
            enddo
         enddo
         
         ! Convert to eV/Angstrom^4
         do i=1,nsta
            do j=1,ncoo
               quart(i,j)=quart(i,j)*27.2116d0/(0.529177249d0**4)
            enddo
         enddo
         
!         ! Mass-weighting the quartic terms
!         do i=1,nsta
!            do j=1,ncoo
!               quart(i,j)=quart(i,j)/sqrt(mass(j)**4)
!            enddo
!         enddo
      endif

!------------------------------------------------------------------
! Read the Hessian files
!------------------------------------------------------------------
      if (lhess) then
         write(6,'(/,3x,a)') 'Reading the Hessian files...'
         do i=1,nhess
            call rdhess(hessfile(i),nsta,ncoo,pospos,negneg,&
                 posneg,negpos,maxdim,qcxcoo,e0,diff,pos,neg)
         enddo
         
         ! Calculate the elements of the Hessians
         ! (i) Bilinear terms
         do i=1,nsta
            do j=1,ncoo-1
               do k=j+1,ncoo
                  hess(i,j,k)=pospos(j,k,i)+negneg(j,k,i)&
                       -posneg(j,k,i)-negpos(j,k,i)
                  hess(i,j,k)=hess(i,j,k)/(4*diff*diff)
                  hess(i,k,j)=hess(i,j,k)
               enddo
            enddo
         enddo
         
         ! (ii) Quadratic terms
         do i=1,nsta
            do j=1,ncoo
               hess(i,j,j)=pos(j,i)+neg(j,i)-2*(e0(i)-e0(1))
               hess(i,j,j)=hess(i,j,j)/(diff*diff)
            enddo
         enddo
         
         ! Convert to eV/Angstrom^2
         do i=1,nsta
            do j=1,ncoo
               do k=1,ncoo
                  hess(i,j,k)=hess(i,j,k)*27.2116d0/(0.529177249d0**2)
               enddo
            enddo
         enddo
         
!         ! Mass-weighting the Hessians
!         do i=1,nsta
!            do j=1,ncoo
!               do k=1,ncoo
!                  hess(i,j,k)=hess(i,j,k)/sqrt(mass(j)*mass(k))
!               enddo
!            enddo
!         enddo
      endif

      return
    end subroutine rdqc

!##################################################################

    subroutine rdgrad(filename,nsta,ncoo,grad,maxdim)

      implicit none

      integer*8                         :: unit,i,j,k,ncalc,&
                                           nsta,ncoo,maxdim
      real*8, dimension(nsta)           :: en
      real*8                            :: ftmp
      real*8, dimension(maxdim,maxdim)  :: grad
      character(len=80)                 :: filename,string
      character(len=4), dimension(nsta) :: asta
      character(len=4)                  :: atmp

!------------------------------------------------------------------
! Set inital values
!------------------------------------------------------------------
      do i=1,nsta
         en(i)=0.0d0
         asta(i)=''
      enddo

!------------------------------------------------------------------
! Open file
!------------------------------------------------------------------
      unit=20
      open(unit,file=filename,form='formatted',status='old')

!------------------------------------------------------------------
! Read to final MCSCF calculation
!------------------------------------------------------------------
      ! Determine number of MCSCF calculations
      ncalc=0
10    continue
      read(unit,'(a)',end=900) string
      if (string(2:25).ne.'Variable memory released') then
         if (string(12:48).eq.'MULTI (Direct Multiconfiguration SCF)') then
            ncalc=ncalc+1
         endif
         goto 10
      endif
900   continue

      ! Read to the final calculation
      rewind(unit)
      k=0
20    continue
      read(unit,'(a)') string
      if (string(12:48).eq.'MULTI (Direct Multiconfiguration SCF)') k=k+1
      if (k.ne.ncalc) goto 20

!------------------------------------------------------------------
! Read state energies and labels
!------------------------------------------------------------------
      k=0
30    continue
      read(unit,'(a)') string
      if (string(2:15).ne.'State-averaged') then
         if (string(19:24).eq.'Energy') then            
            k=k+1
            read(string,'(36x,F18.12)') en(k)
            read(string,'(13x,a4)') asta(k)
         endif

         if (k.ne.nsta) goto 30
      endif

!------------------------------------------------------------------
! Sort state energies and labels
!------------------------------------------------------------------
      do i=1,nsta-1
         do j=i+1,nsta
            if (en(i).gt.en(j)) then
               ftmp=en(i)
               en(i)=en(j)
               en(j)=ftmp
               atmp=asta(i)
               asta(i)=asta(j)
               asta(j)=atmp
            endif
         enddo
      enddo

!------------------------------------------------------------------
! Read the gradient vector
!------------------------------------------------------------------
40    continue
      read(unit,'(a)') string
      if (string(2:15).ne.'SA-MC GRADIENT') goto 40

      ! Determine which state we are dealing with
      read(string,'(25x,a4)') atmp
      do i=1,nsta
         if (atmp.eq.asta(i)) k=i
      enddo

      ! Skip to the gradient vector
      do i=1,3
         read(unit,*)
      enddo

      ! Read the gradient vector
      do i=1,ncoo/3
         read(unit,'(4x,3(8x,F12.9))') grad(k,i*3-2),&
              grad(k,i*3-1),grad(k,i*3)         
      enddo

      ! Convert to eV/Angstrom
      do i=1,ncoo
         grad(k,i)=grad(k,i)/0.529177249d0
         grad(k,i)=grad(k,i)*27.2116d0
      enddo

!------------------------------------------------------------------
! Close file
!------------------------------------------------------------------
      close(unit)

      return
    end subroutine rdgrad

!##################################################################
    
    subroutine rdnact(filename,nsta,ncoo,nact,maxdim,qcxcoo)

      implicit none

      integer*8                               :: nsta,ncoo,maxdim,&
                                                 i,j,k,ngeom,unit,&
                                                 itmp,s,s1
      real*8, dimension(maxdim,maxdim,maxdim) :: nact
      real*8, dimension(maxdim)               :: qcxcoo,tmpxcoo
      real*8, parameter                       :: tol=10E-7
      character(len=80)                       :: filename,string

!------------------------------------------------------------------
! Open file
!------------------------------------------------------------------
      unit=20
      open(unit,file=filename,form='formatted',status='old')

!------------------------------------------------------------------
! Read to final geometry
!------------------------------------------------------------------
      ngeom=0
10    continue
      read(unit,'(a)',end=100) string
      if (string(2:19).eq.'ATOMIC COORDINATES') ngeom=ngeom+1
      goto 10

100   continue
      rewind(unit)

      k=0
20    continue
      read(unit,'(a)') string
      if (string(2:19).eq.'ATOMIC COORDINATES') k=k+1
      if (k.ne.ngeom) goto 20

!------------------------------------------------------------------
! Read the displaced geometry
!------------------------------------------------------------------
      ! Skip to the coordinates
      do i=1,3
         read(unit,*)
      enddo

      ! Read the coordinates
      do i=1,ncoo,3
         read(unit,'(18x,3(3x,F12.9))') (tmpxcoo(j), j=i,i+2)
      enddo

      ! Convert to Angstrom for comparison with the reference
      ! geometry
      do i=1,ncoo
         tmpxcoo(i)=tmpxcoo(i)*0.529177249d0
      enddo

!------------------------------------------------------------------
! Determine which atom has been displaced
!------------------------------------------------------------------
      do i=1,ncoo
         if (abs(tmpxcoo(i)-qcxcoo(i)).gt.tol) itmp=i
      enddo

!------------------------------------------------------------------
! Read the NACTs
!------------------------------------------------------------------
30    continue
      read(unit,'(a)',end=200) string
      if (string(2:24).eq.'Construct non-adiabatic') then
       
         ! Skip to state info
         do i=1,13
            read(unit,*)
         enddo
         
         ! Read state numbers
         read(unit,'(56x,i1,5x,i1)') s,s1

         ! Skip to NACT value
         do i=1,6
            read(unit,*)
         enddo

         ! Read NACT
         read(unit,'(30x,F14.8)') nact(s,s1,itmp)

!------------------------------------------------------------------
! No need to do this: NACTS were calculated in units of Angstrom^-1
! using the DDR routine.
! 
! N.B. should really introduce a keyword to specify which units
!      have been used with the DDR routine...
!------------------------------------------------------------------
!         ! Convert to Angstrom^-1
!         nact(s,s1,itmp)=nact(s,s1,itmp)/0.529177249d0
!------------------------------------------------------------------

         ! Fill in the upper triangle of the NACT matrix
         nact(s1,s,itmp)=-nact(s,s1,itmp)

      endif
      goto 30

200   continue

!------------------------------------------------------------------
! Close file
!------------------------------------------------------------------
      close(unit)

      return
    end subroutine rdnact

!##################################################################

    subroutine rdhess(filename,nsta,ncoo,pospos,negneg,posneg,&
         negpos,maxdim,qcxcoo,e0,diff,pos,neg)

      implicit none
      
      integer*8                               :: nsta,ncoo,maxdim,&
                                                 unit,ngeom,i,j,k,&
                                                 ilbl,jlbl
      real*8, dimension(maxdim)               :: qcxcoo,tmpxcoo,&
                                                 en,e0
      real*8, dimension(ncoo,ncoo,nsta)       :: pospos,negneg,&
                                                 posneg,negpos
      real*8, dimension(ncoo,nsta)            :: pos,neg
      real*8                                  :: diff
      real*8, parameter                       :: tol=10E-7
      character(len=80)                       :: filename,string

!------------------------------------------------------------------
! Open file
!------------------------------------------------------------------
      unit=20
      open(unit,file=filename,form='formatted',status='old')

!------------------------------------------------------------------
! Read to the final geometry
!------------------------------------------------------------------
      ngeom=0
10    continue
      read(unit,'(a)',end=100) string
      if (string(2:19).eq.'ATOMIC COORDINATES') ngeom=ngeom+1
      goto 10

100   continue
      rewind(unit)

      k=0
20    continue
      read(unit,'(a)') string
      if (string(2:19).eq.'ATOMIC COORDINATES') k=k+1
      if (k.ne.ngeom) goto 20

!------------------------------------------------------------------
! Read the displaced geometry
!------------------------------------------------------------------
      ! Skip to the coordinates
      do i=1,3
         read(unit,*)
      enddo

      ! Read the coordinates
      do i=1,ncoo,3
         read(unit,'(18x,3(3x,F12.9))') (tmpxcoo(j), j=i,i+2)
      enddo

      ! Convert to Angstrom for comparison with the reference
      ! geometry
      do i=1,ncoo
         tmpxcoo(i)=tmpxcoo(i)*0.529177249d0
      enddo

!------------------------------------------------------------------
! Determine which DOFs have been displaced
!------------------------------------------------------------------
      ilbl=0
      jlbl=0
      do i=1,ncoo
         if (abs(tmpxcoo(i)-qcxcoo(i)).gt.tol) then
            if (ilbl.eq.0) then
               diff=abs(tmpxcoo(i)-qcxcoo(i))/0.529177249d0
               if (tmpxcoo(i)-qcxcoo(i).lt.0) then
                  ilbl=-i
               else
                  ilbl=i
               endif
            else
               if (tmpxcoo(i)-qcxcoo(i).lt.0) then
                  jlbl=-i
               else
                  jlbl=i
               endif
            endif
         endif
      enddo

!------------------------------------------------------------------
! Read the energies
!------------------------------------------------------------------
30    continue
      read(unit,'(a)') string
      if (string(3:13).ne.'MCSCF expec') then
         if (string(3:13).eq.'MCSCF STATE'&
              .and.string(19:24).eq.'Energy') then
            read(string(15:15),'(i1)') k
            read(string,'(36x,F18.12)') en(k)
            en(k)=(en(k)-e0(1))
         endif
         goto 30
      endif

!------------------------------------------------------------------
! Determine which array to save the energies to
!------------------------------------------------------------------
      ! pos,pos
      if (ilbl.gt.0.and.jlbl.gt.0) then
         do i=1,nsta
            pospos(ilbl,jlbl,i)=en(i)
         enddo
      endif

      ! negneg
      if (ilbl.lt.0.and.jlbl.lt.0) then
         do i=1,nsta
            negneg(abs(ilbl),abs(jlbl),i)=en(i)
         enddo
      endif

      ! posneg
      if (ilbl.gt.0.and.jlbl.lt.0) then
         do i=1,nsta
            posneg(ilbl,abs(jlbl),i)=en(i)
         enddo
      endif
      if (ilbl.lt.0.and.jlbl.gt.0) then
         do i=1,nsta
            posneg(jlbl,abs(ilbl),i)=en(i)
         enddo
      endif

      ! negpos
      if (ilbl.lt.0.and.jlbl.gt.0) then
         do i=1,nsta
            negpos(abs(ilbl),jlbl,i)=en(i)
         enddo
      endif
      if (ilbl.gt.0.and.jlbl.lt.0) then
         do i=1,nsta
            negpos(abs(jlbl),ilbl,i)=en(i)
         enddo
      endif

      ! pos
      if (jlbl.eq.0.and.ilbl.gt.0) then
         do i=1,nsta
            pos(ilbl,i)=en(i)
         enddo
      endif

      ! neg
      if (jlbl.eq.0.and.ilbl.lt.0) then
         do i=1,nsta
            neg(abs(ilbl),i)=en(i)
         enddo
      endif

!------------------------------------------------------------------
! Close file
!------------------------------------------------------------------
      close(unit)

      return
    end subroutine rdhess

!##################################################################

    subroutine rdquart(filename,nsta,ncoo,maxdim,qcxcoo,qdiff,&
         qdiff2,e0,pos,neg,pos2,neg2)

      implicit none

      integer*8                    :: nsta,ncoo,maxdim,unit,&
                                      ngeom,i,j,k,ilbl
      real*8, dimension(maxdim)    :: qcxcoo,tmpxcoo,en,e0
      real*8, dimension(ncoo,nsta) :: pos,neg,pos2,neg2
      real*8                       :: qdiff,qdiff2,ftmp
      real*8, parameter            :: tol=10E-7
      character(len=80)            :: filename,string
      logical(kind=4)              :: l1dx,l2dx

!------------------------------------------------------------------
! Open file
!------------------------------------------------------------------
      unit=20
      open(unit,file=filename,form='formatted',status='old')

!------------------------------------------------------------------
! Read to the final geometry
!------------------------------------------------------------------
      ngeom=0
10    continue
      read(unit,'(a)',end=100) string
      if (string(2:19).eq.'ATOMIC COORDINATES') ngeom=ngeom+1
      goto 10

100   continue
      rewind(unit)

      k=0
20    continue
      read(unit,'(a)') string
      if (string(2:19).eq.'ATOMIC COORDINATES') k=k+1
      if (k.ne.ngeom) goto 20

!------------------------------------------------------------------
! Read the displaced geometry
!------------------------------------------------------------------
      ! Skip to the coordinates
      do i=1,3
         read(unit,*)
      enddo
      
      ! Read the coordinates
      do i=1,ncoo,3
         read(unit,'(18x,3(3x,F12.9))') (tmpxcoo(j), j=i,i+2)
      enddo

      ! Convert to Angstrom for comparison with the reference
      ! geometry
      do i=1,ncoo
         tmpxcoo(i)=tmpxcoo(i)*0.529177249d0
      enddo

!------------------------------------------------------------------
! Determine which DOF has been displaced and which stepsize has
! been used
!------------------------------------------------------------------
      ilbl=0
      l1dx=.false.
      l2dx=.false.
      do i=1,ncoo
         ftmp=abs(tmpxcoo(i)-qcxcoo(i))
         if (ftmp.gt.tol) then            
            if (abs(ftmp/0.529177249d0-qdiff).lt.tol) then
               l1dx=.true.
               if (tmpxcoo(i)-qcxcoo(i).lt.0) then
                  ilbl=-i
               else
                  ilbl=i
               endif
            else if (abs(ftmp/0.529177249d0-qdiff2).lt.tol) then
               l2dx=.true.
               if (tmpxcoo(i)-qcxcoo(i).lt.0) then
                  ilbl=-i
               else
                  ilbl=i
               endif
            endif
         endif
      enddo

!------------------------------------------------------------------
! Read the energies
!------------------------------------------------------------------
30    continue
      read(unit,'(a)') string
      if (string(3:13).ne.'MCSCF expec') then
         if (string(3:13).eq.'MCSCF STATE'&
              .and.string(19:24).eq.'Energy') then
            read(string(15:15),'(i1)') k
            read(string,'(36x,F18.12)') en(k)
            en(k)=(en(k)-e0(1))
         endif
         goto 30
      endif 

!------------------------------------------------------------------
! Determine which array to save the energies to
!------------------------------------------------------------------
      if (l1dx) then
         ! x0 + dx
         if (ilbl.gt.0) then
            ! pos
            do i=1,nsta
               pos(ilbl,i)=en(i)
            enddo
         else
            ! neg
            do i=1,nsta
               neg(abs(ilbl),i)=en(i)
            enddo
         endif
      else if (l2dx) then
         ! x0 + 2dx
         if (ilbl.gt.0) then
            ! pos2
            do i=1,nsta
               pos2(ilbl,i)=en(i)
            enddo
         else
            ! neg2
            do i=1,nsta
               neg2(abs(ilbl),i)=en(i)
            enddo
         endif
      endif

!------------------------------------------------------------------
! Close file
!------------------------------------------------------------------
      close(unit)

      return
    end subroutine rdquart

!##################################################################

    subroutine rotate(grad,nact,hess,maxdim,ncoo,xcoo0,qcxcoo,&
         mass,nsta,quart,lgrad,lnact,lhess,lquart)

      implicit none

      integer*8                               :: maxdim,ncoo,i,j,&
                                                 k,ilbl,jlbl,m,n,&
                                                 nsta
      real*8, dimension(maxdim)               :: xcoo0,qcxcoo,&
                                                 mass
      real*8, dimension(maxdim,maxdim)        :: grad,quart,tmp2
      real*8, dimension(maxdim,maxdim,maxdim) :: nact,hess,tmp3
      real*8                                  :: tmass
      real*8, dimension(3)                    :: com
      real*8, dimension(ncoo)                 :: tmpxcoo,tmp
      real*8, dimension(ncoo,ncoo)            :: tmat
      real*8, dimension(3,3)                  :: itensor

      integer*4                               :: error,e2
      real*8, dimension(3)                    :: ieig
      real*8, dimension(9)                    :: work

      logical(kind=4)                         :: lgrad,lnact,&
                                                 lhess,lquart

!------------------------------------------------------------------
! Translate to the centre of mass
!------------------------------------------------------------------
      tmass=0.0d0
      com=0.0d0
      tmpxcoo=0.0d0

      do i=1,ncoo/3
         tmass=tmass+mass(i*3)
         do j=1,3
            com(j)=com(j)+qcxcoo(i*3-3+j)*mass(i*3)
         enddo
      enddo

      do i=1,3
         com(i)=com(i)/tmass
      enddo

      do i=1,ncoo/3
         do j=1,3
            tmpxcoo(i*3-3+j)=qcxcoo(i*3-3+j)-com(j)
         enddo
      enddo

!------------------------------------------------------------------
! Calculate the moment of inertia tensor
!------------------------------------------------------------------
      itensor=0.0d0

      do i=1,ncoo/3
         itensor(1,1)=itensor(1,1)+mass(i*3)*(tmpxcoo(i*3-1)**2+tmpxcoo(i*3)**2)
         itensor(2,2)=itensor(2,2)+mass(i*3)*(tmpxcoo(i*3-2)**2+tmpxcoo(i*3)**2)
         itensor(3,3)=itensor(3,3)+mass(i*3)*(tmpxcoo(i*3-2)**2+tmpxcoo(i*3-1)**2)
         itensor(1,2)=itensor(1,2)-mass(i*3)*(tmpxcoo(i*3-2)*tmpxcoo(i*3-1))
         itensor(1,3)=itensor(1,3)-mass(i*3)*(tmpxcoo(i*3-2)*tmpxcoo(i*3))
         itensor(2,3)=itensor(2,3)-mass(i*3)*(tmpxcoo(i*3-1)*tmpxcoo(i*3))
      enddo
      itensor(2,1)=itensor(1,2)
      itensor(3,1)=itensor(1,3)
      itensor(3,2)=itensor(2,3)
      
!------------------------------------------------------------------
! Diagonalise the moment of inertia tensor
!------------------------------------------------------------------
      e2=9
      call dsyev('V','U',3,itensor,3,ieig,work,e2,error)
      
      if (error.ne.0) then
         write(6,'(/,a,2/,a,x,i20,/)') &
        ' Diagonalisation of the moment of inertia tensor failed.'&
        ,' dsyev error:',error
         STOP
      endif

!------------------------------------------------------------------
! Create the transformation matrix
!------------------------------------------------------------------
      tmat=0.0d0
      
      do i=1,ncoo/3
         k=i*3-2
         ilbl=0
         do m=k,k+2
            ilbl=ilbl+1
            jlbl=0
            do n=k,k+2
               jlbl=jlbl+1
               tmat(m,n)=itensor(ilbl,jlbl)
            enddo
         enddo
      enddo

!------------------------------------------------------------------
! Rotate the coordinates using the eigenvectors of the moment of
! inertia tensor
!------------------------------------------------------------------
      do i=1,ncoo
         tmp(i)=0.0d0
         do j=1,ncoo
            tmp(i)=tmp(i)+tmat(j,i)*tmpxcoo(j)
         enddo
      enddo

      do i=1,ncoo
         qcxcoo(i)=tmp(i)
      enddo

!------------------------------------------------------------------
! Rotate the gradient vectors using the eigenvectors of the 
! moment of inertia tensor
!------------------------------------------------------------------
      if (lgrad) then
         do i=1,nsta
            do j=1,ncoo
               tmp2(i,j)=0.0d0
               do k=1,ncoo
                  tmp2(i,j)=tmp2(i,j)+tmat(k,j)*grad(i,k)
               enddo
            enddo
         enddo
         
         do i=1,nsta
            do j=1,ncoo
               grad(i,j)=tmp2(i,j)
            enddo
         enddo
      endif

!------------------------------------------------------------------
! Rotate the NACT vectors using the eigenvectors of the moment of 
! inertia tensor
!------------------------------------------------------------------
      if (lnact) then
         do i=1,nsta
            do j=1,nsta
               do k=1,ncoo
                  tmp3(i,j,k)=0.0d0
                  do m=1,ncoo
                     tmp3(i,j,k)=tmp3(i,j,k)+tmat(m,k)*nact(i,j,m)
                  enddo
               enddo
            enddo
         enddo
         
         do i=1,nsta
            do j=1,nsta
               do k=1,ncoo
                  nact(i,j,k)=tmp3(i,j,k)
               enddo
            enddo
         enddo
      endif

!------------------------------------------------------------------
! Rotate the Hessians using the eigenvectors of the moment of 
! inertia tensor
!------------------------------------------------------------------
      if (lhess) then
         do i=1,nsta
            do j=1,ncoo
               do k=1,ncoo
                  tmp3(i,j,k)=0.0d0
                  do m=1,ncoo
                     do n=1,ncoo
                        tmp3(i,j,k)=tmp3(i,j,k)&
                             +tmat(m,j)*hess(i,m,n)*tmat(n,k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         
         do i=1,nsta
            do j=1,ncoo
               do k=1,ncoo
                  hess(i,j,k)=tmp3(i,j,k)
               enddo
            enddo
         enddo
      endif

!------------------------------------------------------------------
! Rotate the quartic terms using the eigenvectors of the 
! moment of inertia tensor
!------------------------------------------------------------------
      if (lquart) then
         do i=1,nsta
            do j=1,ncoo
               tmp2(i,j)=0.0d0
               do k=1,ncoo
                  tmp2(i,j)=tmp2(i,j)+(tmat(k,j)**4)*quart(i,k)
               enddo
            enddo
         enddo
      
         do i=1,nsta
            do j=1,ncoo
               quart(i,j)=tmp2(i,j)
            enddo
         enddo
      endif

      return
    end subroutine rotate

!##################################################################

    subroutine transform(grad,nact,hess,maxdim,ncoo,nsta,pvec,&
         nzero,quart,lgrad,lnact,lhess,lquart)

      implicit none

      integer*8                               :: maxdim,ncoo,nsta,&
                                                 i,j,k,m,n,nzero
      real*8, dimension(maxdim,maxdim)        :: grad,pvec,tmp2,&
                                                 quart
      real*8, dimension(maxdim,maxdim,maxdim) :: nact,hess,tmp3
      logical(kind=4)                         :: lgrad,lnact,&
                                                 lhess,lquart
      write(6,'(/,3x,a)') &
           'Calculating the vibronic coupling coefficients...'

!------------------------------------------------------------------
! Transform the gradient vectors
!------------------------------------------------------------------
      if (lgrad) then
         do i=1,nsta
            do j=1,ncoo-nzero
               tmp2(i,j)=0.0d0
               do k=1,ncoo
                  tmp2(i,j)=tmp2(i,j)+pvec(j,k)*grad(i,k)
               enddo
            enddo
         enddo
         
         do i=1,nsta
            do j=1,ncoo-nzero
               grad(i,j)=tmp2(i,j)
            enddo
         enddo
      endif

!------------------------------------------------------------------
! Transform the NACT vectors
!------------------------------------------------------------------
      if (lnact) then
         do i=1,nsta
            do j=1,nsta
               do k=1,ncoo-nzero
                  tmp3(i,j,k)=0.0d0
                  do m=1,ncoo
                     tmp3(i,j,k)=tmp3(i,j,k)+pvec(m,k)*nact(i,j,m)
                  enddo
               enddo
            enddo
         enddo
         
         do i=1,nsta
            do j=1,nsta
               do k=1,ncoo-nzero
                  nact(i,j,k)=tmp3(i,j,k)
               enddo
            enddo
         enddo
      endif

!------------------------------------------------------------------
! Transform the Hessians
!------------------------------------------------------------------
      if (lhess) then
         do i=1,nsta
            do j=1,ncoo-nzero
               do k=1,ncoo-nzero
                  tmp3(i,j,k)=0.0d0
                  do m=1,ncoo
                     do n=1,ncoo
                        tmp3(i,j,k)=tmp3(i,j,k)&
                             +pvec(m,j)*hess(i,m,n)*pvec(n,k)
                     enddo
                  enddo
               enddo
            enddo
         enddo

         do i=1,nsta
            do j=1,ncoo-nzero
               do k=1,ncoo-nzero
                  hess(i,j,k)=tmp3(i,j,k)
               enddo
            enddo
         enddo
      endif

!------------------------------------------------------------------
! Transform the quartic terms
!------------------------------------------------------------------
      if (lquart) then
         do i=1,nsta
            do j=1,ncoo-nzero
               tmp2(i,j)=0.0d0
               do k=1,ncoo
                  tmp2(i,j)=tmp2(i,j)+(pvec(k,j)**4)*quart(i,k)
               enddo
            enddo
         enddo
         
         do i=1,nsta
            do j=1,ncoo-nzero
               quart(i,j)=tmp2(i,j)
            enddo
         enddo
      endif

      return
    end subroutine transform

!##################################################################

    subroutine outcoeffs(grad,nact,hess,e0,maxdim,ncoo,nsta,freq,&
         peig,iimag,nzero,aoper,enact,quart,lgrad,lnact,lhess,&
         lquart)

      implicit none

      integer*8                               :: maxdim,ncoo,nsta,&
                                                 nzero
      integer*8, dimension(maxdim)            :: iimag
      real*8, dimension(maxdim)               :: e0,freq,peig,enact
      real*8, dimension(maxdim,maxdim)        :: grad,quart
      real*8, dimension(maxdim,maxdim,maxdim) :: nact,hess
      character(len=80)                       :: aoper
      logical(kind=4)                         :: lgrad,lnact,&
                                                 lhess,lquart
      
!------------------------------------------------------------------
! Write the operator file
!------------------------------------------------------------------
      call wrop(grad,nact,hess,e0,maxdim,ncoo,nsta,freq,peig,&
           iimag,nzero,aoper,enact,quart,lgrad,lnact,lhess,lquart)

!------------------------------------------------------------------
! Write parameter data file
!------------------------------------------------------------------
      call wrdat(grad,nact,hess,e0,maxdim,ncoo,nsta,freq,peig,&
           iimag,nzero,aoper,enact,quart,lgrad,lnact,lhess,lquart)

      return
    end subroutine outcoeffs

!##################################################################

    subroutine rmblank(string,length)

      implicit none

      integer*4             :: length
      integer*8             :: i,k
      character(len=length) :: string,string1

      string1=''
      k=0
      do i=1,length
         if (string(i:i).ne.' ') then
            k=k+1
            write(string1(k:k),'(a1)') string(i:i)
         endif
      enddo
      string=string1

      return
    end subroutine rmblank

!##################################################################

    subroutine wrop(grad,nact,hess,e0,maxdim,ncoo,nsta,freq,peig,&
         iimag,nzero,aoper,enact,quart,lgrad,lnact,lhess,lquart)

      implicit none


      integer*8                               :: maxdim,ncoo,nsta,&
                                                 i,j,k,ilbl,unit,&
                                                 l,c,d,e,f,g,nzero
      integer*8, dimension(maxdim)            :: iimag
      real*8                                  :: ftmp
      real*8, parameter                       :: tol=0.005d0
      real*8, dimension(maxdim)               :: e0,freq,peig,enact
      real*8, dimension(maxdim,maxdim)        :: grad,quart
      real*8, dimension(maxdim,maxdim,maxdim) :: nact,hess
      character(len=80)                       :: atmp,aoper
      character(len=80)                       :: mout1,mout2,&
                                                 mout3,mout4,mout5
      logical(kind=4)                         :: lgrad,lnact,&
                                                 lhess,lquart

      real*8, dimension(nsta,ncoo) :: efreq

!------------------------------------------------------------------
! Effective-frequency weighted parameters
!------------------------------------------------------------------
      do i=1,nsta
         do j=1,ncoo-nzero
            ftmp=hess(i,j,j)-(freq(j)*0.6373641d0)
            efreq(i,j)=sqrt(abs(freq(j)*0.6373641d0*freq(j)*0.6373641d0 &
                 + freq(j)*0.6373641d0*ftmp))
         enddo
      enddo

!      do j=1,ncoo-nzero
!         do i=1,nsta      
!            if (abs(grad(i,j)).lt.tol) cycle
!            write(6,'(i2,x,i2,x,F7.4)'),j,i,grad(i,j)/efreq(i,j)
!         enddo
!      enddo


!      do k=1,ncoo-nzero
!         do i=1,nsta-1
!            do j=i+1,nsta
!
!               ftmp=nact(i,j,k)*(enact(j)-enact(i))*27.2116d0
!               if (abs(ftmp).lt.tol) cycle
!
!               write(6,'(3(i2,x),F7.4)'),k,i,j,ftmp/efreq(i,j)
!
!            enddo
!         enddo
!      enddo

!------------------------------------------------------------------
! Open file
!------------------------------------------------------------------
      write(6,'(/,3x,a,/)') 'Writing the operator file...'
      unit=20
      open(unit,file=aoper,form='formatted',status='unknown')

!------------------------------------------------------------------
! Write parameter section
!------------------------------------------------------------------
      write(unit,'(a)') 'parameter-section'

! VEEs
      write(unit,*)
      do i=1,nsta
         write(unit,'(a1,i1,1x,a1,1x,F10.7,1x,a4)') &
              'E',i,'=',(e0(i)-e0(1))*27.2116d0,', ev'
      enddo

! Frequencies
      write(unit,*)
      do i=1,ncoo-nzero
         atmp=''
         write(atmp,'(a6,i2)') 'omega_',i
         call rmblank(atmp,len(atmp))
         ilbl=len_trim(atmp)
         write(unit,'(2(a,x),F10.7,x,a)') &
              atmp(1:ilbl),'=',freq(i)*0.6373641d0,', ev'
      enddo

! Kappa values
      if (lgrad) then
         write(unit,*)
         do i=1,nsta
            do j=1,ncoo-nzero
               if (abs(grad(i,j)).lt.tol) cycle
               atmp=''
               write(atmp,'(a5,2(a1,i2))') &
                    'kappa','_',i,'_',j
               call rmblank(atmp,len(atmp))
               ilbl=len_trim(atmp)
               write(unit,'(2(a,x),F10.7,x,a)') &
                    atmp(1:ilbl),'=',grad(i,j),', ev'
            enddo
         enddo
      endif

! Lambda values
      if (lnact) then
         write(unit,*)
         do i=1,nsta-1
            do j=i+1,nsta
               do k=1,ncoo-nzero
                  ftmp=nact(i,j,k)*(enact(j)-enact(i))*27.2116d0
                  if (abs(ftmp).lt.tol) cycle
                  atmp=''
                  write(atmp,'(a6,3(a1,i2))') &
                       'lambda','_',i,'_',j,'_',k
                  call rmblank(atmp,len(atmp))
                  ilbl=len_trim(atmp)
                  write(unit,'(2(a,x),F10.7,x,a)') &
                       atmp(1:ilbl),'=',ftmp,', ev'
               enddo
            enddo
         enddo
      endif

! Gamma values
      if (lhess) then
         ! Quadratic
         write(unit,*)
         do i=1,nsta
            do j=1,ncoo-nzero
               !            ftmp=abs(hess(i,j,j))
               ftmp=hess(i,j,j)
               if (iimag(j).eq.1) then
                  ftmp=ftmp+(freq(j)*0.6373641d0)
               else
                  ftmp=ftmp-(freq(j)*0.6373641d0)
               endif
               if (abs(ftmp).lt.tol) cycle
               atmp=''
               write(atmp,'(a5,3(a1,i2))') &
                    'gamma','_',i,'_',j,'_',j
               call rmblank(atmp,len(atmp))
               ilbl=len_trim(atmp)
               write(unit,'(2(a,x),F10.7,x,a)') &
                    atmp(1:ilbl),'=',ftmp,', ev'
            enddo
         enddo

         ! Bilinear
         write(unit,*)
         do i=1,nsta
            do j=1,ncoo-1-nzero
               do k=j+1,ncoo-nzero
                  ftmp=hess(i,j,k)
                  if (abs(ftmp).lt.tol) cycle
                  atmp=''
                  write(atmp,'(a5,3(a1,i2))') &
                       'gamma','_',i,'_',j,'_',k
                  call rmblank(atmp,len(atmp))
                  ilbl=len_trim(atmp)
                  write(unit,'(2(a,x),F10.7,x,a)') &
                       atmp(1:ilbl),'=',ftmp,', ev'
               enddo
            enddo
         enddo
      endif

! Quartic terms
      if (lquart) then
         write(unit,*)
         do i=1,nsta
            do j=1,ncoo-nzero
               if (abs(quart(i,j)).lt.tol) cycle
               atmp=''
               write(atmp,'(a7,2(a1,i2))') &
                    'epsilon','_',i,'_',j
               call rmblank(atmp,len(atmp))
               ilbl=len_trim(atmp)
               write(unit,'(2(a,x),F10.7,x,a)') &
                    atmp(1:ilbl),'=',quart(i,j),', ev'
            enddo
         enddo
      endif

! End parameter section
      write(unit,'(a)') 'end-parameter-section'

!------------------------------------------------------------------
! Write Hamiltonian section
!------------------------------------------------------------------
      write(unit,'(/,a)') 'hamiltonian-section'

! Write mode labels
      write(unit,'(65a)') ('-', i=1,65)
      c=7
      d=7
      e=7
      f=7
      g=7
      mout1 = ' '
      mout2 = ' '
      mout3 = ' '
      mout4 = ' '
      mout5 = ' '

      l=ncoo+1
      do i=1,ncoo-nzero
         if (i .lt. 11) then
            write(mout1(1:6),112) 'modes|'
            if (i.lt. 10) then
               write(mout1(c:c+5), 110)  ' v',i,'  |'
            else
               write(mout1(c:c+6), 111)  ' v',i,' |'
            endif
            c=c+6
         else if (i .lt. 21 .and. i .gt. 10) then
            write(mout2(1:6),112) 'modes|'
            write(mout2(d:d+6), 111)  ' v',i,' |'
            d=d+6
         else if (i .lt. 31 .and. i .gt. 20) then
            write(mout3(1:6),112) 'modes|'
           write(mout3(e:e+6), 111)  ' v',i,' |'
           e=e+6
        else if (i .lt. 41 .and. i .gt. 30) then
           write(mout4(1:6),112) 'modes|'
           write(mout4(f:f+6), 111)  ' v',i,' |'
           f=f+6
        else if (i .lt. 51 .and. i .gt. 40) then
           write(mout5(1:6),112) 'modes|'
           write(mout5(g:g+6), 111)  ' v',i,' |'
           g=g+6
        endif

        if (i .eq. ncoo-nzero) then
           if (i .lt. 11) then
              write(mout1(c:c+3),84) ' el'
           else if ((i .lt. 21) .and. (i .gt. 10)) then
              write(mout2(d:d+3),84) ' el'
           else if ((i .lt. 31) .and. (i .gt. 20)) then
              write(mout3(e:e+3),84) ' el'
           else if ((i .lt. 41) .and. (i .gt. 30)) then
              write(mout4(f:f+3),84) ' el'
           else if ((i .lt. 51) .and. (i .gt. 40)) then
              write(mout5(g:g+3),84) ' el'
           endif
        endif
     enddo

     if (mout1.ne.' ') write(unit,20) mout1(1:len_trim(mout1))
     if (mout2.ne.' ') write(unit,20) mout2(1:len_trim(mout2))
     if (mout3.ne.' ') write(unit,20) mout3(1:len_trim(mout3))
     if (mout4.ne.' ') write(unit,20) mout4(1:len_trim(mout4))
     if (mout5.ne.' ') write(unit,20) mout5(1:len_trim(mout5))

     write(unit,'(65a)') ('-', i=1,65)

20   format(a)
110  format(a2,i1,a3)
111  format(a2,i2,a2)
112  format(a6)
84   format(a3)

! Kinetic energy
     write(unit,*)
     do i=1,ncoo-nzero
        if (i.lt.10) then
           write(unit,'(a11,i1,4x,a1,i1,2x,a4)') &
                '-0.5*omega_',i,'|',i,'dq^2'
        else
           write(unit,'(a11,i2,3x,a1,i2,1x,a4)') &
                '-0.5*omega_',i,'|',i,'dq^2'
        endif
     enddo

! VEEs
     write(unit,*)
     do i=1,nsta
        write(unit,'(a1,i1,4x,a1,i2,4x,a1,i1,a1,i1)') &
             'E',i,'|',ncoo-nzero+1,'S',i,'&',i
     enddo

! Zero-order potential
     write(unit,*)
     do i=1,ncoo-nzero
        if (i.lt.10) then
           write(unit,'(a10,i1,4x,a1,i1,3x,a3)') &
                '0.5*omega_',i,'|',i,'q^2'
        else
           write(unit,'(a10,i2,3x,a1,i2,2x,a3)') &
                '0.5*omega_',i,'|',i,'q^2'
        endif
     enddo

! Kappa terms
     if (lgrad) then
        write(unit,*)
        do i=1,nsta
           do j=1,ncoo-nzero
              if (abs(grad(i,j)).lt.tol) cycle
              if (j.lt.10) then
                 write(unit,'(a6,i1,a1,i1,4x,a1,i1,4x,a1,4x,a1,i2,2x,a1,i1,a1,i1)')&
                      'kappa_',i,'_',j,'|',j,'q','|',ncoo-nzero+1,'S',i,'&',i
              else
                 write(unit,'(a6,i1,a1,i2,3x,a1,i2,3x,a1,4x,a1,i2,2x,a1,i1,a1,i1)')&
                      'kappa_',i,'_',j,'|',j,'q','|',ncoo-nzero+1,'S',i,'&',i
              endif
           enddo
        enddo
     endif

! Lambda terms
     if (lnact) then
        write(unit,*)
        do i=1,nsta-1
           do j=i+1,nsta
              do k=1,ncoo-nzero
                 ftmp=nact(i,j,k)*(enact(j)-enact(i))*27.2116d0
                 if (abs(ftmp).lt.tol) cycle
                 if (k.lt.10) then
                    write(unit,'(a7,i1,a1,i1,a1,i1,4x,a1,i1,4x,a1,4x,a1,i2,2x,a1,i1,a1,i1)')&
                         'lambda_',i,'_',j,'_',k,'|',j,'q','|',ncoo-nzero+1,'S',i,'&',j
                 else
                    write(unit,'(a7,i1,a1,i1,a1,i2,3x,a1,i2,3x,a1,4x,a1,i2,2x,a1,i1,a1,i1)')&
                         'lambda_',i,'_',j,'_',k,'|',k,'q','|',ncoo-nzero+1,'S',i,'&',j
                 endif
              enddo
           enddo
        enddo
     endif

! Gamma terms
     if (lhess) then
        ! Quadratic terms
        write(unit,*)
        do i=1,nsta
           do j=1,ncoo-nzero
              ftmp=abs(hess(i,j,j))
              if (iimag(j).eq.1) then
                 ftmp=ftmp+(freq(j)*0.6373641d0)
              else
                 ftmp=ftmp-(freq(j)*0.6373641d0)
              endif
              if (abs(ftmp).lt.tol) cycle
              if (j.lt.10) then
                 write(unit,'(a10,i1,a1,i1,a1,i1,4x,a1,i1,4x,a3,4x,a1,i2,2x,a1,i1,a1,i1)') &
                      '0.5*gamma_',i,'_',j,'_',j,'|',j,'q^2','|',ncoo-nzero+1,'S',i,'&',i
              else
                 write(unit,'(a10,i1,a1,i2,a1,i2,2x,a1,i2,3x,a3,4x,a1,i2,2x,a1,i1,a1,i1)') &
                      '0.5*gamma_',i,'_',j,'_',j,'|',j,'q^2','|',ncoo-nzero+1,'S',i,'&',i
              endif
           enddo
        enddo
        
        ! Bilinear terms
        ! N.B. do NOT multiply these by 0.5 as we only use the pair
        ! alpha,beta and not beta,alpha (alpha < beta)
        write(unit,*)
        do i=1,nsta
           do j=1,ncoo-nzero-1
              do k=j+1,ncoo-nzero
                 ftmp=hess(i,j,k)
                 if (abs(ftmp).lt.tol) cycle
                 if (j.lt.10.and.k.lt.10) then
                    write(unit,'(a6,i1,a1,i1,a1,i1,4x,a1,i1,4x,a1,4x,a1,i1,4x,a1,4x,a1,i2,2x,a1,i1,a1,i1)') &
                         'gamma_',i,'_',j,'_',k,'|',j,'q','|',k,'q','|',ncoo-nzero+1,'S',i,'&',i
                 else if (j.lt.10.and.k.ge.10) then
                    write(unit,'(a6,i1,a1,i1,a1,i2,3x,a1,i1,4x,a1,4x,a1,i2,3x,a1,4x,a1,i2,2x,a1,i1,a1,i1)') &
                         'gamma_',i,'_',j,'_',k,'|',j,'q','|',k,'q','|',ncoo-nzero+1,'S',i,'&',i
                 else ! j>10, k>10
                    write(unit,'(a6,i1,a1,i2,a1,i2,2x,a1,i2,3x,a1,4x,a1,i2,3x,a1,4x,a1,i2,2x,a1,i1,a1,i1)') &
                         'gamma_',i,'_',j,'_',k,'|',j,'q','|',k,'q','|',ncoo-nzero+1,'S',i,'&',i
                 endif
              enddo
           enddo
        enddo
     endif

! Quartic terms
     if (lquart) then
        write(unit,*)
        do i=1,nsta
           do j=1,ncoo-nzero
              if (abs(quart(i,j)).lt.tol) cycle
              if (j.lt.10) then
                 write(unit,'(a16,i1,a1,i1,4x,a1,i1,4x,a3,4x,a1,i2,2x,a1,i1,a1,i1)')&
                      '0.04167*epsilon_',i,'_',j,'|',j,'q^4','|',ncoo-nzero+1,'S',i,'&',i
              else
                 write(unit,'(a16,i1,a1,i2,3x,a1,i2,3x,a3,4x,a1,i2,2x,a1,i1,a1,i1)')&
                      '0.04167*epsilon_',i,'_',j,'|',j,'q^4','|',ncoo-nzero+1,'S',i,'&',i
              endif
           enddo
        enddo
     endif

! End Hamiltonian section
     write(unit,'(/,a)') 'end-hamiltonian-section'

!------------------------------------------------------------------
! Write excitation operators
!------------------------------------------------------------------
   do j=2,nsta
      write(unit,'(/,a26,i1)') 'hamiltonian-section_excite',j

      ! Write mode labels
      write(unit,'(65a)') ('-', i=1,65)
      c=7
      d=7
      e=7
      f=7
      g=7
      mout1 = ' '
      mout2 = ' '
      mout3 = ' '
      mout4 = ' '
      mout5 = ' '

      l=ncoo+1
      do i=1,ncoo-nzero
         if (i .lt. 11) then
            write(mout1(1:6),112) 'modes|'
            if (i.lt. 10) then
               write(mout1(c:c+5), 110)  ' v',i,'  |'
            else
               write(mout1(c:c+6), 111)  ' v',i,' |'
            endif
            c=c+6
         else if (i .lt. 21 .and. i .gt. 10) then
            write(mout2(1:6),112) 'modes|'
            write(mout2(d:d+6), 111)  ' v',i,' |'
            d=d+6
         else if (i .lt. 31 .and. i .gt. 20) then
            write(mout3(1:6),112) 'modes|'
           write(mout3(e:e+6), 111)  ' v',i,' |'
           e=e+6
        else if (i .lt. 41 .and. i .gt. 30) then
           write(mout4(1:6),112) 'modes|'
           write(mout4(f:f+6), 111)  ' v',i,' |'
           f=f+6
        else if (i .lt. 51 .and. i .gt. 40) then
           write(mout5(1:6),112) 'modes|'
           write(mout5(g:g+6), 111)  ' v',i,' |'
           g=g+6
        endif

        if (i .eq. ncoo-nzero) then
           if (i .lt. 11) then
              write(mout1(c:c+3),84) ' el'
           else if ((i .lt. 21) .and. (i .gt. 10)) then
              write(mout2(d:d+3),84) ' el'
           else if ((i .lt. 31) .and. (i .gt. 20)) then
              write(mout3(e:e+3),84) ' el'
           else if ((i .lt. 41) .and. (i .gt. 30)) then
              write(mout4(f:f+3),84) ' el'
           else if ((i .lt. 51) .and. (i .gt. 40)) then
              write(mout5(g:g+3),84) ' el'
           endif
        endif
     enddo

     if (mout1.ne.' ') write(unit,20) mout1(1:len_trim(mout1))
     if (mout2.ne.' ') write(unit,20) mout2(1:len_trim(mout2))
     if (mout3.ne.' ') write(unit,20) mout3(1:len_trim(mout3))
     if (mout4.ne.' ') write(unit,20) mout4(1:len_trim(mout4))
     if (mout5.ne.' ') write(unit,20) mout5(1:len_trim(mout5))

     write(unit,'(65a)') ('-', i=1,65)

     ! Write operator
     write(unit,'(a8,i2,4x,a3,i1)') &
          '1.0    |',ncoo-nzero+1,'S1&',j

     write(unit,'(a)') 'end-hamiltonian-section'
   enddo

!------------------------------------------------------------------
! End operator
!------------------------------------------------------------------
   write(unit,'(/,a)')'end-operator'

!------------------------------------------------------------------
! Close file
!------------------------------------------------------------------
      close(unit)

      return
    end subroutine wrop

!##################################################################

    subroutine wrdat(grad,nact,hess,e0,maxdim,ncoo,nsta,freq,peig,&
         iimag,nzero,aoper,enact,quart,lgrad,lnact,lhess,lquart)

      implicit none

      integer*8                               :: maxdim,ncoo,nsta,&
                                                 nzero,i,j,k,unit
      integer*8, dimension(maxdim)            :: iimag
      real*8, dimension(maxdim)               :: e0,freq,peig,enact
      real*8, dimension(maxdim,maxdim)        :: grad,quart
      real*8, dimension(maxdim,maxdim,maxdim) :: nact,hess
      real*8                                  :: ftmp
      character(len=80)                       :: aoper,filename
      logical(kind=4)                         :: lgrad,lnact,&
                                                 lhess,lquart

!------------------------------------------------------------------
! Open file
!------------------------------------------------------------------
      write(6,'(3x,a,/)') 'Writing the data file...'
      k=index(aoper,'.')
      filename=''
      write(filename(1:k),'(a)') aoper(1:k)
      write(filename(k+1:k+3),'(a)') 'dat'

      unit=20
      open(unit,file=filename,form='unformatted',status='unknown')

!------------------------------------------------------------------
! Output system information
!------------------------------------------------------------------
      write(unit) nsta
      write(unit) ncoo-nzero

!------------------------------------------------------------------
! Ouput parameter values
!------------------------------------------------------------------
      ! VEEs
      do i=1,nsta
         write(unit) (enact(i)-enact(1))*27.211d0
      enddo

      ! Frequencies
      do i=1,ncoo-nzero
         write(unit) freq(i)*0.6373641d0
      enddo

      ! kappa
      do i=1,nsta
         do j=1,ncoo-nzero
            write(unit) grad(i,j)
         enddo
      enddo

      ! lambda
      do i=1,nsta
         do j=1,nsta
            do k=1,ncoo-nzero
               write(unit) &
                    nact(i,j,k)*(enact(j)-enact(i))*27.2116d0
            enddo
         enddo
      enddo

      ! gamma
      do i=1,nsta
         do j=1,ncoo-nzero
            do k=1,ncoo-nzero
               if (j.eq.k) then
                  ! bilinear
                  ftmp=hess(i,j,j)
                  if (iimag(j).eq.1) then
                     ftmp=ftmp+(freq(j)*0.6373641d0)
                  else
                     ftmp=ftmp-(freq(j)*0.6373641d0)
                  endif
               else
                  ! quadratic
                  ftmp=hess(i,j,k)
               endif
               write(unit) ftmp
            enddo
         enddo
      enddo

!------------------------------------------------------------------
! Close file
!------------------------------------------------------------------
      close(unit)

      return
    end subroutine wrdat

!##################################################################
  
  end program vibcp
