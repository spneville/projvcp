  module rddatamod

    implicit none
    
    contains

!##################################################################

      subroutine rddat(adat,ncoo,pvec,peig,iimag,xcoo0,maxdim,mass)

      implicit none

      integer*8                        :: ncoo,i,j,maxdim,iin
      integer*8, dimension(maxdim)     :: iimag
      real*8, dimension(maxdim)        :: peig,xcoo0,mass
      real*8, dimension(maxdim,maxdim) :: pvec
      character(len=80)                :: adat

!------------------------------------------------------------------
! Open the data file
!------------------------------------------------------------------
      iin=20
      open(iin,file=adat,form='unformatted',status='old')

!------------------------------------------------------------------
! Read data
!------------------------------------------------------------------
! System size
      read(iin) ncoo

! Reference coordinates (Cartesian)
      xcoo0=0.0d0
      do i=1,ncoo
         read(iin) xcoo0(i)
      enddo

! Masses
      do i=1,ncoo
         read(iin) mass(i)
      enddo

! Projected eigenvectors
      do i=1,ncoo
         do j=1,ncoo
            read(iin) pvec(i,j)
         enddo
      enddo

! Projected eigenvalues
      do i=1,ncoo
         read(iin) peig(i)
      enddo

! Imaginary frequency information
      do i=1,ncoo
         read(iin) iimag(i)
      enddo

! Close data file
      close(iin)

      return
    end subroutine rddat 

!##################################################################
  end module rddatamod
