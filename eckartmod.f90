  module eckartmod

    implicit none
    
    contains

!##################################################################

      subroutine ecktransmat(xcoo,ncoo,mass,tmat)

        implicit none
        
        integer*8                    :: ncoo,i,j,k,m,n,ilbl,jlbl
        real*8, dimension(ncoo)      :: xcoo,mass,tmpxcoo
        real*8, dimension(3)         :: com
        real*8                       :: tmass
        real*8, dimension(3,3)       :: itensor
        integer*4                    :: error,e2
        real*8, dimension(3)         :: ieig
        real*8, dimension(9)         :: work
        real*8, dimension(ncoo,ncoo) :: tmat

!------------------------------------------------------------------
! Translate to the centre of mass
!------------------------------------------------------------------
      tmass=0.0d0
      com=0.0d0
      tmpxcoo=0.0d0

      do i=1,ncoo/3
         tmass=tmass+mass(i*3)
         do j=1,3
            com(j)=com(j)+xcoo(i*3-3+j)*mass(i*3)
         enddo
      enddo

      do i=1,3
         com(i)=com(i)/tmass
      enddo

      do i=1,ncoo/3
         do j=1,3
            tmpxcoo(i*3-3+j)=xcoo(i*3-3+j)-com(j)
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

        return
      end subroutine ecktransmat

!##################################################################

    end module eckartmod
