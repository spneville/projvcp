!##################################################################
! I/O routines taken from MCTDH
!##################################################################

  module iomod

    implicit none

  contains

!******************************************************************
!
!        SUBROUTINE RDINPF
!
!     read input file : reads a line of the input-file and returns
!                       the separate words in character strings
!
! on input:
! iz:  no. of lines previously read from input file
! iin: channel of file to be read
! maxkey: maximum no. of keywords on line. (Constant)
! ierr = 0  : Error flag.
!
! on return:
! ic: no. of words on line 
! lc: length of words
! iz: line no. in input file read
! keyword: separate keywords read from line inptit
!  IF:
!   ierr=0 no error
!   inptit: line read (line no. iz in file)
!  ELSE IF:
!   ierr=1 error in reading iin file. 
!   ierr=2 end of file found
!   inptit: Error message
!
!    03/97  HDM
!
! If you change maxkeylen, also change it in the include files:
! - global.f90
! - vcglobal.f90
!
!******************************************************************

    subroutine rdinpf(iin,keyword,keyorig,inptit,lc,ic,iz,ierr,&
                      maxkey)

      implicit none
      
      integer*8                                   :: maxkey,maxeof
      integer*8, parameter                        :: maxbuf=360
      integer*8, parameter                        :: c3=64
      integer*8, parameter                        :: c5=240
      integer*8, parameter                        :: maxkeylen=200
      integer*8                                   :: i,j,k,i0,iz,&
                                                     ic,laenge,&
                                                     iin,ierr
      integer*8, dimension(maxkey)                :: lc
      character(len=maxbuf)                       :: buffer,ctmp
      character(len=c5)                           :: inptit
      character(len=maxkeylen), dimension(maxkey) :: keyword,&
                                                     keyorig
      character(len=4)                            :: rest
      data               maxeof /0/

      ierr=0
      ic = 0
      do  i = 1, maxkey
         lc(i) = 0
         keyword(i) = '  '
         keyorig(i) = '  ' 
      end do
      do i=1,maxbuf
         buffer(i:i)= ' '
      enddo

!------------------------------------------------------------------
! Read one line of the input-file. Ignore lines which begin
! with '#' or '-----' or '_____' .
!------------------------------------------------------------------
30    read(iin,'(2a)',end=11,err=12) inptit,rest
      iz = iz + 1
      if( inptit(1:1) .eq. '#' )     goto 30
      if( inptit(1:2) .eq. '/*' )    goto 30
      if( inptit(1:5) .eq. '-----' ) goto 30
      if( inptit(1:5) .eq. '_____' ) goto 30

!------------------------------------------------------------------
! Remove tabs and blanks at the beginning. ( ichar(\tab)=9 )
! Skip lines that start with 32 blanks.
!------------------------------------------------------------------
      do k=1,32
         if(inptit(k:k).ne.' ' .and. ichar(inptit(k:k)).ne.9) then
            if(inptit(k:k)   .eq. '#')      goto 30
            if(inptit(k:k+1) .eq. '/*')     goto 30
            if(inptit(k:k+4) .eq. '-----' ) goto 30
            if(inptit(k:k+4) .eq. '_____' ) goto 30
            buffer(1:1) = inptit(k:k)
            i0 = k + 1
            goto 15
         endif
      enddo
      goto  30                  ! Ignore lines that start with 16 blanks

!------------------------------------------------------------------
! The four last characters must be blank. Otherwise input-line too
! long.
!------------------------------------------------------------------
 15   If( rest .ne. '    ' ) then
         write(6,'(a)') inptit
         write(6,'(2a)') ' rest : ', rest
         ierr=1
         inptit = '****  Input line too long! *****   line no.:   '
         write(inptit(46:50),'(i5)') iz
      endif

!------------------------------------------------------------------
! transfer input line to buffer, with modifications
! (i)   remove semicolons and brackets so they are separators. 
! (ii)  ensure blanks are either side of equals or comma signs so
!       they are keywords. 
! (iii) add end marker '#'
!------------------------------------------------------------------
      j=1
      do i = i0,c5
         if ( j+3 .gt. maxbuf ) goto 16
         j=j+1
         if (inptit(i-1:i) .eq. '  ' ) then
            j=j-1               ! Shrink multiple blanks to a 
                                ! single blank. 
         else if (inptit(i:i) .eq. ';') then
            buffer(j:j) = ' '
         else if (ichar(inptit(i:i)) .eq. 9) then  ! Tab->blank
            buffer(j:j) = ' '
         else if (inptit(i:i) .eq. '(') then
            buffer(j:j) = ' '
         else if (inptit(i:i) .eq. ')') then
            buffer(j:j) = ' '
         else if (inptit(i:i) .eq. '=') then
            buffer(j:j+2) = ' = '
            j=j+2
         else if (inptit(i:i) .eq. '<') then
            buffer(j:j+2) = ' < '
            j=j+2
         else if (inptit(i:i) .eq. '>') then
            buffer(j:j+2) = ' > '
            j=j+2
         else if (inptit(i:i) .eq. ',') then
            buffer(j:j+2) = ' , '
            j=j+2
         else if (inptit(i:i) .eq. '[') then
            buffer(j:j+2) = ' [ '
            j=j+2
         else if (inptit(i:i) .eq. ']') then
            buffer(j:j+2) = ' ] '
            j=j+2
         else if (inptit(i:i) .eq. '{') then
            buffer(j:j+2) = ' { '
            j=j+2
         else if (inptit(i:i) .eq. '}') then
            buffer(j:j+2) = ' } '
            j=j+2
         else if (inptit(i:i) .eq. '|') then
            buffer(j:j+1) = ' |'
            j=j+1
         else if (inptit(i:i) .eq. '#') then
            buffer(j:j+1) = ' #'
            goto 33
         else
            buffer(j:j) = inptit(i:i)
         endif
      end do

 16   if (j+3 .gt. maxbuf) then
         ierr=1
         inptit = '*****  Buffer overflow! *****   line no.:'
         write(inptit(44:48),'(i5)') iz
         return
      endif
      buffer(j+1:j+2) = ' #'

!------------------------------------------------------------------
! remove leading blank and check if at end of line 
!------------------------------------------------------------------
 33   continue
      if (buffer(1:1) .eq. '#') then
         if (ic .eq. 0) then
            goto 30
         else
            return
         endif
      else if (buffer(1:1) .eq. ' ' ) then 
         ctmp = buffer(2:maxbuf)
         buffer=ctmp
         goto 33  
      else if (buffer(1:1) .eq. '!' ) then 
         laenge = index(buffer,' ')
         ctmp = buffer(laenge+1:maxbuf)
         buffer=ctmp
         goto 33  
      endif
      
!------------------------------------------------------------------
! store the next word into keyword(ic)
!------------------------------------------------------------------
      laenge = index(buffer,' ')-1
      ic = ic + 1
      if (ic .ge. maxkey ) then
         ierr=1
         inptit = '*****  Too many keywords! *****   line no.:'
         write(inptit(44:48),'(i5)') iz
         return
      endif
      if (laenge .gt. maxkeylen ) then
         ierr=1
         inptit = '*****  keyword too long! *****   line no.:'
         write(inptit(43:47),'(i5)') iz
         inptit=inptit(1:47)//' '//buffer
         return
      endif

      if(buffer(1:1) .eq. '!' ) then
         ic = ic -1
      else
         lc(ic) = laenge
         keyword(ic) = buffer(1:laenge)
         keyorig(ic) = buffer(1:laenge)
         call keylowc(keyword(ic))
      end if
      ctmp = buffer(laenge+2:maxbuf)
      buffer=ctmp
      goto 33

 11   ierr=1
      inptit = '*****  End of file found! *****   last line no.:'
      write(inptit(49:53),'(i5)') iz
      maxeof = maxeof +1
      if (maxeof .gt. 9 ) goto 20
      return

 12   ierr=1
      inptit = '*****  Error reading file! *****   last line no.:'
      write(inptit(50:54),'(i5)') iz
      maxeof = maxeof +1
      if (maxeof .gt. 9 ) goto 20
      return

 20   write(6,*) 
      write(6,'(a)') &
           '################################################'
      write(6,'(a,i3)') ' Too many errors while reading channel:',&
           iin
      write(6,'(a)')    ' Program is canceled by "rdinpf" '
      write(6,'(a)') &
           '################################################'
      stop 1

    end subroutine rdinpf

!##################################################################

!******************************************************************
! 
!          Subroutine KEYLOWC
!
!   Replaces upper case letters [A-Z] by lower case letters [a-z].
!   All other symbols remain unchanged.
!
!   11/97 AJ
!
!******************************************************************
      
    subroutine keylowc (string)
      
      implicit none
      
      integer*8     ::  i,j,length
      character(*) :: string
      
      length=len(string)
      
      do i=1,length
         do j=65,90
            if (ichar(string(i:i)).eq.j) string(i:i)=char(j+32)
         enddo
      enddo
      
      return
    end subroutine keylowc

!##################################################################
    
  end module iomod
