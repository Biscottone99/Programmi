program scala
integer::i,j,k,nev
real*8::a,b,c, l, m
real*8,allocatable::sz0(:), sz1(:), fin(:), S(:), t(:)
character*1,allocatable::state(:)


open(1,file='eigenvalues_sz0_6.dat')
open(2,file='eigenvalues_sz1_6.dat')
open(3,file='stati.dat')
!open(7,file='tab.dat')

 
nev=30

allocate(sz0(nev), sz1(nev), fin(nev), state(nev), S(5), t(5))

do i=1,nev
   read(1,*) k, a
   sz0(i)=a
   read(2,*) k, b
   sz1(i)=b
enddo
close(1)
close(2)
state='S'
do i=1,nev

   fin(i)=sz0(i)-sz0(1)
   do j=1,nev
      c=dabs(sz0(i)-sz1(j))
      if(c.le.1d-10)then
         state(i)='T'
         
         goto 10
      endif
   enddo
10 continue
enddo
do i=1,nev
   write(3,*) i, fin(i), state(i)
enddo

open(4,file='aa.dat')
open(5,file='bb.dat')

do i=1,20
   if(state(i).eq.'S')then
      write(4,*) fin(i)
   else
      write(5,*) fin(i)
   endif
enddo
close(4)
close(5)
open(4,file='aa.dat')
open(5,file='bb.dat')
do i=1,5
   !read(4,*) s(i)
   read(5,*) t(i)
  
enddo
do i=1,4
   read(4,*) s(i)
  ! read(5,*) t(i)
  
enddo
write(3,*) 'S0-S1=',s(2)-s(1),'ST=', s(2)-t(1)
!write(7,*) 'ST=', s(2)-t(1)
close(3)
end program
