program transizioni
implicit none
real*8, allocatable:: autovalori(:), mux(:), muy(:)
real*8:: energia, mu_x, mu_y, intensita, pi, omegastart, omegaend, domega, omega, sigma
integer:: i, nev, a, j, nomega
pi =dacos(-1.d0)
!write(*,*) 'numero autovalori'
nev=30
open(1,file='eigenvalues_sz0_6.dat')
open(2,file='momenti_di_dipolo_sz0_6.dat')
open(3,file='transizioni.dat')
allocate(autovalori(nev), mux(nev), muy(nev))
do i=1,nev
   read(1,*) a, energia
   autovalori(i)=energia
   read(2,*) a, mu_x, mu_y
   mux(i)=mu_x
   muy(i)=mu_y
enddo
close(1)
close(2)

do i=2,nev
  ! if(i==2)write(3,*)i,autovalori(i)-autovalori(1)
   if((dabs(mux(i)).ge.1d-6).or.(dabs(muy(i)).ge.1d-6))then
      write(3,*) i, autovalori(i)-autovalori(1), (autovalori(i)-autovalori(1))*1239.813383
   endif
enddo
end program transizioni
