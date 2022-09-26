program spettro
implicit none
real*8, allocatable:: autovalori(:), mux(:), muy(:)
real*8:: energia, mu_x, mu_y, intensita, pi, omegastart, omegaend, domega, omega, sigma
integer:: i, nev, a, j, nomega
pi =dacos(-1.d0)
write(*,*) 'numero autovalori'
read(*,*) nev
open(1,file='eigenvalues_sz0_7.dat')
open(2,file='momenti_di_dipolo_sz0_7.dat')
open(3,file='spettro.dat')
allocate(autovalori(nev), mux(nev), muy(nev))
do i=1,nev
   read(1,*) a, energia
   autovalori(i)=energia
enddo
close(1)
sigma=0.1d0
do i=1,nev
   read(2,*) a, mu_x, mu_y
   mux(i)=mu_x
   muy(i)=mu_y
enddo
close(2)
nomega=1000000
omegastart=0.d0
omegaend=8.0d0
domega=(omegaend-omegastart)/nomega
omega=omegastart-domega
do j=1,nomega
   omega=omega+domega
   intensita=0.d0
   do i= 2, nev
      intensita=intensita+(autovalori(i)-autovalori(1)) * (mux(i)**2+muy(i)**2) *((2*pi)**(-0.5)*sigma)* dexp(-(omega-(autovalori(i)-autovalori(1)))**2/(2*sigma**2))
   enddo
   write(3,*) omega, intensita
enddo
end program spettro
