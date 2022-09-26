program oscillator
implicit none
real*8::strenght,hbar,hbarev,echarge,emass
real*8,allocatable::energy(:),mux(:),muy(:)
integer::i, nev, aa, bb
open(1,file="momenti_di_dipolo.dat")
open(2,file="eigenvalues.dat")

nev=4


hbar=1.0545718d-34 ![J*s]
hbarev=6.5821195d-16 ![eV]
echarge=1.6021766d-19 ![C]
emass=9.1093837d-32 ![kg]

allocate(energy(nev),mux(nev),muy(nev))


mux=0d0
muy=0d0
energy=0d0

do i=1,nev
   read(1,*) bb, mux(i), muy(i)
   !write(*,*) mux(i), muy(i)
   read(2,*) aa, energy(i)
  !write(*,*) energy(i)
enddo

strenght=0.d0

do i=2,nev
   strenght=strenght+(4.*emass)/((echarge**2.)*hbar*3.)*((energy(i)-energy(1))/hbarev)*((mux(i)**2+muy(i)**2)*10.d-20*echarge**2)
enddo

write(*,*) 'forza di oscillatore:',strenght
write(*,*) energy(nev)-energy(1)

endprogram oscillator
