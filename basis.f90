!!!! trovare configurazione massima e minima in bit representation
program basis
  implicit none
  integer:: nsiti,i,max,min, n, max1,max2, min1, min2,nso
  double precision::o
  write(*,*) 'doppie occupazioni'
  read(*,*) n
  open(1,file='basis.dat')
  open(2,file='numero_doppie.dat')
  nsiti=14 !siti di disposizione elettroni
  nso=14*2
  max2=0
  do i=2*n,nsiti-1
     max2=max2+2**i
  enddo
  max1=0
  do i=nso-1,nso-2*n,-1
     max1=max1+2**i
  enddo
  max=max1+max2
  write(1,*)max
  write(*,*) max
  min1=0
  do i=0,nsiti-1
     min1=min1+2**i
     if(i.ne.(i/2)*2)stop
  enddo
  min2=0

  min=min1+min2
  write(1,*)min
  write(*,*)min

  write(2,*) n
endprogram
