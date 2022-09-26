program eccitoni
  implicit none
  integer,allocatable::numero(:),config(:)
  real*8,allocatable:: peso(:)
  character::array(26)
  integer::n, i, l, s, b,c,dim,kk,k
  real*8::a, ecc, pesomax
  logical::bool, bool1
!!$  write(*,*) 'dammi n'
!!$  read(*,*) n
  n=10000
  pesomax=0.0006
  allocate(numero(n),peso(n))
  open(1,file='eigenvector_6_sz0_7.dat')
  open(2,file='output.dat')
  open(3,file='configurazioni.dat')
  open(4,file='dim.dat')
  read(4,*) dim
  allocate(config(dim))
  do i=1,n
     read(1,*,iostat=kk) l,a
     if(is_iostat_end(kk))then
        k=i-1
        exit
     else
        k=i
     endif
    ! write(*,*) a
     numero(i)=l
     peso(i)=a**2
  enddo
  do i=1,dim
     read(3,*,iostat=kk) l
     if(is_iostat_end(kk))then
        k=i-1
        exit
     else
        k=i
     endif
     config(i)=l
  enddo
  

!!$  do i=1,n
!!$     write(2,*) funzioni(i), peso(i)
!!$  enddo
!!$
  write(2,*) 'configurazione       ','peso     ', 'eccitoni    '
  do i=1,n
     l=config(numero(i))
     ecc=0
     if(peso(i).ge.pesomax)then
        write(55,*) '------------------'
        write(55,*) l, peso(i)
     endif
     do s=0,24,2
        bool=btest(l,s)
        bool1=btest(l,s+1)
        if(bool)then
           b=1
           array(s+1)='1'
        else
           b=0
           array(s+1)='0'
        endif
        if(bool1)then
           c=1
           array(s+1)='1'
        else
           c=0
           array(s+1)='0'
        endif
        if(b+c.eq.2)then
           ecc=ecc+1
           if(peso(i).ge.pesomax)then
              write(55,*) '-',(s+2)/2
           endif
        endif
        if(b+c.eq.0)then
           ecc=ecc+1
           if(peso(i).ge.pesomax)then
              write(55,*) '+',(s+2)/2
           endif
        endif
     enddo
     if(peso(i).ge.pesomax)then
        !  write(2,*)(array(s),s=25,1,-1)
        
        write(2,*) config(numero(i)),peso(i),ecc
      !  write(2,*) '---------------------------------'
     endif
     
  enddo
     
endprogram eccitoni
