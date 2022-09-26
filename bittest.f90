!!bit test sul set di base scelto
program bittest
  implicit none 
  integer:: i,n, max, min,nsiti,nso,nelettroni, count, a,b,count2,prod, d, e, countd, counte, nfunctions, config,  check, g, h, l, p, coppie
  character:: array(26)
  double precision::c, sz, spin
  logical::bool
  open(1,file='bittest.dat')
  open(2,file='dim.dat')
  open(3, file='configurazioni.dat')
  open(4, file='basis.dat')
  open(5,file='numero_doppie.dat')
  read(5,*) coppie
  close(5)
  write(*,*)'dammi il numero di siti'

  read(*,*) nsiti
  write(*,*) 'dammi lo spazio di spin'
  
  read(*,*) spin
 ! write(*,*) 'dammi il numero di coppie'
  !read(*,*) coppie
  nso=nsiti*2
  nelettroni=nsiti+1

  read(4,*) max
  read(4,*)min
  close(4)
 
  nfunctions=0
  do n=min,max
     count=0
     a=0; b=0; d=0; e=0
     count2=0
     countd=0
     counte=0
     prod=1
     sz=0.d0
     config=0
!!$     g=0
!!$     l=0
     h=0
     p=0
     check=0

     do i=0,nso-1
        bool=btest(n,i)
        ! write(*,*)i, bool
        count2=count2+1
        if(bool)then
           array(i+1)='1'
           count=count+1
           if(i/2*2.eq.i)then
              a=a+1
              !l=l+1
              countd=countd+1
           endif

           if(i/2*2.ne.i)then
              b=b+1
            !  g=g+1
              counte=counte+1
           endif
          
         
       
           ! write(1,*)'1'
        else
           array(i+1)='0'
           !write(1,*)'0'
           p=0
        endif
        prod=prod*(a+b)
       ! if(prod.eq.0) write(*,*) prod
        h=a+b
        if(h.eq.2)then
           check=check+1
         !  write(*,*) check
        endif
        if(count2.eq.2)then
           a=0
           b=0
           d=0
           e=0
!!$           l=0
!!$           g=0
           
           p=0
           count2=0
        endif
     enddo


     sz=(countd-counte)*0.5d0
     if((count.eq.nelettroni).and.((prod.eq.0).or.(prod.eq.2)).and.(check.le.coppie).and.(dabs(sz-spin).lt.1d-8))then
        write(1,*)(array(i),i=nso,1,-1)
        nfunctions=nfunctions+1

        do i=0,nso-1
           if(array(i+1).eq.'1')then
              config=config+2**i
           endif
        enddo
        write(3,*) config
     endif
     !write(*,*) sz

  enddo
  write(*,*)'numero di funzioni di base uguale', nfunctions
  write(2,*) nfunctions
endprogram bittest
