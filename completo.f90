program trianguleni
  !Questo programma permette di diagonalizzare completamente i trianguleni, dando in input solamente il numero di eccitoni permessi e lo spazio di spin
  implicit none
  integer::i,max,min, n, max1,max2, min1, min2,nsiti,nso,nelettroni, a,b,count2,prod, d, e, countd, counte, nfunctions, config,  check, h, coppie, countelem, ab, cd,p, k, count, countelet,z, np,np1,q, r, g, dim2, contasito, f, temp, m, binarysearch, ef, s, nomega, ecc
  integer,allocatable:: vecconfig(:), vecconfig_all(:), vecrow(:), veccol(:), doppia(:), doppian(:), nvec(:), occupazioni(:), doppian2(:), numero(:), peso(:),occ(:), nz(:)
  real*8,allocatable::hamval(:), xvec(:), yvec(:), mux(:), muy(:), mmu_x(:), mmu_y(:), autovalori(:),es(:), U(:)
  real*8::  energia, t, carica4, carica8, carica12, caricaN, caricax, caricay, l, l1, PP1, PP, V2, V1, esn2, PP2, Un2, mu_x, mu_y, carica1, carica2, carica3, carica5, carica6, carica7, carica9, carica10, carica11, y, x, pi, alpha, pesomax,omegastart, omegaend, domega, omega, sigma2,int, aaa, energia1
  character:: array(26), array1(26)
  double precision::c, sz, spin
  logical::bool, bool1, bool2, bool3, bool4, bool5, bool6, bool7, bool8,bool9, n1, n2
    external binarysearch

!!!===================ARPACK VARIABLES=============================
  !     %-----------------------------%
  !     | Define leading dimensions   |
  !     | for all arrays.             |
  !     | MAXN:   Maximum dimension   |
  !     |         of the A allowed.   |
  !     | MAXNEV: Maximum NEV allowed |
  !     | MAXNCV: Maximum NCV allowed |
  !     %-----------------------------%
  !
  integer          maxn, maxnev, maxncv, ldv
  !parameter       (maxn=400000, maxnev=100, maxncv=100000,& 
  !                ldv=maxn )
  parameter       (maxn=400000, maxnev=10000, maxncv=100000)
  !     %--------------%
  !     | Local Arrays |
  !     %--------------%
  !
  !Double precision&
  !                v(ldv,maxncv), workl(maxncv*(maxncv+8)),&
  !                workd(3*maxn), dd(maxncv,2), resid(maxn),&
  !                ax(maxn)
  double precision, allocatable :: v(:,:),workl(:),workd(:),dd(:,:),resid(:),ax(:)
  !logical          select(maxncv)
  logical,allocatable :: select(:)
  integer          iparam(11), ipntr(11)
  !
  !     %---------------%
  !     | Local Scalars |
  !     %---------------%
  !
  character        bmat*1, which*2
  integer          ido, nn, nev, ncv, lworkl, info, ierr, j, &
       nx, nconv, maxitr, mode, ishfts
  logical          rvec
  Double precision   &   
       tol, sigma
  !
  !     %------------%
  !     | Parameters |
  !     %------------%
  !
  Double precision&
       zero
  parameter        (zero = 0.0D+0)
  !  
  !     %-----------------------------%
  !     | BLAS & LAPACK routines used |
  !     %-----------------------------%
  !
  Double precision &          
       dnrm2
  external         dnrm2, daxpy
  !
  !     %--------------------%
  !     | Intrinsic function |
  !     %--------------------%
  !
  !      intrinsic        abs
  !
  character*1:: uplo
!!!==================================================
  !file da aprire
  open(1,file='bittest.dat')
  open(2, file='configurazioni.dat')
  open(3, file='geometria.dat')
  open(5,file='eigenvectors.dat')
  open(4,file='eigenvalues.dat')
  open(6, file='momenti_di_dipolo.dat')
  open(7,file='cariche.dat')
  open(8,file='spettro.dat')
  open(9,file='post-processing.dat')


  !=========================SET DI BASE=========================   
  write(*,*) 'doppie occupazioni'
  read(*,*) coppie

  max2=0
  do i=25-2*coppie,2*coppie-1,-2
     max2=max2+2**i
     if((i/2)*2.eq.i)stop
  enddo
  max1=0
  do i=25,26-2*coppie,-1
     max1=max1+2**i
  enddo
  max=max1+max2
  min1=0
  do i=2*coppie,26-2*coppie,2
     min1=min1+2**i
     if(i.ne.(i/2)*2)stop
  enddo
  min2=0
  do i=0,2*coppie-1
     min2=min2+2**i
  enddo
  min=min1+min2
  
  !=========================CONFIGURAZIONI====================
  nsiti=13
  write(*,*) 'dammi lo spazio di spin'
  read(*,*) spin
  nso=nsiti*2
  nelettroni=nsiti+1
  nfunctions=0
  do s=min,max
     count=0
     a=0; b=0; d=0; e=0
     count2=0
     countd=0
     counte=0
     prod=1
     sz=0.d0
     config=0
     h=0
     p=0
     check=0

     do i=0,nso-1
        bool=btest(s,i)
        count2=count2+1
        if(bool)then
           array(i+1)='1'
           count=count+1
           if(i/2*2.eq.i)then
              a=a+1
              countd=countd+1
           endif

           if(i/2*2.ne.i)then
              b=b+1
              counte=counte+1
           endif
        else
           array(i+1)='0'
           p=0
        endif
        prod=prod*(a+b)
        h=a+b
        if(h.eq.2)then
           check=check+1
        endif
        if(count2.eq.2)then
           a=0
           b=0
           d=0
           e=0
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
        write(2,*) config
     endif
  enddo
  close(1)
  close(2)
  write(*,*)'numero di funzioni di base uguale', nfunctions
  !=========================GEOMETRIA MOLECOLARE=========================
  l=1.50
  l1=l
  pi=dacos(-1.d0)
  alpha=pi/6
  allocate(xvec(nsiti), yvec(nsiti))
  write(*,*) dsin(alpha), dcos(alpha)
  open(1,file='geometrie.dat')
  do i=1,13
     if(i==1)then
        x=-l*dcos(alpha)
        y=l1+l*dsin(alpha)
     endif
     if(i==2)then
        x=-2*l*dcos(alpha)
        y=l
     endif
     if(i==3)then
        x=-(l+l1)*dcos(alpha)
        y=0
     endif
     if(i==4)then
        x=-l1*dcos(alpha)
        y=-l1*dsin(alpha)
     endif
     if(i==5)then
        x=-l*dcos(alpha)
        y=-(l1*dsin(alpha)+l)
     endif
     if(i==6)then
        x=0
        y=-(l1*dsin(alpha)+l+l*dsin(alpha))
     endif
     if(i==7)then
        x=l*dcos(alpha)
        y=-(l1*dsin(alpha)+l)
     endif
     if(i==8)then
        x=l1*dcos(alpha)
        y=-(l1*dsin(alpha))
     endif
     if(i==9)then
        x=(l+l1)*dcos(alpha)
        y=0
     endif
     if(i==10)then
        x=2*l*dcos(alpha)
        y=l
     endif
     if(i==11)then
        x=(l)*dcos(alpha)
        y=(l1+l*dsin(alpha))
     endif
     if(i==12)then
        x=0
        y=(l1)
     endif
     if(i==13)then
        x=0
        y=0
     endif
     xvec(i)=x
     yvec(i)=y
     write(3,*)i, x, y
  enddo
  close(3)
  !=========================SCRITTURA E DIAGONALIZZAZIONE=========================
  allocate(U(nsiti),occ(nsiti),es(nsiti),nz(i))
  open(10,file='input.dat')
  do i=1,nsiti
     read(10,*) aaa
     U(i)=aaa
  enddo
  do i=1,nsiti
     read(10,*) aaa
     es(i)=aaa
  enddo
  do i=1,nsiti
     read(10,*) a
     nz(i)=a
  enddo
  z=coppie
  open(2,file='configurazioni.dat')

  dim2=nfunctions

  allocate(vecconfig_all(dim2),hamval(10*dim2), vecrow(dim2+1), veccol(10*dim2), vecconfig(dim2), doppia(dim2),doppian(dim2),occupazioni(nsiti),doppian2(dim2))

  do n=1,dim2
     read(2,*) vecconfig_all(n)
  enddo
  close(2)

  hamval=0
  t=2.5d0
  countelem=0
  doppian=0
  doppia=0
  doppian2=0

!!! Ohno parametrization as reported in jpca (2013) 117, 7804

  do n=1,dim2 
     !  ifunc=binarysearch(1,dim2,vecconfig_all,vecconfig_all(n))

     !   diagonale     =================================================== 

     do p=0,2*nsiti-2,2
        bool=btest(vecconfig_all(n),p)
        bool1=btest(vecconfig_all(n),p+1)
        if(bool.and.bool1)then
           doppia((p+2)/2)=1
        else
           doppia((p+2)/2)=0
        endif

     enddo
     a=0
     b=0
     do i=0,2*nsiti-2,2
        bool=btest(vecconfig_all(n),i)
        bool1=btest(vecconfig_all(n),i+1)
        if(bool)then
           a=1
        else
           a=0
        endif
        if(bool1)then
           b=1
        else
           b=0
        endif
        occ((i+2)/2)=a+b
     enddo
     !==================================================
     contasito=0
     do i=0,2*nsiti-2,2
        contasito=contasito+1
        bool=btest(vecconfig_all(n),i)
        bool1=btest(vecconfig_all(n),i+1)
        if(bool)then
           a=1
        else
           a=0
        endif
        if(bool1)then
           b=1
        else
           b=0
        endif
        occupazioni(contasito)=a+b
     enddo


     PP=0.d0
     do i=1,nsiti
        do p=1,nsiti
           if(i.ne.p)then
              V2=(14.397)/((28.794/(U(i)+U(p)))**2+((xvec(i)-xvec(p))**2+(yvec(i)-yvec(p))**2))**0.5
              PP=PP+V2*(nz(i)-occupazioni(i))*(nz(p)-occupazioni(p))
           endif
        enddo
     enddo
     !==================================================
     energia1=0.d0
     do i=1,nsiti
        energia1=energia1 +U(i)*doppia(i) + es(i)*occ(i)
     enddo
     energia=energia1+0.5*PP
!!$!  write(*,*) conta

     countelem=countelem+1
     hamval(countelem)=energia
     vecrow(n)=countelem
     veccol(countelem)=n

!!$!  hopping anello  ==================================================

     a=0
     b=0
     c=0
     d=0
     ab=0
     cd=0
     doppia=0

     do p=0,2*nsiti-2,2
        bool=btest(vecconfig_all(n),p)
        bool1=btest(vecconfig_all(n),p+1)
        if(bool.and.bool1)doppia(n)=doppia(n)+1

     enddo

     do i=0,2*nsiti-2,2

        if(ab.le.3.and.i/2*2.eq.i)then
           ab=0
        endif

        if(cd.le.3.and.i/2*2.eq.i)then
           cd=0
        endif

        if(i/2*2.eq.i)then
           a=0;b=0; c=0; d=0
        endif

        bool=btest(vecconfig_all(n),i)
        if(bool)then
           a=1
        else
           a=0
        endif

        bool1=btest(vecconfig_all(n),i+1)
        if(bool1)then
           b=2
        else
           b=0
        endif

        ab= a+b               !deifinisco ab

        bool2=btest(vecconfig_all(n),i+2)
        if(bool2)then
           c=1
        else
           c=0
        endif

        bool3=btest(vecconfig_all(n),i+3)
        if(bool3)then
           d=2
        else
           d=0
        endif

        cd= c+d               !definisco cd

        if((ab.eq.1).and.(cd.eq.0))then
           temp=ibclr(vecconfig_all(n),i)
           temp=ibset(temp,i+2)
           countelem=countelem+1
           m=binarysearch(1,dim2,vecconfig_all,temp)
           veccol(countelem)=m
           hamval(countelem)=-t
        endif

        if((ab.eq.1).and.(cd.eq.2).and.(doppia(n).lt.z))then
           temp=ibclr(vecconfig_all(n),i)
           temp=ibset(temp,i+2)
           countelem=countelem+1
           m=binarysearch(1,dim2,vecconfig_all,temp)
           veccol(countelem)=m
           hamval(countelem)=-t
        endif

        if((ab.eq.2).and.(cd.eq.0))then
           temp=ibclr(vecconfig_all(n),i+1)
           temp=ibset(temp,i+3)
           countelem=countelem+1
           m=binarysearch(1,dim2,vecconfig_all,temp)
           veccol(countelem)=m
           hamval(countelem)=-t
        endif

        if((ab.eq.2).and.(cd.eq.1).and.(doppia(n).lt.z))then
           temp=ibclr(vecconfig_all(n),i+1)
           temp=ibset(temp,i+3)
           countelem=countelem+1
           m=binarysearch(1,dim2,vecconfig_all,temp)
           veccol(countelem)=m
           hamval(countelem)=+t
        endif

        if((ab.eq.3).and.(cd.eq.0))then
           temp=ibclr(vecconfig_all(n),i)
           temp=ibset(temp,i+2)
           countelem=countelem+1
           m=binarysearch(1,dim2,vecconfig_all,temp)
           veccol(countelem)=m
           hamval(countelem)=+t
        endif

        if((ab.eq.3).and.(cd.eq.2))then
           temp=ibclr(vecconfig_all(n),i)
           temp=ibset(temp,i+2)
           countelem=countelem+1
           m=binarysearch(1,dim2,vecconfig_all,temp)
           veccol(countelem)=m
           hamval(countelem)=+t
        endif

        if((ab.eq.3).and.(cd.eq.0))then
           temp=ibclr(vecconfig_all(n),i+1)
           temp=ibset(temp,i+3)
           countelem=countelem+1
           m=binarysearch(1,dim2,vecconfig_all,temp)
           veccol(countelem)=m
           hamval(countelem)=-t
        endif

        if((ab.eq.3).and.(cd.eq.1))then
           temp=ibclr(vecconfig_all(n),i+1)
           temp=ibset(temp,i+3)
           countelem=countelem+1
           m=binarysearch(1,dim2,vecconfig_all,temp)
           veccol(countelem)=m
           hamval(countelem)=+t
        endif
!!!!===================PERIODIC BOUNDARY CONDITION=========================================
        if(i.eq.0)then
           countelet=0

           do k=2,21
              bool4=btest(vecconfig_all(n),k)
              if(bool4)countelet=countelet+1
           enddo

           ab=0; cd=0
           a=0;b=0; c=0; d=0; 
           bool=btest(vecconfig_all(n),0)
           if(bool)then
              a=1
           else
              a=0
           endif

           bool1=btest(vecconfig_all(n),1)
           if(bool1)then
              b=2
           else
              b=0
           endif
           ab= a+b               !deifinisco ab

           bool2=btest(vecconfig_all(n),22)
           if(bool2)then
              c=1
           else
              c=0
           endif

           bool3=btest(vecconfig_all(n),23)
           if(bool3)then
              d=2
           else
              d=0
           endif
           cd= c+d               !definisco cd

           if((ab.eq.1).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),0)
              temp=ibset(temp,22)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.1).and.(cd.eq.2).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),0)
              temp=ibset(temp,22)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.2).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),1)
              temp=ibset(temp,23)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.2).and.(cd.eq.1).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),1)
              temp=ibset(temp,23)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),0)
              temp=ibset(temp,22)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.2))then
              temp=ibclr(vecconfig_all(n),0)
              temp=ibset(temp,22)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),1)
              temp=ibset(temp,23)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.1))then
              temp=ibclr(vecconfig_all(n),1)
              temp=ibset(temp,23)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif
        endif
!!!================HOPPING 4-13============================================
        if(i==6)then
           ab=0; cd=0
           a=0;b=0; c=0; d=0; 
           bool=btest(vecconfig_all(n),6)
           if(bool)then
              a=1
           else
              a=0
           endif
           bool1=btest(vecconfig_all(n),7)
           if(bool1)then
              b=2
           else
              b=0
           endif

           ab= a+b               !deifinisco ab
           !    ef=ab

           bool2=btest(vecconfig_all(n),24)
           if(bool2)then
              c=1
           else
              c=0
           endif

           bool3=btest(vecconfig_all(n),25)
           if(bool3)then
              d=2
           else
              d=0
           endif

           cd= c+d               !definisco cd
           countelet=0
           do e=8,23
              bool4=btest(vecconfig_all(n),e)
              if(bool4)countelet=countelet+1
           enddo
           if((ab.eq.1).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),6)
              temp=ibset(temp,24)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.1).and.(cd.eq.2).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),6)
              temp=ibset(temp,24)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.2).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),7)
              temp=ibset(temp,25)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.2).and.(cd.eq.1).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),7)
              temp=ibset(temp,25)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),6)
              temp=ibset(temp,24)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),7)
              temp=ibset(temp,25)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.2))then
              temp=ibclr(vecconfig_all(n),6)
              temp=ibset(temp,24)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.1))then
              temp=ibclr(vecconfig_all(n),7)
              temp=ibset(temp,25)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif
        endif

!!!====================HOPPING 8-13========================================
        if(i==8)then
           a=0; b=0; c=0; d=0; cd=0; ab=0
           bool=btest(vecconfig_all(n),14)
           if(bool)then
              a=1
           else
              a=0
           endif

           bool1=btest(vecconfig_all(n),15)
           if(bool1)then
              b=2
           else
              b=0
           endif
           ab= a+b               !deifinisco ab

           bool2=btest(vecconfig_all(n),24)
           if(bool2)then
              c=1
           else
              c=0
           endif

           bool3=btest(vecconfig_all(n),25)
           if(bool3)then
              d=2
           else
              d=0
           endif

           cd= c+d
           countelet=0
           do e=16,23
              bool4=btest(vecconfig_all(n),e)
              if(bool4)countelet=countelet+1
           enddo

           if((ab.eq.1).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),14)
              temp=ibset(temp,24)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=t
              endif
           endif

           if((ab.eq.1).and.(cd.eq.2).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),14)
              temp=ibset(temp,24)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.2).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),15)
              temp=ibset(temp,25)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.2).and.(cd.eq.1).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),15)
              temp=ibset(temp,25)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=t
              else
                 hamval(countelem)=-t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),14)
              temp=ibset(temp,24)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),15)
              temp=ibset(temp,25)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.2))then
              temp=ibclr(vecconfig_all(n),14)
              temp=ibset(temp,24)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif

           if((ab.eq.3).and.(cd.eq.1))then
              temp=ibclr(vecconfig_all(n),15)
              temp=ibset(temp,25)
              countelem=countelem+1
              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
           endif
        endif
     enddo

  enddo
  vecrow(dim2+1)=countelem+1
  write(*,*) 'STARTING DIAGONALIZATION'
!!!=====================ARPACK STARTING=============================
  uplo='U'

  !     %-----------------------%
  !     | Executable Statements |
  !     %-----------------------%
  !
  !     %----------------------------------------------------%
  !     | The number NX is the number of interior points     |
  !     | in the discretization of the 2-dimensional         |
  !     | Laplacian on the unit square with zero Dirichlet   |
  !     | boundary condition.  The number NN(=NX*NX) is the  |
  !     | dimension of the matrix.  A standard eigenvalue    |
  !     | problem is solved (BMAT = 'I'). NEV is the number  |
  !     | of eigenvalues to be approximated.  The user can   |
  !     | modify NEV, NCV, WHICH to solve problems of        |
  !     | different sizes, and to get different parts of the |
  !     | spectrum.  However, The following conditions must  |
  !     | be satisfied:                                      |
  !     |                  NN <= MAXN,                       | 
  !     |                 NEV <= MAXNEV,                     |
  !     |             NEV + 1 <= NCV <= MAXNCV               | 
  !     %----------------------------------------------------% 
  !
  nx = dim2
  nn = maxn !nx*nx
  nev =  10 
  ncv =  2*nev +2!22 
  ldv=nx

  allocate(v(nx,ncv),workl(ncv*(ncv+8)),workd(3*nx),dd(nev,2),resid(nx),ax(nx),&
       select(ncv))
  v=0.d0; workl=0.d0; workd=0.d0; d=0.d0; resid=0d0; ax=0.d0; select=0.d0

  if ( nn .gt. maxn ) then
     print *, ' ERROR with _SDRV1: N is greater than MAXN '
     go to 9000
  else if ( nev .gt. maxnev ) then
     print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
     go to 9000
  else if ( ncv .gt. maxncv ) then
     print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
     go to 9000
  end if
  bmat = 'I'
  which = 'SA'
  !
  !     %--------------------------------------------------%
  !     | The work array WORKL is used in DSAUPD as        |
  !     | workspace.  Its dimension LWORKL is set as       |
  !     | illustrated below.  The parameter TOL determines |
  !     | the stopping criterion.  If TOL<=0, machine      |
  !     | precision is used.  The variable IDO is used for |
  !     | reverse communication and is initially set to 0. |
  !     | Setting INFO=0 indicates that a random vector is |
  !     | generated in DSAUPD to start the Arnoldi         |
  !     | iteration.                                       |
  !     %--------------------------------------------------%
  !
  lworkl = ncv*(ncv+8)
  tol = zero 
  info = 0
  ido = 0
  !
  !     %---------------------------------------------------%
  !     | This program uses exact shifts with respect to    |
  !     | the current Hessenberg matrix (IPARAM(1) = 1).    |
  !     | IPARAM(3) specifies the maximum number of Arnoldi |
  !     | iterations allowed.  Mode 1 of DSAUPD is used     |
  !     | (IPARAM(7) = 1).  All these options may be        |
  !     | changed by the user. For details, see the         |
  !     | documentation in DSAUPD.                          |
  !     %---------------------------------------------------%
  !
  ishfts = 1
  maxitr = 300000
  mode   = 1
  !      
  iparam(1) = ishfts 
  iparam(3) = maxitr 
  iparam(7) = mode


  !     %-------------------------------------------%
  !     | M A I N   L O O P (Reverse communication) |
  !     %-------------------------------------------%

10 continue

  !        %---------------------------------------------%
  !        | Repeatedly call the routine DSAUPD and take | 
  !        | actions indicated by parameter IDO until    |
  !        | either convergence is indicated or maxitr   |
  !        | has been exceeded.                          |
  !        %---------------------------------------------%

  call dsaupd ( ido, bmat,nx, which, nev, tol, resid, &
       ncv, v, ldv, iparam, ipntr, workd, workl,&
       lworkl, info )

  if (ido .eq. -1 .or. ido .eq. 1) then

     !           %--------------------------------------%
     !           | Perform matrix vector multiplication |
     !           |              y <--- OP*x             |
     !           | The user should supply his/her own   |
     !           | matrix vector multiplication routine |
     !           | here that takes workd(ipntr(1)) as   |
     !           | the input, and return the result to  |
     !           | workd(ipntr(2)).                     |
     !           %--------------------------------------%

     !            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
     call mkl_dcsrsymv(uplo, dim2, hamval, vecrow, veccol, workd(ipntr(1)), workd(ipntr(2)))
     !           %-----------------------------------------%
     !           | L O O P   B A C K to call DSAUPD again. |
     !           %-----------------------------------------%

     go to 10

  end if

  !     %----------------------------------------%
  !     | Either we have convergence or there is |
  !     | an error.                              |
  !     %----------------------------------------%

  if ( info .lt. 0 ) then

     !        %--------------------------%
     !        | Error message. Check the |
     !        | documentation in DSAUPD. |
     !        %--------------------------%

     print *, ' '
     print *, ' Error with _saupd, info = ', info
     print *, ' Check documentation in _saupd '
     print *, ' '

  else 

     !        %-------------------------------------------%
     !        | No fatal errors occurred.                 |
     !        | Post-Process using DSEUPD.                |
     !        |                                           |
     !        | Computed eigenvalues may be extracted.    |  
     !        |                                           |
     !        | Eigenvectors may also be computed now if  |
     !        | desired.  (indicated by rvec = .true.)    | 
     !        %-------------------------------------------%

     rvec = .true.

     call dseupd ( rvec, 'All', select, dd, v, ldv, sigma, &
          bmat,nx, which, nev, tol, resid, ncv, v, ldv, &
          iparam, ipntr, workd, workl, lworkl, ierr )

     !        %----------------------------------------------%
     !        | Eigenvalues are returned in the first column |
     !        | of the two dimensional array D and the       |
     !        | corresponding eigenvectors are returned in   |
     !        | the first NEV columns of the two dimensional |
     !        | array V if requested.  Otherwise, an         |
     !        | orthogonal basis for the invariant subspace  |
     !        | corresponding to the eigenvalues in D is     |
     !        | returned in V.                               |
     !        %----------------------------------------------%
     !
     if ( ierr .ne. 0) then

        !            %------------------------------------%
        !            | Error condition:                   |
        !            | Check the documentation of DSEUPD. |
        !            %------------------------------------%

        print *, ' '
        print *, ' Error with _seupd, info = ', ierr
        print *, ' Check the documentation of _seupd. '
        print *, ' '

     else

        nconv =  iparam(5)
        allocate(autovalori(nev))
        do i=1,nev
           write (4,*) i,dd(i,1)
           autovalori(i)=dd(i,1)
           if(i.le.7)write(5,*)'EIGENVECTORS OF EIGENSTATE', i
           do j= 1,dim2
              if((i.le.7).and.(dabs(v(j,i)).gt.0.001))write(5,*) j, v(j,i)
           enddo
        enddo
        !write (700,*) (dd(i,1),i=1,nev)

!!$             do 20 j=1, nconv
!!$
!!$!               %---------------------------%
!!$!               | Compute the residual norm |
!!$!               |                           |
!!$!               |   ||  A*x - lambda*x ||   |
!!$!               |                           |
!!$!               | for the NCONV accurately  |
!!$!               | computed eigenvalues and  |
!!$!               | eigenvectors.  (iparam(5) |
!!$!               | indicates how many are    |
!!$!               | accurate to the requested |
!!$!               | tolerance)                |
!!$!               %---------------------------%
!!$
!!$                call av(nx, v(1,j), ax)
!!$                !call daxpy(nn, -dd(j,1), v(1,j), 1, ax, 1)
!!$                dd(j,2) = dnrm2(n, ax, 1)
!!$                dd(j,2) = dd(j,2) / abs(dd(j,1))
!!$
!!$ 20          continue

        !            %-------------------------------%
        !            | Display computed residuals    |
        !            %-------------------------------%
        !
        !call dmout(6, nconv, 2, dd, maxncv, -6,&
        !     'Ritz values and relative residuals')
     end if

     !        %------------------------------------------%
     !        | Print additional convergence information |
     !        %------------------------------------------%

     if ( info .eq. 1) then
        print *, ' '
        print *, ' Maximum number of iterations reached.'
        print *, ' '
     else if ( info .eq. 3) then
        print *, ' ' 
        print *, ' No shifts could be applied during implicit',&
             ' Arnoldi update, try increasing NCV.'
        print *, ' '
     end if

     print *, ' '
     print *, ' _SDRV1 '
     print *, ' ====== '
     print *, ' '
     print *, ' Size of the matrix is ', nn
     print *, ' The number of Ritz values requested is ', nev
     print *, ' The number of Arnoldi vectors generated',&
          ' (NCV) is ', ncv
     print *, ' What portion of the spectrum: ', which
     print *, ' The number of converged Ritz values is ', &
          nconv 
     print *, ' The number of Implicit Arnoldi update',&
          ' iterations taken is ', iparam(3)
     print *, ' The number of OP*x is ', iparam(9)
     print *, ' The convergence criterion is ', tol
     print *, ' '
     !
  end if
  !
  !     %---------------------------%
  !     | Done with program dsdrv1. |
  !     %---------------------------%
  !
9000 continue
  !
!!!=====================ARPACK ENDING=============================        
  !=========================Valore attesa carica=========================
  do i=1,nev
     carica4=0.d0
     carica8=0.d0
     carica12=0.d0
     caricaN=0.d0
     carica1=0.d0
     carica2=0.d0
     carica3=0.d0
     carica5=0.d0
     carica6=0.d0
     carica7=0.d0
     carica9=0.d0
     carica10=0.d0
     carica11=0.d0
     carica12=0.d0
     
     do n=1,dim2
        !_________________________carbonio1________________________
        bool1=btest(vecconfig_all(n),0)
        bool2=btest(vecconfig_all(n),1)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica1=carica1+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio2________________________
        bool1=btest(vecconfig_all(n),2)
        bool2=btest(vecconfig_all(n),3)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica2=carica2+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio3________________________
        bool1=btest(vecconfig_all(n),4)
        bool2=btest(vecconfig_all(n),5)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica3=carica3+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio4________________________
        bool1=btest(vecconfig_all(n),6)
        bool2=btest(vecconfig_all(n),7)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica4=carica4+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio5________________________
        bool1=btest(vecconfig_all(n),8)
        bool2=btest(vecconfig_all(n),9)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica5=carica5+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio6________________________
        bool1=btest(vecconfig_all(n),10)
        bool2=btest(vecconfig_all(n),11)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica6=carica6+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio7________________________
        bool1=btest(vecconfig_all(n),12)
        bool2=btest(vecconfig_all(n),13)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica7=carica7+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio8________________________
        bool1=btest(vecconfig_all(n),14)
        bool2=btest(vecconfig_all(n),15)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica8=carica8+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio9________________________
        bool1=btest(vecconfig_all(n),16)
        bool2=btest(vecconfig_all(n),17)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica9=carica9+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio10________________________
        bool1=btest(vecconfig_all(n),18)
        bool2=btest(vecconfig_all(n),19)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica10=carica10+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio11________________________
        bool1=btest(vecconfig_all(n),20)
        bool2=btest(vecconfig_all(n),21)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica11=carica11+(1-a-b)*(V(n,i)**2)
        !_________________________carbonio12________________________
        bool1=btest(vecconfig_all(n),22)
        bool2=btest(vecconfig_all(n),23)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        carica12=carica12+(1-a-b)*(V(n,i)**2)
        !_________________________azoto________________________
        bool1=btest(vecconfig_all(n),24)
        bool2=btest(vecconfig_all(n),25)
        if(bool1)then
           a=1
        else
           a=0
        endif
        if(bool2)then
           b=1
        else
           b=0
        endif
        caricaN=caricaN+(2-a-b)*(V(n,i)**2)
     enddo
     write(7,'(2x,i3,14(2x,f20.13))') i,carica1, carica2, carica3, carica4,carica5, carica6, carica7, carica8,carica9, carica10, carica11, carica12, caricaN

  enddo
  close(7)
  !=========================operatore momento di dipolo========================
  allocate(mux(dim2),muy(dim2))
  mux=0.d0
  muy=0.d0
  do n=1,dim2
     contasito=0
     do i=0,2*nsiti-2,2
        contasito=contasito+1
        bool=btest(vecconfig_all(n),i)
        bool1=btest(vecconfig_all(n),i+1)
        if(bool)then
           a=1
        else
           a=0
        endif
        if(bool1)then
           b=1
        else
           b=0
        endif
        occupazioni(contasito)=a+b
     enddo


     do i=1,11,2 !!!azoti aza
        mux(n)=mux(n)+xvec(i)*(1-occupazioni(i))
     enddo

     do i=2,12,2
        mux(n)=mux(n)+xvec(i)*(1-occupazioni(i))
     enddo

     do i=1,11,2 !!carboni 
        muy(n)=muy(n)+yvec(i)*(1-occupazioni(i))
     enddo

     do i=2,12,2
        muy(n)=muy(n)+yvec(i)*(1-occupazioni(i))
     enddo

  enddo

  !=========================momenti di dipolo=========================
allocate(mmu_x(nev), mmu_y(nev))
  do i=1,nev
     mu_x=0.d0
     mu_y=0.d0
     do j=1,dim2
        mu_x=mu_x+v(j,i)*v(j,1)*mux(j)
        mu_y=mu_y+v(j,i)*v(j,1)*muy(j)
     enddo
     mmu_x(i)=mu_x
     mmu_y(i)=mu_y
     write(6,*)i, mu_x, mu_y
  enddo
  close(6)
  !=========================ECCITONI================================
  n=10000
  pesomax=0.0006

  do i=1,dim2
     j=vecconfig_all(i)
     ecc=0
     if((v(i,1))**2.ge.pesomax)then
        write(9,*)'configurazione=',j

        write(9,*)'peso=',(v(i,1))**2
        do s=0,22,2
           bool=btest(j,s)
           bool1=btest(j,s+1)
           if(bool)then
              b=1
              array1(s+1)='1'
           else
              b=0
              array1(s+1)='0'
           endif
           if(bool1)then
              c=1
              array1(s+2)='1'
           else
              c=0
              array1(s+2)='0'
           endif
           if(b+c.eq.2)then
              ecc=ecc+1
              write(9,*) '-',(s+2)/2
           endif
           if(b+c.eq.0)then
              write(9,*) '+',(s+2)/2
           endif
        enddo
        do s=24,24,2
           bool=btest(j,s)
           bool1=btest(j,s+1)
           if(bool)then
              b=1
              array1(25)='1'
           else
              b=0
              array1(25)='0'
           endif
           if(bool1)then
              c=1
              array1(26)='1'
           else
              c=0
              array1(26)='0'
           endif
           if(b+c.eq.1)then
              write(9,*) '+',(s+2)/2
           endif
           if(b+c.eq.0)then
              write(9,*) '2+',(s+2)/2
           endif
        enddo
        write(9,*) 'NUMERO DI ECCITONI=',ecc
        write(9,*)(array1(s),s=25,1,-1)
        write(9,*) '---------------------------------'
     endif
  enddo
  !========================SPETTRO=============================
  sigma2=0.1d0
  nomega=1000000
  omegastart=0.d0
  omegaend=8.0d0
  domega=(omegaend-omegastart)/nomega
  omega=omegastart-domega

  do j=1,nomega
   omega=omega+domega
   int=0.d0
   do i= 2, nev
      int=int+(autovalori(i)-autovalori(1)) * (mmu_x(i)**2+mmu_y(i)**2) *((2*pi)**(-0.5)*sigma2)* dexp(-(omega-(autovalori(i)-autovalori(1)))**2/(2*sigma2**2))
   enddo
   write(8,*) omega, int
enddo
! ==========================TRANSIZIONI==========================
write(9,*)'TRANSIZIONI'
WRITE(9,*) 'In ordine si hanno: stato di partenza, stato di arrivo ed energia della transizione in eV'
do i=1,nev
   do j=i,nev
      if((dabs(mmu_x(i)-mmu_x(j)).ge.0.001).and.(dabs(mmu_y(i)-mmu_y(j)).ge.0.001))then
         write(9,*) i,j,dd(j,1)-dd(i,1)
      endif
   enddo
enddo

endprogram trianguleni





!!!!!!!!!!!!!!
 subroutine av (nx, v, w)
      integer           nx, j, lo, n2
      Double precision &
                       v(nx*nx), w(nx*nx), one, h2
      parameter         ( one = 1.0D+0 ) 

      call tv(nx,v(1),w(1))
      call daxpy(nx, -one, v(nx+1), 1, w(1), 1)

      do 10 j = 2, nx-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call daxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
         call daxpy(nx, -one, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 

      lo = (nx-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call daxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
      
      !     Scale the vector w by (1/h^2), where h is the mesh size

      n2 = nx*nx
      h2 = one / dble((nx+1)*(nx+1))
      call dscal(n2, one/h2, w, 1) 
      return
      end
!
!-------------------------------------------------------------------
      subroutine tv (nx, x, y)
!
      integer           nx, j 
      Double precision &
                       x(nx), y(nx), dd, dl, du
!
      Double precision &
                      one
      parameter        (one = 1.0D+0 )
!
!     Compute the matrix vector multiplication y<---T*x
!     where T is a nx by nx tridiagonal matrix with DD on the 
!     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
!     
!
      dd  = 4.0D+0
      dl  = -one 
      du  = -one
! 
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1) 
 10   continue 
      y(nx) =  dl*x(nx-1) + dd*x(nx) 
      return
      end


!!!!!!!!!!!!!!



integer function binarysearch(i, length, array, val)
  ! Given an array and a value, returns the index of the element that
  ! is closest to, but less than, the given value.
  ! Uses a binary search algorithm.
  ! "delta" is the tolerance used to determine if two values are equal
  ! if ( abs(x1 - x2) <= delta) then
  ! assume x1 = x2
  ! endif

  implicit none

  integer, intent(in) :: length, i
  integer, dimension(length), intent(in) :: array
  integer, intent(in) :: val

  !integer :: binarysearch

  integer :: left, middle, right

  left = i
  right = length
  binarysearch=0

  
  if (val.lt.array(left) .or. val.gt.array(right)) go to 10
  
  do

     if (left .gt. right) then
        exit
        !write(*,*) 'ERRORE!!!'
     endif

     !divisione=((left+right) / 2.0)
     !middle = jnint(divisione)
     middle=(left+right)/2
     
     if ( array(middle) .eq. val ) then
        binarySearch = middle
        return
     else 
        if (array(middle) .gt. val) then
           right = middle - 1
        else
           left = middle + 1
        end if
     end if
  end do

  binarysearch = right
10 continue
end function binarysearch


