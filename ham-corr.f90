program hoppingtotal

  implicit none

  integer::nsiti, dim, dim2, i, countelem, n, a, b, c, d, e, f, ab, cd, ef, p, countelet, k, h, count, countelet1,z, np,np1,q, r, g, y, q1, r1,npn, np2, t2, y2, np3, np4, np5,y3, t3, r3 ,q3, y1, g1, r2, q2, t1,jj
  integer::  temp, m, binarysearch, contasito
  integer,allocatable:: vecconfig(:), vecconfig_all(:), vecrow(:), veccol(:), doppia(:), doppian(:), nvec(:), occupazioni(:), doppian2(:), hop(:),nz(:), occ(:,:)
  real*8,allocatable::hamval(:), xvec(:), yvec(:), mux(:), muy(:), pot(:), onsite(:), siten(:),u(:),evec(:),  mu_x(:), mu_y(:),zumba(:), valzer(:), tiptap(:),carica(:,:), nnn(:,:), spinden(:,:), hiphop(:,:)
  real*8,allocatable:: umat(:,:), emat(:,:), corrint(:,:), correx(:,:), vmat(:,:,:),double(:), double2(:,:), dsite(:,:), ssite(:)
  real*8:: esc, esn, conta, Uc, Un, energia, t, Valoreattesa, l, l1, PP1, PP, V2, V1, esn2, PP2, Un2,att,toll, ss
  logical:: bool, bool1, bool2, bool3, bool4, bool5, bool6, bool7, bool8,bool9, n1, n2
  real*8::strenght,hbar,hbarev,echarge,emass,dpg, bubu, balu
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
     
  nsiti=13
  esc=0.d0
  esn=-12.0d0
  esn2=-2.d0
  Uc=11.26d0
  Un=15.d0
  Un2=Uc
  toll=0.0001d0
 
 ! open(17,file='input.dat')
 ! read(17,*) Un
 ! close(17)
  open(2,file='numero_doppie.dat')
  read(2,*)z
  close(2)
 ! open(4,file='ham-unico.dat')
  open(5,file='configurazioni.dat')
  open(6, file='dim.dat')
 
  open(700,file='eigenvalues.dat')
  open(800,file='geometrie.dat')
  open(1000, file='momenti_di_dipolo.dat')
  read(6,*) dim2
  write(*,*) 'dim2=', dim2
!!$  dim2=228228
!  read(2,*) dim
!  write(*,*)'dim=',dim
  
  allocate(vecconfig_all(dim2),hamval(10*dim2), vecrow(dim2+1), veccol(10*dim2), vecconfig(dim2), doppia(dim2),doppian(dim2),occupazioni(nsiti), xvec(nsiti), yvec(nsiti),doppian2(dim2),hop(dim2),pot(dim2),onsite(dim2),siten(dim2),nz(nsiti),u(nsiti),evec(nsiti))

  allocate(correx(nsiti,nsiti),vmat(nsiti,nsiti,dim2),double(dim2),double2(nsiti,dim2),occ(nsiti,dim2),corrint(nsiti,nsiti),dsite(nsiti,nsiti), ssite(nsiti),spinden(dim2,nsiti))

   
 
  do n=1,dim2
     read(5,*) vecconfig_all(n)
  !   write(101,*) vecconfig_all(n)
  enddo
!!$  do i=1,dim2
!!$     read(5,*) vecconfig(i)
!!$  enddo
  close(1)
  close(2)
  do i=1,nsiti
     read(800,*) xvec(i), yvec(i)
  enddo
  close(800)
  hamval=0
  t=2.5d0
  countelem=0
  doppian=0
  doppia=0
  doppian2=0
  hop=0
  pot=0
  onsite=0
  siten=0
  vmat=0
  spinden=0

  do i=1,13
     if(i.ne.13)then
        nz(i)=1
        u(i)=uc
        evec(i)=esc
     else
        nz(i)=2
        u(i)=un
        evec(i)=esn 
     endif
  enddo

!!! Ohno parametrization as reported in jpca (2013) 117, 7804
  
  do n=1,dim2 
     !  ifunc=binarysearch(1,dim2,vecconfig_all,vecconfig_all(n))
!!!!!!!!!!!!!!!!!

     !   diagonale     =================================================== 
     ss=0.0d0
     do p=0,2*nsiti-2,2
        bool=btest(vecconfig_all(n),p)
        bool1=btest(vecconfig_all(n),p+1)
        if(bool.and.bool1)then
           ss=ss+u((p+2)/2)
        endif
        if(bool)then
           bubu=0.5
        else
           bubu=0
        endif
        if(bool1)then
           balu=-0.5
        else
           balu=0
        endif
        spinden(n,(p+2)/2)=bubu+balu
     enddo
       
      
     conta=0.0d0
  
     do p=0,2*nsiti-1
        bool=btest(vecconfig_all(n),p)
        if(bool)then
           conta=conta+evec(p/2+1)
        endif
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
              V2=(14.397)/((28.794/(u(i)+U(p)))**2+((xvec(i)-xvec(p))**2+(yvec(i)-yvec(p))**2))**0.5
              PP=PP+V2*(nz(i)-occupazioni(i))*(nz(p)-occupazioni(p))
              vmat(i,p,n)=V2*(nz(i)-occupazioni(i))*(nz(p)-occupazioni(p))
           endif
        enddo
     enddo
     


          
        
 !==================================================    
     energia=conta+ss+0.5d0*PP
     
     
     pot(n)=0.5*PP
     countelem=countelem+1
     hamval(countelem)=energia
     vecrow(n)=countelem
     veccol(countelem)=n
     onsite(n)=ss
     siten(n)=conta
!!$!  hopping anello  ==================================================
    
     a=0
     b=0
     c=0
     d=0
     e=0
     f=0
     ab=0
     cd=0
     ef=0
     doppia=0
     
     do p=0,2*nsiti-2,2
        bool=btest(vecconfig_all(n),p)
        bool1=btest(vecconfig_all(n),p+1)
        if(bool.and.bool1)doppia(n)=doppia(n)+1
        if(bool.and.bool1)double(n)=double(n)+1
        if(bool.and.bool1)double2((p+2)/2,n)=double2((p+2)/2,n)+1
        if(bool)then
           e=1
        else
           e=0
        endif
        if(bool1)then
           f=1
        else
           f=0
        endif
        occ((p+2)/2,n)=e+f
     enddo
     ! write(200,*) n, doppia(n)
     
     do i=0,2*nsiti-3,2

        if(ab.le.3.and.i/2*2.eq.i)then
           ab=0
        endif

        if(cd.le.3.and.i/2*2.eq.i)then
           cd=0
        endif
        if(ef.le.3.and.i/2*2.eq.i)then
           ef=0
        endif
                
        if(i/2*2.eq.i)then
           a=0;b=0; c=0; d=0; e=0; f=0
        endif
      
        bool=btest(vecconfig_all(n),i)
        if(bool)then
           a=1
           !countelet=countelet+1
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
        
!!$        if(ef.eq.3)then
!!$           doppia=doppia+1    !conta doppie occupazioni    
!!$        endif

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
           !count=count+1

           m=binarysearch(1,dim2,vecconfig_all,temp)
           if(m.ne.0)hop(n)=hop(n)+1
           veccol(countelem)=m
          ! if(m.eq.0)write(300,*)'10=',temp
          
           hamval(countelem)=-t
          ! if(count.eq.1)then
           !   vecrow(n)=countelem
          ! endif

        endif

        if((ab.eq.1).and.(cd.eq.2).and.(doppia(n).lt.z))then
           temp=ibclr(vecconfig_all(n),i)
           temp=ibset(temp,i+2)

           countelem=countelem+1
          ! count=count+1

           m=binarysearch(1,dim2,vecconfig_all,temp)
           if(m.ne.0)hop(n)=hop(n)+1

           veccol(countelem)=m
           ! if(m.eq.0)write(300,*)'12=',temp
           hamval(countelem)=-t
!!$           if(count.eq.1)then
!!$              vecrow(n)=countelem
!!$           endif

        endif
        
        if((ab.eq.2).and.(cd.eq.0))then
           temp=ibclr(vecconfig_all(n),i+1)
           temp=ibset(temp,i+3)

           countelem=countelem+1
           !count=count+1

           m=binarysearch(1,dim2,vecconfig_all,temp)
           veccol(countelem)=m
           if(m.ne.0)hop(n)=hop(n)+1
           ! if(m.eq.0)write(300,*)'20=',temp
           hamval(countelem)=-t
!!$           if(count.eq.1)then
!!$              vecrow(n)=countelem
!!$           endif

        endif

        if((ab.eq.2).and.(cd.eq.1).and.(doppia(n).lt.z))then
           temp=ibclr(vecconfig_all(n),i+1)
           temp=ibset(temp,i+3)

           countelem=countelem+1
         !  count=count+1

           m=binarysearch(1,dim2,vecconfig_all,temp)
           if(m.ne.0)hop(n)=hop(n)+1

           veccol(countelem)=m
           ! if(m.eq.0)write(300,*)'21=',temp
           hamval(countelem)=+t
!!$           if(count.eq.1)then
!!$              vecrow(n)=countelem
!!$           endif

        endif
        
        if((ab.eq.3).and.(cd.eq.0))then
           temp=ibclr(vecconfig_all(n),i)
           temp=ibset(temp,i+2)

           countelem=countelem+1
          ! count=count+1

           m=binarysearch(1,dim2,vecconfig_all,temp)
           if(m.ne.0)hop(n)=hop(n)+1
                      
           veccol(countelem)=m
           ! if(m.eq.0)write(300,*)'30up=',temp
           hamval(countelem)=+t
!!$           if(count.eq.1)then
!!$              vecrow(n)=countelem
!!$           endif
        endif

        if((ab.eq.3).and.(cd.eq.2))then
           temp=ibclr(vecconfig_all(n),i)
           temp=ibset(temp,i+2)

           countelem=countelem+1
         !  count=count+1

           m=binarysearch(1,dim2,vecconfig_all,temp)
           if(m.ne.0)hop(n)=hop(n)+1

           veccol(countelem)=m
           ! if(m.eq.0)write(300,*)'32=',temp
           hamval(countelem)=+t
!!$           if(count.eq.1)then
!!$              vecrow(n)=countelem
!!$           endif

        endif

        if((ab.eq.3).and.(cd.eq.0))then
           temp=ibclr(vecconfig_all(n),i+1)
           temp=ibset(temp,i+3)

           countelem=countelem+1
         !  count=count+1

           m=binarysearch(1,dim2,vecconfig_all,temp)
           if(m.ne.0)hop(n)=hop(n)+1

           veccol(countelem)=m
           ! if(m.eq.0)write(300,*)'30down=',temp
           hamval(countelem)=-t
!!$           if(count.eq.1)then
!!$              vecrow(n)=countelem
!!$           endif
           
        endif

        if((ab.eq.3).and.(cd.eq.1))then
           temp=ibclr(vecconfig_all(n),i+1)
           temp=ibset(temp,i+3)
        
           countelem=countelem+1
       !    count=count+1
           
           m=binarysearch(1,dim2,vecconfig_all,temp)
           if(m.ne.0)hop(n)=hop(n)+1

           veccol(countelem)=m
           ! if(m.eq.0)write(300,*)'31=',temp
           hamval(countelem)=+t
!!$           if(count.eq.1)then
!!$              vecrow(n)=countelem
!!$           endif
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


!!$        if(ef.eq.3)then
!!$           doppia=doppia+1    !conta doppie occupazioni    
!!$        endif

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
            !  count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'10=',temp

              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.1).and.(cd.eq.2).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),0)
              temp=ibset(temp,22)

              countelem=countelem+1
             ! count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'12=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
             else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.2).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),1)
              temp=ibset(temp,23)

              countelem=countelem+1
            !  count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'20=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.2).and.(cd.eq.1).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),1)
              temp=ibset(temp,23)

              countelem=countelem+1
           !   count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'21=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),0)
              temp=ibset(temp,22)

              countelem=countelem+1
             ! count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'30up=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif
           endif

           if((ab.eq.3).and.(cd.eq.2))then
              temp=ibclr(vecconfig_all(n),0)
              temp=ibset(temp,22)

              countelem=countelem+1
            !  count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'32=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),1)
              temp=ibset(temp,23)

              countelem=countelem+1
           !   count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'30down=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.3).and.(cd.eq.1))then
              temp=ibclr(vecconfig_all(n),1)
              temp=ibset(temp,23)

              countelem=countelem+1
             ! count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'31=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif
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

!!$        if(ef.eq.3)then
!!$           doppia=doppia+1    !conta doppie occupazioni    
!!$        endif

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
            !!!!!!!!!!!  count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(m.ne.0)hop(n)=hop(n)+1
              ! if(m.eq.0)write(300,*)'10=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.1).and.(cd.eq.2).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),6)
              temp=ibset(temp,24)

              countelem=countelem+1
            !  count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'12=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.2).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),7)
              temp=ibset(temp,25)

              countelem=countelem+1
           !   count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              veccol(countelem)=m
              if(m.ne.0)hop(n)=hop(n)+1
              ! if(m.eq.0)write(300,*)'20=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.2).and.(cd.eq.1).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),7)
              temp=ibset(temp,25)

              countelem=countelem+1
             ! count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'21=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),6)
              temp=ibset(temp,24)

              countelem=countelem+1
           !   count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'30up=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif
           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),7)
              temp=ibset(temp,25)

              countelem=countelem+1
           !   count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'30down=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif
           endif

           if((ab.eq.3).and.(cd.eq.2))then
              temp=ibclr(vecconfig_all(n),6)
              temp=ibset(temp,24)

              countelem=countelem+1
          !    count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'32=',temp

              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif

!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif



           if((ab.eq.3).and.(cd.eq.1))then
              temp=ibclr(vecconfig_all(n),7)
              temp=ibset(temp,25)

              countelem=countelem+1
           !   count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'31=',temp

              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif

!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif
           endif
        endif

!!!====================HOPPING 8-13========================================
        if(i==14)then
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
           !    ef=ab

!!$        if(ef.eq.3)then
!!$           doppia=doppia+1    !conta doppie occupazioni    
!!$        endif

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
           !   count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'10=',temp

              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.1).and.(cd.eq.2).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),14)
              temp=ibset(temp,24)

              countelem=countelem+1
             ! count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'12=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.2).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),15)
              temp=ibset(temp,25)

              countelem=countelem+1
         !     count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'20=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.2).and.(cd.eq.1).and.(doppia(n).lt.z))then
              temp=ibclr(vecconfig_all(n),15)
              temp=ibset(temp,25)

              countelem=countelem+1
            !  count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'21=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),14)
              temp=ibset(temp,24)

              countelem=countelem+1
           !   count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'30up=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif
           endif

           if((ab.eq.3).and.(cd.eq.0))then
              temp=ibclr(vecconfig_all(n),15)
              temp=ibset(temp,25)

              countelem=countelem+1
           !   count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'30down=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=-t
              else
                 hamval(countelem)=+t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif

           if((ab.eq.3).and.(cd.eq.2))then
              temp=ibclr(vecconfig_all(n),14)
              temp=ibset(temp,24)

              countelem=countelem+1
            !  count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'32=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif

           endif



           if((ab.eq.3).and.(cd.eq.1))then
              temp=ibclr(vecconfig_all(n),15)
              temp=ibset(temp,25)

              countelem=countelem+1
          !    count=count+1

              m=binarysearch(1,dim2,vecconfig_all,temp)
              if(m.ne.0)hop(n)=hop(n)+1
              veccol(countelem)=m
              ! if(m.eq.0)write(300,*)'31=',temp
              if(countelet/2*2.eq.countelet)then
                 hamval(countelem)=+t
              else
                 hamval(countelem)=-t
              endif
!!$              if(count.eq.1)then
!!$                 vecrow(n)=countelem
!!$              endif
           endif
        endif
  
     enddo
  enddo
  vecrow(dim2+1)=countelem+1

write(*,*)'=================================='
write(*,*)'reached ARPACK'
write(*,*)'=================================='

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
      nev =  15 
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

 10   continue

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

            do i=1,nev
               write (700,*) i,dd(i,1)
               do j= 1,dim2
!!$                  if(i==1.and.(dabs(v(j,i)).gt.toll))write(900,*) j, v(j,i)
!!$                  if(i==2.and.(dabs(v(j,i)).gt.toll))write(901,*) j, v(j,i)
!!$                  if(i==3.and.(dabs(v(j,i)).gt.toll))write(902,*) j, v(j,i)
!!$                  if(i==4.and.(dabs(v(j,i)).gt.toll))write(903,*) j, v(j,i)
!!$                  if(i==5.and.(dabs(v(j,i)).gt.toll))write(904,*) j, v(j,i)
!!$                  if(i==6.and.(dabs(v(j,i)).gt.toll))write(905,*) j, v(j,i)
               enddo
            enddo
            
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
9000  continue
      !
!!!=====================ARPACK ENDING=============================
!===========Valore attesa carica=========================
open(115,file='cariche.dat')
allocate(carica(nev,nsiti),nnn(nev,nsiti))
carica=0
do i=1,nev
   do p=0,2*nsiti-2,2
      do n=1,dim2
         bool1=btest(vecconfig_all(n),p)
         bool2=btest(vecconfig_all(n),p+1)
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
         carica(i,(p+2)/2)=carica(i,(p+2)/2)+(nz((p+2)/2)-a-b)*(V(n,i)**2)
      enddo
   enddo
   nnn=carica
enddo
!!$do i=1,nev-1
!!$   if((dabs(dd(i,1)-dd(i+1,1)).le.1d-10))then
!!$      do p=1,nsiti
!!$         carica(i,p)=(nnn(i,p)+nnn(i+1,p))/2
!!$         carica(i+1,p)=carica(i,p)
!!$      enddo
!!$   endif
!!$enddo


do p=1,nev
   write(115,*) 'STATO',p
   write(115,'(<nsiti>(2x,f10.5))')(carica(i,p),i=1,nsiti)
   write(115,*)
enddo
allocate(zumba(nev), valzer(nev), tiptap(nev),hiphop(nev,nsiti))
!=========================spin-density=========================
open(1234,file='spin-density.dat')

do i=1,nev
   do j=1,nsiti
      do n=1,dim2
         hiphop(i,j)=hiphop(i,j)+spinden(n,j)*v(n,i)**2
         
      enddo
   enddo
enddo
do i=1,nev
   write(1234,*) 'STATO', i
   write(1234,*)(hiphop(i,j),j=1,nsiti)
   write(1234,*)
enddo


!=========================Valore attesa potenziale=========================
open(1113,file='attending.dat')
open(1114,file='config-pot.dat')
zumba=0
do i=1,nev
 
   write(1113,*) 'AUTOVALORE=',i
   write(1114,*) 'AUTOVALORE=',i
   do n=1,dim2
      zumba(i)=zumba(i)+v(n,i)**2*pot(n)
      if(v(n,i)**2.ge.toll) write(1114,*) vecconfig_all(n), pot(n), v(n,i)**2
   enddo
  ! write(1113,*) i, zumba(i)
enddo
!=========================Valore attesa onsite=========================
!open(1115,file='onsite.dat')
open(1116,file='config-onsite.dat')
valzer=0
do i=1,nev
 
   write(1115,*) 'AUTOVALORE=',i
   write(1116,*) 'AUTOVALORE=',i
   do n=1,dim2
      valzer(i)=valzer(i)+v(n,i)**2*onsite(n)
      if(v(n,i)**2.ge.toll) write(1116,*) vecconfig_all(n), onsite(n), v(n,i)**2
   enddo
  ! write(1115,*) i, valzer(i)
enddo
!=========================Valore attesa sitenergy=========================
!open(1117,file='site-energy.dat')
open(1118,file='config-site-energy.dat')
tiptap=0
do i=1,nev
  
   write(1117,*) 'AUTOVALORE=',i
   write(1118,*) 'AUTOVALORE=',i
   do n=1,dim2
      tiptap(i)=tiptap(i)+v(n,i)**2*siten(n)
      if(v(n,i)**2.ge.toll) write(1118,*) vecconfig_all(n), siten(n), v(n,i)**2
   enddo
 !  write(1117,*) i, tiptap(i)
enddo



!=========================Valore attesa hopping=========================
open(1111,file='hopping.dat')
write(1113,*)'stato, potenziale, onsite, site-energy, cinetica'
do i=1,nev
   att=dd(i,1)-tiptap(i)-valzer(i)-zumba(i)
   write(1113,*) i, zumba(i), valzer(i), tiptap(i), att
enddo
!=========================operatore momento di dipolo========================
allocate(mux(dim2),muy(dim2), mu_x(nev), mu_y(nev))
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


   do i=1,13 
      mux(n)=mux(n)+xvec(i)*(nz(i)-occupazioni(i))
      muy(n)=muy(n)+yvec(i)*(nz(i)-occupazioni(i))
   enddo
enddo

!=========================momenti di dipolo=========================
mu_x=0.d0
mu_y=0.d0
do i=1,nev
   do j=1,dim2
      mu_x(i)=mu_x(i)+v(j,i)*v(j,1)*mux(j)
      mu_y(i)=mu_y(i)+v(j,i)*v(j,1)*muy(j)
   enddo
   write(1000,*)i, mu_x(i), mu_y(i)
enddo

!========================VALORE ATTESA DOPPIA=========================
do jj=1,nev
   write(111,*) 'STATO', jj
   dpg=0.d0
   do i=1,dim2

      dpg=dpg + (v(i,jj)**2)*double(i)
   enddo
   write(111,*) 'doppie media',dpg

   !=========================VALORE ATTESA CORRELAZIONE=========================


   do i=1,nsiti
      do p=1,nsiti
         do n=1,dim2
            if(i.ne.p)then
               correx(i,p)=correx(i,p)+vmat(i,p,n)*(v(n,jj)**2)
            else
               correx(i,p)=correx(i,p)+u(i)*double2(i,n)*(v(n,jj)**2)
            endif
         enddo
      enddo
   enddo

   do i=1,nsiti
      write(111,'(<nsiti>(2x,f10.5))')(correx(i,j),j=1,nsiti)
   enddo

   !=========================CORRELAZIONE=========================
   corrint=0
   ssite=0
   do i=1,nsiti
      do n=1,dim2
         ssite(i)=ssite(i)+(v(n,jj)**2)*occ(i,n)
      enddo
   enddo
   dsite=0
   do i=1,nsiti
      do p=1,nsiti
         do n=1,dim2
            dsite(i,p)=dsite(i,p)+(v(n,jj)**2)*(occ(i,n)*occ(p,n))
         enddo
      enddo
   enddo

   do i=1,nsiti
      do p=1,nsiti
         corrint(i,p)=dsite(i,p)-ssite(i)*ssite(p)
      enddo
   enddo


   write(111,*)
   do i=1,nsiti
      write(111,'(<nsiti>(2x,f10.5))')(corrint(i,j),j=1,nsiti)
   enddo
enddo

!=========================FORZA DI OSCILLATORE=========================
open(77,file='forza-oscillatore.dat')
hbar=1.0545718d-34 ![J*s]
hbarev=6.5821195d-16 ![eV]
echarge=1.6021766d-19 ![C]
emass=9.1093837d-32 ![kg]
do i=2,nev
   strenght=(4.*emass)/((echarge**2.)*hbar*3.)*((dd(1,i)-dd(1,1))/hbarev)*((mu_x(i)**2+mu_y(i)**2)*10.d-20*echarge**2)
   write(77,*) strenght, dd(i,1)-dd(1,1)
enddo

endprogram hoppingtotal



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


