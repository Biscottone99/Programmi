program check
  implicit none
  integer::i,a,j,k,n
  integer,allocatable:: config(:)
  character:: array(26)
  logical::bool
  write(*,*) 'dammi n'
  read(*,*) n
  allocate(config(n))!,array(26))
  open(1,file='eigenvectors_sz0_7.dat')
  do i=1,n
     read(1,*) a,k
     config(i)=a
  enddo
  do i=1,n
     do j=0,25
        bool=btest(config(i),j)
        if(bool)then
           array(j+1)='1'
        else
           array(j+1)='0'
        endif
     enddo
     write(2,*) config(i)
     write(2,*)(array(j),j=25,1,-1)
     write(2,*) '---------'
    ! enddo
  enddo
           
end program check
