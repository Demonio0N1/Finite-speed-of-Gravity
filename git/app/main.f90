PROGRAM ver
  use parametros
  
  IMPLICIT NONE
 
    
  integer,parameter :: k15 = selected_int_kind(32)
  integer(kind=k15) :: step,nt


  call cpu_time(start)

  nt=3650000q0 ! ##numero de pasos termpora/

  read(*,'(a)') datin
  dat = len_trim(datin)
r=0.0q0
v=0.0q0
a=0.0q0
!root=0.0q0

  !m(1)=333000
  ! m(2)=0.0553
  ! m(3)=0.815
  ! m(4)=1
  ! m(5)=0.1074
  ! m(6)=317.8
  ! m(7)=95.2
  ! m(8)=14.54
  ! m(9)=17.15

     m(1)=1.0q0
      m(2)=0.0123q0
  

! if(N.ge.3)then
!   do i=1,N
!     call random_seed()
!         call random_number(az)
!     m(i)=az
!     !write(*,*) m(i)/10
!   end do
! end if
!  m(1)=10000
! v(:,1,:)=0.d0
!  m(11)=10000
! v(:,11,:)=0.d0

    ! if(N.ge.3)then
    !   do i=1,N
    !     m(i)=1
    !   end do
    ! end if
    ! m(1)=333000

!vp=999999999
!vp=100000000
!vp=10000
!vp=1
!vp=0.0001
vp=0.00001
!vp=10000000
!vp=10000
!vp=100
!vp=10
!vp=1
!vp=0.1q0
!vp=0.01
!vp=0.0017q0
!vp=0.001q0
!vp=0.0005
!vp=0.00032
!vp=0.0003
!vp=0.0001
!vp=0.00001
!vp=0.000001
!vp=0.0000001
!vp=0.00000001
!vp=0.000000001
!vp=0.0000000001
!vp=0.00000000001
!vp=0.000000000001
!vp=0.0000000000001
!vp=0.00000000000001
!vp=0.0000000000000000001
!vp=0

vp=vp*63242

 
  open(unit=1,file=datin(1:dat),status='unknown')
     do l=1,N
      read (1,*) r(:,l,1),v(:,l,1),a(:,l,1)
     end do

     coef=0
      a=0
      aux=0
      step=1
      
      open(file='rac.dat',unit=21,status='unknown')
      open(file='dacc.dat',unit=14,status='unknown')
    open(file='ener.dat',unit=77,status='unknown')
    open(file='dat.dat',unit=7,status='unknown')
    open(file='orb.dat',unit=9,status='unknown')
    open(file='dist.dat',unit=13,status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        

    do  i = 1,N-1
      do concurrent (j = i+1:N)
  
    
        
        coef(2)= dot_product(v(:,j,1),v(:,j,1))-vp**2
        coef(1)= 2*dot_product(r(:,i,1)-r(:,j,1),v(:,j,1))
        coef(0)= dot_product(r(:,i,1)-r(:,j,1),r(:,i,1)-r(:,j,1))       
        ret1 = (-2*coef(1)-sqrt(2*coef(1)**2 - 4*(coef(2)*coef(0))))/(2*coef(2))

        write(*,*) ret1
        coef(2)= dot_product(v(:,i,1),v(:,i,1))-vp**2
        coef(1)= 2*dot_product(r(:,j,1)-r(:,i,1),v(:,i,1))
        coef(0)= dot_product(r(:,j,1)-r(:,i,1),r(:,j,1)-r(:,i,1))
        ret2 = (-2*coef(1)-sqrt(2*coef(1)**2 - 4*(coef(2)*coef(0))))/(2*coef(2))


       
        rret(:,i,1) = r(:,i,1) - (ret1)*v(:,i,1) 
        rret(:,j,1) = r(:,j,1) - (ret2)*v(:,j,1) 

        

        Rij = r(:,i,1) - rret(:,j,1)
        Rji = r(:,j,1) - rret(:,i,1)

        Rsqij= dot_product(Rij ,Rij)
        Rsqji= dot_product(Rji ,Rji)

        phij  = -(m(j)*CG*m(i))/sqrt(Rsqij)
        dphij = -(m(j)*CG)/((sqrt(Rsqij))**3)

        phji  = -(m(j)*CG*m(i))/sqrt(Rsqji)
        dphji = -(CG*m(i))/((sqrt(Rsqji))**3)

        p(i,1) = p(i,1) + 0.5q0*phij
        p(j,1) = p(j,1) + 0.5q0*phji

        a(:,i,1) = a(:,i,1) + dphij*Rij
        a(:,j,1) = a(:,j,1) + dphji*Rji
  
      enddo
    enddo

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        do step=1,nt
    
          rret1=r
          v(:,:,1) = v(:,:,1) + 0.5q0*dt*a(:,:,1)
          r(:,:,1) = r(:,:,1) + dt*v(:,:,1) 
            a(:,:,2)=0.0q0
            eka=0
            epa=0
            p=0
            eta=0
            
            do  i = 1,N-1
              do concurrent (j = i+1:N)

                coef(2)= dot_product(v(:,j,1),v(:,j,1))-vp**2
                coef(1)= 2*dot_product(r(:,i,1)-r(:,j,1),v(:,j,1))
                coef(0)= dot_product(r(:,i,1)-r(:,j,1),r(:,i,1)-r(:,j,1))
                ret1 = (-2*coef(1)-sqrt(2*coef(1)**2 - 4*(coef(2)*coef(0))))/(2*coef(2))
                               
                coef(2)= dot_product(v(:,i,1),v(:,i,1))-vp**2
                coef(1)= 2*dot_product(r(:,j,1)-r(:,i,1),v(:,i,1))
                coef(0)= dot_product(r(:,j,1)-r(:,i,1),r(:,j,1)-r(:,i,1))
                ret2 = (-2*coef(1)-sqrt(2*coef(1)**2 - 4*(coef(2)*coef(0))))/(2*coef(2))
                          
                rret(:,i,1) = r(:,i,1)-ret1*v(:,i,1) 
                rret(:,j,1) = r(:,j,1)-ret2*v(:,j,1) 

                Rij = r(:,i,1) - rret(:,j,1)
                Rji = r(:,j,1) - rret(:,i,1)

            

                Rsqij= dot_product(Rij ,Rij)
                Rsqji= dot_product(Rji ,Rji)

                phij  = -(m(j)*CG*m(i))/sqrt(Rsqij)
                dphij = -(m(j)*CG)/((sqrt(Rsqij))**3)

                phji  = -(m(j)*CG*m(i))/sqrt(Rsqji)
                dphji = -(CG*m(i))/((sqrt(Rsqji))**3)

                p(i,1) = p(i,1) + 0.5q0*phij
                p(j,1) = p(j,1) + 0.5q0*phji

                a(:,i,2) = a(:,i,2) + dphij*Rij
                a(:,j,2) = a(:,j,2) + dphji*Rji
          
              enddo
            enddo

            v(:,:,1) = v(:,:,1) + 0.5q0*dt*(a(:,:,2))
           
            do i=1,N
              real_vel = v(:,i,1)
              eka = eka+0.5q0*m(i)*dot_product(real_vel,real_vel)
              epa=epa+p(i,1)
            enddo
  
            eta = eka+epa
          
            aux=aux+1

          if(aux.eq.1000)then !!!!!!aqui se cambia el numero datos que se imprimen
          
            
            

            do nn=1,N
              write(7,*) r(1,nn,1),r(2,nn,1),r(3,nn,1)
            end do
            write(7,*) ' '
            write(7,*) ' ' 

            do i=1,N
              dist(i)=sqrt(dot_product(r(:,i,1)-r(:,1,1), r(:,i,1)-r(:,1,1)))
              orb(:,i)= r(:,i,1)-r(:,1,1)
              acc1(1,i)=sqrt(dot_product((a(:,i,2)-a(:,i,1)),(a(:,i,2)-a(:,i,1)))) / sqrt(dot_product(a(:,i,1),a(:,i,1)))
              acc1(2,i)=sqrt(dot_product((v(:,i,2)-v(:,i,1)),(v(:,i,2)-v(:,i,1)))) / sqrt(dot_product(v(:,i,1),v(:,i,1)))
            enddo

            write(13,*) dist
            write(9,*)  orb
            write(77,*) step*dt, eka,epa,eta
            write(14,*) step*dt, acc1
            aux=0
          endif

          a(:,:,1)=a(:,:,2)
          v(:,:,2)=v(:,:,1)
          
        enddo
    close(13)
    close (9)
    close (7)
    close (77)
    close(14)
    close(21)
  close(1)
 

  call cpu_time(finish)
  write(*,*) finish-start
    

  END PROGRAM ver