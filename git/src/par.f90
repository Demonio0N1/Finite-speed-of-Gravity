module parametros

    integer,parameter::ikind=selected_real_kind(16)
    integer, parameter :: DIM=3
    !real(kind=ikind), parameter :: CG=0.0002959867255922384281755837575197813572649020894500715428756906
    !real(kind=ikind), parameter :: CG=6988.191152154 !tierra luna
    real(kind=ikind), parameter :: CG=0.0001184462068986038q0 !unidades astronomicas-years
    real :: start, finish
    integer, parameter :: N=2
    
    
        
    !integer, parameter :: nt=1 !itearion number
    !integer, parameter :: nt=365000000*5
    !integer, parameter :: nt=3650
        
        
    real(kind=ikind), parameter :: Rc=0.03270138 !radio de los planetas
    real(kind=ikind) :: vp
    
        
    real(kind=ikind), parameter :: dt=0.000000273972602739726Q0 !itearation step
    !real(kind=ikind), parameter :: dt=0.0000000317098 !itearation step
    real(kind=ikind), dimension(3,N,2) :: r,v,a,rret,rret1 !!posicion,velocidad,aceleracion
    real(kind=ikind), dimension(N) :: m !!masa
    real(kind=ikind) :: az,min,acc2
    real(kind=ikind) :: rtx1, rtx2,x
    
    character*80 :: datin
    integer :: dat
    
    logical :: es_real   
    integer :: l,i,j,o,vpg,jj,stp,ir,ic,jc
    real(kind=ikind) :: sp,sk,st !!energia potencial, cinetica
    real(kind=ikind), dimension (N,3):: p,k
    integer :: nn,aux
    real(kind=ikind) :: eka,epa,eta,AC,BC,CC,DC,EC,ret3,rqt,rcb,rqd,ret1,ret2
    real(kind=ikind), dimension(DIM) :: Sij,Rij,Rji,Dij
    real(kind=ikind) :: Rsqij,phij,dphij,Rsqji,phji,dphji,rtc,Rsq
    real(kind=ikind), dimension(DIM) :: real_vel


    real(kind=ikind), dimension(3,N) :: orb
    real(kind=ikind), dimension(N) :: dist
    real(kind=ikind), dimension(2,N) :: acc1


    REAL (kind=ikind)    :: coef(0:2)
    !COMPLEX (kind=ikind) :: root(4)
    
    end module parametros