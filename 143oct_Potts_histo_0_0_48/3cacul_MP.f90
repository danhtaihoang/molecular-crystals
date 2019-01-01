      
      PROGRAM test_open
      IMPLICIT NONE
      
      INTEGER (KIND=4), PARAMETER :: nF=5
      
      CHARACTER (LEN=150) :: tamp

      INTEGER (KIND=4) :: natx,naty,natz,i,j,nE
      REAL    (KIND=8) :: T

      REAL    (KIND=8) :: Cv,Ksi
      REAL    (KIND=8) :: T2,E_MP,M_MP,Cv_MP,Ksi_MP,E_2_MP,M_2_MP
          
      REAL    (KIND=8),DIMENSION(nF) :: E_moy,M_moy,E_2_moy,M_2_moy
      
      REAL    (KIND=8),DIMENSION(:,:),ALLOCATABLE  :: E,P
      REAL    (KIND=8),DIMENSION(:),ALLOCATABLE :: P_MP
      

 
!!!=================================================================================================
!!!=================================================================================================
!!!=================================================================================================

      
      OPEN(11,file='1parameter.in')
      
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I5))')    tamp, natx
      READ(11, '(A30,(I5))')    tamp, naty
      READ(11, '(A30,(I5))')    tamp, natz
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp,nE
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 


      ALLOCATE(E(nF,nE))
      ALLOCATE(P(nF,nE))
      ALLOCATE(P_MP(nE))
      


!!!=================================================================================================
!!!=================================================================================================




!!! Doc gia tri P tu cac file       
      DO i=1,nF

            IF (i==1) THEN
                  OPEN(unit=12,file='1/value_P_at_To.dat')
                  OPEN(unit=13,file='1/average_thermal.dat')
            END IF

            IF (i==2) THEN
                  OPEN(unit=12,file='2/value_P_at_To.dat')
                  OPEN(unit=13,file='2/average_thermal.dat')
            END IF
                             
            IF (i==3) THEN
                  OPEN(unit=12,file='3/value_P_at_To.dat')
                  OPEN(unit=13,file='3/average_thermal.dat')
            END IF 
            
             IF (i==4) THEN
                  OPEN(unit=12,file='4/value_P_at_To.dat')
                  OPEN(unit=13,file='4/average_thermal.dat')
            END IF

            IF (i==5) THEN
                  OPEN(unit=12,file='5/value_P_at_To.dat')
                  OPEN(unit=13,file='5/average_thermal.dat')
            END IF
                             
            IF (i==6) THEN
                  OPEN(unit=12,file='6/value_P_at_To.dat')
                  OPEN(unit=13,file='6/average_thermal.dat')
            END IF 
            
             IF (i==7) THEN
                  OPEN(unit=12,file='7/value_P_at_To.dat')
                  OPEN(unit=13,file='7/average_thermal.dat')
            END IF

            IF (i==8) THEN
                  OPEN(unit=12,file='8/value_P_at_To.dat')
                  OPEN(unit=13,file='8/average_thermal.dat')
            END IF
                             
            IF (i==9) THEN
                  OPEN(unit=12,file='9/value_P_at_To.dat')
                  OPEN(unit=13,file='9/average_thermal.dat')
            END IF 
            
            IF (i==10) THEN
                  OPEN(unit=12,file='10/value_P_at_To.dat')
                  OPEN(unit=13,file='10/average_thermal.dat')
            END IF 
            
            
            DO j=1,nE
            READ(12,*)T,E(i,j),P(i,j)
            !WRITE(*,*)T,E(i,j),P(i,j)
            END DO

            
            READ(13,*)T2,E_moy(i),M_moy(i),Cv,Ksi,E_2_moy(i),M_2_moy(i)

            
      END DO 

      CLOSE(12) 
                 

!!! Tinh gia tri trung binh P_MP

      OPEN(unit=21,file='P_MP.dat')  
       
      P_MP(:)=0.
          
      DO j=1,nE
            DO i=1,nF
            P_MP(j)=P_MP(j)+P(i,j)
            
            END DO
            
            P_MP(j)=P_MP(j)/real(nF)
            
            WRITE(21,*)T,E(1,j),P_MP(j)
      END DO      

      CLOSE(21) 

!!! Tinh gia tri trung binh E,M,Cv,Ksi

      OPEN(unit=22,file='average_thermal.dat')      
       
      E_MP=0.
      M_MP=0.
      E_2_MP=0.
      M_2_MP=0.
 
      DO i=1,nF
            E_MP=E_MP+E_moy(i)
            M_MP=M_MP+M_moy(i)
            E_2_MP=E_2_MP+E_2_moy(i)
            M_2_MP=M_2_MP+M_2_moy(i)         
      
      END DO
      E_MP=E_MP/real(nF)
      M_MP=M_MP/real(nF)
      E_2_MP=E_2_MP/real(nF)
      M_2_MP=M_2_MP/real(nF)
      Cv_MP=real(natx*naty*natz)*(E_2_MP-E_MP**2.)/(T2**2.)
      Ksi_MP=real(natx*naty*natz)*(M_2_MP-M_MP**2.)/T2

      WRITE(22,*)T2,E_MP,M_MP,Cv_MP,Ksi_MP

      CLOSE(22) 


!!! ================================================================================================



      
      END PROGRAM


