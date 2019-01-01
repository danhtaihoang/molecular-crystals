!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Laboratoire de Physique Théorique et Modélisation 
!     UMR 8089 CNRS-Université de Cergy-Pontoise
!     2, Avenue Adolphe Chauvin, 95032 Cergy-Pontoise Cedex, France
!----------------------------------------------------------------------------------------------------!
!     PROGRAMME: MONTE CARLO TRANSPORT OF SPIN ON MODERN POTTS
!     4/5/2011: sua lai, tinh r 1 lan de chay nhanh hon
!     CHUONG TRINH NAY DA CHUAN, DA TEST, THOI GIAN CHAY LA NHANH NHAT (chieu toi ngay 4/5/2011)
!     13/5/2011: sua lai GS cho he N=24, sua lai cach tinh Mz
!!    17/05/2011: sua lai config_GS cho toan bo gia tri r0 va d. (d=1,2,3,4,... )
!!    19/05/2011: sua thanh phan thu 2 cua dipolar, chi su dung cac lien ket song song
!!    01/06/2011: histogram
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  

      PROGRAM main_transport_modern_potts
      IMPLICIT NONE

      CHARACTER (LEN=150):: CONFIG_INI_DIREC_X
      CHARACTER (LEN=50) :: name
      CHARACTER (LEN=15) :: tmp
      CHARACTER (LEN=3)  :: spinAt
      CHARACTER (256)    :: Ligne20,Ligne126

      REAL    (KIND=8),PARAMETER :: aaa=5.
      REAL    (KIND=8),PARAMETER :: nul=0.

      INTEGER (KIND=4) :: iT,nT,i_case,d,i_layer,n_layer,m
      INTEGER (KIND=4) :: i,j,k,ip,jp,kp,im,jm,km,natx,naty,natz,natx_p,naty_p,natz_p,number

      INTEGER (KIND=4) :: i_loop,i_loop2,n_equi_reseau1,n_equi_reseau2,n_average_thermal,i_times,n_times

      INTEGER (KIND=4) :: deli0,delj0,delk0,r01,iv,n_vois,i1,j1,k1,k2,na


      REAL    (KIND=8) :: T,delT,Tmax,Tmin,rdn_etat,rdn_mtp,tab_tmp,rdn_config,n_atom
      REAL    (KIND=8) :: x,y,z,EN_tmp1,EN_tmp2
      REAL    (KIND=8) :: energy,energy_2,Mz,Mz_2
      REAL    (KIND=8) :: EN_moy,EN_2_moy,Mz_moy,Mz_2_moy,Cv,Ksi
      REAL    (KIND=8) :: r0,r1,D_tmp,EN_D_tmpx,EN_D_tmpy,EN_D_tmpz
      REAL    (KIND =8):: M12,M13,M23,M12x,M12y,M12z,M13x,M13y,M13z,M23x,M23y,M23z
      
         
      INTEGER (KIND=8),DIMENSION(3)                  :: clock
      
      REAL    (KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE  :: spin


      REAL    (KIND=8),DIMENSION(:),      ALLOCATABLE :: deli,delj,delk
      INTEGER (KIND=4),DIMENSION(:),      ALLOCATABLE :: int_deli,int_delj,int_delk
      REAL    (KIND=8),DIMENSION(:),      ALLOCATABLE :: r

       !!! Khai bao phan Histogram
      INTEGER (KIND=4) :: i_histo,n_EN_histo
      REAL    (KIND=8) :: EN_min_histo,EN_max_histo,del_EN_histo

      REAL    (KIND=8),DIMENSION(:),ALLOCATABLE :: EN_histo,H_histo,P_histo


!!!==================================================================================================
!!!==================================================================================================
      CALL system('rm -r config_ini_3D')
      CALL system('mkdir config_ini_3D')
      CALL system('rm -r config_3D')
      CALL system('mkdir config_3D')
      CALL system('rm *.dat*')


      CALL ini_rdm_number()

      CALL read_input_parameter_file()

      natx_p=natx+1
      naty_p=naty+1
      natz_p=natz+1

      n_atom=real(natx*naty*natz)
      na=natx*naty*natz
      r01=int(r0)

      n_layer=natz/(2*d)

      ALLOCATE(spin(natx_p,naty_p,natz_p,3))
      ALLOCATE(EN_histo(n_EN_histo))
      ALLOCATE(H_histo(n_EN_histo))
      ALLOCATE(P_histo(n_EN_histo))

      IF (nT==1) THEN
            delT=0.
      ELSE
            delT=(Tmax-Tmin)/real(nT-1)
      END IF


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM =====
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      OPEN(unit=20,file='average_thermal.dat')

      CALL load_config_ini()  

      CALL number_voisin()
      
      CALL energy_gs()      

      CALL write_config_ini_3D()


      CALL calcul_time_run()
!!!======================================  
!!!======================================  
      DO iT=1,nT
            WRITE(*,*)'iT = ', iT                        

            T=Tmin+delT*real(iT-1)


            DO i_times=1,n_times
                  CALL equi_reseau1()
                  CALL average_thermal()
            END DO

            !CALL value_thermal()
      
            CALL write_config_3D()


      END DO

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

      CLOSE(20)
     

      DEALLOCATE(spin)
     
      CONTAINS


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE init_rdm_number()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ini_rdm_number()
      IMPLICIT NONE

      INTEGER (KIND=4) :: i_time,i

      CALL ITIME(clock)
      i_time=(clock(1)+1)*(clock(2)+1)*(clock(3)+1)
        
      DO i=1,i_time
            CALL random_number(tab_tmp)
      ENDDO 

      END SUBROUTINE ini_rdm_number

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE read_input_parameter_file() 
!!! OPEN the parameter from file "parameter.in"
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE read_input_parameter_file()
      IMPLICIT NONE

      CHARACTER (LEN=150) :: tamp
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
      READ(11, '(A30,(A10))')   tamp, CONFIG_INI_DIREC_X
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
      READ(11, '(A30,(I5))')    tamp, d
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F7.4))')  tamp, D_tmp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F7.4))')  tamp, r0
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp, nT
      READ(11, '(A30,(F7.4))')  tamp, Tmin
      READ(11, '(A30,(F7.4))')  tamp, Tmax
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp,n_equi_reseau1
      READ(11, '(A30,(I8))')    tamp,n_equi_reseau2
      READ(11, '(A30,(I8))')    tamp,n_average_thermal
      READ(11, '(A30,(I8))')    tamp,n_times
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp,n_EN_histo
      READ(11, '(A30,(F7.4))')  tamp,EN_min_histo
      READ(11, '(A30,(F7.4))')  tamp,EN_max_histo
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 

      END SUBROUTINE read_input_parameter_file

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE calcul_time_run
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE calcul_time_run()

      IMPLICIT NONE      
      
      CHARACTER (8)  :: date
      CHARACTER (10) :: time
      CHARACTER (5)  :: zone
      INTEGER,DIMENSION(8) :: value_time1,value_time2,value_time3,value_time4

      INTEGER (KIND=4) :: n_date_run,n_hour_run,n_minute_run,i
      REAL    (KIND=8) :: n_time_run_total

      CALL system('rm time_run.dat')

      OPEN (90,file='time_run.dat')

!!! Tinh thoi gian moi lan thermalisation ==============================

      CALL date_and_time(date,time,zone,value_time1)
      WRITE(90,'(5I6)')value_time1(5),value_time1(6),value_time1(7),value_time1(8)
      
      DO i=1,1
      
      CALL equi_reseau()
      
      END DO
      
      CALL date_and_time(date,time,zone,value_time2)
      WRITE(90,'(5I6)')value_time2(5),value_time2(6),value_time2(7),value_time2(8)
!!! ===================================================================

!!! Tinh thoi gian moi lan tinh gia tri ==============================

      CALL date_and_time(date,time,zone,value_time3)
      WRITE(90,'(5I6)')value_time3(5),value_time3(6),value_time3(7),value_time3(8)
      
      DO i=1,1
      
      CALL value_thermal()

      END DO
      
      CALL date_and_time(date,time,zone,value_time4)
      WRITE(90,'(5I6)')value_time4(5),value_time4(6),value_time4(7),value_time4(8)
!!! =================================================================

      n_time_run_total = real(nT*(n_equi_reseau1+n_equi_reseau2*n_average_thermal))&
           *(real(value_time2(5)-value_time1(5))*60.+real(value_time2(6)-value_time1(6))&
            +real(value_time2(7)-value_time1(7))/60.+real(value_time2(8)-value_time1(8))/60000.)&
 
            +real(nT*n_average_thermal)&
           *(real(value_time4(5)-value_time3(5))*60.+real(value_time4(6)-value_time3(6))&
            +real(value_time4(7)-value_time3(7))/60.+real(value_time4(8)-value_time3(8))/60000.)

      n_date_run=int(n_time_run_total/1440.)
      n_hour_run=int(n_time_run_total/60.-n_date_run*24.)
      n_minute_run=int(n_time_run_total-n_date_run*1440.-n_hour_run*60.)

      WRITE(90,*)'n_time_run_total:',n_time_run_total,'mm'
      WRITE(90,*)'n_date_run      :',n_date_run
      WRITE(90,*)'n_hour_run      :',n_hour_run
      WRITE(90,*)'n_minute_run    :',n_minute_run

      CLOSE(90)

      END SUBROUTINE calcul_time_run

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Initial position configuration
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE load_config_ini()
      IMPLICIT NONE      

!!!====================================================  
!!! Hinh dang ban dau theo truc x

      IF (CONFIG_INI_DIREC_X == 'YES') THEN

      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx
                  spin(i,j,k,1)=1.
                  spin(i,j,k,2)=0.
                  spin(i,j,k,3)=0.

                  ENDDO
            ENDDO
      ENDDO

      END IF

!!!====================================================  
!!! Hinh dang ban dau bat ki

      IF (CONFIG_INI_DIREC_X == 'NO') THEN

      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx

                  CALL random_number(rdn_config)

                  IF (rdn_config<1./3.) THEN     
                        spin(i,j,k,1)=1.
                        spin(i,j,k,2)=0.
                        spin(i,j,k,3)=0.
                  ELSE
                        IF (rdn_config<2./3.) THEN
                        spin(i,j,k,1)=0.
                        spin(i,j,k,2)=1.
                        spin(i,j,k,3)=0.
                        
                        ELSE
                        spin(i,j,k,1)=0.
                        spin(i,j,k,2)=0.
                        spin(i,j,k,3)=1.
                        END IF

                  END IF

                  ENDDO
            ENDDO
      ENDDO

      END IF

!!!====================================================  
!!! Hinh dang ban dau theo GS
      IF (CONFIG_INI_DIREC_X == 'GS') THEN

      spin(:,:,:,:)=0.
      
      DO i_layer=1,n_layer
            DO m=1,d
                  k1 = (i_layer-1)*2*d+m
                  k2= (i_layer-1)*2*d+d+m
                  spin(:,:,k1,1)=1.
                  spin(:,:,k2,2)=1.

            END DO
      END DO

      END IF

   
      END SUBROUTINE load_config_ini


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_ini_3D()
      IMPLICIT NONE

          
      OPEN(unit=12,file='config_ini_3D/config_ini_3D_POTTS.pdb')


      DO k=1,natz 
            DO j=1,naty
                  DO i=1,natx

                        x=real(i-1)*aaa
                        y=real(j-1)*aaa
                        z=real(k-1)*aaa


                        IF (int(spin(i,j,k,1))==1) THEN 
                              spinAt='Au'
                              WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul

                        ELSE
                              
                              IF (int(spin(i,j,k,2))==1) THEN 
                              spinAt='Cu'
                              WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                  spinAt,x,y,z,nul

                              ELSE
                              
                              spinAt='Mn'
                              WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul

                              END IF
                        END IF

                       
                  END DO
            END DO
      END DO

      CLOSE(12)

      END SUBROUTINE write_config_ini_3D


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_3D()
      IMPLICIT NONE

      number=iT+10000000
      WRITE(tmp,'(I8)') number

      name='_POTTS_'//TRIM(tmp)
      
      OPEN(unit=13,file='config_3D/config_3D'//trim(name)//'.pdb')

      DO k=1,natz 
            DO j=1,naty
                  DO i=1,natx

                        x=real(i-1)*aaa
                        y=real(j-1)*aaa
                        z=real(k-1)*aaa

                        
                       IF (int(spin(i,j,k,1))==1) THEN 
                              spinAt='Au'
                              WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul

                        ELSE
                              
                              IF (int(spin(i,j,k,2))==1) THEN 
                              spinAt='Cu'
                              WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                  spinAt,x,y,z,nul

                              ELSE
                              
                              spinAt='Mn'
                              WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul

                              END IF
                        END IF

                        
                  END DO
            END DO
      END DO


      CLOSE(13)

      END SUBROUTINE write_config_3D


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_reseau1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_reseau1()
      IMPLICIT NONE

      DO i_loop=1,n_equi_reseau1
            CALL equi_reseau()
      END DO

      END SUBROUTINE equi_reseau1

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_reseau()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_reseau()
      IMPLICIT NONE

      DO k=1,natz
            kp=k+1-(k/natz)*natz
            km=k-1+(1/k)*natz
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                  jm=j-1+(1/j)*naty
                  DO i=1,natx
                        ip=i+1-(i/natx)*natx
                        im=i-1+(1/i)*natx
                       

                        i_case=int(spin(i,j,k,1))+2*int(spin(i,j,k,2))+3*int(spin(i,j,k,3))
                              
                        SELECT CASE (i_case)
!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                        !!! Truong hop spin nam doc theo truc x
                        !!=================================================================================
                        CASE(1)
                        
                        CALL random_number(rdn_etat)
                        IF (rdn_etat<0.5) THEN

                        CALL value_EN_D_tmp12()                        

                        EN_tmp1=-(spin(ip,j,k,1)+spin(im,j,k,1)+spin(i,jp,k,1)&
                                 +spin(i,jm,k,1)+spin(i,j,kp,1)+spin(i,j,km,1))+D_tmp*EN_D_tmpx
                       
                        EN_tmp2=-(spin(ip,j,k,2)+spin(im,j,k,2)+spin(i,jp,k,2)&
                                  +spin(i,jm,k,2)+spin(i,j,kp,2)+spin(i,j,km,2))+D_tmp*EN_D_tmpy

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,1)=0.
                                    spin(i,j,k,2)=1.
                                   
                              END IF

                        ELSE

                        CALL value_EN_D_tmp13()

                        EN_tmp1=-(spin(ip,j,k,1)+spin(im,j,k,1)+spin(i,jp,k,1)&
                                 +spin(i,jm,k,1)+spin(i,j,kp,1)+spin(i,j,km,1))+D_tmp*EN_D_tmpx        
                        
                        EN_tmp2=-(spin(ip,j,k,3)+spin(im,j,k,3)+spin(i,jp,k,3)&
                                 +spin(i,jm,k,3)+spin(i,j,kp,3)+spin(i,j,km,3))+D_tmp*EN_D_tmpz

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,1)=0.
                                    spin(i,j,k,3)=1.
                                    
                              END IF

                        END IF

                        !!! Truong hop spin nam doc theo truc y
                        !!=================================================================================
                        CASE(2)
                        CALL random_number(rdn_etat)

                        IF (rdn_etat<0.5) THEN
                        CALL value_EN_D_tmp12()

                        EN_tmp1=-(spin(ip,j,k,2)+spin(im,j,k,2)+spin(i,jp,k,2)&
                                 +spin(i,jm,k,2)+spin(i,j,kp,2)+spin(i,j,km,2))+D_tmp*EN_D_tmpy
                                               
                        EN_tmp2=-(spin(ip,j,k,1)+spin(im,j,k,1)+spin(i,jp,k,1)&
                                 +spin(i,jm,k,1)+spin(i,j,kp,1)+spin(i,j,km,1))+D_tmp*EN_D_tmpx

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,2)=0.
                                    spin(i,j,k,1)=1.
                                    
                              END IF

                        ELSE

                        CALL value_EN_D_tmp23()
                        EN_tmp1=-(spin(ip,j,k,2)+spin(im,j,k,2)+spin(i,jp,k,2)&
                                 +spin(i,jm,k,2)+spin(i,j,kp,2)+spin(i,j,km,2))+D_tmp*EN_D_tmpy

                        EN_tmp2=-(spin(ip,j,k,3)+spin(im,j,k,3)+spin(i,jp,k,3)&
                                 +spin(i,jm,k,3)+spin(i,j,kp,3)+spin(i,j,km,3))+D_tmp*EN_D_tmpz

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,2)=0.
                                    spin(i,j,k,3)=1.
                                   
                              END IF

                        END IF

                        !!! Truong hop spin nam doc theo truc z
                        !!=================================================================================
                        CASE(3)
                        CALL random_number(rdn_etat)
                        IF (rdn_etat<0.5) THEN
                        CALL value_EN_D_tmp13()

                        EN_tmp1=-(spin(ip,j,k,3)+spin(im,j,k,3)+spin(i,jp,k,3)&
                                 +spin(i,jm,k,3)+spin(i,j,kp,3)+spin(i,j,km,3))+D_tmp*EN_D_tmpz
                                              
                        EN_tmp2=-(spin(ip,j,k,1)+spin(im,j,k,1)+spin(i,jp,k,1)&
                                 +spin(i,jm,k,1)+spin(i,j,kp,1)+spin(i,j,km,1))+D_tmp*EN_D_tmpx

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,3)=0.
                                    spin(i,j,k,1)=1.
                                  
                              END IF

                        ELSE

                        CALL value_EN_D_tmp23()
                        EN_tmp1=-(spin(ip,j,k,3)+spin(im,j,k,3)+spin(i,jp,k,3)&
                                 +spin(i,jm,k,3)+spin(i,j,kp,3)+spin(i,j,km,3))+D_tmp*EN_D_tmpz
                        EN_tmp2=-(spin(ip,j,k,2)+spin(im,j,k,2)+spin(i,jp,k,2)&
                                 +spin(i,jm,k,2)+spin(i,j,kp,2)+spin(i,j,km,2))+D_tmp*EN_D_tmpy
                        

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,3)=0.
                                    spin(i,j,k,2)=1.
                                   
                              END IF

                        END IF

                        !!=================================================================================

                        END SELECT

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                          
                   ENDDO
             ENDDO
      ENDDO

        
      END SUBROUTINE equi_reseau

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_thermal()

      IMPLICIT NONE
          
      energy=0.
    
      DO k=1,natz
            kp=k+1-(k/natz)*natz
            km=k-1+(1/k)*natz
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                  jm=j-1+(1/j)*naty
                  DO i=1,natx
                        ip=i+1-(i/natx)*natx
                        im=i-1+(1/i)*natx

                        

                        i_case=int(spin(i,j,k,1))+2*int(spin(i,j,k,2))+3*int(spin(i,j,k,3))
                              
                        SELECT CASE (i_case)
!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                        !!! Truong hop spin nam doc theo truc x
                        !!=================================================================================
                        CASE(1)
                        
                        CALL random_number(rdn_etat)
                        IF (rdn_etat<0.5) THEN

                        CALL value_EN_D_tmp12()                        

                        EN_tmp1=-(spin(ip,j,k,1)+spin(im,j,k,1)+spin(i,jp,k,1)&
                                 +spin(i,jm,k,1)+spin(i,j,kp,1)+spin(i,j,km,1))+D_tmp*EN_D_tmpx
                       
                        EN_tmp2=-(spin(ip,j,k,2)+spin(im,j,k,2)+spin(i,jp,k,2)&
                                  +spin(i,jm,k,2)+spin(i,j,kp,2)+spin(i,j,km,2))+D_tmp*EN_D_tmpy

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,1)=0.
                                    spin(i,j,k,2)=1.
                                    EN_tmp1=EN_tmp2
                              END IF

                        ELSE

                        CALL value_EN_D_tmp13()

                        EN_tmp1=-(spin(ip,j,k,1)+spin(im,j,k,1)+spin(i,jp,k,1)&
                                 +spin(i,jm,k,1)+spin(i,j,kp,1)+spin(i,j,km,1))+D_tmp*EN_D_tmpx        
                        
                        EN_tmp2=-(spin(ip,j,k,3)+spin(im,j,k,3)+spin(i,jp,k,3)&
                                 +spin(i,jm,k,3)+spin(i,j,kp,3)+spin(i,j,km,3))+D_tmp*EN_D_tmpz

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,1)=0.
                                    spin(i,j,k,3)=1.
                                    EN_tmp1=EN_tmp2
                              END IF

                        END IF

                        !!! Truong hop spin nam doc theo truc y
                        !!=================================================================================
                        CASE(2)
                        CALL random_number(rdn_etat)

                        IF (rdn_etat<0.5) THEN
                        CALL value_EN_D_tmp12()

                        EN_tmp1=-(spin(ip,j,k,2)+spin(im,j,k,2)+spin(i,jp,k,2)&
                                 +spin(i,jm,k,2)+spin(i,j,kp,2)+spin(i,j,km,2))+D_tmp*EN_D_tmpy
                                               
                        EN_tmp2=-(spin(ip,j,k,1)+spin(im,j,k,1)+spin(i,jp,k,1)&
                                 +spin(i,jm,k,1)+spin(i,j,kp,1)+spin(i,j,km,1))+D_tmp*EN_D_tmpx

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,2)=0.
                                    spin(i,j,k,1)=1.
                                    EN_tmp1=EN_tmp2
                              END IF

                        ELSE

                        CALL value_EN_D_tmp23()
                        EN_tmp1=-(spin(ip,j,k,2)+spin(im,j,k,2)+spin(i,jp,k,2)&
                                 +spin(i,jm,k,2)+spin(i,j,kp,2)+spin(i,j,km,2))+D_tmp*EN_D_tmpy

                        EN_tmp2=-(spin(ip,j,k,3)+spin(im,j,k,3)+spin(i,jp,k,3)&
                                 +spin(i,jm,k,3)+spin(i,j,kp,3)+spin(i,j,km,3))+D_tmp*EN_D_tmpz

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,2)=0.
                                    spin(i,j,k,3)=1.
                                    EN_tmp1=EN_tmp2
                              END IF

                        END IF

                        !!! Truong hop spin nam doc theo truc z
                        !!=================================================================================
                        CASE(3)
                        CALL random_number(rdn_etat)
                        IF (rdn_etat<0.5) THEN
                        CALL value_EN_D_tmp13()

                        EN_tmp1=-(spin(ip,j,k,3)+spin(im,j,k,3)+spin(i,jp,k,3)&
                                 +spin(i,jm,k,3)+spin(i,j,kp,3)+spin(i,j,km,3))+D_tmp*EN_D_tmpz
                                              
                        EN_tmp2=-(spin(ip,j,k,1)+spin(im,j,k,1)+spin(i,jp,k,1)&
                                 +spin(i,jm,k,1)+spin(i,j,kp,1)+spin(i,j,km,1))+D_tmp*EN_D_tmpx

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,3)=0.
                                    spin(i,j,k,1)=1.
                                    EN_tmp1=EN_tmp2
                              END IF

                        ELSE

                        CALL value_EN_D_tmp23()
                        EN_tmp1=-(spin(ip,j,k,3)+spin(im,j,k,3)+spin(i,jp,k,3)&
                                 +spin(i,jm,k,3)+spin(i,j,kp,3)+spin(i,j,km,3))+D_tmp*EN_D_tmpz
                        EN_tmp2=-(spin(ip,j,k,2)+spin(im,j,k,2)+spin(i,jp,k,2)&
                                 +spin(i,jm,k,2)+spin(i,j,kp,2)+spin(i,j,km,2))+D_tmp*EN_D_tmpy
                        

                              CALL random_number(rdn_mtp)

                              IF (exp(-(EN_tmp2-EN_tmp1)/T) > rdn_mtp) THEN
                                    spin(i,j,k,3)=0.
                                    spin(i,j,k,2)=1.
                                    EN_tmp1=EN_tmp2
                              END IF

                        END IF

                        !!=================================================================================

                        END SELECT

                        energy=energy+EN_tmp1

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                          
                   ENDDO
             ENDDO
      ENDDO

      CALL value_Mz()

      energy=energy/(2.*n_atom)
      energy_2=energy**2.
      Mz_2=Mz**2.
     


      !WRITE(*,*)'EN=',energy,'Mz=',Mz,'Mz1=',Mz1,'Mz2=',Mz2,'Mz3=',Mz3,'Mz1+Mz2+Mz3=',Mz1+Mz2+Mz3

      END SUBROUTINE value_thermal

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_Mz()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_Mz()
      IMPLICIT NONE

      !!! Tinh theo cac mat phang xy
      M12=0.
      DO k= 1,natz

            M12x=0.; M12y=0.; M12z=0.

            DO j=1,naty
                  DO i=1,natx
                  M12x=M12x+spin(i,j,k,1)
                  M12y=M12y+spin(i,j,k,2)
                  M12z=M12z+spin(i,j,k,3)
                  END DO
            END DO

            M12=M12+abs(3.*max(M12x,M12y,M12z)/real(natx*naty)-1.)/2.
      END DO

      M12=M12/real(natz)

      !!! Tinh theo cac mat phang xz
      M13=0.
      DO j= 1,naty

            M13x=0.; M13y=0.; M13z=0.

            DO k=1,natz
                  DO i=1,natx
                  M13x=M13x+spin(i,j,k,1)
                  M13y=M13y+spin(i,j,k,2)
                  M13z=M13z+spin(i,j,k,3)
                  END DO
            END DO

            M13=M13+abs(3.*max(M13x,M13y,M13z)/real(natx*natz)-1.)/2.
      END DO

      M13=M13/real(naty)

      !!! Tinh theo cac mat phang yz
      M23=0.
      DO i= 1,natx

            M23x=0.; M23y=0.; M23z=0.

            DO k=1,natz
                  DO j=1,naty
                  M23x=M23x+spin(i,j,k,1)
                  M23y=M23y+spin(i,j,k,2)
                  M23z=M23z+spin(i,j,k,3)
                  END DO
            END DO

            M23=M23+abs(3.*max(M23x,M23y,M23z)/real(naty*natz)-1.)/2.
      END DO

      M23=M23/real(natx)

      Mz=max(M12,M13,M23)      

      END SUBROUTINE value_Mz
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE average_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE average_thermal()

      IMPLICIT NONE

      EN_moy=0.
      EN_2_moy=0.
      Mz_moy=0.
      Mz_2_moy=0.

      !!!------------------------------------------
      OPEN(unit=126,file='value_P_at_To.dat')

      H_histo(:)=0.
      EN_histo(:)=0.

      del_EN_histo = (EN_max_histo-EN_min_histo)/real(n_EN_histo-1)

      !! END khai bao bien cho histogram


      DO i_loop=1,n_average_thermal

            DO i_loop2=1,n_equi_reseau2
                  CALL equi_reseau()
            END DO

            CALL value_thermal()        
            
            EN_moy=EN_moy+energy
            EN_2_moy=EN_2_moy+energy_2
            Mz_moy=Mz_moy+Mz
            Mz_2_moy=Mz_2_moy+Mz_2

            IF ((EN_min_histo <= energy) .and. (energy <= EN_max_histo)) THEN              
                  i_histo=int((energy-EN_min_histo)/del_EN_histo)+1
                  H_histo(i_histo)=H_histo(i_histo)+1.            

            END IF      
                       
      END DO

      
      DO i_histo=1,n_EN_histo
            P_histo(i_histo)=H_histo(i_histo)/real(n_average_thermal)
            EN_histo(i_histo)=EN_min_histo+(real(i_histo)-0.5)*del_EN_histo

            WRITE(Ligne126,*)T,EN_histo(i_histo),P_histo(i_histo)
            WRITE(126,'(a)') trim(Ligne126)
      ENDDO     
 
      CLOSE(126)


      EN_moy=EN_moy/real(n_average_thermal)
      EN_2_moy=EN_2_moy/real(n_average_thermal)
      Mz_moy=Mz_moy/real(n_average_thermal)
      Mz_2_moy=Mz_2_moy/real(n_average_thermal)

      Cv = n_atom*(EN_2_moy-EN_moy**2.)/(T**2.)
      Ksi= n_atom*(Mz_2_moy-Mz_moy**2.)/T

      WRITE(Ligne20,*) T,EN_moy,Mz_moy,Cv,Ksi,EN_2_moy,Mz_2_moy

      WRITE(20,'(a)') trim(Ligne20)

      WRITE(*,*)'i_times=',i_times,'EN_moy=',EN_moy,'Mz_moy=',Mz_moy,'Cv=',Cv,'Ksi=',Ksi

      END SUBROUTINE average_thermal

      
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_EN_D_tmp12()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE value_EN_D_tmp12()
      IMPLICIT NONE

      EN_D_tmpx=0.
      EN_D_tmpy=0.

      DO iv=1,n_vois

            i1=i+int_deli(iv)
            j1=j+int_delj(iv)
            k1=k+int_delk(iv)

            IF (i1<1)    i1=natx+i1
            IF (i1>natx) i1=i1-natx
            IF (j1<1)    j1=naty+j1                        
            IF (j1>naty) j1=j1-naty                        
            IF (k1<1)    k1=natz+k1                       
            IF (k1>natz) k1=k1-natz

            EN_D_tmpx=EN_D_tmpx+spin(i1,j1,k1,1)*(1./r(iv)**3.-3.*(deli(iv)**2.)/(r(iv)**5.))   
            EN_D_tmpy=EN_D_tmpy+spin(i1,j1,k1,2)*(1./r(iv)**3.-3.*(delj(iv)**2.)/(r(iv)**5.))             


      END DO

      END SUBROUTINE value_EN_D_tmp12

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_EN_D_tmp13()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE value_EN_D_tmp13()
      IMPLICIT NONE


      EN_D_tmpx=0.
      EN_D_tmpz=0.

      DO iv=1,n_vois

            i1=i+int_deli(iv)
            j1=j+int_delj(iv)
            k1=k+int_delk(iv)

            IF (i1<1)    i1=natx+i1
            IF (i1>natx) i1=i1-natx
            IF (j1<1)    j1=naty+j1                        
            IF (j1>naty) j1=j1-naty                        
            IF (k1<1)    k1=natz+k1                       
            IF (k1>natz) k1=k1-natz

            EN_D_tmpx=EN_D_tmpx+spin(i1,j1,k1,1)*(1./r(iv)**3.-3.*(deli(iv)**2.)/(r(iv)**5.))
            EN_D_tmpz=EN_D_tmpz+spin(i1,j1,k1,3)*(1./r(iv)**3.-3.*(delk(iv)**2.)/(r(iv)**5.))

      END DO

      END SUBROUTINE value_EN_D_tmp13


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_EN_D_tmp23()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE value_EN_D_tmp23()
      IMPLICIT NONE


      EN_D_tmpy=0.
      EN_D_tmpz=0.

      DO iv=1,n_vois

            i1=i+int_deli(iv)
            j1=j+int_delj(iv)
            k1=k+int_delk(iv)

            IF (i1<1)    i1=natx+i1
            IF (i1>natx) i1=i1-natx
            IF (j1<1)    j1=naty+j1                        
            IF (j1>naty) j1=j1-naty                        
            IF (k1<1)    k1=natz+k1                       
            IF (k1>natz) k1=k1-natz

            EN_D_tmpy=EN_D_tmpy+spin(i1,j1,k1,2)*(1./r(iv)**3.-3.*(delj(iv)**2.)/(r(iv)**5.))  
            EN_D_tmpz=EN_D_tmpz+spin(i1,j1,k1,3)*(1./r(iv)**3.-3.*(delk(iv)**2.)/(r(iv)**5.))

      END DO


      !WRITE(*,*)'ENDx=',EN_D_tmpx,'ENDy=',EN_D_tmpy,'ENDz=',EN_D_tmpz


      END SUBROUTINE value_EN_D_tmp23

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_EN_D_tmp()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE value_EN_D_tmp()
      IMPLICIT NONE


      EN_D_tmpx=0.
      EN_D_tmpy=0.
      EN_D_tmpz=0.

      DO iv=1,n_vois

            i1=i+int_deli(iv)
            j1=j+int_delj(iv)
            k1=k+int_delk(iv)

            IF (i1<1)    i1=natx+i1
            IF (i1>natx) i1=i1-natx
            IF (j1<1)    j1=naty+j1                        
            IF (j1>naty) j1=j1-naty                        
            IF (k1<1)    k1=natz+k1                       
            IF (k1>natz) k1=k1-natz

            EN_D_tmpx=EN_D_tmpx+spin(i1,j1,k1,1)*(1./r(iv)**3.-3.*(deli(iv)**2.)/(r(iv)**5.))
            EN_D_tmpy=EN_D_tmpy+spin(i1,j1,k1,2)*(1./r(iv)**3.-3.*(delj(iv)**2.)/(r(iv)**5.))  
            EN_D_tmpz=EN_D_tmpz+spin(i1,j1,k1,3)*(1./r(iv)**3.-3.*(delk(iv)**2.)/(r(iv)**5.))

      END DO

      END SUBROUTINE value_EN_D_tmp


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE number_voisin()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE number_voisin()

      IMPLICIT NONE
!!!! Tinh so voisin (n_vois) cua moi atom, de su dung cho khai bao cac tab
      iv=0
                
      DO deli0=-r01,r01
      DO delj0=-r01,r01
      DO delk0=-r01,r01

            r1=sqrt(real(deli0)**2.+real(delj0)**2.+real(delk0)**2.)
                        
            IF ((0.<r1).and.(r1<=r0)) THEN

                  iv=iv+1
            END IF

      END DO
      END DO
      END DO

      n_vois=iv
  

!!!===============================================================================

!!! Tinh deli, delj, delk, ban kinh r (co gia tri bang nhau cho cac nguyen tu)

      iv=0
      ALLOCATE(deli(n_vois),delj(n_vois),delk(n_vois))
      ALLOCATE(int_deli(n_vois),int_delj(n_vois),int_delk(n_vois))
      ALLOCATE(r(n_vois))

      OPEN(14,file='voisin.dat')
      
      DO deli0=-r01,r01
      DO delj0=-r01,r01
      DO delk0=-r01,r01

            r1=sqrt(real(deli0)**2.+real(delj0)**2.+real(delk0)**2.)
                        
            IF ((0.<r1).and.(r1<=r0)) THEN

                  iv=iv+1

                  deli(iv)=real(deli0)
                  delj(iv)=real(delj0)
                  delk(iv)=real(delk0)
                  r(iv)=r1

                  int_deli(iv)=deli0
                  int_delj(iv)=delj0
                  int_delk(iv)=delk0
      
            WRITE(14,*)iv,deli(iv),delj(iv),delk(iv),r(iv)

            END IF

      END DO
      END DO
      END DO

      WRITE(*,*)'Xong buoc tinh del,r,i1,j1,k1'
      
      CLOSE(14)

      END SUBROUTINE number_voisin

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE energy_gs()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE energy_gs()
      IMPLICIT NONE

      REAL (KIND=8) :: EN_gs
               
      DO k=1,natz
            kp=k+1-(k/natz)*natz
            km=k-1+(1/k)*natz
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                  jm=j-1+(1/j)*naty
                  DO i=1,natx
                        ip=i+1-(i/natx)*natx
                        im=i-1+(1/i)*natx

                        CALL value_EN_D_tmp12()

                        i_case=int(spin(i,j,k,1))+2*int(spin(i,j,k,2))
                              
                        SELECT CASE (i_case)

                        CASE(1)                   
                             
                        EN_tmp1=EN_tmp1+spin(i,j,k,1)*&
                                    (-(spin(ip,j,k,1)+spin(im,j,k,1)+spin(i,jp,k,1)&
                                 +spin(i,jm,k,1)+spin(i,j,kp,1)+spin(i,j,km,1))+D_tmp*EN_D_tmpx)
                       
                        CASE(2)

                        EN_tmp2=EN_tmp2+spin(i,j,k,2)*&
                                    (-(spin(ip,j,k,2)+spin(im,j,k,2)+spin(i,jp,k,2)&
                                  +spin(i,jm,k,2)+spin(i,j,kp,2)+spin(i,j,km,2))+D_tmp*EN_D_tmpy)

                        END SELECT                   

                  END DO
            END DO

      END DO

      EN_gs = (EN_tmp1+EN_tmp2)/(2.*n_atom)

      WRITE(*,*)'EN_gs=',EN_gs

      END SUBROUTINE energy_gs


      END PROGRAM main_transport_modern_potts
      
