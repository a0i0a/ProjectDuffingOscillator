!Program details:
!Numerically simulates the Duffing Oscillator
!Outputs: 

program ProjectDuffingOscillator
    implicit none !all variables must be declared
    
    !Precision constants. Note some may not exist on different compilers/systems/platforms
    integer, parameter :: SGL = Selected_real_kind(p=6)    !Single Precision:    kind =  4 , Precision =  6 , Range =   37
    integer, parameter :: DGL = Selected_real_kind(p=15)   !Double Precision:    kind =  8 , Precision = 15 , Range =  307
    integer, parameter :: EGL = Selected_real_kind(p=18)   !Extended_D Precison: kind = 10 , Precision = 18 , Range = 4931
    integer, parameter :: QGL = Selected_real_kind(p=33)   !Quadriple Precision: kind = 16 , Precision = 33 , Range = 4931
    
    !Set precision to be used for REAL numbers throughout program
    integer, parameter :: PRECS = DGL; !Double Precison for all real variables
    
    !Simulation type
    integer, parameter  ::  sTY = 2  !1 for constant time step, 2 for adaptive time step
    
    !Data for Time Series, Phase Plot and Poincare sections, i.e [time: t, displacement: x1 and velocity: x2]
    logical, parameter  ::  gTXdata = .false. !if true, then get data [t, x1, x2] at the points in variable 'pT'
    logical, parameter  ::  gXdata = .false. !if true, then get data [x1, x2] at the points in variable 'pT' for poincare section
    integer, dimension(6) :: pT = (/1,1, 1,1, 2,2/)!Points (Inital,Final, I,F, I,F) at which to get [t, x1, x2] data
    !max values for pT == (/ Damping constant (1:101) | Forcing Amplitude (1:101) | RungeKutta Version (1:55) /)
    
    !Required Data (Initial conditions)
    real (KIND=PRECS), parameter ::  t_i = 0.0_PRECS  !Initial Time
    real (KIND=PRECS), parameter :: x2_i = 0.0_PRECS  !Initial Velocity
    
    real (KIND=PRECS) :: x1_i! = 0.0_PRECS  !Initial Displacement
    
    !Duffing equation parameters: (d^2(x)/dt^2) + ((dC)*d(x)/dt) + (a1*x) + (b1*(x^3)) = (Po*sin(w*t))
    real (KIND=PRECS), parameter  :: Po = 0.1_PRECS         !Amplitude of the periodic driving force (Forcing amplitude)
    real (KIND=PRECS), parameter  :: dC = 0.0168_PRECS      !Damping constant
    real (KIND=PRECS), parameter  :: a1 = -0.5_PRECS        !amount of linear stiffness
    real (KIND=PRECS), parameter  :: bt = 0.5_PRECS         !amount of non-linearity in the restoring force
    real (KIND=PRECS), parameter  :: Po_inc = 0.0011_PRECS, dC_inc = 0.001512_PRECS !increments of Po and dC
    integer                       :: numOfP_d = 101         !equal number of points for Po and dC (101  points), from 0 to 100
    integer                          idC, numOfd            !starting point and ending point of dC, goes from 1 to numOfP_d
    
    !Other required data
    real (KIND=PRECS), parameter  ::  nT = 500.0_PRECS        !number of excitation periods
    real (KIND=PRECS), parameter  ::  Eq = 1.0E-06_PRECS      !Required accuracy/tolerance
    real (KIND=PRECS), parameter  ::  Pi = 3.141592654_PRECS  !Value of Pi
    
    real (KIND=PRECS) ::   w! = 1.0_PRECS          !Excitation frequency (angular)
    real (KIND=PRECS) ::  T1! = 2.0_PRECS*Pi/w     !time for one excitation, i.e one period
    real (KIND=PRECS) :: T_t! = nT*T1              !Total simulation time
    real (KIND=PRECS) :: h_i! = T1/500.0_PRECS     !initial length of interval/time step
    real (KIND=PRECS) :: uST! = 50_PRECS*T1        !unsteady period of oscillation. Any period after is steady
    integer           :: nTS! = int(T_t/h_i)       !Total number of steps for constant time
    
    !Runge Kutta Data
    integer, parameter                     :: nOfV = 55      !Number of Runge Kutta Versions
    real (KIND=PRECS), dimension (nOfV,13) :: rC = 0         !Runge Kutta constants [13 constants, 55 versions]
    
    !Variables for Results      !No of Steady(successful,failed), No of Total Steps(successful,failed)
    integer, dimension (nOfV) :: nOfSs = 0, nOfSf = 0,            nOfTSs = 0, nOfTSf = 0   !delcare and initialize all elements
    !In a successful step, the error is less than the tolerance. While in failed step, the error is greater than the tolerance
    !The failed steps exists only in step halving.
    
    integer i,j,k, m,n, p,q !to be used by loops in the program
    real (KIND=PRECS) dC_1,Po_1,t,t_,t_2,x1,x1_,x1_1,x1_2,x2,x2_1,x2_2,h,h_,h_e,tt,xt,M1,M2,M3,M4,K1,K2,K3,K4,Er !Others
    logical cnd, cnd0, cnd1, cnd2 !others
    CHARACTER(len = 128) :: cmd, SN;
    
    !read data from Command line
    !Note that w, x1, idC, numOfd must be set from the command line    
    CALL get_command_argument(1, cmd); read (cmd, *)    w; !read the excitation frequency
    CALL get_command_argument(2, cmd); read (cmd, *) x1_i; !read the initial displacement
    CALL get_command_argument(3, cmd); read (cmd, *)  idC; !read the starting damping constant
    CALL get_command_argument(4, cmd); read (cmd, *) numOfd; !read the ending damping constant
    T1 = 2.0_PRECS*Pi/w;  T_t = nT*T1;  h_i = T1/500.0_PRECS;   uST = 50_PRECS*T1;  nTS = int(T_t/h_i); !calculate others
    
    !Input files
    open ( unit=20, file='input.txt' )                     !input file for runge kutta constants [13 constants, 55 versions]
    
    !Output files
    write (SN, fmt = "('w_',F4.1,'x_',F4.1,'s_',I1,'d_',I3)") w, x1_i, sTY, idC;
    open ( unit=21, file= 'results\otDatFl_'//trim(SN)//'.csv') !output file for no of steps for each runge kutta version [55 versions]
    if(gTXdata) open ( unit=22, file= 'results\otTXXFl_'//trim(SN)//'.csv') !output file for t, x1 and x2  [optional] for Time Series
    if(gXdata)  open ( unit=23, file= 'results\otXXFl_'//trim(SN)//'.csv')  !output file for t, x1 and x2  [optional] for Poincare Sections
    
    !1. Heading for file containing number of steps (unit=21)
    write (21, fmt = "(1X, A5, ',', 12(A5,',',E15.6,','))") 'Data:','t_i:',t_i,'x1_i:',x1_i,'x2_i:',x2_i,'w:',w,'T1:',T1,'nT1:', &
                    nT,'T_t:',T_t,'h:',h_i,'Sim:',float(sTY),'UnSt:',uST/T1,'Eq:',Eq,'PRES',float(PRECS);
    write (21, fmt = "( 1X, /, ', ,' ,55('nTSs',','), 55('nTSf',','), 55('nSs',','), 55('nSf',','), /, 1X, 57(A5,',') )")&
        'dC','Po','RK_1','RK_2','RK_3','RK_4','RK_5','RK_6','RK_7','RK_8','RK_9','RK_10','RK_11','RK_12','RK_13','RK_14','RK_15',&
         'RK_16','RK_17','RK_18','RK_19','RK_20','RK_21','RK_22','RK_23','RK_24','RK_25','RK_26','RK_27','RK_28','RK_29','RK_30',&
         'RK_31','RK_32','RK_33','RK_34','RK_35','RK_36','RK_37','RK_38','RK_39','RK_40','RK_41','RK_42','RK_43','RK_44','RK_45',&
         'RK_46','RK_47','RK_48','RK_49','RK_50','RK_51','RK_52','RK_53','RK_54','RK_55';
    
    !2. Heading for file containing [t,x1,x2] (unit=22)
    if(gTXdata) write (22, fmt = "(1X, A5, ',', 12(A5,',',E15.6,','))") 'ALL:','t_i:',t_i,'x1_i:',x1_i,'x2_i:',x2_i, &
                    'w:',w,'T1:',T1,'n_T1:',nT,'T_t:',T_t,'h:',h_i,'Sim:',float(sTY),'UnSt:',uST/T1,'Eq:',Eq,'PRES',float(PRECS);
    if(gXdata)  write (23, fmt = "(1X, A5, ',', 12(A5,',',E15.6,','))") 'ALL:','t_i:',t_i,'x1_i:',x1_i,'x2_i:',x2_i, &
                    'w:',w,'T1:',T1,'n_T1:',nT,'T_t:',T_t,'h:',h_i,'Sim:',float(sTY),'UnSt:',uST/T1,'Eq:',Eq,'PRES',float(PRECS);
    
    !Read the Runge Kutta constants from Input file
    do m = 1,nOfV
        read(20,*)(rC(m,n),n=1,13) !read 13 constants [using implied Do loop]
    end do !end do loop
    close (20); !since constants have been read, close the file.
    
    !Calculation starts
    print *,'Program starts...';
    
    DodC: do i = idC,numOfd !Run for each damping constant value, min(1) : max(101)
        dC_1 = dC + ((i-1)*dC_inc) !calculate for the current iteration
        
        DoPo: do j = 1,numOfP_d !Run for each amplitude value
            Po_1 = Po + ((j-1)*Po_inc); !calculate for the current iteration
            
            DoRK: do k= 1,1 !,nOfV !Run for each Runge Kutta version
                !print *,'Runge Kutta Version: ', k
                
                t = t_i; x1 = x1_i; x2 = x2_i; h = h_i; cnd1 = .false.; cnd2 = .false.; !initialize initial conditions and others
                nOfSf(k) = 0; nOfSs(k) = 0; nOfTSf(k) = 0; nOfTSs(k) = 0;!initialize steady and total steps counters.
                
                !2.1 Write heading for file (unit=22), only if all conditions are true
                cnd = (gTXdata .and. (i>=pT(1) .and. i<=pT(2)) .and. (j>=pT(3) .and. j<=pT(4)) .and. (k>=pT(5) .and. k<=pT(6)))
                if (cnd) write (22, fmt = "(1X,/, A5, ',', 2(A5,',',E15.6,',')   , (A5,I2,','), /, 5(A3,','))") &
                                                'Data',     'dC:',dC_1,'Po:',Po_1, 'RK_',k,        'S/N','D/N','t','x1','x2';
                
                !2.2 Write heading for file (unit=23), only if all conditions are true
                cnd0 = (gXdata .and. (i>=pT(1) .and. i<=pT(2)) .and. (j>=pT(3) .and. j<=pT(4)) .and. (k>=pT(5) .and. k<=pT(6)))
                if (cnd0) write (23, fmt = "(1X,/, A5, ',', 2(A5,',',E15.6,',')   , (A5,I2,','), /, 3(A3,','))") &
                                                'Data',     'dC:',dC_1,'Po:',Po_1, 'RK_',k,        'S/N','x1','x2';
                
                Endless: do ! [endless do] Run while total simulation time is greater than current time (t <= T_t) 
                    !print *, 'Step:', nOfTSs(k), '.....................';
                    tt =  t;                       xt = x1;                                                            !Step 1
                    M1 = x2;                                   K1 = Po_1*sin(w*tt)-(dC_1*M1)-(bt*(xt**3))-(a1*xt);
                    !print *, 't x1', tt, xt;   !print *, 'M1 K1', M1, K1;
                    
                    tt =  t + rC(k,1)*h;           xt = x1 + rC(k,4)*M1*h;                                             !Step 2
                    M2 = x2 + rC(k,4)*K1*h;                    K2 = Po_1*sin(w*tt)-(dC_1*M2)-(bt*(xt**3))-(a1*xt);
                    !print *, 't x1', tt, xt;   !print *, 'M2 K2', M2, K2;
                    
                    tt =  t + rC(k,2)*h;           xt = x1 + rC(k,5)*M1*h + rC(k,6)*M2*h;                              !Step 3
                    M3 = x2 + rC(k,5)*K1*h + rC(k,6)*K2*h;     K3 = Po_1*sin(w*tt)-(dC_1*M3)-(bt*(xt**3))-(a1*xt);
                    !print *, 't x1', tt, xt;   print *, 'M3 K3', M3, K3;
                    
                    tt =  t + rC(k,3)*h;           xt = x1 + rC(k,7)*M1*h + rC(k,8)*M2*h + rC(k,9)*M3*h;               !Step 4
                    M4 = x2 + rC(k,7)*K1*h + rC(k,8)*K2*h + rC(k,9)*K3*h;  
                                                               K4 = Po_1*sin(w*tt)-(dC_1*M4)-(bt*(xt**3))-(a1*xt);
                    !print *, 't x1', tt, xt;   print *, 'M4 K4', M4, K4;
                    
                    !compute new (final) values
                    x1_1 = x1 + ((M1*rC(k,10)) + (M2*rC(k,11)) + (M3*rC(k,12)) + (M4*rC(k,13)))*h
                    !print *, 'x1_1', x1_1;
                    
                    x2_1 = x2 + ((K1*rC(k,10)) + (K2*rC(k,11)) + (K3*rC(k,12)) + (K4*rC(k,13)))*h
                    !print *, 'x2_1', x2_1;
                    
                    if (sTY == 2) then !only if adaptive step size is to be used
                        t_ = t; x1_ = x1; M1 = x2; h_ = h/2 !initialize initial conditions
                        
                        StepHalving: do p = 1,2 !since we are using step halving, run twice.
                            !print *, 't x1', t_, x1_; print *, 'M1_1/2 K1_1/2', M1, K1;
                            !Step 1 is not needed here. Simply reuse values from the previous loop
                            
                            tt =  t_ + rC(k,1)*h_;           xt = x1_ + rC(k,4)*M1*h_;                                    !Step 2
                            M2 = M1 + rC(k,4)*K1*h_;                  K2 = Po_1*sin(w*tt)-(dC_1*M2)-(bt*(xt**3))-(a1*xt);
                            !print *, 'M2_1/2 K2_1/2', M2, K2;
                            
                            tt =  t_ + rC(k,2)*h_;           xt = x1_ + rC(k,5)*M1*h_ + rC(k,6)*M2*h_;                    !Step 3
                            M3 = M1 + rC(k,5)*K1*h_ + rC(k,6)*K2*h_;  K3 = Po_1*sin(w*tt)-(dC_1*M3)-(bt*(xt**3))-(a1*xt);
                            !print *, 'M3_1/2 K3_1/2', M3, K3;
                            
                            tt =  t_ + rC(k,3)*h_;           xt = x1_ + rC(k,7)*M1*h_ + rC(k,8)*M2*h_ + rC(k,9)*M3*h_;    !Step 4
                            M4 = M1 + rC(k,7)*K1*h_ + rC(k,8)*K2*h_ + rC(k,9)*K3*h_;
                                                                      K4 = Po_1*sin(w*tt)-(dC_1*M4)-(bt*(xt**3))-(a1*xt);
                            !print *, 'M4_1/2 K4_1/2', M4, K4;
                            
                            !compute new (final) values
                            t_2  =  t_ + h_;   !print *, 't_1/2', t_2;
                            
                            x1_2 = x1_ + ((M1*rC(k,10)) + (M2*rC(k,11)) + (M3*rC(k,12)) + (M4*rC(k,13)))*h_;
                            !print *, 'x1_1/2', x1_2;
                            
                            x2_2 = M1 + ((K1*rC(k,10)) + (K2*rC(k,11)) + (K3*rC(k,12)) + (K4*rC(k,13)))*h_;
                            !print *, 'x2_1/2', x2_2;
                            
                            !final values become initial values of next loop (needed once)
                            if(p==1) t_ = t_2; x1_ = x1_2; M1 = x2_2; K1 = Po_1*sin(w*t_)-(dC_1*M1)-(bt*(x1_**3))-(a1*x1_);!Step 1
                        end do StepHalving !end do for step halving
                        
                        !Determine the error and compute the next interval step size (adaptive step size only)
                        Er = abs(x2_2 - x2_1)
                        
                        if (Er > Eq) then !Error is greater than tolerance [reduce h and restart]
                            h = 0.95 * h * ((Eq/Er)**(0.25)) !decrease interval size
                            nOfTSf(k) = nOfTSf(k) + 1;             !increment the number of (failed) total steps
                            if (t > uST) nOfSf(k) = nOfSf(k) + 1;  !increment the number of (failed) steady steps
                            !if (cnd) write (22, fmt = "(1X, 2(I16, ','), 5(E15.6, ','))") 0,nOfTSf(k), t+h, x1_1, x2_1,x2_2, h;
                            cycle Endless !repeat/restart step with adjusted interval/step size
                        else !error is less than tolerance, [we continue]
                            h_e = 0.95 * h * ((Eq/Er)**(0.20)) !increase interval size, and continue
                        end if
                    end if !end if for adaptive steps only
                    
                    !increment steps
                    nOfTSs(k) = nOfTSs(k) +  1;             !increment the number of (successful) total steps
                    if (t > uST) nOfSs(k) = nOfSs(k) +  1;  !increment the number of (successful) steady steps
                    
                    !3. write data [t, x1, x2] to output file                      S/N, D/N, t, x1, x2 
                    if (cnd)  write (22, fmt = "(1X, 2(I16, ','), 3(E15.6, ','))") nOfTSs(k),nOfTSf(k), t, x1, x2;
                    !3.1. write data [t, x1, x2] to output file                      S/N, D/N, x1, x2 [poincare section]
                    !NOTE: Also, must use constant time step; Starts from steady steps; Also, if mod(A,B) = x, then 0 <= x < B;
                    if (sTY == 1 .and. (t > uST) .and. cnd0 .and. (mod(nOfSs(k),int(T1/h_i)) == 0.0)) & !pick one point per period
                                    write (23, fmt = "(1X, 1(I16, ','), 2(E15.6, ','))") nOfTSs(k),x1, x2;
                    
                    !final values become initial values of next loop
                    t = t + h; x1 = x1_1; x2 = x2_1;
                    if(sTY == 2) h = h_e; !step size for the next loop, if adaptive
                    !if (cnd) write (22, fmt = "(1X, 2(I16, ','), 3(E15.6, ','))") nOfTSs(k),nOfTSf(k), t, x1, x2; !...to remove
                    
                    !Tests
                    !1. Test for end of simulation: If [current time] is greater than [total simulation time], then exit (true)
                    !Hence the loop continues while total simulation time is greater than or equal to current time
                    if (t >= T_t) cnd1 = .true.; !Note that this would SOMETIMES result in an extra step
                    
                    !2. Test for divergence: If simulation is adaptive and If [current number of successful steps] is greater
                    !than [the number of steps it would have taken a constant time step program to run], then exit (true)
                    if (sTY == 2 .and. nOfTSs(k) > nTS) cnd2 = .true.;
                    
                    !If the result of any of the above tests are true, then exit.
                    if (cnd1 .or. cnd2) exit Endless;
                    
                end do Endless !end do (endless do)
                    
            end do DoRK!end do for Runge Kutta
            
            !write data to output files
            !4. No of steps for each Runge Kutta version: 'dC','Po','RK_1','RK_2', ... ,'RK_55'
            !OUTPUT FORMAT:     [dC, Po]  ,[RK_nOfSs 1-55],| |,[RK_nOfSd 1-55],| |,[RK_nOfTSs 1-55],| |,[RK_nOfTSd 1-55]
            write (21, fmt = "(1X, 2(E15.6,','), 220(I16,',') )" )  &    !, 57(I8,','), 55(I8,','), 55(I8,',')
                         dC_1,Po_1,                                 &    !Damping constant and Forcing amplitude
                        (nOfTSs(p),p=1,nOfV), (nOfTSf(q),q=1,nOfV), &    !total number of steps, successful and failed
                        (nOfSs(m) ,m=1,nOfV), (nOfSf(n) ,n=1,nOfV);      !number of steady steps, successful and failed
        end do DoPo!end do for Amplitude value
        
        !write to screen the completed percentage of the program numOfP_d
        write (*, fmt = "(/,5X,'..............',F7.3,'% complete')") (float(i)/float(numOfd))*100.0;
        
    end do DodC!end do for Damping constant
    
    close (21); close (22); !close opened files
    
    print *,'Program End';
    !print *,'Program End, Press enter to exit?'; read (*,*);
end program
