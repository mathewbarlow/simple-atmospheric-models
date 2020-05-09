CC
CC  Coded by Michelle Weatherwax
CC
C     05-12-08
C     Fortran 77 program
C     Baroclinic Instablility 2-layer model
C     Solve Equations 8.16 & 8.17 in Holton
C     Status: Model runs for 24 hours
C             Advects to west when no mean wind
C             Advects to east when U1 = U3
C             Runs when U1 does not equal U3



C      The model is global in longitude
C      Standard SI units
C      Ignore curvature terms
C
C      Use Center Difference in Time w/ single Forward Difference step
C          every 10 time steps


C       Definitions
        REAL rpi, re, clat, dt, dx, dP
        INTEGER nx, nt
        PARAMETER (rpi = 3.14159) ! Set Pi
        PARAMETER (re = 6378100 ) ! Radius of the Earth (meters)
        PARAMETER (clat = 45) ! Latitude for the model (degrees)
        PARAMETER (nx = 144) ! Number of grid spaces
        PARAMETER (nk = nx/2 - 1)
        PARAMETER (dt = 60*60) ! Time Step of 3600 seconds = 1 hr
        PARAMETER (nt = 24*31) ! Time of run is 24 hours * X days
        PARAMETER (dP = 500) ! Pressure difference b/w 2 layers in mb
        
        REAL f, Beta, lambda, sigma, delta_P
        REAL U1(nx), U3(nx), Um(nx), Ut(nx)
        REAL vort_1(nx, nt + 1), vort_3(nx, nt + 1)
        REAL vort_m(nx, nt + 1), vort_t(nx, nt + 1)
        REAL psi_1(nx, nt + 1), psi_3(nx, nt + 1), w_2(nx, nt + 1)
        REAL psi_m(nx, nt + 1), psi_t(nx, nt + 1)
        REAL v_m(nx, nt + 1), v_t(nx, nt + 1)
        REAL advect_mean(nx), advect_combo(nx)
        REAL combo(nx, nt + 1)
        REAL w2(nx, nt + 1)
C
C
C      Set up input and output files
       OPEN(11, file = "two_layer.dat", status = "unknown"
     &  , access = "direct", recl = 4)
C
C
C      Set initial conditions
       it = 0 ! Initial time step
       dx = ((2*rpi*re) * COS(clat*rpi/180))/nx ! Grid spacing in m
       f = 2.0 * 7.292E-5 * SIN(clat*rpi/180.0) ! Coriolis term
       Beta = (2.0 * 7.292E-5 * COS(clat*rpi/180.0))/re ! df/dy term
       sigma = 2.0E-6 ! Potential Temp, Given pg. 153 for midtroposhere, units m^2 Pa^-2 s^-2
       delta_P = dP * 100 ! Convert mb to Pa
       lambda = f / ((sigma**0.5) * delta_P) ! Equ. given on pg. 234 Holton
	
       DO ix = 1, nx
          U1(ix) = 15 ! Assume mean wind of 20 m/s around circ. at level 1
          U3(ix) = 5 ! Assume mean wind of 10 m/s around circ. at level 3
C
C  random barotropic initial conditions
C
C         psi_1(ix,it+1)=rand(0)-0.5
C         psi_3(ix,it+1)=psi_1(ix,it+1)
C
C  random baroclinic initital conditions
C
         psi_1(ix,it+1)=rand(0)-0.5
         psi_3(ix,it+1)=rand(0)-0.5
C
C  gaussian initial condition at upper level
C
C         psi_1(ix,it+1)=exp(-real(ix-72)**2/144.0)
C	 psi_3(ix,it+1)=0

       ENDDO
C       CALL Smooth_random_psi(nx, nk, rpi, it, nt, re, clat, psi_1)
C       CALL Smooth_random_psi(nx, nk, rpi, it, nt, re, clat, psi_3)
C      Calculate wind field
       CALL calc_Um(nx, U1, U3, Um)
       CALL calc_Ut(nx, U1, U3, Ut)
       write(6,*) 'U_m',Um(1),' U_t',Ut(1)
C      Calculate vorticity field
       CALL Fourier_calc_vort(nx, nk, rpi, it, nt, f, psi_1, re,
     &                                       clat, vort_1)
       CALL Fourier_calc_vort(nx, nk, rpi, it, nt, f, psi_3, re,
     &                                       clat, vort_3)
C      Calculate mean and thermal height and vorticity fields
       CALL calc_psi_m(nx, nt, it, psi_1, psi_3, psi_m)
       CALL calc_psi_t(nx, nt, it, psi_1, psi_3, psi_t)

       CALL Fourier_calc_vort(nx, nk, rpi, it, nt, f, psi_m, re,
     &                                       clat, vort_m)
       CALL Fourier_calc_vort(nx, nk, rpi, it, nt, f, psi_t, re,
     &                                       clat, vort_t)
C      Calculate velocity v
       CALL calc_v(it, nt, f, nx, psi_m, dx, v_m)
       CALL calc_v(it, nt, f, nx, psi_t, dx, v_t)
       
C      Initialize "combo"
       DO ix = 1, nx
          combo(ix, it+1) = vort_t(ix,it+1)-
     &                     ((2*(lambda**2))*psi_t(ix,it+1))
       ENDDO

       write(6,*) 'done with time step',it
C      **All initial conditions are now set and/or calculated**


       
C      Now start to project forward in time
C      First time step, and every 10th step thereafter uses FORWARD DIFFERENCING
C      All other time steps use CENTER DIFFERENCING
C
       DO it = 1, nt ! Run for nt hours in increments of 1 hour
C         Put all non-time-dependent quantities from Eqs. 8.16 and 8.17
C             on the RHS and solve for the value
          CALL calc_advect_mean(nx, dx, it, nt, dt, Beta, Um, Ut,
     &              vort_m, vort_t, psi_m, psi_t, v_m, v_t, advect_mean)
          CALL calc_advect_combo(nx, dx, it, nt, dt, Beta, lambda, Um,
     &         Ut, vort_m, vort_t, psi_m, psi_t, v_m, v_t, advect_combo)
          IF (mod(it,10) .NE. 1) THEN ! CD for time is NOT a mult. of 10
             DO ix = 1, nx ! Advect the quanities solved for above forward in time
                vort_m(ix, it+1) = vort_m(ix, it-1) +
     &                                         (2*dt*advect_mean(ix))
                combo(ix, it+1) = combo(ix, it-1) +
     &                                         (2*dt*advect_combo(ix))
             ENDDO
          ELSE ! FD for time IS a mult. of 10
             DO ix = 1, nx ! Advect the quanities solved for above forward in time
                vort_m(ix, it+1) = vort_m(ix, it) +
     &                                         (dt*advect_mean(ix))
                combo(ix, it+1) = combo(ix, it) +
     &                                         (dt*advect_combo(ix))
             ENDDO
          ENDIF
       CALL Smooth_random_psi(nx, nk, rpi, it, nt, re, clat, vort_m)
       CALL Smooth_random_psi(nx, nk, rpi, it, nt, re, clat, combo)

C         Calculate the height field
          CALL Fourier_psi_combo(nx, nk, rpi, it, nt, f, combo, re,
     &                                       lambda, clat, psi_t)
          CALL Fourier_calc_psi(nx, nk, rpi, it, nt, f, vort_m, re,
     &                                       clat, psi_m)
          CALL Fourier_calc_vort(nx, nk, rpi, it, nt, f, psi_t, re,
     &                                       clat, vort_t)
          CALL calc_psi_1(nx, nt, it, psi_m, psi_t, psi_1)
          CALL calc_psi_3(nx, nt, it, psi_m, psi_t, psi_3)
          CALL Fourier_calc_vort(nx, nk, rpi, it, nt, f, psi_1, re,
     &                                       clat, vort_1)
          CALL Fourier_calc_vort(nx, nk, rpi, it, nt, f, psi_3, re,
     &                                       clat, vort_3)
C         Once we have phi, we can calculate v (use center diff method)
          CALL calc_v(it, nt, f, nx, psi_m, dx, v_m)
          CALL calc_v(it, nt, f, nx, psi_t, dx, v_t)
          CALL calc_w2(nx, it, nt, f, sigma, delta_P, psi_m, psi_t,
     &                                       dt, dx, Um, Ut, w2)

C          write(6,*) 'psi_1',psi_1(6,it+1),' psi_3',psi_3(6,it+1)
C          write(6,*) 'vort_1',vort_1(6,it+1),' vort_3',vort_3(6,it+1)
          write(6,*) 'done with time step',it
C
       ENDDO

C      Write out forecasts to file
       irec = 1
       DO it = 1, nt + 1
          DO ix = 1, nx
             write(11, rec = irec)vort_1(ix,it)
             irec = irec + 1
          ENDDO
          DO ix = 1, nx
             write(11, rec = irec)vort_3(ix,it)
             irec = irec + 1
          ENDDO
          DO ix = 1, nx
             write(11, rec = irec)psi_1(ix,it)
             irec = irec + 1
          ENDDO
          DO ix = 1, nx
             write(11, rec = irec)psi_3(ix,it)
             irec = irec + 1
          ENDDO
          DO ix = 1, nx
             write(11, rec = irec)w2(ix,it)
             irec = irec + 1
          ENDDO
       ENDDO

       STOP
       END
        
        
C      *************************************************************

       SUBROUTINE calc_Um(nx, U1, U3, Um)
C      Calculated Um (Mean Wind) based on U1 and U3
C      Variable "Um" is returned; no other variables modified
       INTEGER nx
       REAL Um(nx), U1(nx), U3(nx)
       
       DO ix = 1, nx
          Um(ix) = (U1(ix) + U3(ix))/2.0
       ENDDO
       
       RETURN
       END
       
       
C      *************************************************************

       SUBROUTINE calc_Ut(nx, U1, U3, Ut)
C      Calculated Ut (Thermal Wind) based on U1 and U3
C      Variable "Ut" is returned; no other variables modified
       INTEGER nx
       REAL Ut(nx), U1(nx), U3(nx)

       DO ix = 1, nx
          Ut(ix) = (U1(ix) - U3(ix))/2.0
       ENDDO

       RETURN
       END


C      *************************************************************

       SUBROUTINE Fourier_calc_psi(nx, nk, rpi, it, nt, f, vort, re,
     &                                       clat, psi)
C      Use the Fourier Transform to calculate the height field
C      Variable "psi" is returned; no other variables modified
       INTEGER nx, it, nt, nk
       REAL vort(nx, nt + 1), psi(nx, nt + 1), f, re, clat, rpi

       REAL k_const
       REAL a(nk+1), b(nk), p(nk+1), q(nk)
       REAL cos_term, sin_term
C      Define constants
       k_const = re * COS(clat*rpi/180.0) ! re * cos

C      Calculate the a and b coefficients for vorticity field
       DO ik = 1, nk
          a(ik) = 0 !Reset coefficient
          b(ik) = 0 !Reset coefficient
          DO ix = 1, nx
             a(ik) = a(ik) + vort(ix,it+1)*COS(2*rpi*real(ix)*real(ik)
     &           / real(nx))
             b(ik) = b(ik) + vort(ix,it+1)*SIN(2*rpi*real(ix)*real(ik)
     &           / real(nx))
          ENDDO
          a(ik) = a(ik)*(2/real(nx))
          b(ik) = b(ik)*(2/real(nx))
       ENDDO
C      Calculate the a(0) and a(nx/2) coefficients
       f_ave = 0
       a(nk+1) = 0
       DO ix = 1, nx
C             f_ave = f_ave + (vort(ix,it+1)/real(nx))
             a(nk+1) = a(nk+1) + vort(ix,it+1)*COS(2*rpi*real(ix)
     &           *real(nk+1) / real(nx))
       ENDDO
       a(nk + 1) = (a(nk + 1)/2)*(2/real(nx))

C      Calculate the p and q coefficients for the height field
       DO ik = 1, nk
          p(ik) = -(1.0/(real(ik)**2)) * ((k_const)**2) * a(ik)
          q(ik) = -(1.0/(real(ik)**2)) * ((k_const)**2) * b(ik)
       ENDDO
C      For a(nx/2) coefficient
       p(nk+1) = -(1.0/(real(nk+1)**2)) * ((k_const)**2) * a(nk+1)

C      Calculate Height Field by using inverse Fourier Transform
       DO ix = 1, nx
          cos_term = 0
          sin_term = 0
          DO ik = 1, nk
             cos_term = cos_term + p(ik)*cos(2*rpi*real(ik)*real(ix)
     &                 / real(nx))
             sin_term = sin_term + q(ik)*sin(2*rpi*real(ik)*real(ix)
     &                 / real(nx))
          ENDDO
          psi(ix, it+1) = f_ave + cos_term
     &                 + (p(nk+1)*cos(rpi*real(ix))) + sin_term
       ENDDO

       RETURN
       END



C      *************************************************************

       SUBROUTINE calc_psi_m(nx, nt, it, psi_1, psi_3, psi_m)
C      Calculate psi_m (Barotropic pertabation) based on psi_1 and psi_3
C      Variable "psi_m" is returned; no other variables modified
       INTEGER nx, nt, it
       REAL psi_1(nx, nt + 1), psi_3(nx, nt + 1)
       REAL psi_m(nx, nt + 1)

       DO ix = 1, nx
          psi_m(ix, it + 1) = (psi_1(ix, it + 1)+ psi_3(ix, it + 1))/2.0
       ENDDO

       RETURN
       END


C      *************************************************************

       SUBROUTINE calc_psi_t(nx, nt, it, psi_1, psi_3, psi_t)
C      Calculate psi_t (Baroclinic pertabation) based on psi_1 and psi_3
C      Variable "psi_t" is returned; no other variables modified
       INTEGER nx, nt, it
       REAL psi_1(nx, nt + 1), psi_3(nx, nt + 1)
       REAL psi_t(nx, nt + 1)

       DO ix = 1, nx
          psi_t(ix, it + 1) = (psi_1(ix, it + 1)- psi_3(ix, it + 1))/2.0
       ENDDO

       RETURN
       END



C      *************************************************************

       SUBROUTINE calc_v(it, nt, f, nx, psi, dx, v)
C      Once we have psi, we can calculate v (use center diff method)
C      Variable "v" is returned; no other variables modified
       INTEGER nx, it, nt
       REAL v(nx, nt + 1), psi(nx, nt + 1), dx, f
C      Calculate Endpoints
       v(1,it + 1) = ((psi(2, it + 1) - psi(nx, it + 1)) / (2*dx))
       v(nx,it + 1) = ((psi(1, it + 1) - psi(nx-1, it + 1)) / (2*dx))
C      Calculate Interior points
       DO ix = 2, nx-1
          v(ix,it + 1) = ((psi(ix+1,it + 1) - psi(ix-1,it + 1)) /(2*dx))
       ENDDO

       RETURN
       END
       
       
C      *************************************************************

       SUBROUTINE calc_advect_mean(nx, dx, it, nt, dt, Beta, Um, Ut,
     &              vort_m, vort_t, psi_m, psi_t, v_m, v_t, advect_mean)
C      Calculates all terms in Equation 8.16 that will be multiplied by
C                 the "dt" term
C      Results used to calculate new vort_m term
C      Works for both forward and center differencing
       REAL dt, dx, Beta
       INTEGER nx, nt, it
       REAL Um(nx), Ut(nx)
       REAL vort_m(nx, nt + 1), vort_t(nx, nt + 1)
       REAL psi_m(nx, nt + 1), psi_t(nx, nt + 1)
       REAL v_m(nx, nt + 1), v_t(nx, nt + 1)
       REAL advect_mean(nx)
C      Solve for all "interior" points
        DO ix = 2, nx - 1 ! Solve for each element at every time
           advect_mean(ix) = -(Um(ix)*(vort_m(ix+1,it)-vort_m(ix-1,it))/
     &                                                          (2*dx))
     &                -(Ut(ix)*(vort_t(ix+1,it)-vort_t(ix-1,it))/(2*dx))
     &                -(Beta*v_m(ix,it))
        ENDDO
C      Solve for the "end" points
         advect_mean(1) = -(Um(1)*(vort_m(2,it)-vort_m(nx,it))/(2*dx))
     &                   -(Ut(1)*(vort_t(2,it)-vort_t(nx,it))/(2*dx))
     &                   -(Beta*v_m(1,it))
         advect_mean(nx) = -(Um(nx)*(vort_m(1,it)-vort_m(nx-1,it))/
     &                                                           (2*dx))
     &                   -(Ut(nx)*(vort_t(1,it)-vort_t(nx-1,it))/(2*dx))
     &                   -(Beta*v_m(nx,it))

       RETURN
       END
       
       
C      *************************************************************

       SUBROUTINE calc_advect_combo(nx, dx, it, nt, dt, Beta, lambda,Um,
     &         Ut, vort_m, vort_t, psi_m, psi_t, v_m, v_t, advect_combo)
C      Calculates all terms in Equation 8.17 that will be multiplied by
C                 the "dt" term
C      Results used to calculate new vort_m term
C      Works for both forward and center differencing
       REAL dt, dx, Beta, lambda
       INTEGER nx, nt, it
       REAL Um(nx), Ut(nx)
       REAL vort_m(nx, nt + 1), vort_t(nx, nt + 1)
       REAL psi_m(nx, nt + 1), psi_t(nx, nt + 1)
       REAL v_m(nx, nt + 1), v_t(nx, nt + 1)
       REAL advect_combo(nx)
C      Solve for all "interior" points
        DO ix = 2, nx - 1 ! Solve for each element at every time
           advect_combo(ix) = -(Um(ix)*(vort_t(ix+1,it)-
     &                                          vort_t(ix-1,it))/(2*dx))
     &    +(Um(ix)*2*(lambda**2)*(psi_t(ix+1,it)-psi_t(ix-1,it))/(2*dx))
     &                -(Beta*v_t(ix,it))
     &                -(Ut(ix)*(vort_m(ix+1,it)-vort_m(ix-1,it))/(2*dx))
     &    -(Ut(ix)*2*(lambda**2)*(psi_m(ix+1,it)-psi_m(ix-1,it))/(2*dx))

        ENDDO
C      Solve for the "end" points
         advect_combo(1) = -(Um(1)*(vort_t(2,it)-vort_t(nx,it))/(2*dx))
     &          +(Um(1)*2*(lambda**2)*(psi_t(2,it)-psi_t(nx,it))/(2*dx))
     &                  -(Beta*v_t(1,it))
     &                  -(Ut(1)*(vort_m(2,it)-vort_m(nx,it))/(2*dx))
     &          -(Ut(1)*2*(lambda**2)*(psi_m(2,it)-psi_m(nx,it))/(2*dx))
         advect_combo(nx) = -(Um(nx)*(vort_t(1,it)-vort_t(nx-1,it))/
     &                                                          (2*dx))
     &       +(Um(nx)*2*(lambda**2)*(psi_t(1,it)-psi_t(nx-1,it))/(2*dx))
     &                  -(Beta*v_t(nx,it))
     &                  -(Ut(nx)*(vort_m(1,it)-vort_m(nx-1,it))/(2*dx))
     &       -(Ut(nx)*2*(lambda**2)*(psi_m(1,it)-psi_m(nx-1,it))/(2*dx))
     
       RETURN
       END
       
       
C      *************************************************************

       SUBROUTINE calc_psi_1(nx, nt, it, psi_m, psi_t, psi_1)
C      Calculate psi_1 based on psi_m and psi_t
C      Variable "psi_1" is returned; no other variables modified
       INTEGER nx, nt, it
       REAL psi_m(nx, nt + 1), psi_t(nx, nt + 1)
       REAL psi_1(nx, nt + 1)

       DO ix = 1, nx
          psi_1(ix, it + 1) = psi_m(ix, it + 1)+ psi_t(ix, it + 1)
       ENDDO

       RETURN
       END
       
       
C      *************************************************************

       SUBROUTINE calc_psi_3(nx, nt, it, psi_m, psi_t, psi_3)
C      Calculate psi_3 based on psi_m and psi_t
C      Variable "psi_3" is returned; no other variables modified
       INTEGER nx, nt, it
       REAL psi_m(nx, nt + 1), psi_t(nx, nt + 1)
       REAL psi_3(nx, nt + 1)

       DO ix = 1, nx
          psi_3(ix, it + 1) = psi_m(ix, it + 1)- psi_t(ix, it + 1)
       ENDDO

       RETURN
       END
       
       
C      *************************************************************

       SUBROUTINE Fourier_calc_vort(nx, nk, rpi, it, nt, f, psi, re,
     &                                       clat, vort)
C      Use the Fourier Transform to calculate the vorticity field
C      Variable "vort" is returned; no other variables modified
       INTEGER nx, it, nt, nk
       REAL vort(nx, nt + 1), psi(nx, nt + 1), f, re, clat, rpi

       REAL k_const
       REAL a(nk+1), b(nk), p(nk+1), q(nk)
       REAL cos_term, sin_term
C      Define constants
       k_const = re * COS(clat*rpi/180.0) ! re * cos

C      Calculate the a and b coefficients for height field
       DO ik = 1, nk
          a(ik) = 0 !Reset coefficient
          b(ik) = 0 !Reset coefficient
          DO ix = 1, nx
             a(ik) = a(ik) + psi(ix,it+1)*COS(2*rpi*real(ix)*real(ik)
     &           / real(nx))
             b(ik) = b(ik) + psi(ix,it+1)*SIN(2*rpi*real(ix)*real(ik)
     &           / real(nx))
          ENDDO
          a(ik) = a(ik)*(2/real(nx))
          b(ik) = b(ik)*(2/real(nx))
       ENDDO
C      Calculate the a(0) and a(nx/2) coefficients
       f_ave = 0
       a(nk+1) = 0
       DO ix = 1, nx
             f_ave = f_ave + (psi(ix,it+1)/real(nx))
             a(nk+1) = a(nk+1) + psi(ix,it+1)*COS(2*rpi*real(ix)
     &           *real(nk+1) / real(nx))
       ENDDO
       a(nk + 1) = (a(nk + 1)/2)*(2/real(nx))

C      Calculate the p and q coefficients for the vorticity field
       DO ik = 1, nk
          p(ik) = -((real(ik)**2) / ((k_const)**2)) * a(ik)
          q(ik) = -((real(ik)**2) / ((k_const)**2)) * b(ik)
       ENDDO
C      For a(nx/2) coefficient
       p(nk+1) = -((real(nk+1)**2) / ((k_const)**2)) * a(nk+1)

C      Calculate Vorticity Field by using inverse Fourier Transform
       DO ix = 1, nx
          cos_term = 0
          sin_term = 0
          DO ik = 1, nk
             cos_term = cos_term + p(ik)*cos(2*rpi*real(ik)*real(ix)
     &                 / real(nx))
             sin_term = sin_term + q(ik)*sin(2*rpi*real(ik)*real(ix)
     &                 / real(nx))
          ENDDO
          vort(ix, it+1)= cos_term +(p(nk+1)*cos(rpi*real(ix)))
     &                 + sin_term
       ENDDO

       RETURN
       END
       
       
C      *************************************************************

       SUBROUTINE Fourier_psi_combo(nx, nk, rpi, it, nt, f, combo, re,
     &                                       lambda, clat, psi_t)
C      Use the Fourier Transform to calculate the height field from "combo" variable
C      Variable "psi_t" is returned; no other variables modified
       INTEGER nx, it, nt, nk
       REAL combo(nx, nt + 1), psi_t(nx, nt + 1), f, re, clat, rpi

       REAL k_const, lambda
       REAL a(nk+1), b(nk), p(nk+1), q(nk)
       REAL cos_term, sin_term
C      Define constants
       k_const = re * COS(clat*rpi/180.0) ! re * cos

C      Calculate the a and b coefficients for "combo" field
       DO ik = 1, nk
          a(ik) = 0 !Reset coefficient
          b(ik) = 0 !Reset coefficient
          DO ix = 1, nx
             a(ik) = a(ik) + combo(ix,it+1)*COS(2*rpi*real(ix)*real(ik)
     &           / real(nx))
             b(ik) = b(ik) + combo(ix,it+1)*SIN(2*rpi*real(ix)*real(ik)
     &           / real(nx))
          ENDDO
          a(ik) = a(ik)*(2/real(nx))
          b(ik) = b(ik)*(2/real(nx))
       ENDDO
C      Calculate the a(0) and a(nx/2) coefficients
       f_ave = 0
       a(nk+1) = 0
       DO ix = 1, nx
             f_ave = f_ave + (combo(ix,it+1)/real(nx))
             a(nk+1) = a(nk+1) + combo(ix,it+1)*COS(2*rpi*real(ix)
     &           *real(nk+1) / real(nx))
       ENDDO
       a(nk + 1) = (a(nk + 1)/2)*(2/real(nx))

C      Calculate the p and q coefficients for the thermal height field
       DO ik = 1, nk
          p(ik) = -(1.0/((real(ik)**2)/((k_const)**2)+
     &       (2*(lambda**2)))) * a(ik)
          q(ik) = -(1.0/((real(ik)**2)/((k_const)**2)+
     &       (2*(lambda**2)))) * b(ik)
       ENDDO
C      For a(nx/2) coefficient
       p(nk+1) = -(1.0/((real(ik)**2)/((k_const)**2)+
     &       (2*(lambda**2)))) * a(nk+1)

C      Calculate Thermal Height Field by using inverse Fourier
C        Transform
       DO ix = 1, nx
          cos_term = 0
          sin_term = 0
          DO ik = 1, nk
             cos_term = cos_term + p(ik)*cos(2*rpi*real(ik)*real(ix)
     &                 / real(nx))
             sin_term = sin_term + q(ik)*sin(2*rpi*real(ik)*real(ix)
     &                 / real(nx))
          ENDDO
          psi_t(ix, it + 1) = cos_term + (p(nk+1)*cos(rpi*real(ix)))
     &                + sin_term
       ENDDO

       RETURN
       END
       
       
C      *************************************************************

       SUBROUTINE calc_w2(nx, it, nt, f, sigma, delta_P, psi_m, psi_t,
     &                                      dt, dx, Um, Ut, w2)
C      Use Holton eq. 8.12 to calculate upward motion, w2
C      Variable "w2" is returned; no other variables modified
       INTEGER nx, it, nt
       REAL psi_m(nx, nt + 1), psi_t(nx, nt + 1), f, sigma, delta_P
       REAL Um(nx), Ut(nx), dt, dx
       REAL w2(nx, nt + 1)

C      Solve for all "interior" points
       DO ix = 2, nx-1
          w2(ix, it+1) = ((2*(psi_t(ix, it+1)-psi_t(ix,it))/dt)
     &                  + (Um(ix)*(psi_t(ix+1, it)-psi_t(ix-1,it))/dx)
     &                  - (Ut(ix)*(psi_m(ix+1, it)-psi_m(ix-1,it))/dx))
     &                    * (f / (sigma*delta_P))
       ENDDO
C      Solve for endpoints
       w2(1, it+1) = ((2*(psi_t(1, it+1)-psi_t(1,it))/dt)
     &                  + (Um(1)*(psi_t(2, it)-psi_t(nx,it))/dx)
     &                  - (Ut(1)*(psi_m(2, it)-psi_m(nx,it))/dx))
     &                    * (f / (sigma*delta_P))
       w2(nx, it+1) = ((2*(psi_t(nx, it+1)-psi_t(nx,it))/dt)
     &                  + (Um(nx)*(psi_t(1, it)-psi_t(nx-1,it))/dx)
     &                  - (Ut(nx)*(psi_m(1, it)-psi_m(nx-1,it))/dx))
     &                    * (f / (sigma*delta_P))
     
       RETURN
       END


C      *************************************************************

       SUBROUTINE Smooth_random_psi(nx, nk, rpi, it, nt, re,
     &                                       clat, psi)
C      Take random initial psi field and remove higher order waves
C      Variable "psi" is returned; no other variables modified
       INTEGER nx, it, nt, nk
       REAL psi(nx, nt + 1), re, clat, rpi

       REAL k_const
       REAL a(nk+1), b(nk), p(nk+1), q(nk)
       REAL cos_term, sin_term
       INTEGER smooth
C      Define constants
       k_const = re * COS(clat*rpi/180.0) ! re * cos
C       smooth = 37 ! Remove the highest "smooth" frequencies from psi
	smooth = 10

C      Calculate the a and b coefficients for height field
       DO ik = 1, nk
          a(ik) = 0 !Reset coefficient
          b(ik) = 0 !Reset coefficient
          DO ix = 1, nx
             a(ik) = a(ik) + psi(ix,it+1)*COS(2*rpi*real(ix)*real(ik)
     &           / real(nx))
             b(ik) = b(ik) + psi(ix,it+1)*SIN(2*rpi*real(ix)*real(ik)
     &           / real(nx))
          ENDDO
          a(ik) = a(ik)*(2/real(nx))
          b(ik) = b(ik)*(2/real(nx))
       ENDDO
C      Calculate the a(0) and a(nx/2) coefficients
       f_ave = 0
       a(nk+1) = 0
       DO ix = 1, nx
             f_ave = f_ave + (psi(ix,it+1)/real(nx))
             a(nk+1) = a(nk+1) + psi(ix,it+1)*COS(2*rpi*real(ix)
     &           *real(nk+1) / real(nx))
       ENDDO
       a(nk + 1) = (a(nk + 1)/2)*(2/real(nx))

C      Rebuild Height Field by using using nk - smooth terms
       DO ix = 1, nx
          cos_term = 0
          sin_term = 0
          DO ik = 1, nk-smooth
             cos_term = cos_term + a(ik)*cos(2*rpi*real(ik)*real(ix)
     &                 / real(nx))
             sin_term = sin_term + b(ik)*sin(2*rpi*real(ik)*real(ix)
     &                 / real(nx))
          ENDDO
          psi(ix, it+1) = f_ave + cos_term + sin_term
       ENDDO

       RETURN
       END



