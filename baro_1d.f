CC
CC  Fortran 77 program to integrate the linearized barotropic
CC  vorticity equation at 500hPa for a single latitude, only 
CC  considering advection by the mean wind and ignoring all
CC  variation in y:
CC
CC    d/dt vort = -U*d/dx vort - beta*v
CC
CC  The model is global in longitude
CC 
C
	parameter(nx=144) ! same as reanalysis grid
        parameter(dt=60*15) ! time step (seconds)
	parameter(nt=21*24*4) ! number of time steps
	parameter(ntout=4) ! how many time steps to skip on output
	parameter(clat=45.0) ! latitude 
	parameter(re=6.378E6) ! radius of Earth (meters)
	parameter(rpi=3.1415926) ! Pi
        parameter(omega=7.292E-5) ! angular velocity of Earth
	parameter(s=4e6) ! spatial scale of initial disturbance
C
	real vortold(nx),vortnew(nx),vort(nx)
	real u,temp(nx)
	real phi(nx),ddxphi(nx),ddxvort(nx),v(nx)
C
C  set up input and output files and data
C
	open(13,file="baro_1d_output.dat",status="unknown",
     &    access="direct",recl=4*nx)
C
C  calculate grid spacing from number of grid points and latitude
C
	dx=cos(clat*rpi/180.0)*2.0*rpi*re/real(nx)
        f0=2.0*omega*sin(45.0*3.14159/180.0)
	beta=2.0*omega*cos(45.0*3.14159/180.0)/re
C
	irec=1
C
C  specify value for u and for initial phi
C
	u=15.0
C
	do ix=1,nx
           x=(ix-nx/2)*dx
	   gauss=exp(-x*x/(s*s))
	   phi(ix)=500.0*gauss*sin(10.0*2.0*rpi*real(ix-1)/real(nx))
	enddo
	phiav=0.0
	do ix=1,nx
	  phiav=phiav+phi(ix)/real(nx)
	enddo
	do ix=1,nx
	  phi(ix)=phi(ix)-phiav
	enddo
	write(13,rec=irec)phi
	irec=irec+1
C
C get initial vor from initial phi
C
  	call lap(phi,vort,dx)
	do ix=1,nx
	  vort(ix)=(1.0/f0)*vort(ix)
	enddo	
C
C  time step through solution
C
	do it=1,nt
C
           call getphi(vort,phi)
           call ddx(phi,ddxphi,dx)
           call ddx(vort,ddxvort,dx)
C
           do ix=1,nx
              v(ix)=(1.0/f0)*ddxphi(ix)
           enddo
C
	   if(mod(it-1,10).eq.0) then
C
C  forward time step initially and every 10 steps
C
	   	do ix=1,nx
	     	   vortnew(ix)=vort(ix)-
     &             dt*(u*ddxvort(ix) + beta*v(ix))
	   	enddo
C
	    else
C
C  else centered time step
C
                do ix=1,nx
                   vortnew(ix)=vortold(ix)-
     &             2.0*dt*(u*ddxvort(ix) + beta*v(ix))
                enddo
C
	    endif
C
C
           if(mod(it,ntout).eq.0)then
                write(13,rec=irec)phi
                irec=irec+1
           endif
C
C update variables
C
	   do ix=1,nx
             vortold(ix)=vort(ix)
	     vort(ix)=vortnew(ix) 
	   enddo
C
	   write(6,*)'done with time step',it
C
	enddo ! end loop over time integration
C
	write(6,*)'program done'
	write(6,*)'number of time steps written:',irec-1
C
C
C
	stop 
	end
C
C
C
	subroutine ddx(var,ddxvar,dx)
CC
CC  center difference assuming cyclic b.c.
CC
	parameter(nx=144)
C
	real var(nx),ddxvar(nx),dx
C
	do ix=2,nx-1
	   ddxvar(ix)=(var(ix+1)-var(ix-1))/(2.0*dx)
	enddo
	ddxvar(1)=(var(2)-var(nx))/(2.0*dx)
        ddxvar(nx)=(var(1)-var(nx-1))/(2.0*dx)
C
	return
	end
C
C
C
        subroutine lap(var,lapvar,dx)
CC
CC  1D Laplacian assuming cyclic b.c.
CC
        parameter(nx=144)
C
        real var(nx),lapvar(nx),dx
C
        do ix=2,nx-1
           lapvar(ix)=(var(ix+1)-2.0*var(ix)+var(ix-1))/(dx*dx)
        enddo
        lapvar(1)=(var(2)-2.0*var(1)+var(nx))/(dx*dx)
        lapvar(nx)=(var(1)-2.0*var(nx)+var(nx-1))/(dx*dx)
C
        return
        end
C
C
C
	subroutine getphi(vort,phi)
CC
CC  fortran 77 program to compute a fourier spectrum --
CC  a "slow fourier transform" or sft
CC
	parameter(nx=144)  ! length of series
	parameter(nf=72)  ! number of frequencies
	parameter(rpi=3.14159265)
C
        parameter(omega=7.292E-5)
        parameter(re=6.37E6)
C
C
	real r(nx) ! r is the series of interest
	real rout(nx) ! reconstructed time series
	real f(2,nf) ! fourier coefficienxs
	real all(nf,nx) ! all sines and cosines
	real both(2,nx)
	real temp(nx)
	real phi(nx),vort(nx)
C
	real vort1(nx)
C
        f0=2.0*omega*sin(45.0*3.14159/180.0)
C
	do if=1,nf
	   f(1,if)=0.0
	   f(2,if)=0.0
	enddo
C
	do it=1,nx
	   rout(it)=0.0
	enddo
C
	do it=1,nx
	do if=1,nf
	   all(if,it)=0.0
	enddo
	enddo
C
C  read in data
C
	do ix=1,nx
	   r(ix)=vort(ix)
	enddo
C
	rav=0.0
	do it=1,nx
	   rav=rav+r(it)/real(nx)
	enddo
	do it=1,nx
	   r(it)=r(it)-rav
	enddo

C  calculate fourier spectrum the slow way
C
C
	do if=1,nf
	do it=1,nx
C
	   f(1,if)=f(1,if)+
     &  2.0*r(it)*cos(real(it-1)*real(if)*2.0*rpi/real(nx))/real(nx)
	   f(2,if)=f(2,if)+
     &  2.0*r(it)*sin(real(it-1)*real(if)*2.0*rpi/real(nx))/real(nx)
C
	enddo
	enddo
C
	f(1,nf)=f(1,nf)/2.0 ! nyquist freq has differenx scaling
	f(2,nf)=0.0  ! correct any rounding errors
        irec=1
C
C  change to phi coefficients
C
	do if=1,nf
          cs=cos(45.0*rpi/180.0)
	  fac=-(f0/(real(if)*real(if)))*(re*re*cs*cs)
 	  f(1,if)=fac*f(1,if)
	  f(2,if)=fac*f(2,if)
	enddo
C
C
C
C  now reconstruct series
C
	do if=1,nf
	do it=1,nx
	   rout(it)=rout(it)+
     &       f(1,if)*cos(real(it-1)*real(if)*2.0*rpi/real(nx))
     &      +f(2,if)*sin(real(it-1)*real(if)*2.0*rpi/real(nx))
	enddo
	enddo
C
	do ix=1,nx
	  phi(ix)=rout(ix) 
	enddo
C
C
	return
	end
C
C
C
