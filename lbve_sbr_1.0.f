CC
CC  Fortran 77 program to solve the linearized barotropic
CC  vorticity equation on the sphere, for idealized ENSO forcing and
CC  a basic state of solid body rotation, following:
CC
CC  Branstator, G., 1985: Analysis of General Circulation Model
CC  Sea-Surface Temperature Anomaly Simulations Using a Linear 
CC  Model. Part I: Forced Solutions. J. Atmos. Sci., 42, 2225–2241.
CC
CC  Coded by Matt Barlow, in my youth, to learn how spectral models
CC  work.  If you find any errors, please send me
CC  an email at Mathew_Barlow@uml.edu.
CC
CC  CAUTION: has not been extensively tested but could be compared
CC  to Fig. 5a in Branstator (1985). 
CC
CC  NOTES:
CC
CC  global spectral barotropic model with triangular truncation
CC
CC  Maximum wave number and grid dimensions set separately, 
CC  up to the user to ensure the correct relationship
CC
CC  Variables ending in 's' are spherical harmonic (spectral) 
CC  coefficients and variables ending in 'g' are gridded values
CC
CC  For time differencing, vorticity at the current time is 'vort',
CC  vorticity at the previous time is 'oldvort', and vorticity
CC  at the next time is 'newvort'
CC
CC  Single index, ik, for coefficients for convenience, related to
CC  the indices of the SHs by ik=1+(l+1)*l-m
CC
CC  Keeps more waves than the maximum triangular truncation, for
CC  calculation of the meridional derivative.  For convenience, all
CC  coefficients have these extra waves, although they're only used
CC  for the calcuation of the meridional derivative
CC
CC  The standard transform to U=u*cos(lat) and V=v*cos(lat) is not 
CC  made here, to keep the code simple.  See ddx and ddy subroutines
CC  for more discussion.
CC
CC  Throughout, "SH" = spherical harmonic.
CC
	implicit none
        integer lmax,nx,ny,nk,nk2,nt,fwdskip,outskip
        real pi,re,dt,omega,r
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
	parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
	parameter(r=-1.0/(60.0*60.0*24.0*7.0)) ! drag coefficient
        parameter(pi=3.1415927)
	parameter(re=6.371e6)
	parameter(omega=7.292E-5)
C
	parameter(dt=60.0*60.0) ! time step
	parameter(nt=24*21)        ! number of times
C
	parameter(fwdskip=10) ! how often to use a forward time step
	parameter(outskip=1)  ! how often to write output
C
        real lats(ny),w(ny)
	real vortg(nx,ny),vorts(nk2),psig(nx,ny),psis(nk2)
	real varg(nx,ny),vars(nk2)
	real ddyg(nx,ny),ddxg(nx,ny)
	real fac,eta,lat,etalp1,rlp1,rm,rl,rlp2,rlm1
	real ddxvortg(nx,ny),ddyvortg(nx,ny)
	real advecg(nx,ny),advecs(nk2)
	real newvorts(nk2),oldvorts(nk2)
	real ug(nx,ny),vg(nx,ny),beta(ny)
	real alp(nk2,ny),diffs(nk2)
	real hgtg(nx,ny)
C
	real vortbarg(nx,ny),ubarg(nx,ny),vbarg(nx,ny)
	real vortbars(nk2),psibars(nk2),psibarg(nx,ny)
	real ddxvortbarg(nx,ny),ddyvortbarg(nx,ny)
	real ubars(nk2)
C
	real forceg(nx,ny),forces(nk2)
C
	integer l,m,ik,ik2,ix,iy,iklp1,iklm1,it
	integer irec
C
C open input and output files
C
C
C NB: this file can be used to read in the basic state
C for cases other than solid body rotation
C
C        open(13,file='vortT42.dat',access='direct',
C     &    recl=nx*ny*4,status='old')
C
C output
C
	open(14,file='forecast.dat',access='direct',
     &    recl=nx*ny*4,status='unknown')
C
C summary of output
C
	open(15,file='summary.dat',access='direct',
     &    recl=nx*ny*4,status='unknown')
C
C calculate Gaussian latitudes and weights, and 
C Associated Legendre Polynomials (ALPs) ahead of time
C
        call makegauss(w,lats,alp)
C
C  set forcing
C
        call makeforce(forceg,forces,alp,w,lats)
C
C
C  set up basic state: solid body rotation
C
        do ix=1,nx
        do iy=1,ny
           lat=lats(iy)*pi/180.0
           ubarg(ix,iy)=14.0*cos(lat)
        enddo
        enddo
C
        call grid2specfft(ubarg,ubars,alp,w)
        call ddyfft(-ubars,vortbarg,alp,w,lats)
        do ix=1,nx
        do iy=1,ny
           lat=lats(iy)*pi/180.0
           vortbarg(ix,iy)=vortbarg(ix,iy)+ubarg(ix,iy)*tan(lat)/re
        enddo
        enddo
        call grid2specfft(vortbarg,vortbars,alp,w)
        call vort2psi(vortbars,psibars)
        call spec2gridfft(psibars,psibarg,alp)
        call ddxfft(vortbars,ddxvortbarg,alp,w,lats)
        call ddyfft(vortbars,ddyvortbarg,alp,w,lats)
C
C  set intial conditions (spectral coefficients of vorticity)
C
        do ix=1,nx
        do iy=1,ny
           vortg(ix,iy)=0.0
        enddo
        enddo
        call grid2specfft(vortg,vorts,alp,w)
C
C
C  write out forcing and initial conditions to summary file
C
	write(15,rec=1)forceg
        write(15,rec=2)vortbarg
        write(15,rec=3)psibarg
	write(15,rec=4)ubarg
	write(15,rec=5)vbarg
        call psi2hgt(psis,hgtg,alp,w,lats)
        write(15,rec=8)hgtg
C
C  set up beta (df/dy)
C
	do iy=1,ny
            lat=lats(iy)*pi/180.0
	    beta(iy)=2.0*omega*cos(lat)/re
	enddo
C
	irec=1
C
C START TIME-STEPPING
C
	do it=1,nt
	   write(6,'(A9,I4,A4,I4)')'timestep ',it,' of ',nt
C
C  current vorticity is spectral:  vorts
C
C
C  calculate diffusion
C
	   call makediff(vorts,diffs)
C
C  calculate advection, starting with spectral values and ending with
C  spectral values, but doing multiplication on the grid
C
           call vort2psi(vorts,psis)
	   call ddyfft(psis,ug,alp,w,lats)
	   do ix=1,nx
	   do iy=1,ny
	      ug(ix,iy)=-ug(ix,iy)
	   enddo
	   enddo
	   call ddxfft(psis,vg,alp,w,lats)
           call ddxfft(vorts,ddxvortg,alp,w,lats)
	   call ddyfft(vorts,ddyvortg,alp,w,lats)
	   do ix=1,nx
	   do iy=1,ny
	      advecg(ix,iy)=
     &            -ubarg(ix,iy)*ddxvortg(ix,iy)
     &            -vbarg(ix,iy)*ddyvortg(ix,iy)
     &            -ug(ix,iy)*ddxvortbarg(ix,iy)
     &            -vg(ix,iy)*ddyvortbarg(ix,iy)
     &            -beta(iy)*vg(ix,iy)
	   enddo
	   enddo
	   call grid2specfft(advecg,advecs,alp,w)

C
	   if(mod(it-1,fwdskip).eq.0)then
C
C  forward difference time step, spectral
C
              do ik=1,nk2
                newvorts(ik)=
     &              vorts(ik)+dt*advecs(ik)
     &            + dt*(r*vorts(ik)+diffs(ik))
     &            + dt*forces(ik)

	      enddo
C
	   else
C
C  center difference time step, spectral
C
              do ik=1,nk2
                newvorts(ik)=
     &              oldvorts(ik)+2.0*dt*advecs(ik)
     &            + 2.0*dt*(r*vorts(ik)+diffs(ik))
     &            + 2.0*dt*forces(ik)

              enddo

	   endif
C
C  write out values along the way
C
	   if(mod(it-1,outskip).eq.0.or.it.eq.nt)then
	       call spec2gridfft(vorts,vortg,alp)
	       write(14,rec=irec)vortg
	       irec=irec+1
               call vort2psi(vorts,psis)
	       call spec2gridfft(psis,psig,alp)
	       write(14,rec=irec)psig
	       irec=irec+1
	   endif
C
C  update variables (in terms of spectral coefficients)
C
	   do ik=1,nk2
	     oldvorts(ik)=vorts(ik)
	     vorts(ik)=newvorts(ik)
	   enddo
C
	enddo
C
C  TIME-STEPPING DONE
C
C
C  write out to summary file
C
        call spec2gridfft(newvorts,vortg,alp)
	write(15,rec=6)vortg
        call vort2psi(vorts,psis)
	call spec2gridfft(psis,psig,alp)
	write(15,rec=7)psig
        call psi2hgt(psis,hgtg,alp,w,lats)
        write(15,rec=8)hgtg
C
C
	stop
	end
C
C
C
        subroutine makegauss(w,lats,alp)
C
        integer lmax,nx,ny,nk,nk2
        real pi
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
C
        real lats(ny),lat,mu,norm,plgndr,alp(nk2,ny)
        integer ix,iy,l,m,fact,ik
        real w(ny),x(ny)
C
C  calculate Gaussian latitudes and weights
C
        call gauleg(-1.0,1.0,x,w,ny)
C
        do iy=1,ny
           lats(iy)=asin(x(iy))*180.0/pi
        enddo
C
        ik=1
        do l=0,lmax+1
        do m=-l,l
           do iy=1,ny
              lat=lats(iy)*pi/180.0
              mu=sin(lat)
              alp(ik,iy)=plgndr(l,abs(m),mu)
           enddo
           ik=ik+1
        enddo
        enddo

C
	return
	end
C
C
C
        subroutine grid2specfft(fieldg,fields,alp,w)
C
	implicit none
        integer lmax,nx,ny,nk,nk2
        real pi
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
C
        real w(ny),fieldg(nx,ny),fields(nk2)
	real alp(nk2,ny),temp(nx),norm,tempg(nx,ny)
        integer ix,iy,ik,l,m
C
        do iy=1,ny
          do ix=1,nx
            temp(ix)=fieldg(ix,iy)
          enddo
          call realft(temp,nx,1)
          do ix=1,nx
            tempg(ix,iy)=temp(ix)
          enddo
        enddo
C
        ik=1
        do l=0,lmax
        do m=-l,l
           norm=2.0*((-1.0)**m)*sqrt(pi/real(ny)) 
           if(m.ne.0)then
             norm=norm/sqrt(2.0)
           else
             norm=norm/2.0
           endif
           fields(ik)=0.0
C
C  have to relate fourier wavenumber m to how the fft subroutine
C  stores output
C
           if(m.gt.0.and.m.lt.nx)then
             ix=2*m+1
           endif
           if(m.lt.0)then
             ix=-2*m+2
           endif
           if(m.eq.0)then
             ix=1
           endif
           if(m.eq.nx)then
             ix=2
           endif 
           do iy=1,ny
             fields(ik)=fields(ik)+
     &         norm*tempg(ix,iy)*w(iy)*alp(ik,iy)
           enddo
           ik=ik+1
        enddo
        enddo
C
C  zero out extra waves
C
        do ik=nk+1,nk2
           fields(ik)=0.0
        enddo
C
        return
        end
C
C
C
        subroutine spec2gridfft(fields,fieldg,alp)
C
C  This subroutine takes the spectral coefficients of a field and returns
C  the gridded field
C
	implicit none
        integer lmax,nx,ny,nk,nk2
        real pi
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
C
        real fieldg(nx,ny),fields(nk2)
	real norm,temp(nx),alp(nk2,ny)
        integer ix,iy,ik,l,m
C
        do ix=1,nx
        do iy=1,ny
          fieldg(ix,iy)=0.0
        enddo
        enddo
C
        ik=1
        do l=0,lmax+1
        do m=-l,l
           norm=2.0*((-1.0)**m)*sqrt(pi/real(ny))
           if(m.ne.0)then
             norm=norm/sqrt(2.0)
           else
             norm=norm/2.0
           endif
C
C  have to relate fourier wavenumber m to how the fft subroutine
C  stores output
C
           if(m.gt.0.and.m.lt.nx)then
             ix=2*m+1
           endif
           if(m.lt.0)then
             ix=-2*m+2
           endif
           if(m.eq.0)then
             ix=1
           endif
           if(m.eq.nx)then
             ix=2
           endif
           do iy=1,ny
               fieldg(ix,iy)=fieldg(ix,iy)+
     &           fields(ik)*alp(ik,iy)/norm
           enddo
           ik=ik+1
        enddo
        enddo
C
        do iy=1,ny
          do ix=1,nx
            temp(ix)=fieldg(ix,iy)
          enddo
          call realft(temp,nx,-1)
          do ix=1,nx
            fieldg(ix,iy)=temp(ix)*2.0*pi*2.0/real(nx)
          enddo
        enddo
C
	return
	end
C
C
C
      	subroutine realft(data,n,isign)
CC
CC from Numerical Recipes
CC
      	integer isign,n
      	real data(n)
C    	uses four1
C   
C        Calculates the Fourier tranform of a set of n real-valued
C        data points.  Replaces this data (which is stored in the array
C        data(1:n)) by the positive frequency half of its complex Fourier
C        transform.  The real-valued first and last components of the 
C        complex transform are returned as elements data(1) and data(2),
C        respectively.  n must be a power of 2.  This routine also
C        calculates the inverse transform of a complex data array if it
C        is the transform of real data.  (Result in this case must be
C        multiplied by 2/n).  Use isign=1 for forward transform.

      	integer i,i1,i2,i3,i4,n2p3
      	real c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      	double precision theta,wi,wpi,wpr,wr,wtemp
      	theta=3.141592653589793d0/dble(n/2)
      	c1=0.5
      	if (isign.eq.1) then
           c2=-0.5
           call four1(data,n/2,+1)
      	else
           c2=0.5
           theta=-theta
      	endif
      	wpr=-2.0d0*sin(0.5d0*theta)**2
      	wpi=sin(theta)
     	wr=1.0d0+wpr
      	wi=wpi
      	n2p3=n+3
      	do i=2,n/4
           i1=2*i-1
           i2=i1+1
           i3=n2p3-i2
           i4=i3+1
           wrs=sngl(wr)
           wis=sngl(wi)
           h1r=c1*(data(i1)+data(i3))
           h1i=c1*(data(i2)-data(i4))
           h2r=-c2*(data(i2)+data(i4))
           h2i=c2*(data(i1)-data(i3))
           data(i1)=h1r+wrs*h2r-wis*h2i
           data(i2)=h1i+wrs*h2i+wis*h2r
           data(i3)=h1r-wrs*h2r+wis*h2i
           data(i4)=-h1i+wrs*h2i+wis*h2r
           wtemp=wr
           wr=wr*wpr-wi*wpi+wr
           wi=wi*wpr+wtemp*wpi+wi
	enddo
      	if (isign.eq.1) then
           h1r=data(1)
           data(1)=h1r+data(2)
           data(2)=h1r-data(2)
      	else
           h1r=data(1)
           data(1)=c1*(h1r+data(2))
           data(2)=c1*(h1r-data(2))
           call four1(data,n/2,-1)
      	endif
      	return
      	end
C
C
C
	subroutine four1(data,nn,isign)
C
C from Numerical Recipes
C
      	integer isign,nn
      	real data(2*nn)
      	integer i,istep,j,m,mmax,n
      	real tempi,tempr
      	double precision theta,wi,wpi,wpr,wr,wtemp
      	n=2*nn
      	j=1
      	do i=1,n,2
           if(j.gt.i)then
          	tempr=data(j)
          	tempi=data(j+1)
          	data(j)=data(i)
          	data(j+1)=data(i+1)
          	data(i)=tempr
          	data(i+1)=tempi
           endif
           m=n/2
1          if ((m.ge.2).and.(j.gt.m)) then
          	j=j-m
          	m=m/2
        	goto 1
           endif
           j=j+m
	enddo
      	mmax=2
2     	if (n.gt.mmax) then
           istep=2*mmax
           theta=6.28318530717959d0/(isign*mmax)
           wpr=-2.d0*sin(0.5d0*theta)**2
           wpi=sin(theta)
           wr=1.d0
           wi=0.d0
           do m=1,mmax,2
              do i=m,n,istep
                 j=i+mmax
                 tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
                 tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
                 data(j)=data(i)-tempr
                 data(j+1)=data(i+1)-tempi
                 data(i)=data(i)+tempr
                 data(i+1)=data(i+1)+tempi
	      enddo
              wtemp=wr
              wr=wr*wpr-wi*wpi+wr
              wi=wi*wpr+wtemp*wpi+wi
	   enddo
           mmax=istep
      	   goto 2
      	endif
C
      	return
	end
C
C
C
	function plgndr(l,m,x)
CC
CC  Numerical Recipes f77 function to calculate the associated Legendre
CC  polynomials
CC
CC  from NR website (ww.nr.com/forum) - improvement of 
CC  of code in book, to better handle large l
CC
CC  These ALP are normalized
CC
	INTEGER l,m
	real plgndr,x
	INTEGER i,ll
	real fact,oldfact,pll,pmm,pmmp1,omx2,PI
	PARAMETER (PI=3.14159265)
C
	if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)then
	   write(6,*) 'bad arguments in plgndr'
	   stop
	endif
	pmm=1.
	if(m.gt.0) then
	   omx2=(1.-x)*(1.+x)
	   fact=1.
	   do i=1,m
	      pmm=pmm*omx2*fact/(fact+1.)
	      fact=fact+2.
	   enddo
	endif
	pmm=sqrt((2*m+1)*pmm/(4.*PI))
	if(mod(m,2).eq.1)pmm=-pmm
	if(l.eq.m) then
	   plgndr=pmm
	else
	   pmmp1=x*sqrt(2.*m+3.)*pmm
	   if(l.eq.m+1) then
	      plgndr=pmmp1
	   else
	      oldfact=sqrt(2.*m+3.)
	      do ll=m+2,l
	         fact=
     &             sqrt((4.*ll**2-1.)/(ll**2-m**2))
	         pll=(x*pmmp1-pmm/oldfact)*fact
	         oldfact=fact
	         pmm=pmmp1
                 pmmp1=pll
	      enddo
	      plgndr=pll
	   endif
	endif
C
	return
	END
C
C
C
  	subroutine gauleg(x1,x2,x,w,n)
CC
CC  Based on Numerical Recipes subroutine gauleg, which
CC  calculates the abscissas and weights for Gaussian integration
CC
    	implicit none
    	integer i, j, m, n
	real x1,x2,x(n),w(n)
	double precision eps
    	parameter (eps=3.D-14)
    	double precision p1,p2,p3,pp,xl,xm,z,z1

    	m=(n+1)/2
    	xm=0.5d0*(x2+x1)
    	xl=0.5d0*(x2-x1)
    	do i=1,m
       	   z1=0.
       	   z=COS(3.141592654d0*(i-.25d0)/(n+.5d0))
       	   do while ( ABS(z-z1) > EPS)
              p1=1.
              p2=0.
              do j=1,n
                p3=p2
                p2=p1
                p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j
              enddo
              pp=n*(z*p1-p2)/(z*z-1.)
              z1=z
              z=z-p1/pp
           enddo
           x(i)=xm-xl*z
           x(n+1-i)=xm+xl*z
           w(i)=2.*xl/((1.-z*z)*pp*pp)
           w(n+1-i)=w(i)
	enddo
C
	return
	end
C
C
C
	subroutine ddyfft(vars,ddyg,alp,w,lats)
CC
CC  ddy cacluates the meridional derivative using SHs, starting with the spectral form of 
CC  the variable and returning the derivative on the grid (because of the 1/cos(lat) factor)
CC
CC  calls spec2gridfft
CC
        integer lmax,nx,ny,nk,nk2
        real pi,re
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
	parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
	parameter(re=6.371e6)
C
        real lats(ny),w(ny)
	real varg(nx,ny),vars(nk2)
	real ddyg(nx,ny),ddys(nk2)
	real fac,eta,lat,etalp1,rlp1,rm,rl,rlp2,rlm1
	real alp(nk2,ny)
C
	integer l,m,ik,ik2,ix,iy,iklp1,iklm1
C
C
C
        do l=0,lmax+1   ! extra meridional mode for derivative
        do m=-l,l
C
	   rl=real(l)
	   rm=real(m)
	   eta=sqrt((rl*rl-rm*rm)/(4.0*rl*rl-1.0))
	   rlp1=rl+1.0
	   rlm1=rl-1.0
	   rlp2=rl+2.0
	   etalp1=sqrt((rlp1*rlp1-rm*rm)/(4.0*rlp1*rlp1-1.0))
C
           ik=1+(l+1)*l-m         ! index as function of l and m
	   iklp1=1+(l+2)*(l+1)-m  ! index for l+1
	   iklm1=1+l*(l-1)-m      ! index for l-1
C
	   if(l.le.lmax)then
	      if(abs(m).lt.l)then
	         ddys(ik)=rlp2*etalp1*vars(iklp1)-rlm1*eta*vars(iklm1)
	      else
                 ddys(ik)=rlp2*etalp1*vars(iklp1) ! no l-1 at edges
	      endif
	   else
              if(abs(m).lt.l)then
                 ddys(ik)=-rlm1*eta*vars(iklm1)
              else
                 ddys(ik)=0.0
              endif
	   endif
	
C
	enddo
	enddo
C
	call spec2gridfft(ddys,ddyg,alp)
C
C  still need to multiply by 1/cos(theta) -- note that this could be a problem
C  if we included the poles in our grid and also that this prevents a straightforward
C  spectral representation of the meridional derivative -- two reasons why the 
C  equations are often recast in terms of U=u*cos(theta) and V=v*cos(theta).  
C  It is also convenient
C  that Gaussian integration does not require the endpoints.
C
	do ix=1,nx
	do iy=1,ny
	  lat=lats(iy)*pi/180.0
	  ddyg(ix,iy) = ddyg(ix,iy)/(re*cos(lat))
	enddo
	enddo
C
	return
	end
C
C
C
	subroutine ddxfft(vars,ddxg,alp,w,lats)
CC
CC  ddx cacluates the zonal derivative using SHs -- input is the variable in
CC  terms of its spectral coefficients, output is on the grid (because the
CC  cos(lat) factor prevents a simple spectral representation)
CC
        integer lmax,nx,ny,nk,nk2
        real pi,re
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
	parameter(re=6.371e6)
C
        real alp(nk2,ny),lats(ny),w(ny)
	real varg(nx,ny),vars(nk2)
	real ddxg(nx,ny),ddxs(nk2)
	real fac,lat
C
	integer l,m,ik,ik2,ix,iy
C
        do l=0,lmax+1  ! including extra waves for consistency, but should = 0.0
        do m=-l,l
           if(m.ne.0.0)then
	      fac=-real(m)
           else
	      fac=0.0
	   endif
C
C  In addition to the factor of m that comes from d/dx, there's also
C  a switch from cos to sin and vice versa, which is effected by flipping
C  the sign of m in the referencing index, ik2.  Also note that the derivative
C  of both the sine and cosine have a minus since we defined our sine terms
C  using a -m in the original calculation of the SHs.
C
	   ik=1+(l+1)*l-m
	   ik2=1+(l+1)*l+m
 	   ddxs(ik)=fac*vars(ik2)
	enddo
	enddo
C

	call spec2gridfft(ddxs,ddxg,alp)
C
C
C  still need to multiply by 1/cos(theta) -- note that this could be a problem
C  if we included the poles in our grid and also that this prevents a straightforward
C  spectral representation of the meridional derivative -- two reasons why the
C  equations are often recast in terms of U=u*cos(theta) and V=v*cos(theta).  
C  It is also convenient
C  that Gaussian integration does not require the endpoints.
C
        do ix=1,nx
        do iy=1,ny
          lat=lats(iy)*pi/180.0
          ddxg(ix,iy) = ddxg(ix,iy)/(re*cos(lat))
        enddo
        enddo

C
	return
	end
C
C
C
	subroutine vort2psi(vorts,psis)
CC
CC  This subroutine takes the spectral coefficients of vorticity
CC  and returns the spectral coefficients of streamfunction
CC
        integer lmax,nx,ny,nk,nk2
        real pi,re
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
	parameter(re=6.371E6)
C
        real lats(ny),w(ny)
	real vortg(nx,ny),vorts(nk2),outg(nx,ny)
	real psig(nx,ny),psis(nk2)
	real fac
C
	integer l,m,ik
C
C  inverting the Laplacian is cake in spectral space
C
        ik=1
        do l=0,lmax+1 ! extra modes for consistency, but should be 0.0
        do m=-l,l
           if(l.gt.0)then
	      fac=-re*re/(real(l)*real(l+1))
           else
	      fac=0.0
	   endif
 	   psis(ik)=fac*vorts(ik)
	   ik=ik+1
	enddo
	enddo
C
	return
	end
C
C
C
        subroutine invlap(fields,invlaps)
CC
CC  This subroutine takes the spectral coefficients of a field
CC  and returns the spectral coefficients of the inverse laplacian
CC  of the field
CC
        integer lmax,nx,ny,nk,nk2
        real pi,re
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
        parameter(re=6.371E6)
C
        real lats(ny),w(ny)
        real fields(nk2),invlaps(nk2)
        real fac
C
        integer l,m,ik
C
C  inverting the Laplacian is cake in spectral space
C
        ik=1
        do l=0,lmax+1 ! extra modes for consistency, but should be 0.0
        do m=-l,l
           if(l.gt.0)then
              fac=-re*re/(real(l)*real(l+1))
           else
              fac=0.0
           endif
           invlaps(ik)=fac*fields(ik)
           ik=ik+1
        enddo
        enddo
C
        return
        end
C

C
C
C
        subroutine makediff(vorts,diffs)
CC
CC  This subroutine takes the spectral coefficients of vorticity
CC  and returns the spectral coefficients of a diffusion term 
CC
        integer lmax,nx,ny,nk,nk2
        real pi,re
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
        parameter(re=6.371E6)
C
        real lats(ny),w(ny)
        real vortg(nx,ny),vorts(nk2),outg(nx,ny)
        real diffs(nk2)
        real fac,nu
C
        integer l,m,ik
C
	nu=2.5E5
C
        ik=1
        do l=0,lmax+1 ! extra modes for consistency, but should be 0.0
        do m=-l,l
           if(l.gt.0)then
	      fac=-real(l)*real(l+1)/(re*re)
           else
              fac=0.0
           endif
	   if(l-m.gt.15)then
             diffs(ik)=nu*fac*vorts(ik)
	   else
	     diffs(ik)=0.0
	   endif
           ik=ik+1
        enddo
        enddo
C
        return
        end
C
C
C
        subroutine makeforce(forceg,forces,alp,w,lats)
CC
CC  This subroutine sets the forcing, returning both the gridded form and
CC  the spectral coefficients
CC
CC  Note that the zonal mean of the forcing is removed, following 
CC  Branstator (1985)
CC
        integer lmax,nx,ny,nk,nk2
        real pi,re,nu,ndiff
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(ndiff=1.0) ! order of diffusion
        parameter(nu=2.5E5) ! diffusion constant
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
        parameter(re=6.371E6)
        parameter(omega=7.292E-5)
C
        real lats(ny),w(ny)
        real forceg(nx,ny),forces(nk2),zav(ny)
        real diffs(nk2)
        real dl,dp,lon,lon0,lat,lat0,dll,dpp
	real f,d0,dx
C
        integer l,m,ik,ix,iy
C
	dl=30.0
	dp=25.0
	d0=3E-6
	lon0=200.0
	lat0=0.0
	dx=360.0/real(nx)
C
	do ix=1,nx
	do iy=1,ny
	  lon=real(ix-1)*dx+dx/2.0
	  lat=lats(iy)
          f=2.0*omega*sin(lat*pi/180.0)
	  dll=(lon-lon0)*(lon-lon0)/(dl*dl)
	  dpp=(lat-lat0)*(lat-lat0)/(dp*dp)
	  if(dll+dpp.le.1.0)then
	    forceg(ix,iy)=-f*d0*(1-sqrt(dll+dpp))
	  else
	    forceg(ix,iy)=0.0
	  endif
	enddo
	enddo
C
	do iy=1,ny
	  zav(iy)=0.0
	  do ix=1,nx
	     zav(iy)=zav(iy)+forceg(ix,iy)/real(nx)
	  enddo
	enddo
C
	do ix=1,nx
	do iy=1,ny
	   forceg(ix,iy)=forceg(ix,iy)-zav(iy)
	enddo
	enddo
C
       call grid2specfft(forceg,forces,alp,w)
C
        return
        end
C
C
C
        subroutine makediv(us,vs,divg,alp,w,lats)
CC
CC  This subroutine takes the spectral coefficients of two variables 
CC  and returns the gridded divergence
CC
        integer lmax,nx,ny,nk,nk2
        real pi,re
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
        parameter(re=6.371E6)
C
        real lats(ny),w(ny),lat
	real us(nk2),vs(nk2),divg(nx,ny)
	real ddxug(nx,ny),ddyvg(nx,ny),vg(nx,ny)
C
	integer ix,iy
C
        call spec2gridfft(vs,vg,alp)
C
        call ddxfft(us,ddxug,alp,w,lats)
        call ddyfft(vs,ddyvg,alp,w,lats)
C
	do ix=1,nx
	do iy=1,ny
	  lat=lats(iy)*pi/180.0
	  divg(ix,iy)=ddxug(ix,iy)+ddyvg(ix,iy)+vg(ix,iy)*tan(lat)/re
	enddo
	enddo	
C
        return
        end
C
C
C
C
C
C
        subroutine psi2hgt(psis,hgtg,alp,w,lats)
CC
CC  This subroutine takes the spectral coefficients of streamfunction
CC  and returns the height field, using the linearized balance equation:
CC    g*lap(h)=div(f*grad(psi))
CC
        integer lmax,nx,ny,nk,nk2
        real pi,re,omega,g
C
        parameter(lmax=21) ! max degree of SH, triang. trunc.
        parameter(nk=lmax*(lmax+2)+1) ! number of waves for lmax triang. trunc.
        parameter(nk2=lmax*lmax+4*lmax+4) ! extra modes for merid. deriv.
        parameter(nx=128)
        parameter(ny=64)
        parameter(pi=3.1415927)
        parameter(re=6.371E6)
        parameter(omega=7.292E-5)
	parameter(g=9.8)
C
	real psis(nk2),hgtg(nx,ny)
        real fac,f(ny),lat,lats(ny)
	real psig(nx,ny),ddxpsig(nx,ny),ddypsig(nx,ny)
	real ug(nx,ny),vg(nx,ny),us(nk2),vs(nk2)
	real divg(nx,ny),divs(nk2),hgts(nk2)
C
        integer l,m,ik,ix,iy
C
        do iy=1,ny
            lat=lats(iy)*pi/180.0
            f(iy)=2.0*omega*sin(lat)
        enddo
C
	call ddxfft(psis,ddxpsig,alp,w,lats)	
	call ddyfft(psis,ddypsig,alp,w,lats)
C
	do ix=1,nx
	do iy=1,ny
	   ug(ix,iy)=f(iy)*ddxpsig(ix,iy)
	   vg(ix,iy)=f(iy)*ddypsig(ix,iy)
	enddo
	enddo
C
	call grid2specfft(ug,us,alp,w)
	call grid2specfft(vg,vs,alp,w)
C
	call makediv(us,vs,divg,alp,w,lats)
C
	call grid2specfft(divg,divs,alp,w)
C
	call invlap(divs,hgts)
C
	call spec2gridfft(hgts,hgtg,alp)
C
	do ix=1,nx
	do iy=1,ny
	   hgtg(ix,iy)=hgtg(ix,iy)/g
	enddo
	enddo
C
        return
        end
C
