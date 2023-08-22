program wrapper

! gfortran amodules.f90 super_kerr.f90
  
  implicit none
  integer ne,i,ifl,j,jmax
  parameter (ne=400,jmax=8)
  real Emax,Emin,ear(0:ne),param(1),photar(ne),E,dE
  real Ed,dEd,ratio(ne),t0,t1
  real norm_super_kerr, sk_param(4)
  
! Set energy grid
  Emax  = 50.0
  Emin  = 0.01
  do i = 0,ne
    ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
  end do

! Call super_kerr
  sk_param(1)  = 1.0   !a        " "  spin
  sk_param(2)  = 50.0   !inc      deg  inclination
  sk_param(3)  = 10.0  !M        Msun Black hole mass
  sk_param(4)  = 0.057  !mdot     Edd  Accretion rate as a fraction of Eddington
  
  ! Important note: I am defining the Eddington mass accretion rate to be
  ! \dot M_edd = L_edd/c^2 = 1.26e31/c^2 * (M_bh/M_sun) [kg/s]. 
  ! or explicitly 
  ! \dot M_edd = 0.1402 (M_bh/M_sun) [10^15 kg/s]. 

  call super_kerr(ear,ne,sk_param,ifl,photar)
  norm_super_kerr = 1e-2 !Assuming D = 10 kpc, norm = 1/10^2 = 1e-2
  
! Write out model output
  do i = 1,ne
     E  = 0.5 * ( ear(i) + ear(i-1) )
     dE =         ear(i) - ear(i-1)
     write(99,*)E,norm_super_kerr*E**2*photar(i)/dE
  end do
  
end program wrapper

! Can't use an include for compiling within XSPEC
! include 'amodules.f90'

!=======================================================================
subroutine super_kerr(ear,ne,param,ifl,photar)
! Calculates observed disk spectrum
  use internal_grids
  implicit none
  integer ne,ifl,i,j,k
  real ear(0:ne),param(4),photar(ne)
  double precision a,inc,m,mdot,pi,rin,rout,mu0
  double precision Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale
  double precision chi, psi, c_a, r_b, i_b! super-extremal parameters
  double precision rnmin,rnmax,d,rfunc,disco,mudisk,re,kT_sk
  double precision alpha(nro,nphi),beta(nro,nphi),dOmega(nro,nphi)
  double precision alphan(nro,nphi),betan(nro,nphi),dOmegan(nro,nphi)
  double precision g,dlgfac,fcol,kT,phie
  real mybbody,E,dE,kTcol,dNbydE(nec),Eem,dEem
  logical needtrace
  pi  = acos(-1.d0)
  ifl = 1

! Parameters
  a       = dble( param(1) )                !Spin parameter
  inc     = dble( param(2) ) * pi / 180.d0  !Inclination (degrees)
  m       = dble( param(3) )                !Black hole mass (solar)
  mdot    = dble( param(4) )                !Mass accretion rate (Eddington)
  ! Important note: I am defining the Eddington mass accretion rate to be
  ! \dot M_edd = L_edd/c^2 = 1.26E31/c^2 * (M_bh/M_sun) [kg/s]. 
  ! or explicitly 
  ! \dot M_edd = 0.1402 (M_bh/M_sun) [10^15 kg/s]. 
  
! Derived and hardwired quantities
  rin     = disco(a)
  rout    = 5e3
  mu0     = cos(inc)
  mdot    = mdot * m * 0.1402!Convert mdot from Eddington ratio to units of 10^{15} kg/s.
  ! See above for comment on Eddington mass accretion rate definition.  
  
! Initialize
  if( firstcall )then
     firstcall = .false.
     !Define coarse internal energy grid
     Emax  = 50.0
     Emin  = 0.05
     do i = 0,nec
       earc(i) = Emin * (Emax/Emin)**(real(i)/real(nec))
     end do
     !Assign impossible initial values to previous parameters
     ! Note that only impossible parameter is now cos(inc). 
     aprev   = 10.d0
     mu0prev = 10.d0
  end if

! Initialise temperature calculation
  call Tconstants(a,rin,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,chi,psi,c_a,r_b,i_b) 
  !Assumed Mdot is in units of 10^{15} kg/s
  
! Set up full GR grid
  rnmax = 300d0                            !Sets outer boundary of full GR grid
  rnmin = rfunc(a,mu0)                    !Sets inner boundary of full GR grid
  call impactgrid(rnmin,rnmax,mu0,nro,nphi,alpha,beta,dOmega)
  d     = max( 1.0d4 , 2.0d2 * rnmax**2 ) !Sensible distance to BH  

! Set up `straight lines' grid
  rnmax = rout                            !Sets outer boundary of Newtonian grid
  rnmin = 300.d0                          !Sets inner boundary of Newtonian grid
  call impactgrid(rnmin,rnmax,mu0,nro,nphi,alphan,betan,dOmegan)
  
! Do the ray tracing in full GR
  needtrace = .false.
  if( abs( a - aprev ) .gt. tiny(a) ) needtrace = .true.
  if( abs(mu0 - mu0prev) .gt. tiny(mu0) ) needtrace = .true.
  mudisk = 0.d0       !razor thin disk
  if( needtrace )then
     call dGRtrace(nro,nphi,alpha,beta,mu0,a,rin,rout,mudisk,d,pem1,re1)
  end if
  aprev = a
  mu0prev = mu0


! Loop through inner relativistic grid
  dNbydE = 0.0
  do j = 1,nphi
     do i = 1,nro
       if( pem1(i,j) .gt. 0.d0 )then
          re = re1(i,j)
           if( re .gt. rin .and. re .le. rout )then
              !Calclulate temperature
              kT = kT_sk(re,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,mdot,m,chi,psi,c_a,r_b,i_b)
              !Calculate colour-temperature
              kTcol = kT * fcol(kT)
              !Calculate g-factor
              g = dlgfac( a,mu0,alpha(i,j),re )
              do k = 1,nec
                 E          = 0.5 * ( earc(k) + earc(k-1) )
                 dE         =       ( earc(k) - earc(k-1) )
                 Eem        =  E / g
                 dEem       = dE / g
                 dNbydE(k) = dNbydE(k) + g**3 * kT**4 * mybbody(kTcol,Eem,dEem) * dOmega(i,j) / dE
              end do
           end if
       end if
     end do
  end do

! Loop through outer Newtonian grid
  do j = 1,nphi
     do i = 1,nro
        call drandphithick(alphan(i,j),betan(i,j),mu0,mudisk,re,phie)
        if( re .gt. rin .and. re .le. rout )then
           !Calclulate temperature
           kT = kT_sk(re,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,mdot,m,chi,psi,c_a,r_b,i_b)
           !Calculate colour-temperature
           kTcol = kT * fcol(kT)
           !Calculate g-factor
           g = dlgfac( a,mu0,alphan(i,j),re )
           do k = 1,nec
              E         = 0.5 * ( earc(k) + earc(k-1) )
              dE        = earc(k) - earc(k-1)
              Eem        =  E / g
              dEem       = dE / g
              dNbydE(k) = dNbydE(k) + g**3 * kT**4 * mybbody(kTcol,Eem,dEem) * dOmegan(i,j) / dE
           end do
        end if
     end do
  end do

! Rebin onto input grid
  call myinterp(nec,earc,dNbydE,ne,ear,photar)
!  rebinE(earc,dNbydE,nec,ear,photar,ne)

! Multiply by dE
  do i = 1,ne
     E  = 0.5 * ( ear(i) + ear(i-1) )
     dE = ear(i) - ear(i-1)
     photar(i) = photar(i) * dE
  end do
  
! Re-normalise so that norm = [ 1 / Dkpc ]^2
  photar = 0.467842078 * m**2 * photar
  
  return
end subroutine super_kerr
!=======================================================================




!-----------------------------------------------------------------------
subroutine myinterp(nfx,farx,Gfx,nf,far,Gf)
! Interpolates the function Gfx from the grid farx(0:nfx) to the
! function Gf on the grid far(0:nf)
  implicit none
  integer nfx,nf
  real farx(0:nfx),Gfx(nfx),far(0:nf),Gf(nf)
  integer ix,j
  real fx(nfx),f,fxhi,Gxhi,fxlo,Gxlo
! Define grid of central input frequencies
  do ix = 1,nfx
     fx(ix) = 0.5 * ( farx(ix) + farx(ix-1) )
  end do
! Run through grid of central output frequencies
  ix = 1
  do j = 1,nf
     !Find the input grid frequencies either side of the current
     !output grid frequency
     f = 0.5 * ( far(j) + far(j-1) )
     do while( fx(ix) .lt. f .and. ix .lt. nfx )
        ix = ix + 1
     end do
     ix = max( 2 , ix )
     fxhi = fx(ix)
     Gxhi = Gfx(ix)
     ix = ix - 1
     fxlo = fx(ix)
     Gxlo = Gfx(ix)
     !Interpolate
     Gf(j) = Gxlo + ( Gxhi - Gxlo ) * ( f - fxlo ) / ( fxhi - fxlo )
  end do
  return
end subroutine myinterp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function kT_sk(re,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,mdot,m,chi,psi,c_a,r_b,i_b)
! mdot in units of 10^{15} kg s^{-1}
! m in units of solar
  implicit none
  double precision kT_sk,re,a,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale
  double precision chi,psi,c_a,r_b,i_b
  double precision mdot,m
  double precision x,br1,br2
  x    = sqrt(re)

  if (abs(a) .lt. 1.0) then
    br1 = dble((1.0 - (3.0*a)/(2.0*x)*log(x) + Ca/x*log(x-xa) + Cb/x*log(x-xb) + Cg/x*log(x-xg) - j0*xI/x)**0.25)
    br2 = dble( (1.0/(1.0 - 3.0/(x**2) + 2.0*a/(x**3)))**0.25 )
  end if
  if (abs(a) .gt. 1.0) then
    br1 = dble((1.0 - (3.0*a)/(2.0*x)*log(x) ))
    br1 = dble( br1 + c_a/x * log(x-chi) + r_b/x * log((x+0.5*chi)**2.0+psi**2.0) )
    br1 = dble( br1 +  4.0 * i_b/x * atan(psi/(x+0.5*chi + ((x+0.5*chi)**2.0+psi**2.0)**0.5) ) )
    br1 = dble( br1 - j0*xI/x )
    br1 = dble( br1**0.25 )
    br2 = dble( (1.0/(1.0 - 3.0/(x**2) + 2.0*a/(x**3)))**0.25 )
  end if
  if (a .eq. 1.0) then
    if (x .gt. 1.0) then
      br1 = dble((1.0 - (3.0)/(2.0*x)*log(x) + 3.0/(2.0*x)*log(x+2.0) - j0*xI/x ))
      br1 = dble( br1**0.25 )  
      br2 = dble( (1.0/(1.0 - 3.0/(x**2) + 2.0*a/(x**3)))**0.25 )
    else
      br1 = 0.0
      br2 = 0.0
    end if 
  end if
  if (a .eq. -1.0) then
    br1 = dble((1.0 + (3.0)/(2.0*x)*log(x) - 3.0/(2.0*x)*log(x-2.0) - j0*xI/x ))
    br1 = dble( br1**0.25 )  
    br2 = dble( (1.0/(1.0 - 3.0/(x**2) + 2.0*a/(x**3)))**0.25 )
  end if 

  kT_sk = Tscale * (mdot**0.25) / (m**0.5) * re**(-0.75) * br1 * br2

  return
end function kT_sk  
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine Tconstants(a,rin,Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale,chi,psi,c_a,r_b,i_b)
! Defines physical parameters and derived quantities for relativistic
! temperature calculation
! mdot in units of 10^{15} kg s^{-1}
! m in units of solar
! Input: a,rin
! Output: xa,xb,xg,j0,xI,Tscale,chi,psi,c_a,r_b,i_b
  implicit none
  double precision a,rin
  double precision Ca,Cb,Cg,xa,xb,xg,j0,xI,Tscale
  double precision chi, psi, c_a, r_b, i_b! super-extremal parameters
  double precision rI,disco,pi
  pi = acos(-1.d0)
  Tscale  = 8.31697777 !Tempaerature scale, in keV, for a one solar mass BH and an 10^15 kg/s accretion rate)
  if (abs(a) .lt. 1.0) then
    xa = dble( 2.0*cos(1.0/3.0 * acos(-a) ) )
    xb = dble( 2.0*cos(1.0/3.0 * acos(-a) - 2.0*pi/3.0 ) )
    xg = dble( 2.0*cos(1.0/3.0 * acos(-a) - 4.0*pi/3.0 ) )
    
    Ca = dble( (2.0*xa - a*(1 + xa*xa))/(2.0*(1 - xa*xa)) )
    Cb = dble( (2.0*xb - a*(1 + xb*xb))/(2.0*(1 - xb*xb)) )
    Cg = dble( (2.0*xg - a*(1 + xg*xg))/(2.0*(1 - xg*xg)) )
    
    rI = disco(a) ! Finds the ISCO in gravitational units
    xI = dble( sqrt(rI) )
    
    j0=dble( (1.0 - (3.0*a)/(2.0*xI) * log(xI) + Ca/xI * log(xI-xa) + Cb/xI * log(xI-xb) + Cg/xI * log(xI-xg)) )
    ! Vanishing stress angular momentum constant (in natural units).     
  end if
  if (abs(a) .gt. 1.0) then
    chi = dble(  -2.0 * cosh(1.0/3.0 * acosh(abs(a)) ) )
    psi = dble( (3.0/4.0*chi**2.0 - 3.0)**0.5 )

    c_a = dble(2.0 * chi - a * (1+chi**2.0))/(2.0*(1-chi**2.0))    

    r_b= dble((-chi-a*(1.0+0.25*chi**2.0-psi**2.0)))
    r_b = dble( r_b*(1.0-0.25*chi**2.0+psi**2.0)+2.0*chi*psi**2.0*(1.0+0.5*chi*a) )
    r_b = dble(r_b/(2.0*((1.0 - 0.25*chi**2.0 + psi**2.0)**2.0 + chi**2.0 * psi**2.0)))

    i_b=dble((chi+a*(1.0+0.25*chi**2.0 - psi**2.0))*chi*psi + 2.0*psi*(1.0+0.5*a*chi)*(1-0.25*chi**2.0+psi**2.0))
    i_b = dble(i_b/(2.0*((1.0 - 0.25*chi**2.0 + psi**2.0)**2.0 + chi**2.0 * psi**2.0)))

    rI = disco(a) ! Finds the ISCO in gravitational units
    xI = dble( sqrt(rI) )

    j0 = dble( (1.0 - (3.0*a)/(2.0*xI) * log(xI) ))
    j0 = dble(j0 + (c_a/xI * log(xI-chi) + r_b/xI * log((xI+0.5*chi)**2.0+psi**2.0) ) )
    j0 = dble(j0 + 4.0 * i_b/xI * atan(psi/(xI+0.5*chi + ((xI+0.5*chi)**2.0+psi**2.0)**0.5) ) )
    ! Vanishing stress angular momentum constant (in natural units). 
  end if
  if (a .eq. 1.0) then
    rI = disco(a)
    xI = dble( sqrt(rI))
    j0 = 1.0 + 3.0/(2.0*xI) * log(xI + 2.0) - 3.0/(2.0*xI)*log(xI)
    ! Vanishing stress angular momentum constant (in natural units). 
  end if 
  if (a .eq. -1.0) then
    rI = disco(a)
    xI = dble( sqrt(rI))
    j0 = 1.0 - 3.0/(2.0*xI) * log(xI - 2.0) + 3.0/(2.0*xI)*log(xI)
    ! Vanishing stress angular momentum constant (in natural units). 
  end if 
  return
end subroutine Tconstants
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
function fcol(T)
! fcol = colour-temperature correction
! T    = true temperature in keV
  implicit none
  double precision fcol,T
  if( T .lt. 2.585d-3 )then
     fcol = 1.d0
  else if( T .lt. 8.617e-3 )then
     fcol = ( T / 2.585d-3 )**0.833
  else
     fcol = ( 72.d0 / T )**(1.d0/9.d0)
  end if
  return
end function fcol
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
function mybbody(kT,E,dE)
! Blackbody function in terms of number of photons with energy
! between E-dE/2 and E+dE/2. i.e. This is photar!
! Function is normalized such that the integrated energy flux is 1
! i.e. sum E * photar(E) = 1
  implicit none
  real mybbody,E,dE,kT
  real pi,fac,f
  pi   = acos(-1.0)
  fac  = E/kT
  if(fac .lt. 1e-3)then
    f = E * kT   !Using a Taylor expansion
  else
    if (fac .lt. 70.) then
      f = E**2 / ( exp(fac) - 1.0 ) 
    else
      f = E**2 * exp(-fac)
    end  if
  end if
  mybbody = f * (15.0/pi**4) / kT**4 * dE
  return
end function mybbody
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function dlgfac(a,mu0,alpha,r)
!c Calculates g-factor for a disk in the BH equatorial plane
  implicit none
  double precision dlgfac,a,mu0,alpha,r
  double precision sin0,omega,Delta,Sigma2,gtt,gtp,gpp
  sin0   = sqrt( 1.0 - mu0**2 )
  omega  = 1. / (r**1.5+a)
  Delta  = r**2 - 2*r + a**2
  Sigma2 = (r**2+a**2)**2 - a**2 * Delta
  gtt    = 4*a**2/Sigma2 - r**2*Delta/Sigma2
  gtp    = -2*a/r
  gpp    = Sigma2/r**2
  dlgfac = sqrt( -gtt - 2*omega*gtp - omega**2.*gpp )
  dlgfac = dlgfac / ( 1.+omega*alpha*sin0 )
  return
end function dlgfac
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine dGRtrace(nro,nphi,alpha,beta,mu0,spin,rmin,rout,mudisk,d,pem1,re1)
! Traces rays in the Kerr metric for a camera defined by the impact
! parameters at infinity: alpha(nro,nphi) and beta(nro,nphi).
! Traces back to a disk defined by mudisk = cos(theta_disk), where
! theta_disk is the angle between the vertical and the disk surface.
! i.e. tan( theta_disk ) = 1 / (h/r)
! OUTPUT:
! pem1(nro,nphi)
! pem > 1: there is a solution
! pem = -1 photon goes to infinity without hitting disk surface
! pem = -2 photon falls into horizon without hitting disk surface
! re1(nro,nphi)      radius that the geodesic hits the disc
  use blcoordinate     ! This is a YNOGK module
  implicit none
  integer nro,nphi,i,j
  double precision alpha(nro,nphi),beta(nro,nphi),mu0,spin,rmin,rout,mudisk,d
  double precision pem1(nro,nphi),re1(nro,nphi)
  double precision cos0,sin0,scal,velocity(3),f1234(4),lambda,q
  double precision pem,re,mucros,phie,taudo,sigmacros      
  cos0  = mu0
  sin0  = sqrt(1.0-cos0**2)
  scal     = 1.d0
  velocity = 0.d0
  re1      = 0.0
  do i = 1,nro
    do j = 1,NPHI
      call lambdaq(-alpha(i,j),-beta(i,j),d,sin0,cos0,spin,scal,velocity,f1234,lambda,q)
      pem = Pemdisk(f1234,lambda,q,sin0,cos0,spin,d,scal,mudisk,rout,rmin)  !Can try rin instead of rmin to save an if statement
      pem1(i,j) = pem
      !pem > 1 means there is a solution
      !pem < 1 means there is no solution
      if( pem .gt. 0.0d0 )then
        call YNOGK(pem,f1234,lambda,q,sin0,cos0,spin,d,scal,re,mucros,phie,taudo,sigmacros)
        re1(i,j)    = re
      end if
    end do
  end do
  return
end subroutine dGRtrace
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine impactgrid(rnmin,rnmax,mu0,nro,nphi,alpha,beta,dOmega)
! Calculates a grid of impact parameters
! INPUT:
! rnmin        Sets inner edge of impact parameter grid
! rnmax        Sets outer edge of impact parameter grid
! mu0          Sets `eccentricity' of the grid
! nro          Number of steps in radial impact parameter (b)
! nphi         Number of steps in azimuthal impact parameter (phi)
! OUTPUT:
! alpha(nro,nphi)   Horizontal impact parameter
! beta(nro,nphi)    Vertical impact parameter
! dOmega(nro,nphi)  dalpha*dbeta
  implicit none
  integer nro,nphi,i,j
  double precision rnmin,rnmax,mu0,alpha(nro,nphi),beta(nro,nphi)
  double precision dOmega(nro,nphi),mueff,pi,rar(0:nro),dlogr,rn(nro)
  double precision logr,phin
  pi     = acos(-1.d0)

  mueff = max( mu0 , 0.3d0 )
  
  rar(0) = rnmin
  dlogr  = log10( rnmax/rnmin ) / dble(nro)
  do i = 1,NRO
    logr = log10(rnmin) + dble(i) * dlogr
    rar(i)    = 10.d0**logr
    rn(i)     = 0.5 * ( rar(i) + rar(i-1) )
    do j = 1,nphi
       domega(i,j) = rn(i) * ( rar(i) - rar(i-1) ) * mueff * 2.d0 * pi / dble(nphi)
       phin       = (j-0.5) * 2.d0 * pi / dble(nphi) 
       alpha(i,j) = rn(i)  * sin(phin)
       beta(i,j)  = rn(i) * cos(phin) * mueff
    end do
  end do
  
  return
end subroutine impactgrid
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function dISCO(a)
  !ISCO in Rg 
  implicit none
  double precision a,dISCO,z1,z2
  if(abs(a).ge.1.0)then
    z1 = -(a**2.0 - 1.0)**(1.0/3.0)
    z1 = z1 * ( (1.0+abs(a))**(1.0/3.0)-(abs(a)-1.0)**(1.0/3.0))+1.0
  else
    z1 = ( 1.0 - a**2.0 )**(1.0/3.0)
    z1 = z1 * ( (1.0+a)**(1.0/3.0)+(1.0-a)**(1.0/3.0))+1.0
  end if 
	  
  z2 = sqrt( 3.0 * a**2.0 + z1**2.0 )
  if(a.ge.0.0)then
    dISCO = 3.0 + z2 - sqrt( (3.0-z1) * (3.0 + z1 + 2.0*z2) )
  else
    dISCO = 3.0 + z2 + sqrt( (3.0-z1) * (3.0 + z1 + 2.0*z2) )
  end if
  return
end function dISCO
!-----------------------------------------------------------------------

  

!-----------------------------------------------------------------------
function rfunc(a,mu0)
! Sets minimum rn to use for impact parameter grid depending on mu0
! This is just an analytic function based on empirical calculations:
! I simply set a=0.998, went through the full range of mu0, and then
! calculated the lowest rn value for which there was a disk crossing.
! The function used here makes sure the calculated rnmin is always
! slightly lower than the one required.
  implicit none
  double precision rfunc,mu0,a
  if( a .gt. 0.8 )then
    rfunc = 1.5d0 + 0.5d0 * mu0**5.5d0
    rfunc = min( rfunc , -0.1d0 + 5.6d0*mu0 )
    rfunc = max( 0.1d0 , rfunc )
  else
    rfunc = 3.0d0 + 0.5d0 * mu0**5.5d0
    rfunc = min( rfunc , -0.2d0 + 10.0d0*mu0 )
    rfunc = max( 0.1d0 , rfunc )
  end if
  end function rfunc
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
! Calculates an r-grid that will be used to define impact parameters
  implicit none
  integer nro,nphi,i
  double precision rnmin,rnmax,mueff,rn(nro),domega(nro)
  double precision rar(0:nro),dlogr,logr,pi
  pi     = acos(-1.d0)
  rar(0) = rnmin
  dlogr  = log10( rnmax/rnmin ) / dble(nro)
  do i = 1,NRO
    logr = log10(rnmin) + dble(i) * dlogr
    rar(i)    = 10.d0**logr
    rn(i)     = 0.5 * ( rar(i) + rar(i-1) )
    domega(i) = rn(i) * ( rar(i) - rar(i-1) ) * mueff * 2.d0 * pi / dble(nphi)
  end do
  return
end subroutine getrgrid
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine drandphithick(alpha,beta,cosi,costheta,r,phi)
!
! A disk with an arbitrary thickness
! The angle between the normal to the midplane and the disk surface is theta
! The inclination angle is i
      implicit none
      double precision alpha,beta,cosi,sini,r,phi
      double precision pi,costheta,sintheta,x,a,b,c,det
      double precision mu,sinphi
!      double precision muplus,muminus,ra,rb,rab,xplus1,xminus1,xplus2,xminus2
      pi = acos(-1.d0)
      sintheta = sqrt( 1.d0 - costheta**2 )
      sini     = sqrt( 1.d0 - cosi**2 )
      x        = alpha / beta
      if( abs(alpha) .lt. abs(tiny(alpha)) .and. abs(beta) .lt. abs(tiny(beta))  )then
        mu = 0.d0
        r  = 0.d0
      else if( abs(beta) .lt. abs(tiny(beta)) )then
        mu     = sini*costheta/(cosi*sintheta)
        sinphi = sign( 1.d0 , alpha ) * sqrt( 1.d0 - mu**2 )
        r      = alpha / ( sintheta * sinphi )
      else if( abs(alpha) .lt. abs(tiny(alpha)) )then
        mu     = 1.d0
        sinphi = 0.d0
        r      = beta / ( sini*costheta - cosi*sintheta )
      else
        a      = sintheta**2 + x**2*cosi**2*sintheta**2
        b      = -2*x**2*sini*cosi*sintheta*costheta
        c      = x**2*sini**2*costheta**2-sintheta**2
        det    = b**2 - 4.d0 * a * c
        if( det .lt. 0.d0 ) write(*,*)"determinant <0!!!"
        if( beta .gt. 0.d0 )then
          mu     = ( -b + sqrt( det ) ) / ( 2.d0 * a )
        else
          mu     = ( -b - sqrt( det ) ) / ( 2.d0 * a )
        end if
        sinphi = sign( 1.d0 , alpha ) * sqrt( 1.d0 - mu**2 )
        r      = alpha / ( sintheta * sinphi )
      end if
      phi = atan2( sinphi , mu )
      return
      end subroutine drandphithick
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine rebinE(earx,px,nex,ear,p,ne)
      !General rebinning scheme, should be nice and robust - BUT IT FUCKING ISN'T
      !i,nex,earx,px = input
      !j,ne,ear,p    = output
      implicit none
      integer i,nex,j,ne
      real earx(0:nex),ear(0:ne),px(nex),p(ne),Ehigh,Elow,upper,lower
      real FRAC,Ej,Ei,pi,Ei2,pi2,grad,cons
      logical interp
      i = 1
      do j = 1,ne
        p(j) = 0.0
        Ehigh = ear(j)
        Elow  = ear(j-1)
        do while( earx(i) .le. Elow .and. i .lt. nex )
          i = i + 1
        end do
        interp = .true.
        do while(earx(i-1).lt.Ehigh.and.i.lt.nex)
          lower = MAX( earx(i-1) , Elow  )
          upper = MIN( earx(i)   , Ehigh )
          FRAC  = (upper-lower) / ( Ehigh - Elow )
          p(j)  = p(j) + px(i) * FRAC
          i = i + 1
          interp = .false.
        end do
        if( interp )then
          !Work out if it's ok to interpolate
          if( Elow  .ge. earx(nex) ) interp = .false.
          if( Ehigh .le. earx(0)   ) interp = .false.
          if( i     .ge. nex-1     ) interp = .false.
        end if
        if( interp )then
          !write(*,*)"need to interpolate!"
          !p(j) is interpolation between px(i) and px(i+1)
          !unless i=nex, in which case p(j) = px(i)
          Ej = 0.5 * ( Ehigh + Elow )        !Centre of newbin
          Ei = 0.5 * ( earx(i+1) + earx(i) ) !Centre of one oldbin
          pi = px(i+1)         !Value at bin centre
          if( Ei .eq. Ej )then
             p(j) = pi
          else
            if( Ei .gt. Ej )then
              Ei2 = 0.5 * (earx(i) + earx(i-1) )
              pi2 = px(i)
            else
              Ei2 = 0.5 * (earx(i+2) + earx(i+1) )
              pi2 = px(+2)
            end if
            grad = ( pi - pi2 ) / ( Ei - Ei2 )
            cons = 0.5 * ( pi + pi2 - grad*(Ei+Ei2) )
            p(j) = grad * Ej + cons
          end if
        end if
        if( i .gt. 2 ) i = i - 1
      end do
      RETURN
      END
!-----------------------------------------------------------------------


