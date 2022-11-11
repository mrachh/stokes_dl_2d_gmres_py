c     This program tests the DL representation for Stokes
c     flow for solving interior problem in 2D 
c     for velocity boundary conditions
c     using an iterative solver.

      implicit real *8 (a-h,o-z)
      real *8, allocatable :: src(:,:),dsdt(:),rn(:,:),rkappa(:)
      complex *16, allocatable :: rhs(:),soln(:)

      integer ifbipot

      real *8 a,b,cx,cy,thetai,theta,xpos,ypos,xp,yp
      real *8 cxsource,cysource
      real *8 xvt,yvt,vort,pressure
      real *8 bipot, bipotexact

      real *8 swgt,dwgt,cwgt,pwgt,pi,h2

      real *8 t1,t2,uinf,vinf,w

      complex *16 zdis,zx,ratio,eye
      complex *16 zvel,zvort

      real *8 eps
      integer numit, ngmrec
      complex *16 ima,zn,zu,ztmp,velpert
      real *8, allocatable :: errs(:)

      data ima/(0.0d0,1.0d0)/

      done = 1.0d0
      pi = atan(done)*4.0d0


c     The weight of the double layer in the representation
      dwgt = 2.0d0

c     The weight corresponding to the generalized 1's matrix
      pwgt = 1.0d0


      call prini(6,13)

      a = 1.2d0
      b = 1.3d0
      n = 256
      cxsource = -3.05d0 
      cysource = -0.02d0 

      allocate(src(2,n),dsdt(n),rn(2,n),rkappa(n))
      allocate(rhs(n),soln(n))
      ifbipot = 0
      itype = 0
      velpert = 0
      do i=1,n
        theta = (i-1.0d0)*2*pi/n
        ct = cos(theta)
        st = sin(theta)
        src(1,i) = a*ct
        src(2,i) = b*st

        dxdt = -a*st
        dydt = b*ct
        
        dsdt(i) = sqrt(dxdt**2 + dydt**2)
        rn(1,i) = dydt/dsdt(i)
        rn(2,i) = -dxdt/dsdt(i)
        rkappa(i) = a*b/dsdt(i)**3
        dsdt(i) = dsdt(i)*2*pi/n
        call uexact(itype,src(1,i),src(2,i),cxsource,
     1           cysource,xvel,yvel,vort,
     2           pressure,ifbipot,bipotexact)      
        rhs(i) = -yvel + ima*xvel
        zn = rn(1,i) + ima*rn(2,i)
        velpert = velpert + dimag(conjg(zn)*rhs(i))*dsdt(i)
      enddo
      print *, "velpert=",velpert


      ra = 0
      rk = 0

      do i=1,n
        ra = ra + dsdt(i)
        rk = rk + rkappa(i)
      enddo





c     Set GMRES parameters
      eps = 1.0d-12
      numit = 50
      niter = 0
      allocate(errs(numit+1))
c     Call the solver
      t1 = second()
  
      call stokes_gmres(n,src,dsdt,rn,rkappa,dwgt,pwgt,eps,rhs,
     1  numit,eps,niter,errs,rres,soln)
      
      t2 = second()
      call prin2('Total time to compute solution = *',t2-t1,1)
      call prin2('soln=*',soln,24)
      print *, 'Enter x'
      read *, xp
      print *, 'Enter y'
      read *, yp
      zx = dcmplx(xp,yp)

      zvel = 0


      do i=1,n
        zn = rn(1,i) + ima*rn(2,i)
        zdis = (src(1,i) - xp) + ima*(src(2,i)-yp)
        zu = soln(i)
        rsc = dsdt(i)/(2.0d0*pi)


        zvel = zvel + dwgt*0.5d0*rsc*ima*zn*zu/zdis
        zvel = zvel - dwgt*0.5d0*rsc*ima*dconjg(zn*zu/zdis**2)*
     1                     zdis
        zvel = zvel - dwgt*rsc*dreal(ima*zn*dconjg(zu))/
     1                     dconjg(zdis)

      enddo
      print *, 'zvel = ',zvel
      xvt = 0
      yvt = 0
      vort = 0
      ifw = 0
      call uexact(itype,xp,yp,cxsource,cysource,xvt,yvt,vort,
     1            pressure,ifw,w) 
      print *, 'xvt,yvt',xvt,yvt
      print *, 'err=',dcmplx(xvt,yvt) - zvel

      return
      end
c---------------------------------------------------------------------      
      
      subroutine uexact(itype,xpos,ypos,cx,cy,xvel,yvel,vort,
     1                  pressure,ifw,w)
c      This subroutine evaluates a particular stokes field
c      corresponding to W = x log r where r = sqrt(x**2 + y**2)
c      if itype = 1 and if itype =0, it returns the field
c      phi = 1/zdis and psi = 3/zdis where zdis is the
c      distance from the source located at (cx,cy)
c      ARGUMENTS
c      ------------------------------------------------------
c      itype     IN: integer
c                type of boundary data to be picked up
c
c      xpos      In: real *8
c                x coordinate of target location where field
c                is to be evaluated
c
c      ypos      In: real *8
c                y coordinate of target location where field
c                is to be evaluated
c
c      cx        In: real *8
c                x coordinate of source location
c
c      cy        In: real *8
c                y coordinate of source location
c
c      ifw       In: integer
c                flag to compute the biharmonic potential
c--------------------------------------------------------------
c      OUTPUT
c      xvel      Out: real *8
c                x component of the velocity at the target
c      yvel      Out: real *8
c                y component of the velocity at the target
c      w         Out: real *8
c                biharomnic potential at the target
c -------------------------------------------------------------

      implicit none

      integer itype,ifw
      real *8 xpos,ypos,cx,cy,xvel,yvel,vort,pressure,w
      complex *16 ztarg,zdis,zdis2
      complex *16 phi,phip,psi,gw,chi,bipot

      ztarg = dcmplx(xpos,ypos)
      zdis = ztarg - dcmplx(cx,cy)


c     Option 0: decaying solution at infinity
      if(itype.eq.0) then
         phi = 1.0d0/zdis
         phip = -1.0d0/(zdis*zdis)
         psi = 3.0d0/zdis
         chi = 3.0d0*log(zdis)

         gw = phi + ztarg*dconjg(phip)+dconjg(psi)
         vort = 4*dreal(phip)
         pressure = 4*dimag(phip)

         xvel = dimag(gw)
         yvel = -dreal(gw)
         if(ifw.eq.1)  w = dreal(dconjg(ztarg)*phi + chi)
      else
         yvel = -log(cdabs(ztarg)) - xpos**2/cdabs(ztarg)**2
         xvel = xpos*ypos/cdabs(ztarg)**2
         vort = -2.0d0*xpos/cdabs(ztarg)**2
         if(ifw.eq.1) w = xpos*log(cdabs(ztarg))
      endif
         

      return
      end
c ----------------------------------------------------------------      

      subroutine direval(zx,zvort,zvel,k0,k,nd,xs,ys,rnx,rny,
     1                   rkappa,dsdt,h,wdens,swgt,dwgt,cwgt,pwgt)
c     Given complex density wdens compute the velocity corresponding
c     to the SL+DL representation for Stokes equation. 
c  
c     All the integrals in the representation below are computed
c     using trapezoidal rule.
c
C     The stream function W  is represented
C     by the complex variables formula
C
C           W = Re( CONJG(Z)*PHI + XI)
C     with 
c           W_x + i W_y = PHI + Z*CONJG(PHI') + CONJG(PSI)
c     where u = W_y and v = -W_x are the x and the y components
c     of velocity respectively
c
c     where PSI = XI'. The modified Sherman-Lauricella equation
c     is based on the following representation.
c
c     \phi(z) = -swgt( \int [w(t) log (t-z) ds(t)])/ (8 \pi) 
c               -dwgt( \int [w(t) / (t-z) dt])/ (4 \pi i)
c               + cwgt (\int [w(t) ds(t) ]). (2 \pi) 
c
c     where z,t are complex variables, w(t) is a complex density,
c
c     for any variable x, let xb be dconjg(x)
c     \psi(z) = -swgt(\int [wb log (t-z) ds(t) +
c                    tb*w*ds(t)/(t-z)])/(8\pi)
c               - dwgt(\int [(wb*dt + w*dtb)/(t-z) 
c                      - tb*w(t)/(t-z)**2])/(4 \pi i)
c
c-----------------------------------------------------------------
c     INPUT ARGUMENTS
c     zx       In: complex *16
c              zx = xt + i yt. Location of target in complex
c              notation
C
C     k0,k      in: INTEGER
C
C               k0 = 0 for interior and 1 for exterior problems.
C               k  = number of obstacles.
c 
C     nd(0:k)  in: integer
C
C               ND(J) = number of points in discretization of
C                       Jth boundary component.
C
C**********************************************************************
C     For future reference we define ntot = SUM_{J=K0}^K ND(J)
C**********************************************************************
C
C     The x and y coordinates of the points defining the boundary
C     components are given by
C
C     xs(ntot),ys(NTOT) IN: REAL *8 
C
C               xs(1,...,ND(K0)) xcoordinates of K0th boundary component.
C               ys(1,...,ND(K0)) xcoordinates of K0th boundary component.
C
C               ...
C
C               xs(NTOT-ND(K)+1,...,NTOT) xcoords of Kth boundary.
C               ys(NTOT-ND(K)+1,...,NTOT) xcoords of Kth boundary.
C
C**********************************************************************
C     They are assumed to be equispaced with respect to
C     some parameter theta (see dsdtheta below).
C**********************************************************************
C
C     dsdt(NTOT) IN: REAL *8 
C
C               ds/d(theta) where s denotes arclength and theta is the 
C               parametrization of the curve in which the points are 
C               equispaced (using the same ordering as the X,Y array).
c
C     rnx(NTOT),rny(NTOT) IN: REAL *8 
C
C               x and y components of outward normal at corresponding
C               point in xs,ys array.
C
C     rkappa(NTOT) IN: REAL *8 
C
C               curvature at corresponding point in xs,ys array.
C
C
C     h(0:K)   IN: REAL *8 
C
C               weight of trapezoidal rule, that is h(J) is length of 
C               Jth boundary measured in parameter theta, divided by 
C               nd(J).
C     wdens(1:ntot)   In: complex *16
c               The density for SL + DL in the above representation
c
C     swgt      in: real *8
c               Weight of the single layer in the representation for
c               velocity
c
c     dwgt      in: real *8
c               Weight of the double layer in the representation for
c               velocity
c
c     cwgt      in: real *8
c               Weight of scaled integral of density in the
c               representation for velocity (Check documentation above)
c      
c     pwgt      in: real *8
c               Weight of generalized 1's matrix to compensate
c               for the nullspace in the representation for interior
c               flows
c-------------------------------------------------------------------
c     OUTPUT arguments
c     zvort     out: complex*16
c               zvort = vort + i pressure at the target
c 
c     zvel      out: complex *16
c               complex velocity  -v + i u at the target
c-------------------------------------------------------------------

      implicit none
      integer k0,k,nd(0:k)
      integer i,ii,istart,nbod

      real *8 xs(*),ys(*),dsdt(*),rnx(*),rny(*),rkappa(*)
      real *8 h(0:k),swgt,dwgt,cwgt,pwgt

      real *8 pi

      complex *16 zn,zx,zvel,zvort,zdis,eye,zsc,zu,wdens(*)

      zvel = dcmplx(0.0d0,0.0d0)
      zvort = dcmplx(0.0d0,0.0d0)
      eye = dcmplx(0.0d0,1.0d0)
      pi = 4.0d0*datan(1.0d0)
      istart = 0

      do nbod = k0,k
         do i=1,nd(nbod)
            ii = i+istart
            zn = dcmplx(rnx(ii),rny(ii))
            zdis = dcmplx(xs(ii),ys(ii)) - zx
            zu = wdens(ii)
            zsc = h(nbod)*dsdt(ii)/(2.0d0*pi*eye)


            zvel = zvel - dwgt*0.5d0*zsc*eye*zn*zu/zdis
            zvel = zvel + dwgt*0.5d0*zsc*eye*dconjg(zn*zu/zdis**2)*
     1                     zdis
            zvel = zvel + dwgt*zsc*dreal(eye*zn*dconjg(zu))/
     1                     dconjg(zdis)

         enddo
         istart = istart + nd(nbod)
      enddo
      zvort = zvort*4

cc      call prin2('const=*',wdens(istart+1),2)
cc      call prinf('istart+1=*',istart+1,1)

      if(k0.eq.1) then
         zvel = zvel + wdens(istart+1)
      endif

      return
      end

c------------------------------------------------------------------
