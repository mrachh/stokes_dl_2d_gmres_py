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
      enddo

      call get_bdry_data(n,src,cxsource,cysource,rhs)
      velpert = 0
      do i=1,n
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
      
      call eval_vel(xp,yp,n,src,dsdt,rn,dwgt,soln,zvel)

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
