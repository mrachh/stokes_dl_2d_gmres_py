      subroutine lpcomp_stok_dl_interior(n,src,dsdt,rn,rkappa,dwgt,pwgt,
     1   eps,rdens,rvel)
c
c   This subroutine applies the sherman lauricella version 
c   representation of the stokes double layer potential with an added
c   correction for the 1s matrix corresponding to the normals
c
c   Input arguments:
c     - n: integer
c         number of points
c     - src: real *8 (2,n)
c          x,y coordinates of source locations
c     - dsdt: real *8 (n)
c          quadrature weights for integrating smooth functions
c     - rn: real *8 (2,n)
c         normals
c     - rkappa: real *8 (n)
c         curvature
c     - dwgt: real *8
c          strength of double layer potential (recommended value: 1.0)
c     - pwgt: real *8
c          strength of 1s matrix (recommended value: 1.0)
c     - eps: real *8
c         precision for fmm
c     - dens: complex *16 (n)
c         input density at data points
c
c  Output arguments:
c     - vel: complex *16 (n)
c         velocity corresponding to sherman lauricella representation
c  
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: n
      real *8, intent(in) :: src(2,n),dsdt(n),rn(2,n),rkappa(n),dwgt
      real *8, intent(in) :: pwgt,eps
      real *8, intent(in) :: rdens(2*n)
      real *8, intent(out) :: rvel(2*n)

      complex *16 grad(2),hess(3),veltarg,gradtarg(2),hesstarg(3)
      complex *16 charge
      complex *16, allocatable :: dip(:,:),vel(:),dens(:)
      complex *16 zu,zn,z1,z2,ztmp,ima,zx,zw
      integer nt
      real *8 targ(2)

      data ima/(0.0d0,1.0d0)/

      done = 1.0d0
      pi = atan(done)*4

      ifcharge = 0
      ifdipole = 1

      allocate(dip(2,n),dens(n),vel(n))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,zu,zn)      
      do i=1,n 
        dens(i) = rdens(2*i-1) + ima*rdens(2*i)
        zu = -0.25d0*dens(i)*dsdt(i)/pi
        zn = rn(1,i) + ima*rn(2,i)
        dip(1,i) = -dwgt*zu*zn
        dip(2,i) = dwgt*(zn*dconjg(zu) - zu*dconjg(zn))
      enddo
C$OMP END PARALLEL DO      
      iper = 0
      ifpgh = 1
      ifpghtarg = 0
      nt = 0
      ier = 0
      call bhfmm2d(1,eps,n,src,ifcharge,charge,ifdipole,dip,iper,
     1  ifpgh,vel,grad,hess,nt,targ,ifpghtarg,veltarg,
     2  gradtarg,hesstarg,ier)
      velpert = 0.0d0
      ra = 0.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,zn) REDUCTION(+:ra,velpert)      
      do i=1,n
        zn = rn(1,i) + ima*rn(2,i)
        ra = ra + dsdt(i)
        velpert = velpert + dimag(conjg(zn)*dens(i))*dsdt(i) 
      enddo
C$OMP END PARALLEL DO      

      velpert = velpert/ra
      print *, "velpert=",velpert
c
c  add one's matrix and fix correction
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(zn,zw,zx,i)
      do i=1,n
        zn = rn(1,i) + ima*rn(2,i)
        zw = -0.25d0*dsdt(i)*rkappa(i)/pi
        zx = -0.25d0*dsdt(i)*rkappa(i)*zn*zn/pi
        vel(i) = vel(i) + pwgt*velpert*zn
        vel(i) = vel(i) + dwgt*( (zw+zx)*dreal(dens(i)) +
     1       ima*(zw-zx)*dimag(dens(i)))
        rvel(2*i-1) = dreal(vel(i))
        rvel(2*i) = dimag(vel(i))
      enddo
C$OMP END PARALLEL DO      
      return
      end


      subroutine stokes_gmres(n,src,dsdt,rn,rkappa,dwgt,pwgt,eps,
     1   rhs,numit,eps_gmres,niter,errs,rres,soln)
c
c   This code is a solver for the stokes interior velocity
c   problem using the
c   double layer Sherman lauricella representation 
c
c   Input arguments:
c     - n: integer
c         number of points
c     - src: real *8 (2,n)
c          x,y coordinates of source locations
c     - dsdt: real *8 (n)
c          quadrature weights for integrating smooth functions
c     - rn: real *8 (2,n)
c         normals
c     - rkappa: real *8 (n)
c         curvature
c     - dwgt: real *8
c          strength of double layer potential (recommended value: 1.0)
c     - pwgt: real *8
c          strength of 1s matrix (recommended value: 1.0)
c     - eps: real *8
c         precision for fmm
c     - rhs: complex * 16(n) 
c          boundary data: should be (-vy,vx)
c     - numit: integer
c         max number of itertations
c     - eps_gmres: real *8
c         precision for running gmres until
c
c
c  Output arguments:
c     - niter: integer
c          number of iterations used
c     - errs: real *8 (1:numit+1)
c         errs(1:niter) is the relative residual as a function of
c         iteration number
c     - rres: real *8
c         residual
c     - soln: complex *16 (n)
c         solution of linear system
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: n
      real *8, intent(in) :: src(2,n),dsdt(n),rn(2,n),rkappa(n),dwgt
      real *8, intent(in) :: pwgt,eps
      complex *16, intent(in) :: rhs(n)
      integer, intent(in) :: numit
      real *8, intent(in) :: eps_gmres
      integer, intent(out) :: niter
      real *8, intent(out) :: errs(numit+1),rres
      complex *16, intent(out) :: soln(n)

      real *8 did,ztmp
      real *8 rb,wnrm2
      integer it,iind,it1,k,l
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: vmat(:,:),hmat(:,:)
      real *8, allocatable :: cs(:),sn(:)
      real *8, allocatable :: svec(:),yvec(:),wtmp(:),dsoln(:)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/

      allocate(vmat(2*n,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(2*n),svec(numit+1),yvec(numit+1))
      allocate(dsoln(2*n))

      do i=1,numit+1
        errs(i) = 0
      enddo
c
c
c     start gmres code here
c
c     NOTE: matrix equation should be of the form (z*I + K)x = y
c       the identity scaling (z) is defined via zid below,
c       and K represents the action of the principal value 
c       part of the matvec
c
      did = -0.5d0*dwgt


      niter=0

c
c      compute norm of right hand side and initialize v
c 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo


c
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rb)
      do i=1,n
        rb = rb + abs(rhs(i))**2
      enddo
C$OMP END PARALLEL DO      
      rb = sqrt(rb)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,n
        vmat(2*i-1,1) = dreal(rhs(i))/rb
        vmat(2*i,1) = dimag(rhs(i))/rb
      enddo
C$OMP END PARALLEL DO      

      svec(1) = rb

      print *, "done prepping gmres, about to call lpcomp"
      print *, "rb=",rb

      do it=1,numit
        it1 = it + 1

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c

      call lpcomp_stok_dl_interior(n,src,dsdt,rn,rkappa,dwgt,pwgt,
     1   eps,vmat(1,it),wtmp)

        do k=1,it
          ztmp = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:ztmp)          
          do j=1,2*n
            ztmp = ztmp + wtmp(j)*vmat(j,k)
          enddo
C$OMP END PARALLEL DO          
          hmat(k,it) = ztmp

C$OMP PARALLEL DO DEFAULT(SHARED) 
          do j=1,2*n
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
C$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+did
        wnrm2 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:wnrm2)        
        do j=1,2*n
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
C$OMP END PARALLEL DO        
        wnrm2 = sqrt(wnrm2)

C$OMP PARALLEL DO DEFAULT(SHARED) 
        do j=1,2*n
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo
C$OMP END PARALLEL DO        

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+sn(k)*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        ztmp = wnrm2

        call rotmat_gmres(hmat(it,it),ztmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

c
c            solve the linear system corresponding to
c            upper triangular part of hmat to obtain yvec
c
c            y = triu(H(1:it,1:it))\s(1:it);
c
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo



c
c          estimate x
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
          do j=1,2*n
            dsoln(j) = 0
            do i=1,it
              dsoln(j) = dsoln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo
C$OMP END PARALLEL DO          


          rres = 0
C$OMP PARALLEL DO DEFAULT(SHARED)          
          do i=1,2*n
            wtmp(i) = 0
          enddo
C$OMP END PARALLEL DO          
c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c

      call lpcomp_stok_dl_interior(n,src,dsdt,rn,rkappa,dwgt,pwgt,
     1   eps,dsoln,wtmp)

C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)            
          do i=1,n
            rres = rres + abs(did*dsoln(2*i-1)+wtmp(2*i-1)
     1         -dreal(rhs(i)))**2
            rres = rres + abs(did*dsoln(2*i)+wtmp(2*i)
     1         -dimag(rhs(i)))**2
          enddo
C$OMP END PARALLEL DO          
          rres = sqrt(rres)/rb
          niter = it
          do i=1,n
            soln(i) = dsoln(2*i-1) + ima*dsoln(2*i)
          enddo
          return
        endif
      enddo

c
      return
      end
c
c
c
c

c
c  This file contains the following user-callable routines
c  
c    - rotmat_gmres: 
c        compute the sin and cos such that tan = b/a, with stability
c        close to 0,1
c    - zrotmat_gmres:
c        complex version of rotmat_gmres
c



      subroutine rotmat_gmres(a,b,c,s)
c
c-----------------------
c  Given a,b, compute sin(theta), cos(theta),
c  such that tan(theta) = b/a, note, the routine
c  implements a stabler version of computing
c  b/sqrt(a^2+b^2) and a/sqrt(b^2+a^2)
c
c  Input arguments:
c  
c    - a: double precision
c        cos scaling
c    - b: double precision
c        sin scaling
c
c  Output arguments:
c    - c: double precision
c        cos(theta)
c    - s: double precision
c        sin(theta)
c        
c
c-----------------------
c
      implicit real *8 (a-h,o-z)
      real *8 a,b,c,s

      if(a.eq.0) then
        c = 0
        s = 1
      else
        temp = b/a
        c = 1.0d0/sqrt(1.0d0+temp**2)
        s = temp*c

      endif

      return
      end

          



      subroutine zrotmat_gmres(a,b,c,s)
c-----------------------
c  Given a,b, compute sin(theta), cos(theta),
c  such that tan(theta) = b/a, note, the routine
c  implements a stabler version of computing
c  b/sqrt(a^2+b^2) and a/sqrt(b^2+a^2).
c
c  This routine is the complex version of rotmat_gmres
c
c  Input arguments:
c  
c    - a: double complex
c        cos scaling
c    - b: double complex 
c        sin scaling
c
c  Output arguments:
c    - c: double complex 
c        cos(theta)
c    - s: double complex 
c        sin(theta)
c        
c
c-----------------------
      implicit real *8 (a-h,o-z)
      complex *16 a,b,c,s,temp

      if(a.eq.0) then
        s = 1
        c = 0
      else
        temp = b/a
        c = 1.0d0/sqrt(1.0d0+abs(temp)**2)
        s = temp*c
      endif

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
