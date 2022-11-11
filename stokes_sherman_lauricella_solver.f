      subroutine lpcomp_stok_dl_interior(n,src,dsdt,rn,rkappa,dwgt,pwgt,
     1   eps,dens,vel)
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
      complex *16, intent(in) :: dens(n)
      complex *16, intent(out) :: vel(n)

      complex *16 grad(2),hess(3),veltarg,gradtarg(2),hesstarg(3)
      complex *16 charge
      complex *16, allocatable :: dip(:,:)
      complex *16 zu,zn,z1,z2,ztmp,ima,zx,zw
      integer nt
      real *8 targ(2)

      data ima/(0.0d0,1.0d0)/

      done = 1.0d0
      pi = atan(done)*4

      ifcharge = 0
      ifdipole = 1

      allocate(dip(2,n))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,zu,zn)      
      do i=1,n
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

      complex *16 zid,ztmp
      real *8 rb,wnrm2
      integer it,iind,it1,k,l
      real *8 rmyerr
      complex *16 temp
      complex *16, allocatable :: vmat(:,:),hmat(:,:)
      complex *16, allocatable :: cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:)

      allocate(vmat(n,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(n),svec(numit+1),yvec(numit+1))

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
      zid = -0.5d0*dwgt


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
        vmat(i,1) = rhs(i)/rb
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
          do j=1,n
            ztmp = ztmp + wtmp(j)*conjg(vmat(j,k))
          enddo
C$OMP END PARALLEL DO          
          hmat(k,it) = ztmp

C$OMP PARALLEL DO DEFAULT(SHARED) 
          do j=1,n
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
C$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:wnrm2)        
        do j=1,n
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
C$OMP END PARALLEL DO        
        wnrm2 = sqrt(wnrm2)

C$OMP PARALLEL DO DEFAULT(SHARED) 
        do j=1,n
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo
C$OMP END PARALLEL DO        

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+conjg(sn(k))*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        ztmp = wnrm2

        call zrotmat_gmres(hmat(it,it),ztmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+conjg(sn(it))*wnrm2
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
          do j=1,n
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo
C$OMP END PARALLEL DO          


          rres = 0
C$OMP PARALLEL DO DEFAULT(SHARED)          
          do i=1,n
            wtmp(i) = 0
          enddo
C$OMP END PARALLEL DO          
c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c

      call lpcomp_stok_dl_interior(n,src,dsdt,rn,rkappa,dwgt,pwgt,
     1   eps,soln,wtmp)

C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)            
          do i=1,n
            rres = rres + abs(zid*soln(i) + wtmp(i)-rhs(i))**2
          enddo
C$OMP END PARALLEL DO          
          rres = sqrt(rres)/rb
          niter = it
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

          


