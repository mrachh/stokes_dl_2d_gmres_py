import numpy as np
import numpy.linalg as la
import stokesdl as stk

n = 256
# set major and minor axes of the ellipse
a = 1.2
b = 1.3

# set exterior source location to get boundary data
cx = -3.05
cy = 1.1


# set interior target location
xt = 0.2
yt = 0.3

# Set weight for double layer potential
dwgt = 2.0

# Set weight for 1s matrix
pwgt = 1.0

# Set precision for layer potential
eps = 0.51e-9

# set precision for gmres
eps_gmres = 0.5e-11

# set max number of iterations
numit = 40

# Set source info needed for solvers
thet = np.linspace(0,2*np.pi,n,endpoint=False)
ct = np.cos(thet)
st = np.sin(thet)
src = np.zeros((2,n),order="F")
dsdt = np.zeros((n,),order="F")
rn = np.zeros((2,n),order="F")
rkappa = np.zeros((n,),order="F")

src[0,:] = a*ct
src[1,:] = b*st

dsdt = np.sqrt((a*st)**2 + (b*ct)**2)
rn[0,:] = b*ct/dsdt
rn[1,:] = a*st/dsdt

rkappa = a*b/dsdt**3
dsdt = dsdt*2*np.pi/n

rhs = stk.get_bdry_data(src,cx,cy)
niter,errs,rres,soln = stk.stokes_gmres(src,dsdt,rn,rkappa,dwgt,pwgt,eps,rhs,numit,eps_gmres)

zvel = stk.eval_vel(xt,yt,src,dsdt,rn,dwgt,soln)
targ = [xt,yt]
zvel_ex = stk.get_bdry_data(targ,cx,cy)
# note zvel_ex = -vy + 1j*vx, so need to multiply buy -1j
zvel_ex = -1j*zvel_ex

erra = abs(zvel-zvel_ex)/abs(zvel_ex)
print("error in velocity=",str(erra[0]))
