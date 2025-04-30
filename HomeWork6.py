import basix
from mpi4py import MPI
import dolfinx as dfx
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.io import gmshio
from petsc4py import PETSc
import ufl as ufl
import matplotlib.pyplot as plt
import numpy as np
import time

#-------------------------------------------------------------------------------
# FOR ASSIGNMENT:
# update and implement the integrateFluidStress function
#-------------------------------------------------------------------------------
def integrateFluidStress(a_U,a_P,a_Mu,a_N,a_Mesh,a_GammaP):
    eps = 0.5*(ufl.grad(a_U)+ufl.grad(a_U).T)
    sig = -a_P*ufl.Identity(2)+2.0*a_Mu*eps

    traction = ufl.dot(sig,a_N)

    forceX = traction[0] * a_GammaP
    forceY = traction[1] * a_GammaP

    fXVal = dfx.fem.assemble_scalar(dfx.fem.form(forceX))
    fYVal = dfx.fem.assemble_scalar(dfx.fem.form(forceY))
    return [fXVal, fYVal]

#-------------------------------------------------------------------------------
# FOR ASSIGNMENT:
# Meshing: make sure to change the file names below
#-------------------------------------------------------------------------------
simulationType = 'ellipse'

if simulationType == 'square':
    meshFile = 'squareMesh.msh'
    outFileV = 'square-V.xdmf'
    outFileP = 'square-P.xdmf'
    forceFile = 'square-forces.dat'
    imageFile = 'square-draglift.png'
elif simulationType == 'ellipse':
    meshFile = 'parabollaMesh.msh'
    outFileV = 'ellipse-V.xdmf'
    outFileP = 'ellipse-P.xdmf'
    forceFile = 'ellipse-forces.dat'
    imageFile = 'ellipse-draglift.png'
elif simulationType == 'circle':
    meshFile = 'circleMesh.msh'
    outFileV = 'circle-V.xdmf'
    outFileP = 'circle-P.xdmf'
    forceFile = 'circle-forces.dat'
    imageFile = 'circle-draglift.png'


viscosity = 0.001
density = 1.0
U0 = 1.5
diam = 0.1
boxH = 4.1 * diam
boxL = 22.0 * diam
elemType = 'q1p1'
isStabilize = True

particleFlow = True
injectEvery = 50
particleFile = 'seeds.csv'

dt = 0.005
t_start = 0.0
t_end = 7.0
t_theta = 0.5

#-------------------------------------------------------------------------------
# FOR ASSIGNMENT:
# Meshing: make sure to update/modify the ids as below for your mesh/geometry
#-------------------------------------------------------------------------------
ID_INLET = 12
ID_TOP = 11
ID_OUTLET = 10
ID_BOTTOM = 9
ID_CYL = 13

Reynolds = density*(2.0/3.0)*U0*diam/viscosity
print("The Problem Reynolds Number is:", Reynolds)

def noSlipBC(x):
    return np.stack((np.zeros(x.shape[1]), np.zeros(x.shape[1])))

def pressureBC(x):
    return np.zeros(x.shape[1])

def inletBC(x):
    vals = np.zeros((mesh.geometry.dim, x.shape[1]))
    vals[0] = 4.0 * U0 * x[1] * (boxH - x[1]) / ( boxH * boxH )
    vals[1] = 0.0
    return vals

startTime = time.time()

mesh, cell_markers, facet_markers = gmshio.read_from_msh(meshFile, MPI.COMM_WORLD, gdim=2)

nVec = ufl.FacetNormal(mesh)
tdim = mesh.topology.dim
fdim = tdim - 1

if particleFlow == True:
    #Reading in initial particle locations from csv file
    particles = np.genfromtxt(particleFile,delimiter=',',skip_header=1,dtype=np.float64)
    numParticles = particles.shape[0]
    print("number of particles", numParticles)
    seedParticles = particles.copy()

#setting up mixed function space
PE = basix.ufl.element('CG', mesh.basix_cell(), degree=1)

if elemType == 'q2p1':
    QE = basix.ufl.element('CG', mesh.basix_cell(), degree=2, shape=(mesh.topology.dim,))
elif elemType == 'q1p1':
    QE = basix.ufl.element('CG', mesh.basix_cell(), degree=1, shape=(mesh.topology.dim,))

ME = basix.ufl.mixed_element([QE, PE])
W = dfx.fem.functionspace(mesh, ME)

U_sub, U_submap = W.sub(0).collapse()
P_sub, P_submap = W.sub(1).collapse()

#-------------------------------------------------------------------------------
# FOR ASSIGNMENT:
# Meshing: make sure to update/modify the boundary dofs as below, for your
# mesh/geometry in case if needed
#-------------------------------------------------------------------------------
b_dofs_INLET = dfx.fem.locate_dofs_topological((W.sub(0),U_sub), fdim, facet_markers.find(ID_INLET))
b_dofs_TOP = dfx.fem.locate_dofs_topological((W.sub(0),U_sub), fdim, facet_markers.find(ID_TOP))
b_dofs_OUTLET = dfx.fem.locate_dofs_topological((W.sub(1),P_sub), fdim, facet_markers.find(ID_OUTLET))
b_dofs_BOTTOM = dfx.fem.locate_dofs_topological((W.sub(0),U_sub), fdim, facet_markers.find(ID_BOTTOM))
b_dofs_CYL = dfx.fem.locate_dofs_topological((W.sub(0),U_sub), fdim, facet_markers.find(ID_CYL))

uD_Wall = dfx.fem.Function(U_sub)
uD_Wall.interpolate(noSlipBC)

uD_Inlet = dfx.fem.Function(U_sub)
uD_Inlet.interpolate(inletBC)

uD_Outlet = dfx.fem.Function(P_sub)
uD_Outlet.interpolate(pressureBC)

bc_INLET = dfx.fem.dirichletbc(uD_Inlet, b_dofs_INLET, W.sub(0))
bc_TOP = dfx.fem.dirichletbc(uD_Wall, b_dofs_TOP, W.sub(0))
bc_OUTLET = dfx.fem.dirichletbc(uD_Outlet, b_dofs_OUTLET, W.sub(1))
bc_BOTTOM = dfx.fem.dirichletbc(uD_Wall, b_dofs_BOTTOM, W.sub(0))
bc_CYL = dfx.fem.dirichletbc(uD_Wall, b_dofs_CYL, W.sub(0))

bc = [bc_INLET, bc_TOP, bc_OUTLET, bc_BOTTOM, bc_CYL]

ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_markers)
Gamma_CYL = ds(ID_CYL)

(v,q) = ufl.TestFunctions(W)

mu = dfx.fem.Constant(mesh, dfx.default_scalar_type(viscosity))
rho = dfx.fem.Constant(mesh, dfx.default_scalar_type(density))
idt = dfx.fem.Constant(mesh, dfx.default_scalar_type(1.0/dt))
theta = dfx.fem.Constant(mesh, dfx.default_scalar_type(t_theta))
b = dfx.fem.Constant(mesh, PETSc.ScalarType((0.0,0.0)))
I = ufl.Identity(2)

W1 = dfx.fem.Function(W)
(u,p) = ufl.split(W1)

T1_1 = rho * ufl.inner(v, ufl.grad(u)*u) * ufl.dx
T2_1 = mu * ufl.inner(ufl.grad(v), ufl.grad(u)) * ufl.dx
T3_1 = p * ufl.div(v) * ufl.dx
T4_1 = q * ufl.div(u) * ufl.dx
T5_1 = rho * ufl.dot(v,b) * ufl.dx
L_1  = T1_1 + T2_1 - T3_1 -T4_1 - T5_1

W0 = dfx.fem.Function(W)
(u0,p0) = ufl.split(W0)

T1_0 = rho * ufl.inner(v, ufl.grad(u0)*u0) * ufl.dx
T2_0 = mu * ufl.inner(ufl.grad(v), ufl.grad(u0)) * ufl.dx
T3_0 = p * ufl.div(v) * ufl.dx
T4_0 = q * ufl.div(u0) * ufl.dx
T5_0 = rho * ufl.dot(v,b) * ufl.dx
L_0 = T1_0 + T2_0 - T3_0 -T4_0 - T5_0

F = idt * ufl.inner((u-u0),v) * ufl.dx + (1.0-theta) * L_0 + theta * L_1

uNorm = ufl.sqrt(ufl.inner(u0, u0))
h = ufl.CellDiameter(mesh)
tau = ( (2.0*theta*idt)**2 + (2.0*uNorm/h)**2 + (4.0*mu/h**2)**2 )**(-0.5)

residual = idt*rho*(u - u0) + \
    theta*(rho*ufl.grad(u)*u - mu*ufl.div(ufl.grad(u)) + ufl.grad(p) - rho*b) +\
    (1.0-theta)*(rho*ufl.grad(u0)*u0 - mu*ufl.div(ufl.grad(u0)) + ufl.grad(p) - rho*b)

F_SUPG = tau * ufl.inner(ufl.grad(v)*u, residual) * ufl.dx

if isStabilize == True:
    F_PSPG = - tau * ufl.inner(ufl.grad(q), residual) * ufl.dx

F = F + F_SUPG + F_PSPG

problem = NonlinearProblem(F, W1, bcs=bc)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-7
solver.report = True

ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "gmres"
opts[f"{option_prefix}pc_type"] = "lu"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "umfpack"
ksp.setFromOptions()

t = t_start
tn = 0

vFile = dfx.io.XDMFFile(mesh.comm, outFileV, "w", encoding=dfx.io.XDMFFile.Encoding.ASCII)
vFile.write_mesh(mesh)

pFile = dfx.io.XDMFFile(mesh.comm, outFileP, "w", encoding=dfx.io.XDMFFile.Encoding.ASCII)
pFile.write_mesh(mesh)

time_Arr = []
fX_Arr = []
fY_Arr = []
cD_Arr = []
cL_Arr = []

if particleFlow == True:
    pt = np.zeros(mesh.geometry.dim,dtype=np.float64)
    tree = dfx.geometry.bb_tree(mesh, mesh.topology.dim)

while t < t_end:

    t_in = time.time()
    n, converged = solver.solve(W1)
    assert (converged)
    t_out = time.time()

    print(f"t = {t:.6f}; Number of iterations: {n:d}; compute time: {t_out-t_in:f}")

    uf = W1.split()[0].collapse()
    pf = W1.split()[1].collapse()

    uf.name = 'vel'
    pf.name = 'pres'

    if particleFlow == True:
        vFile.write_function(uf, tn)
        pFile.write_function(pf, tn)
    else:
        vFile.write_function(uf, t)
        pFile.write_function(pf, t)


    #---------------------------------------------------------------------------
    # FOR ASSIGNMENT:
    # implement the call to the stress integration function, and the
    # calculation of drag and lift coefficients
    #---------------------------------------------------------------------------
    [fX,fY] = integrateFluidStress(uf,pf,mu,nVec, mesh, Gamma_CYL)
    cD = fX*2.0/((2.0*U0/3.0)**2*diam)
    cL = fY*2.0/((2.0*U0/3.0)**2*diam)

    time_Arr.append(t)
    fX_Arr.append(fX)
    fY_Arr.append(fY)
    cD_Arr.append(cD)
    cL_Arr.append(cL)

    #-------------------------------------------------------------------------
    #Implementing Lagrangian Particle tracing
    #-------------------------------------------------------------------------
    if particleFlow == True:
        if((t >t_start) and ((tn % injectEvery)==0)):
            particles = np.concatenate((particles,seedParticles),axis=0)
            numParticles = particles.shape[0]
        
        for i in range(numParticles):
            pt = particles[i,:]
            cells = dfx.geometry.compute_collisions_points(tree, pt)
            colliding_cells = dfx.geometry.compute_colliding_cells(mesh,cells,pt)
            cell_candidates = colliding_cells.links(0)

            if len(cell_candidates) > 0:
                cell = cell_candidates[0]
                vel_ploc = uf.eval(pt,cell)
                particles[i,0] = particles[i,0] +vel_ploc[0]*dt
                particles[i,1] = particles[i,1] + vel_ploc[1] *dt
            else: 
                print(f"Point {pt} is not in domain.")
        
        lagrangianFile = particleFile.split('.')[0]+'-'+str(tn)+'.csv'
        np.savetxt(lagrangianFile,particles,delimiter=',', header="x, y, z")


    W0.x.array[:] = W1.x.array
    t += dt
    tn += 1

vFile.close()
pFile.close()

forceData = np.column_stack([time_Arr, fX_Arr, fY_Arr, cD_Arr, cL_Arr])
np.savetxt(forceFile, forceData)

endTime = time.time()

print('Total simulation time:', endTime - startTime)

fig, ax = plt.subplots(1,2)
ax[0].plot(time_Arr[-200:], cD_Arr[-200:], 'r')
ax[1].plot(time_Arr[-200:], cL_Arr[-200:], 'm')
ax[0].set_xlabel('time', fontweight='bold')
ax[1].set_xlabel('time', fontweight='bold')
ax[0].set_ylabel('streamwise force', fontweight='bold')
ax[1].set_ylabel('cross-stream force', fontweight='bold')
plt.savefig(imageFile, bbox_inches='tight',dpi=120)
plt.close()