"""CFD city Flow Final Project"""



from dolfinx.fem.petsc import LinearProblem

import meshio
import basix
from mpi4py import MPI
import dolfinx as dfx
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.io import gmshio
import ufl as ufl
import matplotlib.pyplot as plt
import numpy as np
from petsc4py import PETSc

#------------------------------------------------------------------------------
# definition of problem parameters and inputs

#------------------------------------------------------------------------------
viscosity = 1.0
density = 1.0
outFileV = 'CityFLow2d_V.xdmf'
outFileP = 'CityFlow2d_P.xdmf'
reynolds = 100.0


#------------------------------------------------------------------------------
dt = 0.005
t_start = 0.0
t_end = 2.0
t_theta = 0.5
refLength= 1.0 

#Input wind speed on given wall
#Calculate CFL number fo reference
#-------------------------------------------------------------------------------
#avgWindSpeed =   (reynolds*refLength)/viscosity
avgWindSpeed =   20
#cfl_estimate = ( lidVelocity * dt )/(1.0/Nx)
print('Problem Reynolds Number Input:', reynolds)
#print('Problem Bulk CFL Estimate:', cfl_estimate)

#Define BC types
'''def inletBC(x):
    vals = np.zeros((3, x.shape[1]))
    delta = domainHeight  # boundary layer thickness in Z-direction 
    z = x[2] #x3 = height 
    n = 0.75 # can be changed based on how we want BL curve to look like
    vals[0] = avgWindSpeed * (z / delta)**n  # velocity increases with distance from the wall at z = 0
    vals[1] = 0.0  # no motion in x2-direction
    vals[2] = 0.0  # no motion in x3-direction 
    return vals'''



def noSlipBC(x):
    return np.stack((np.zeros(x.shape[1]), np.zeros(x.shape[1])))

def windBC(x):
    return np.stack((np.full(x.shape[1], avgWindSpeed), np.zeros(x.shape[1])))

#def outletBC(x):
 #   return np.stack((np.zeros(x.shape[1])))


def outletBC(x):
    return np.zeros(x.shape[1])
#ID key surfaces

ID_INLET = 103
ID_OUTLET = 101
ID_EAST = 102
ID_WEST = 104
ID_BUILDINGS = 200

'''
B1_N_id = 33
B1_S_id = 31
B1_E_id = 32
B1_W_id = 34
B1_Sk_id = 35'''

#input mesh file

meshName= '2BLDG_2D_V3.msh'

mesh, cell_markers, facet_markers = gmshio.read_from_msh(meshName, MPI.COMM_WORLD, gdim=2)
#checking BC facet IDS
for name, id_ in zip(["inlet", "outlet", "buildings"], [ID_INLET, ID_OUTLET, ID_BUILDINGS]):
    found = facet_markers.find(id_)
    print(f"{name.upper()} (ID={id_}): Found {len(found)} boundary facets")

#pull fdim and tdim from mesh

tdim = mesh.topology.dim
#mesh2 = meshio.read(meshName)
#print('tdim:', tdim)
#print(mesh2.cells_dict.keys())
fdim = tdim - 1
#mesh.topology.create_connectivity(fdim, tdim)
#Establish Function Space
PE = basix.ufl.element('CG', mesh.basix_cell(), degree=1)
QE = basix.ufl.element('CG', mesh.basix_cell(), degree=1, shape=(mesh.topology.dim,))


ME = basix.ufl.mixed_element([QE, PE])
W = dfx.fem.functionspace(mesh, ME)

#print(len(np.where(facet_tags.values == 1)[0])) 
#-----------------------------------------------------------------------
# extracting the subspace and their mapping to the mixed space is needed
# so that boundary conditions can be appropriately assigned
#-----------------------------------------------------------------------

U_sub, U_submap = W.sub(0).collapse()
P_sub, P_submap = W.sub(1).collapse()

#-------------------------------------------------------------------------------
# Assign IDS to BCs
#We are using cardinal directions (N, S, E, W) to identify the walls. Each building is a single BC
#-------------------------------------------------------------------------------
#b_dofs_In = dfx.fem.locate_dofs_topological((V.sub(0),U_sub), fdim, facet_markers.find(ID_In))

# Use updated facet marker IDs matching Physical Line IDs in your .geo file
ID_INLET = 103   # North
ID_OUTLET = 101    # South
ID_BUILDINGS = 200 # Shared group for both buildings
ID_EAST = 102    # East
ID_WEST = 104    # West

# Locate degrees of freedom (DoFs) on marked facets
b_dofs_inlet = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, facet_markers.find(ID_INLET))
b_dofs_outlet = dfx.fem.locate_dofs_topological((W.sub(1), P_sub), fdim, facet_markers.find(ID_OUTLET))
b_dofs_buildings = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, facet_markers.find(ID_BUILDINGS))
b_dofs_East = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, facet_markers.find(ID_EAST))
b_dofs_West = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, facet_markers.find(ID_WEST))
''''
b_dofs_B1_S = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, ID_B1_S)
b_dofs_B1_E = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, ID_B1_E)
b_dofs_B1_W = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, ID_B1_W)
b_dofs_B1_Sk = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, ID_B1_Sk)
'''
#---------------------------------------------------------------------------
# syntax for defining the 3 different uD values to be used
# recall: u = uD for all x in Gamma_D (this is our template for implementing
# all the Dirichlet conditions)
#---------------------------------------------------------------------------

#We will make the ground and building walls no slip
print(U_sub)
uD_noSlip = dfx.fem.Function(U_sub)
uD_noSlip.interpolate(noSlipBC)

#Inlet will have prescribed velocity
uD_Inlet = dfx.fem.Function(U_sub)
uD_Inlet.interpolate(windBC)

#Outle will have zero pressure
uD_Outlet = dfx.fem.Function(P_sub)
uD_Outlet.interpolate(outletBC)

print("Inlet velocity BC norm:", np.linalg.norm(uD_Inlet.x.array))
print("No-slip BC velocity norm (should be 0):", np.linalg.norm(uD_noSlip.x.array))
print("Outlet pressure BC norm (should be 0):", np.linalg.norm(uD_Outlet.x.array))

#Remainder of surfaces will be "do nothing" BCS


#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
bc_inlet = dfx.fem.dirichletbc(uD_Inlet, b_dofs_inlet, W.sub(0))
bc_outlet = dfx.fem.dirichletbc(uD_Outlet, b_dofs_outlet, W.sub(1))
bc_buildings = dfx.fem.dirichletbc(uD_noSlip, b_dofs_buildings, W.sub(0))
bc_EAST=dfx.fem.dirichletbc(uD_noSlip, b_dofs_East, W.sub(0))
bc_WEST=dfx.fem.dirichletbc(uD_noSlip, b_dofs_West, W.sub(0))


bc = [bc_inlet, bc_outlet, bc_buildings,bc_EAST, bc_WEST]

#---------------------------------------------------------------------
# define the trial and test functions based on the mixedfunctionspace
# NOTE: the call to TrialFunctions and not TrialFunction etc.
#---------------------------------------------------------------------
(v,q) = ufl.TestFunctions(W)

#-------------------------------------------------------------------
# defining all physics properties and discretization constants into
# standard dolfinx constant parameters
#-------------------------------------------------------------------
mu = dfx.fem.Constant(mesh, dfx.default_scalar_type(viscosity))
rho = dfx.fem.Constant(mesh, dfx.default_scalar_type(density))
idt = dfx.fem.Constant(mesh, dfx.default_scalar_type(1.0/dt))
theta = dfx.fem.Constant(mesh, dfx.default_scalar_type(t_theta))
b = dfx.fem.Constant(mesh, PETSc.ScalarType((0.0,0.0)))

#-------------------------------------------------------------------------------
# Theta-Galerkin:
#
# define the variational form terms without time derivative in current timestep
# in theta-Galerkin formulation this is corresponding to the t_n+1
#-------------------------------------------------------------------------------
W1 = dfx.fem.Function(W)
(u,p) = ufl.split(W1)

T1_1 = rho * ufl.inner(v, ufl.grad(u)*u) * ufl.dx
T2_1 = mu * ufl.inner(ufl.grad(v), ufl.grad(u)) * ufl.dx
T3_1 = p * ufl.div(v) * ufl.dx
T4_1 = q * ufl.div(u) * ufl.dx
T5_1 = rho * ufl.dot(v,b) * ufl.dx
L_1  = T1_1 + T2_1 - T3_1 -T4_1 - T5_1

#-------------------------------------------------------------------------------
# Theta-Galerkin:
#
# define the variational form terms without time derivative in previous timestep
# in theta-Galerkin formulation this is corresponding to the t_n
#
# note: however, as a subtle trick, we will retain the pressure to be the same
# as the one defined for t_n+1. this is okay to do considering the fact that
# in reality, the pressure equation is functioning to keep the velocity
# divergence to 0, and the divergence equation does not have a time-derivative
#-------------------------------------------------------------------------------
W0 = dfx.fem.Function(W)
(u0,p0) = ufl.split(W0)

T1_0 = rho * ufl.inner(v, ufl.grad(u0)*u0) * ufl.dx
T2_0 = mu * ufl.inner(ufl.grad(v), ufl.grad(u0)) * ufl.dx
T3_0 = p * ufl.div(v) * ufl.dx
T4_0 = q * ufl.div(u0) * ufl.dx
T5_0 = rho * ufl.dot(v,b) * ufl.dx
L_0 = T1_0 + T2_0 - T3_0 -T4_0 - T5_0

#--------------------------------------------------------------------------
# Theta-Galerkin:
#
# combine variational forms with time derivative as discussed for the
# complete theta-Galerkin formulation with a one-step discretization of the
# time-derivative
#
#  dw/dt + L(t) = 0 is approximated as
#  (w-w0)/dt + (1-theta)*L(t0) + theta*L(t) = 0
#---------------------------------------------------------------------------
F = idt * rho * ufl.inner((u-u0),v) * ufl.dx + (1.0-theta) * L_0 + theta * L_1

#-------------------------------------------------------------------------------
# defining the stabilization parameter for the Petrov-Galerkin stabilization
#-------------------------------------------------------------------------------
uNorm = ufl.sqrt(ufl.inner(u0, u0))
h = ufl.CellDiameter(mesh)
tau = ( (2.0*theta*idt)**2 + (2.0*uNorm/h)**2 + (4.0*mu/h**2)**2 )**(-0.5)

#----------------------------------------------------------------------
# defining the complete residual for the Navier-Stokes momentum balance
#----------------------------------------------------------------------
residual = idt*rho*(u - u0) + \
    theta*(rho*ufl.grad(u)*u - mu*ufl.div(ufl.grad(u)) + ufl.grad(p) - rho*b) +\
    (1.0-theta)*(rho*ufl.grad(u0)*u0 - mu*ufl.div(ufl.grad(u0)) + ufl.grad(p) - rho*b)

#---------------------------------------------------------------------------
# including now the contributions from the SUPG and PSPG stabilization terms
# into the overall theta-Galerkin weak form
#---------------------------------------------------------------------------
F_SUPG = tau * ufl.inner(ufl.grad(v)*u, residual) * ufl.dx
F_PSPG = - tau * ufl.inner(ufl.grad(q), residual) * ufl.dx

F = F + F_SUPG + F_PSPG

#--------------------------------------------------------------------
# define and configure the details of a Newton Solver for the overall
# NonlinearProblem defined in F above, set as F = 0
#--------------------------------------------------------------------
problem = NonlinearProblem(F, W1, bcs=bc)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-7
solver.report = True

ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "gmres"
opts[f"{option_prefix}pc_type"] = "gamg"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "umfpack"
ksp.setFromOptions()

t = t_start
tn = 0

#-------------------------------------------------------------------------------
# set up the output files for velocity and pressure solutions
# if you choose to write into a pvd file, multiple files - one for each time
# step will be generated. if you choose to write into an xdmf file - only one
# combined output file for velocity (and one for pressure) will be generated
#-------------------------------------------------------------------------------
if outFileV.endswith('pvd'):
    vFile = dfx.io.VTKFile(MPI.COMM_WORLD, outFileV, "w")
elif outFileV.endswith('xdmf'):
    vFile = dfx.io.XDMFFile(mesh.comm, outFileV, "w", encoding=dfx.io.XDMFFile.Encoding.ASCII)
    vFile.write_mesh(mesh)

if outFileP.endswith('pvd'):
    pFile = dfx.io.VTKFile(MPI.COMM_WORLD, outFileP, "w")
elif outFileP.endswith('xdmf'):
    pFile = dfx.io.XDMFFile(mesh.comm, outFileP, "w", encoding=dfx.io.XDMFFile.Encoding.ASCII)
    pFile.write_mesh(mesh)

dfx.log.set_log_level(dfx.log.LogLevel.INFO)

#----------------------------------------------------------------------------
# loop through time to solve for the updated version of the solution variable
# it is important to note the structure of the loop as it is implemented:
# - we identify which time step we are at
# - we solve for the updated solution u_n+1, p_n+1
# - note that the u,p update itself is an iteration (Newton method)
# - we set c_n+1 as the c_n for the next time step
# - we update the time step
#
# The INNER LOOP: Newton method iteration
# The OUTER LOOP: time step advancement
#-----------------------------------------------------------------------------
while t < t_end:

    print("t=",t)
    n, converged = solver.solve(W1)
    assert (converged)
    print(f"Number of iterations: {n:d}")

    uf = W1.split()[0].collapse()
    pf = W1.split()[1].collapse()
    print(f"t = {t:.3f} | Velocity norm: {np.linalg.norm(uf.x.array):.6f} | Pressure norm: {np.linalg.norm(pf.x.array):.6f}")

    uf.name = 'vel'
    pf.name = 'pres'

    if outFileV.endswith('pvd'):
        vFile.write_function(uf, tn)
    elif outFileV.endswith('xdmf'):
        vFile.write_function(uf, t)

    if outFileP.endswith('pvd'):
        pFile.write_function(pf, tn)
    elif outFileP.endswith('xdmf'):
        pFile.write_function(pf, t)

    W0.x.array[:] = W1.x.array

    t += dt
    tn += 1

vFile.close()
pFile.close()
print("Simulation finished. Final time =", t)

