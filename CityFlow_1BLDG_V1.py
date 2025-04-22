
"""CFD city Flow Final Project"""
""""""
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
outFileV = 'NSE-cav-V.xdmf'
outFileP = 'NSE-cav-P.xdmf'
reynolds = 100.0


#------------------------------------------------------------------------------
dt = 0.001
t_start = 0.0
t_end = 0.05
t_theta = 0.5

#-------------------------------------------------------------------------------
# we will now compute two entities:
# - first: we set up the velocity of the top wall, or lid using the Reynolds
# number set at the beginning of the problem. this helps us parameterize the
# simulation such that we can vary the Reynolds number and re-run cases
#
# - second: we compute what is an estimate of a bulk CFL number based on the
# computed velocity of the top wall, and a bulk mesh size estimate. this does
# not directly give us an equivalent validation of whether the classic form of
# the CFL condition is met, but it gives us an indicator of how high/low is
# our CFL number estimate for the system overall
#-------------------------------------------------------------------------------
lidVelocity = (reynolds*refLength)/viscosity
cfl_estimate = ( lidVelocity * dt )/(1.0/Nx)
print('Problem Reynolds Number Input:', reynolds)
print('Problem Bulk CFL Estimate:', cfl_estimate)

#-------------------------------------------------------------------------------
# since we are not loading an external mesh with Physical Curve and Surface ids
# marked using Gmsh - here we will have to create functions that help us mark
# all the appropriate boundaries
#
# we will successively define the left, top, right, and bottom boundaries
# of the flow domain comprising the cavity
#-------------------------------------------------------------------------------
def left(x):
    return np.isclose(x[0], 0.0)

def right(x):
    return np.isclose(x[0], 1.0)

def top(x):
    return np.isclose(x[1], 1.0)

def bottom(x):
    return np.isclose(x[1], 0.0)

#-------------------------------------------------------------------------------
# we also identify the left corner of the cavity as a point of interest.
# we will actually use this point to set up a specific condition that the
# pressure is fixed at that point.
# this is a special example of a Dirichlet type condition or a constraint that
# helps set the reference/gauge for the pressure values (since in an
# incompressible flow, the pressure is only determinable up to a constant)
#-------------------------------------------------------------------------------
def leftCorner(x):
    return np.logical_and(np.isclose(x[0], 0.0), np.isclose(x[1], 0.0))

#-------------------------------------------------------------------------------
# define a function that helps set the no slip and no permeation boundary
# conditions at the walls of the cavity
#-------------------------------------------------------------------------------
def noSlipBC(x):
    return np.stack((np.zeros(x.shape[1]), np.zeros(x.shape[1])))

#-------------------------------------------------------------------------------
# define a function that helps set the no slip and no permeation boundary
# conditions specifically at the top wall or lid of the cavity. note that
# this is still a no slip condition - however, because the wall is moving
# the no slip condition results in the flow velocity matching the lid velocity
#-------------------------------------------------------------------------------
def lidBC(x):
    return np.stack((np.full(x.shape[1], lidVelocity), np.zeros(x.shape[1])))

#---------------------------------------------------------------------
# define a function that helps set a zero-pressure boundary condition
#---------------------------------------------------------------------
def pressureBC(x):
    return np.zeros(x.shape[1])

#------------------------------------------------------------------------
# we will now create a mesh using the built-in mesh generator utilities
# by changing the CellType in the last entry in the mesh generation call,
# we can change this into a structured mesh. That is:
# dfx.mesh.CellType.quadrilateral -- for structured
#------------------------------------------------------------------------
if meshType == 'tri':

    mesh = dfx.mesh.create_rectangle(MPI.COMM_WORLD, \
        [np.array([0.0, 0.0]), np.array([1.0, 1.0])], [Nx, Ny], dfx.mesh.CellType.triangle)

elif meshType == 'quad':

    mesh = dfx.mesh.create_rectangle(MPI.COMM_WORLD, \
        [np.array([0.0, 0.0]), np.array([1.0, 1.0])], [Nx, Ny], dfx.mesh.CellType.quadrilateral)

tdim = mesh.topology.dim
fdim = tdim - 1

#-----------------------------------------------------------------------
# the following outlines the syntax for defining a mixed function space
#-----------------------------------------------------------------------
PE = basix.ufl.element('CG', mesh.basix_cell(), degree=1)

if elemType == 'q2p1':
    QE = basix.ufl.element('CG', mesh.basix_cell(), degree=2, shape=(mesh.topology.dim,))
elif elemType == 'q1p1':
    QE = basix.ufl.element('CG', mesh.basix_cell(), degree=1, shape=(mesh.topology.dim,))

ME = basix.ufl.mixed_element([QE, PE])
W = dfx.fem.functionspace(mesh, ME)

#-----------------------------------------------------------------------
# extracting the subspace and their mapping to the mixed space is needed
# so that boundary conditions can be appropriately assigned
#-----------------------------------------------------------------------
U_sub, U_submap = W.sub(0).collapse()
P_sub, P_submap = W.sub(1).collapse()

#-------------------------------------------------------------------------------
# use the boundary definition functions defined previously in the code to
# generate all the ids for the appropriate topological entities for the problem
# - note that the overall domain (Omega) is a 2D box (dim = 2)
# - hence the overall boundary (\delta Omega) is a set of curves (dim = 1)
# - and the lef corner is a point (dim = 0)
# hence we identify the curves as mesh.topology.dim - 1
# and the point as mesh.topology.dim - 2
#-------------------------------------------------------------------------------
ID_L = dfx.mesh.locate_entities_boundary(mesh, fdim, left)
ID_R = dfx.mesh.locate_entities_boundary(mesh, fdim, right)
ID_T = dfx.mesh.locate_entities_boundary(mesh, fdim, top)
ID_B = dfx.mesh.locate_entities_boundary(mesh, fdim, bottom)
ID_C = dfx.mesh.locate_entities_boundary(mesh, fdim-1, leftCorner)

#-------------------------------------------------------------------------------
# syntax to extract the Gamma_D - that is degrees of freedom to be used for
# assigning the Dirichlet boundary conditions
#
# recall that the no slip boundary conditions are all associated with the
# velocity function space (W.sub(0)); while the fixed pressure condition will be
# associated with the pressure function space (W.sub(1))
#-------------------------------------------------------------------------------
b_dofs_L = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, ID_L)
b_dofs_T = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, ID_T)
b_dofs_B = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, ID_B)
b_dofs_R = dfx.fem.locate_dofs_topological((W.sub(0), U_sub), fdim, ID_R)
b_dofs_P = dfx.fem.locate_dofs_topological((W.sub(1), P_sub), fdim-1, ID_C)

#---------------------------------------------------------------------------
# syntax for defining the 3 different uD values to be used
# recall: u = uD for all x in Gamma_D (this is our template for implementing
# all the Dirichlet conditions)
#---------------------------------------------------------------------------
uD_Wall = dfx.fem.Function(U_sub)
uD_Wall.interpolate(noSlipBC)

uD_Lid = dfx.fem.Function(U_sub)
uD_Lid.interpolate(lidBC)

uD_P = dfx.fem.Function(P_sub)
uD_P.interpolate(pressureBC)

#--------------------------------------------------------------------------
# we now assign the actual Dirichlet boundary conditions by identifying the
# correct combinations of uD and Gamma_D
#--------------------------------------------------------------------------
bc_L = dfx.fem.dirichletbc(uD_Wall, b_dofs_L, W.sub(0))
bc_T = dfx.fem.dirichletbc(uD_Lid, b_dofs_T, W.sub(0))
bc_B = dfx.fem.dirichletbc(uD_Wall, b_dofs_B, W.sub(0))
bc_R = dfx.fem.dirichletbc(uD_Wall, b_dofs_R, W.sub(0))
bc_P = dfx.fem.dirichletbc(uD_P, b_dofs_P, W.sub(1))

bc = [bc_L, bc_T, bc_B, bc_R, bc_P]

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
