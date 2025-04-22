"""A hands-on demo code for in-class active learning on implementation of
the finite element solver for a advecton-diffusion mass transport problem

Problem Description:
--------------------
We solve for species concentration in a fluid flow comprising two rotating
vortices inside a box. The box domain is: [0,2] X [0,1]. The flow velocity
field is given as:

    u = -pi*sin(pi*x)*cos(pi*y); v = pi*cos(pi*x)*sin(pi*y)

We assume Dirichlet zero-concentration boundary condition for all walls except
the bottom wall, where a Neumann flux q is specified for the diffusive flux
component.

Disclaimer:
-----------
Developed for computational fluid dynamics class taught at the Paul M Rady
Department of Mechanical Engineering at the University of Colorado Boulder by
Prof. Debanjan Mukherjee.

All inquiries addressed to Prof. Mukherjee directly at debanjan@Colorado.Edu
"""
import basix
from mpi4py import MPI
import dolfinx as dfx
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import gmshio
import ufl as ufl
import matplotlib.pyplot as plt
import numpy as np

#---------------------------------------------------------------
# start the problem by defining a function that programs in the
# analytical expression for the background flow velocity field
#---------------------------------------------------------------
def analyticalVelocity(x):
    vals = np.zeros((mesh.geometry.dim, x.shape[1]))
    p = np.pi
    vals[0] = -p * np.sin(p*x[0]) * np.cos(p*x[1])
    vals[1] = p * np.cos(p*x[0]) * np.sin(p*x[1])
    return vals

#--------------------------------------------------------------------
# problem parameter settings:
#----------------------------
# 1. Enter polynomial order for FEM (pOrder)
# 2. Enter value for constant diffusivity (Dvalue)
# 3. Enter value for constant reaction source/sink (Rvalue)
# 4. Enter any relevant boundary condition values
# 5. Enter the file name to read the mesh (with boundary labels) from
# 6. Enter the file name where the solution will be output
# 7. Enter the file name for storing the background flow
# 8. Enter a boolean choice variable to turn stabilization on/off
# 9. Enter the time-step size
# 10. Enter simulation start-time
# 11. Enter simulation stop time
# 12. Enter the value for theta in the theta-Galerkin method
#---------------------------------------------------------------------
pOrder = 1
Dvalue = 0.01
Rvalue = 0.0
qvalue = 1.0
mshFileName = 'box.msh'
solFileName = 'dg-sol.xdmf'
velFileName = 'dg-vel.pvd'
isStabilize = False
dt = 0.005
t = 0.0
t_end = 2.0
theta = 0.5

#--------------------------------------------------------------------------
# then we recall the IDs that we have used to mark the different boundary
# facets - these must correspond to the values that have been programmed in
# the GMSH geo file
#--------------------------------------------------------------------------
ID_X0 = 5 # ID for the x = 0 boundary
ID_X1 = 7 # ID for the x = 2 boundary
ID_Y0 = 8 # ID for the y = 0 boundary
ID_Y1 = 6 # ID for the y = 1 boundary

#-------------------------------------------------------------------------------
# read the mesh from external mesh file using the built in Gmsh IO functions
# this reads all the mesh data structures into mesh
# then reads all the physical domain markers into cell_markers
# (note: this corresponds to the Physical Surface tag in the Gmsh geo file)
# then reads all the physical boundary markers into facet_markers
# (note: this corresponds to the Physical Curve tags in the Gmsh geo file)
#-------------------------------------------------------------------------------
mesh, cell_markers, facet_markers = gmshio.read_from_msh(mshFileName, MPI.COMM_WORLD, gdim=2)

#--------------------------------------------------------------------
# create a function space on the mesh for the concentration trial and
# test function selection; we need H^1 minimum, so select pOrder >= 1
#--------------------------------------------------------------------
V = dfx.fem.functionspace(mesh, ("Lagrange", pOrder))

tdim = mesh.topology.dim
fdim = tdim - 1

#----------------------------------------------------------
# Defining the essential/Dirichlet boundary conditions
# Step 1: Identify the boundary faces (that is, Gamma_D)
# Rather here we are identifying the degrees of freedom
# located along the Gamma_D
#----------------------------------------------------------
b_dofs_X0 = dfx.fem.locate_dofs_topological(V, fdim, facet_markers.find(ID_X0))
b_dofs_X1 = dfx.fem.locate_dofs_topological(V, fdim, facet_markers.find(ID_X1))
b_dofs_Y0 = dfx.fem.locate_dofs_topological(V, fdim, facet_markers.find(ID_Y0))
b_dofs_Y1 = dfx.fem.locate_dofs_topological(V, fdim, facet_markers.find(ID_Y1))

#----------------------------------------------------------
# Defining the essential/Dirichlet boundary conditions
# Step 2: Specify the boundary value (that is, u_D)
#----------------------------------------------------------
uD_X0 = dfx.fem.Constant(mesh, dfx.default_scalar_type(0))
uD_X1 = dfx.fem.Constant(mesh, dfx.default_scalar_type(0))
uD_Y1 = dfx.fem.Constant(mesh, dfx.default_scalar_type(0))

#------------------------------------------------------------
# Defining the essential/Dirichlet boundary conditions
# Step 3: Apply the Dirichlet Condition (u = u_D at Gamma_D)
#------------------------------------------------------------
bc_X0 = dfx.fem.dirichletbc(uD_X0, b_dofs_X0, V)
bc_X1 = dfx.fem.dirichletbc(uD_X1, b_dofs_X1, V)
bc_Y1 = dfx.fem.dirichletbc(uD_Y1, b_dofs_Y1, V)

#--------------------------------------------------------------------------
# we map the background velocity as a custom expression onto the mesh
# the steps to be followed in general are as follows:
# (a) define a function (Python) that programs the expressions for the
# velocity on the mesh (see analyticalVelocity defined above)
# (b) then create a function space and a function on this space into which
# we can interpolate the expression
# (c) then we can interpolate the custom expression into a FEniCS function
#--------------------------------------------------------------------------
QE = basix.ufl.element('CG', mesh.basix_cell(), pOrder, shape=(mesh.geometry.dim,))
W = dfx.fem.functionspace(mesh, QE)
u0Func = dfx.fem.Function(W)
u0Func.interpolate(lambda x: analyticalVelocity(x))
u0Func.name = 'vel'

#---------------------------------------------------------------------------
# we will write this interpolated background velocity field onto an external
# VTK file to perform flow visualization using an external tool
#---------------------------------------------------------------------------
velFile = dfx.io.VTKFile(MPI.COMM_WORLD, velFileName, "w")
velFile.write_function([u0Func])
velFile.close()

#-------------------------------------------------------
# define the trialfunction and the testfunction objects
#-------------------------------------------------------
c = ufl.TrialFunction(V)
w = ufl.TestFunction(V)

#-------------------------------------------------------------------
# define the solutions c_n and c_n+1 - placeholders for solutions at
# time t_n and t_n+1 respectively
#-------------------------------------------------------------------
c0  = dfx.fem.Function(V)
c1  = dfx.fem.Function(V)

#---------------------------------------------------------------------
# define the diffusion contribution to the weak form and matrix system
# evaluated now at t_n and t_n+1 respectively
#---------------------------------------------------------------------
D = dfx.fem.Constant(mesh, dfx.default_scalar_type(Dvalue))
K0 = D * ufl.inner(ufl.grad(w), ufl.grad(c0)) * ufl.dx
K1 = D * ufl.inner(ufl.grad(w), ufl.grad(c)) * ufl.dx

#---------------------------------------------------------------------
# define the advection contribution to the weak form and matrix system
# evaluated now at t_n and t_n+1 respectively
#---------------------------------------------------------------------
U0 = ufl.inner(w, ufl.dot(u0Func, ufl.grad(c0))) * ufl.dx
U1 = ufl.inner(w, ufl.dot(u0Func, ufl.grad(c))) * ufl.dx

#---------------------------------------------------------------------
# define the reaction contribution to the weak form and matrix system
#---------------------------------------------------------------------
R = dfx.fem.Constant(mesh, dfx.default_scalar_type(Rvalue))
fR = w * R * ufl.dx

#-------------------------------------------------------------------------------
# extract the Neumann boundary portion for the domain
# by extracting out all the boundary facets into the object 'ds' based on the
# facet_marker information, we can now restrict the integral onto only the
# portion of the boundary with id = ID_Y0 by assembling the integral over the
# portion ds(ID_Y0)
#-------------------------------------------------------------------------------
ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_markers)

#-----------------------------------------------------------------------
# define the Neumann BC contribution to the weak form and matrix system
#-----------------------------------------------------------------------
q = dfx.fem.Constant(mesh, dfx.default_scalar_type(qvalue))
fN = w * q * ds(ID_Y0)

#------------------------------------------------------------
# Theta-Galerkin method weak form between t_n, and t_n+1
#------------------------------------------------------------
weakForm    = (1.0/dt) * ufl.inner(w, c) * ufl.dx - (1.0/dt) * ufl.inner(w,c0)*ufl.dx \
            + theta * K1 + theta * U1 \
            + (1.0 - theta) * K0 + (1.0 - theta) * U0 \
            - fN - fR

#-------------------------------------------------------------------------------
# add the stabilization terms for this problem
# uNorm below implements the computation of a norm or magnitude of the flow
# velocity vector; h below represents a measure of the mesh element sizing
# computed locally using a CellDiameter object; and tau below implements
# the definition of the stabilization parameter
#-------------------------------------------------------------------------------
if isStabilize == True:

    uNorm = ufl.sqrt(ufl.inner(u0Func, u0Func))
    h = ufl.CellDiameter(mesh)
    tau = ( (2.0*uNorm/h)**2 + 9.0*(4.0*D/(h*h))**2 )**(-0.5)
    residual = (1.0/dt)*(c1 - c0) + ufl.dot(u0, ufl.grad(c)) - ufl.div(D*ufl.grad(c))
    S = tau * ufl.inner(ufl.dot(u0Func, ufl.grad(w)), residual) * ufl.dx
    weakForm = weakForm + S

#-------------------------------------------------------
# combine all the Dirichlet boundary conditions together
#-------------------------------------------------------
bcSet = [bc_X0, bc_X1, bc_Y1]

#-----------------------------------------------------------------
# assembling the global matrix vector system based on the combined
# theta-Galerkin weak form
#-----------------------------------------------------------------
mat = ufl.lhs(weakForm)
vec = ufl.rhs(weakForm)

problem = LinearProblem(mat, vec, bcs=bcSet, u=c1, \
    petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "umfpack"})

#----------------------------------------------------------------------------
# loop through time to solve for the updated version of the solution variable
# it is important to note the structure of the loop as it is implemented:
# - we identify which time step we are at
# - we solve for the updated solution c_n+1
# - we set c_n+1 as the c_n for the next time step
# - we update the time step
#
# we will also calculate a volume averaged concentration based on the
# computed concentration values from the simulations - and then plot it using
# matplotlib
#
# lastly: we are going to use an XDMF file format - which is more efficient
# in handling these dynamic simulation datasets (and can handle read/write of
# mesh data in parallel)
#----------------------------------------------------------------------------

xdmfOut = dfx.io.XDMFFile(mesh.comm, "unsteady-adr.xdmf", "w", encoding=dfx.io.XDMFFile.Encoding.ASCII)
xdmfOut.write_mesh(mesh)

timeStepList = []
concentrationAvg = []

while t < t_end:

    print('t=', t)
    c1 = problem.solve()
    c0.x.array[:] = c1.x.array

    c1.name = 'conc'
    xdmfOut.write_function(c1, t)

    c_form = dfx.fem.form(c1 * ufl.dx)
    c_total = dfx.fem.assemble_scalar(c_form)
    timeStepList.append(t)
    concentrationAvg.append(c_total)

    t += dt

xdmfOut.close()

plt.plot(timeStepList, concentrationAvg, 'r-')
plt.xlabel('time', fontweight='bold')
plt.ylabel('concentration', fontweight='bold')
plt.title('Computed Average Concentration', fontweight='bold')
plt.savefig('conc-avg.png', dpi=120, bbox_inches='tight')
plt.close()