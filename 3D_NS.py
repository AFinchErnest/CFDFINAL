import fenics as fe
import numpy as np

# writes velocity field to an xdmf so that doesn't have to solve for it each time in the festim simulation

fluid_id = 6
inlet_id = 7
outlet_id = 8
vacuum_id = 9
walls_id = 10


def fluid_dynamics_sim_chorin(
    mesh, volume_markers, surface_markers, id_inlet, id_outlet, id_walls
):
    """
    Solves the Navier-Stokes equations using the Chorin's projection method
    See https://fenicsproject.org/pub/tutorial/html/._ftut1009.html#ftut1:NS for more details

    Args:
        mesh (fenics.Mesh): the mesh
        volume_markers (fenics.MeshFunction): the volume markers (subdomains)
        surface_markers (fenics.MeshFunction): the surface markers (boundaries)
        id_inlet (int): the id of the inlet boundary
        id_outlet (int): the id of the outlet boundary
        id_walls (int or list): the id(s) of the walls

    Returns:
        fenics.Function: the velocity field
        fenics.Function: the pressure field
    """
    V = fe.VectorFunctionSpace(mesh, "CG", 2)
    Q = fe.FunctionSpace(mesh, "CG", 1)

    u = fe.TrialFunction(V)
    p = fe.TrialFunction(Q)
    v = fe.TestFunction(V)
    q = fe.TestFunction(Q)

    u_n = fe.Function(V)
    u_ = fe.Function(V)
    p_n = fe.Function(Q)
    p_ = fe.Function(Q)

    # ##### Boundary conditions ##### #

    inlet_velocity = 2e-03  # units: m s-1
    outlet_pressure = 0  # units: Pa

    inflow = fe.DirichletBC(
        V, fe.Constant((inlet_velocity, 0.0, 0.0)), surface_markers, id_inlet
    )

    # make sure id_walls is a list
    if isinstance(id_walls, int):
        id_walls = [id_walls]

    # iterate through the walls
    walls = []
    for id_wall in id_walls:
        walls.append(
            fe.DirichletBC(V, fe.Constant((0.0, 0.0, 0.0)), surface_markers, id_wall)
        )

    pressure_outlet = fe.DirichletBC(
        Q, fe.Constant(outlet_pressure), surface_markers, id_outlet
    )
    bcu = [inflow] + walls
    bcp = [pressure_outlet]

    # ##### Solver ##### #
    dx = fe.Measure("dx", subdomain_data=volume_markers)
    ds = fe.Measure("ds", subdomain_data=surface_markers)

    dt = 0.1  # Time step size

    k = fe.Constant(dt)
    n = fe.FacetNormal(mesh)
    U = 0.5 * (u_n + u)

    # LiPb
    mu = 1
    rho = 1

    def epsilon(u):
        return fe.sym(fe.nabla_grad(u))

    def sigma(u, p):
        return 2 * mu * epsilon(u) - p * fe.Identity(len(u))

    # ##### Solver ##### #

    # Tentative velocity step
    F1 = rho * fe.dot((u - u_n) / k, v) * dx
    F1 += rho * fe.dot(fe.dot(u_n, fe.nabla_grad(u_n)), v) * dx
    F1 += fe.inner(sigma(U, p_n), epsilon(v)) * dx
    F1 += fe.dot(p_n * n, v) * ds
    F1 -= fe.dot(mu * fe.nabla_grad(U) * n, v) * ds

    a1 = fe.lhs(F1)
    L1 = fe.rhs(F1)

    solver1 = fe.KrylovSolver("bicgstab", "hypre_amg")
    solver1.parameters["absolute_tolerance"] = 1e-08
    solver1.parameters["relative_tolerance"] = 1e-08
    solver1.parameters["maximum_iterations"] = 1000
    solver1.parameters["report"] = False
    solver1.parameters["monitor_convergence"] = False

    # Pressure update
    a2 = fe.dot(fe.nabla_grad(p), fe.nabla_grad(q)) * dx
    L2 = fe.dot(fe.nabla_grad(p_n), fe.nabla_grad(q)) * dx
    L2 -= (1 / k) * fe.div(u_) * q * dx

    solver2 = fe.KrylovSolver("bicgstab", "hypre_amg")
    solver2.parameters["absolute_tolerance"] = 1e-08
    solver2.parameters["relative_tolerance"] = 1e-08
    solver2.parameters["maximum_iterations"] = 1000
    solver2.parameters["report"] = False
    solver2.parameters["monitor_convergence"] = False

    # Velocity update
    a3 = fe.dot(u, v) * dx
    L3 = fe.dot(u_, v) * dx
    L3 -= k * fe.dot(fe.nabla_grad(p_ - p_n), v) * dx

    solver3 = fe.KrylovSolver("cg", "sor")
    solver3.parameters["absolute_tolerance"] = 1e-08
    solver3.parameters["relative_tolerance"] = 1e-08
    solver3.parameters["maximum_iterations"] = 1000
    solver3.parameters["report"] = False
    solver3.parameters["monitor_convergence"] = False

    # Assemble matrices
    A1 = fe.assemble(a1)
    A2 = fe.assemble(a2)
    A3 = fe.assemble(a3)

    [bc.apply(A1) for bc in bcu]
    [bc.apply(A2) for bc in bcp]

    # Compute tentative velocity step
    b1 = fe.assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solver1.solve(A1, u_.vector(), b1)  # <-- culprit!

    # # Pressure correction
    # b2 = fe.assemble(L2)
    # [bc.apply(A2, b2) for bc in bcp]
    # solver2.solve(A2, p_.vector(), b2)

    # # Velocity correction
    # b3 = fe.assemble(L3)
    # [bc.apply(A3, b3) for bc in bcu]
    # solver3.solve(A3, u_.vector(), b3)

    # # Move to next time step
    # u_n.assign(u_)
    # p_n.assign(p_)

    return u_, p_


if __name__ == "__main__":

    volume_file = "mesh_cells.xdmf"
    boundary_file = "mesh_facets.xdmf"
    mesh = fe.Mesh()
    fe.XDMFFile(volume_file).read(mesh)

    # Read tags for volume elements
    volume_markers = fe.MeshFunction("size_t", mesh, mesh.topology().dim())
    fe.XDMFFile(volume_file).read(volume_markers)

    # Read tags for surface elements
    # (can also be used for applying DirichletBC)
    surface_markers = fe.MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
    fe.XDMFFile(boundary_file).read(surface_markers, "f")
    surface_markers = fe.MeshFunction("size_t", mesh, surface_markers)

    print("Succesfully load mesh with " + str(len(volume_markers)) + " cells")

    fluid_dynamics_sim_chorin(
        mesh=mesh,
        volume_markers=volume_markers,
        surface_markers=surface_markers,
        id_inlet=inlet_id,
        id_outlet=outlet_id,
        id_walls=[walls_id, vacuum_id],
        
    )