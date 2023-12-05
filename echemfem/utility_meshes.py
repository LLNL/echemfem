import numpy as np
from pyop2.mpi import COMM_WORLD
from firedrake.cython import dmcommon
from firedrake import mesh

__all__ = ['IntervalBoundaryLayerMesh', 'RectangleBoundaryLayerMesh']


def IntervalBoundaryLayerMesh(ncells, length_or_left, ncells_bdlayer, length_bdlayer,
                              boundary=(2,), right=None, distribution_parameters=None, comm=COMM_WORLD):
    """
    Generate a boundary layer mesh of an interval.

    :arg ncells: The number of the cells over the interval, outside
         the boundary layer.
    :arg length_or_left: The length of the interval (if ``right``
         is not provided) or else the left hand boundary point.
    :arg ncells_bdlayer: The number of the cells in each boundary
         layer.
    :arg length_bdlayer: The length of the boundary layer(s).
    :arg boundary: (optional) boundary markers of the boundary layer(s).
    :arg right: (optional) position of the right
         boundary point (in which case ``length_or_left`` should
         be the left boundary point).
    :kwarg comm: Optional communicator to build the mesh on (defaults to
        COMM_WORLD).

    The left hand boundary point has boundary marker 1,
    while the right hand point has marker 2.
    """
    if right is None:
        left = 0
        right = length_or_left
    else:
        left = length_or_left

    if ncells <= 0 or ncells % 1:
        raise ValueError("Number of cells must be a postive integer")
    right0 = right
    left0 = left
    if 1 in boundary:
        left0 += length_bdlayer
    if 2 in boundary:
        right0 -= length_bdlayer
    length = right0 - left0
    if length < 0:
        raise ValueError("Requested mesh has negative length")
    dx = length / ncells
    dx1 = length_bdlayer / ncells_bdlayer
    # This ensures the rightmost point is actually present.
    coords = np.arange(left0, right0 + 0.01 * dx, dx, dtype=np.double)
    if 1 in boundary:
        coords = np.append(np.arange(left, length_bdlayer - 0.99 * dx1, dx1, dtype=np.double), coords)
    if 2 in boundary:
        coords = np.append(coords, np.arange(right-length_bdlayer + dx1, right + 0.01 * dx1, dx1, dtype=np.double))
    coords = coords.reshape(-1, 1)
    cells = np.dstack((np.arange(0, len(coords) - 1, dtype=np.int32),
                       np.arange(1, len(coords), dtype=np.int32))).reshape(-1, 2)
    plex = mesh.plex_from_cell_list(1, cells, coords, comm)
    # Apply boundary IDs
    plex.createLabel(dmcommon.FACE_SETS_LABEL)
    coordinates = plex.getCoordinates()
    coord_sec = plex.getCoordinateSection()
    vStart, vEnd = plex.getDepthStratum(0)  # vertices
    for v in range(vStart, vEnd):
        vcoord = plex.vecGetClosure(coord_sec, coordinates, v)
        if vcoord[0] == coords[0]:
            plex.setLabelValue(dmcommon.FACE_SETS_LABEL, v, 1)
        if vcoord[0] == coords[-1]:
            plex.setLabelValue(dmcommon.FACE_SETS_LABEL, v, 2)

    return mesh.Mesh(plex, reorder=False, distribution_parameters=distribution_parameters)


def RectangleBoundaryLayerMesh(nx, ny, Lx, Ly, n_bdlayer, L_bdlayer, Ly_bdlayer=None, ny_bdlayer=None,
                               boundary=(1,), reorder=None, distribution_parameters=None, comm=COMM_WORLD):
    """Generate a boundary layer rectangular mesh using quadrilateral elements

    :arg nx: The number of cells in the x direction outside the boundary layers
    :arg ny: The number of cells in the y direction outside the boundary layers
    :arg Lx: The extent in the x direction
    :arg Ly: The extent in the y direction
    :arg n_bdlayer: The number of the cells in each boundary
         layer. If ny_bdlayer is defined, length of the boundary layer
         in the x-direction.
    :arg L_bdlayer: The length of the boundary layer(s). If Ly_bdlayer
                    is defined, length of in the x-direction
    :arg Ly_bdlayer: (optional) Length of the y-direction boundary layer
    :arg ny_bdlayer: (optional) number of cells of the y-direction boundary
         layer
    :arg boundary: (optional) boundary markers of the boundary layer(s).
    :kwarg reorder: (optional), should the mesh be reordered
    :kwarg comm: Optional communicator to build the mesh on (defaults to
        COMM_WORLD).

    The boundary edges in this mesh are numbered as follows:

    * 1: plane x == 0
    * 2: plane x == Lx
    * 3: plane y == 0
    * 4: plane y == Ly
    """

    for n in (nx, ny):
        if n <= 0 or n % 1:
            raise ValueError("Number of cells must be a postive integer")

    x0 = 0.0
    y0 = 0.0
    Lx0 = Lx
    Ly0 = Ly
    if Ly_bdlayer is None:
        Lx_bdlayer = L_bdlayer
        Ly_bdlayer = L_bdlayer
    else:
        Lx_bdlayer = L_bdlayer

    if ny_bdlayer is None:
        nx_bdlayer = n_bdlayer
        ny_bdlayer = n_bdlayer
    else:
        nx_bdlayer = n_bdlayer

    if 1 in boundary:
        x0 += Lx_bdlayer
    if 2 in boundary:
        Lx0 -= Lx_bdlayer
    if 3 in boundary:
        y0 += Ly_bdlayer
    if 4 in boundary:
        Ly0 -= Ly_bdlayer
    xcoords = np.linspace(x0, Lx0, nx + 1, dtype=np.double)
    ycoords = np.linspace(y0, Ly0, ny + 1, dtype=np.double)
    if 1 in boundary:
        xcoords = np.append(np.linspace(0.0, Lx_bdlayer - Lx_bdlayer/nx_bdlayer, nx_bdlayer, dtype=np.double), xcoords)
    if 2 in boundary:
        xcoords = np.append(xcoords, np.linspace(Lx0 + Lx_bdlayer/nx_bdlayer, Lx, nx_bdlayer, dtype=np.double))
    if 3 in boundary:
        ycoords = np.append(np.linspace(0.0, Ly_bdlayer - Ly_bdlayer/ny_bdlayer, ny_bdlayer, dtype=np.double), ycoords)
    if 4 in boundary:
        ycoords = np.append(ycoords, np.linspace(Ly0 + Ly_bdlayer/ny_bdlayer, Ly, ny_bdlayer, dtype=np.double))
    coords = np.asarray(np.meshgrid(xcoords, ycoords)).swapaxes(0, 2).reshape(-1, 2)
    # cell vertices
    if 1 in boundary:
        nx += nx_bdlayer
    if 2 in boundary:
        nx += nx_bdlayer
    if 3 in boundary:
        ny += ny_bdlayer
    if 4 in boundary:
        ny += ny_bdlayer
    i, j = np.meshgrid(np.arange(nx, dtype=np.int32), np.arange(ny, dtype=np.int32))
    cells = [i*(ny+1) + j, i*(ny+1) + j+1, (i+1)*(ny+1) + j+1, (i+1)*(ny+1) + j]
    cells = np.asarray(cells).swapaxes(0, 2).reshape(-1, 4)
    plex = mesh.plex_from_cell_list(2, cells, coords, comm)

    # mark boundary facets
    plex.createLabel(dmcommon.FACE_SETS_LABEL)
    plex.markBoundaryFaces("boundary_faces")
    coords = plex.getCoordinates()
    coord_sec = plex.getCoordinateSection()
    if plex.getStratumSize("boundary_faces", 1) > 0:
        boundary_faces = plex.getStratumIS("boundary_faces", 1).getIndices()
        xtol = Lx/(2*nx)
        xtol_bdlayer = Lx_bdlayer/(2*nx_bdlayer)
        ytol = Ly/(2*ny)
        ytol_bdlayer = Ly_bdlayer/(2*ny_bdlayer)
        if 1 in boundary:
            xtol1 = xtol_bdlayer
        else:
            xtol1 = xtol
        if 2 in boundary:
            xtol2 = xtol_bdlayer
        else:
            xtol2 = xtol
        if 3 in boundary:
            ytol3 = ytol_bdlayer
        else:
            ytol3 = ytol
        if 4 in boundary:
            ytol4 = ytol_bdlayer
        else:
            ytol4 = ytol
        for face in boundary_faces:
            face_coords = plex.vecGetClosure(coord_sec, coords, face)
            if abs(face_coords[0]) < xtol1 and abs(face_coords[2]) < xtol1:
                plex.setLabelValue(dmcommon.FACE_SETS_LABEL, face, 1)
            if abs(face_coords[0] - Lx) < xtol2 and abs(face_coords[2] - Lx) < xtol2:
                plex.setLabelValue(dmcommon.FACE_SETS_LABEL, face, 2)
            if abs(face_coords[1]) < ytol3 and abs(face_coords[3]) < ytol3:
                plex.setLabelValue(dmcommon.FACE_SETS_LABEL, face, 3)
            if abs(face_coords[1] - Ly) < ytol4 and abs(face_coords[3] - Ly) < ytol4:
                plex.setLabelValue(dmcommon.FACE_SETS_LABEL, face, 4)

    return mesh.Mesh(plex, reorder=reorder, distribution_parameters=distribution_parameters, comm=comm)
