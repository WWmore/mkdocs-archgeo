'''meshpy.py: The half-edge mesh data structure'''

__author__ = 'Davide Pellis'


import copy

from io import open

import numpy as np

from scipy.sparse import coo_matrix

from scipy import spatial

# -----------------------------------------------------------------------------

from geometrylab import utilities

from geometrylab.geometry.polyline import Polyline

from geometrylab.geometry.frame import Frame

from geometrylab.geometry.circle import circle_three_points

# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#                               Half-Edge Mesh
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    #                          Build Data Structure
    # -------------------------------------------------------------------------

def make_mesh(self, vertices_list, faces_list):
    """
    Constructs a half-edge mesh data structure from a list of vertices and faces.

    This function initializes the half-edge mesh by creating the necessary data structures,
    including half-edges, vertices, and faces. It ensures that the mesh is manifold and orientable,
    and it reorients the mesh if necessary.

    The main logic involves:

    1. Creating half-edges from the input vertices and faces.

    2. Establishing relationships between half-edges, vertices, and faces.

    3. Handling boundary edges and faces.

    4. Reorienting the mesh if it is not initially orientable.

    Parameters
    ----------
    vertices_list : list of lists
        A list of vertex coordinates, where each vertex is represented as [x, y, z].
    faces_list : list of lists
        A list of faces, where each face is represented as a list of vertex indices.

    Returns
    -------
    None
        The mesh is constructed within the class instance.

    Raises
    ------
    ValueError
        If the mesh is not manifold and orientable.

    Note
    ----
    This function assumes that the input vertices and faces form a valid mesh.
    The mesh must be manifold and orientable for the half-edge structure to be correctly defined.

    See Also
    --------
    read_obj_file : Reads mesh data from an OBJ file.

    """   
    def _make_halfedges(self, vertices_list, faces_list):
        self._V = len(vertices_list)
        self._F = len(faces_list)
        self._vertices = np.array(vertices_list, 'f')
        orig = []
        face = []
        nexx = []
        prev = []
        twin_i = []
        twin_j = []
        h = 0
        for f in range(self.F):
            N = len(faces_list[f])
            orig.append(faces_list[f][0])
            face.append(f)
            nexx.append(h + 1)
            prev.append(h + N - 1)
            twin_i.append(faces_list[f][1])
            twin_j.append(faces_list[f][0])
            for v in range(1, N-1):
                orig.append(faces_list[f][v])
                face.append(f)
                nexx.append(h + v + 1)
                prev.append(h + v - 1)
                twin_i.append(faces_list[f][v+1])
                twin_j.append(faces_list[f][v])
            orig.append(faces_list[f][N-1])
            face.append(f)
            nexx.append(h)
            prev.append(h + N - 2)
            twin_i.append(faces_list[f][0])
            twin_j.append(faces_list[f][N-1])
            h += N
        H = np.zeros((h, 6), 'i')
        H[:,0] = orig
        H[:,1] = face
        H[:,2] = nexx
        H[:,3] = prev
        twin = coo_matrix((np.arange(h) + 1, (twin_i, twin_j)), shape=(h, h))
        twin = twin.tocsc()
        H[:,4] = twin[H[:,0],H[H[:,2],0]] - 1
        b = np.where(H[:,4] == -1)[0]
        boundary = H[b,:]
        boundary[:,0] = H[H[b,2],0]
        boundary[:,1] = -1
        boundary[:,4] = b
        #B = boundary.shape[0]  # test by Davide
        B = len(boundary)
        if B > 0:
            Bh = np.arange(h, h+B)
            H[b,4] = Bh
            zeros = np.zeros(B)
            p = coo_matrix((Bh, (H[b,0], zeros)), shape=(self.V, 1))
            p = p.tocsc()
            # print(boundary[:,0]) # test by Davide
            # print(zeros) # test by Davide
            # print(p) # test by Davide
            boundary[:,3] = p[boundary[:,0],zeros]
            i = boundary[boundary[:,3]-h,0]
            p = coo_matrix((Bh, (i, zeros)), shape=(self.V, 1))
            p = p.tocsc()
            boundary[:,2] = p[boundary[:,0], zeros]
            H = np.vstack((H, boundary))
        K = H[:,(3,4)]
        K[:,0] = np.arange(H.shape[0])
        m = np.amin(K, axis=1)
        u = np.unique(m)
        imap = np.arange(np.max(u) + 1, dtype=int)
        imap[u] = np.arange(u.shape[0])
        H[:,5] = imap[m]
        self._halfedges = H
        self._E = int(H.shape[0] / 2)
        self.topology_update()

    try:
        self._make_halfedges(vertices_list, faces_list)
    except:
        try:
            faces_list = self.orient_faces(vertices_list, faces_list)
            self._make_halfedges(vertices_list, faces_list)
            print('*** Mesh reoriented ***')
        except:
            raise ValueError('The mesh is not manifold and orientable!')
    print(self)

# -------------------------------------------------------------------------
#                                Navigating
# -------------------------------------------------------------------------

def origin(self, halfedge_index=None):
    """
    Retrieves the origin vertex of a specified half-edge.

    This function returns the index of the vertex where the specified half-edge originates.

    Parameters
    ----------
    halfedge_index : int, optional (default=None)
        The index of the half-edge. If None, returns the origin vertices for all half-edges.

    Returns
    -------
    origin_vertex : int or numpy array
        The index of the origin vertex for the specified half-edge, or an array of origin vertices for all half-edges.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The origin vertex is the starting vertex of the half-edge.

    See Also
    --------
    next : Retrieves the next half-edge in the face loop.

    previous : Retrieves the previous half-edge in the face loop.

    twin : Retrieves the twin half-edge.
    """
    H = self.halfedges
    if halfedge_index is None:
        return H[:,0]
    return H[halfedge_index,0]

def face(self, halfedge_index=None):
    """
    Retrieves the face associated with a specified half-edge.

    This function returns the index of the face that the specified half-edge belongs to.

    Parameters
    ----------
    halfedge_index : int, optional (default=None)
        The index of the half-edge. If None, returns the face indices for all half-edges.

    Returns
    -------
    face_index : int or numpy array
        The index of the face associated with the specified half-edge, or an array of face indices for all half-edges.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The face index is -1 if the half-edge is on the boundary.

    See Also
    --------
    origin : Retrieves the origin vertex of a half-edge.

    next : Retrieves the next half-edge in the face loop.

    previous : Retrieves the previous half-edge in the face loop.
    """
    H = self.halfedges
    if halfedge_index is None:
        return H[:,1]
    return H[halfedge_index,1]

def next(self, halfedge_index=None):
    """
    Retrieves the next half-edge in the face loop.

    This function returns the index of the next half-edge in the face loop starting from the specified half-edge.

    Parameters
    ----------
    halfedge_index : int, optional (default=None)
        The index of the half-edge. If None, returns the next half-edges for all half-edges.

    Returns
    -------
    next_halfedge : int or numpy array
        The index of the next half-edge in the face loop, or an array of next half-edges for all half-edges.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The next half-edge is part of the same face as the specified half-edge.

    See Also
    --------
    previous : Retrieves the previous half-edge in the face loop.

    origin : Retrieves the origin vertex of a half-edge.

    face : Retrieves the face associated with a half-edge.
    """
    H = self.halfedges
    if halfedge_index is None:
        return H[:,2]
    return H[halfedge_index,2]

def previous(self, halfedge_index=None):
    """
    Retrieves the previous half-edge in the face loop.

    This function returns the index of the previous half-edge in the face loop starting from the specified half-edge.

    Parameters
    ----------
    halfedge_index : int, optional (default=None)
        The index of the half-edge. If None, returns the previous half-edges for all half-edges.

    Returns
    -------
    previous_halfedge : int or numpy array
        The index of the previous half-edge in the face loop, or an array of previous half-edges for all half-edges.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The previous half-edge is part of the same face as the specified half-edge.

    See Also
    --------
    next : Retrieves the next half-edge in the face loop.

    origin : Retrieves the origin vertex of a half-edge.

    face : Retrieves the face associated with a half-edge.
    """
    H = self.halfedges
    if halfedge_index is None:
        return H[:,3]
    return H[halfedge_index,3]

def twin(self, halfedge_index=None):
    """
    Retrieves the twin half-edge of a specified half-edge.

    This function returns the index of the twin half-edge, which is the counterpart 
    of the specified half-edge across the shared edge.

    Parameters
    ----------
    halfedge_index : int, optional (default=None)
        The index of the half-edge. If None, returns the twin half-edges for all half-edges.

    Returns
    -------
    twin_halfedge : int or numpy array
        The index of the twin half-edge, or an array of twin half-edges for all half-edges.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The twin half-edge is part of the same edge as the specified half-edge but belongs to a different face.

    See Also
    --------
    origin : Retrieves the origin vertex of a half-edge.

    face : Retrieves the face associated with a half-edge.
    """
    H = self.halfedges
    if halfedge_index is None:
        return H[:,4]
    return H[halfedge_index,4]

def edge(self, halfedge_index=None):
    """
    Retrieves the edge associated with a specified half-edge.

    This function returns the index of the edge that the specified half-edge belongs to.

    Parameters
    ----------
    halfedge_index : int, optional (default=None)
        The index of the half-edge. If None, returns the edge indices for all half-edges.

    Returns
    -------
    edge_index : int or numpy array
        The index of the edge associated with the specified half-edge, or an array of edge indices for all half-edges.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The edge index is shared by the specified half-edge and its twin.

    See Also
    --------
    origin : Retrieves the origin vertex of a half-edge.

    face : Retrieves the face associated with a half-edge.

    twin : Retrieves the twin half-edge.
    """
    H = self.halfedges
    if halfedge_index is None:
        return H[:,5]
    return H[halfedge_index,5]

# -------------------------------------------------------------------------
#                               Reading
# -------------------------------------------------------------------------

def read_obj_file(self, file_name):
    """
    Reads mesh data from an OBJ file and constructs the half-edge mesh.

    This function parses the OBJ file to extract vertices and faces, and then calls
    the `make_mesh` method to construct the half-edge mesh data structure.

    The main logic involves:

    1. Reading the OBJ file line by line to extract vertex and face information.
    
    2. Handling different formats of face definitions (e.g., with or without texture coordinates).
    
    3. Calling `make_mesh` to construct the mesh.

    Parameters
    ----------
    file_name : str
        The path to the OBJ file.

    Returns
    -------
    None
        The mesh is constructed within the class instance.

    Note
    ----
    This function assumes that the OBJ file is correctly formatted and contains valid mesh data.
    The mesh must be manifold and orientable for the half-edge structure to be correctly defined.

    See Also
    --------
    make_mesh : Constructs the half-edge mesh from vertices and faces.

    make_obj_file : Writes the mesh data to an OBJ file.
    """
    file_name = str(file_name)
    self.name = file_name.split('.')[0]
    obj_file = open(file_name, encoding='utf-8')
    vertices_list = []
    uv_list = []
    faces_list = []
    for l in obj_file:
        splited_line = l.split(' ')
        if splited_line[0] == 'v':
            split_x = splited_line[1].split('\n')
            x = float(split_x[0])
            split_y = splited_line[2].split('\n')
            y = float(split_y[0])
            split_z = splited_line[3].split('\n')
            try:
                z = float(split_z[0])
            except ValueError:
                print('WARNING: disable line wrap when saving .obj')
            vertices_list.append([x, y ,z])
        elif splited_line[0] == 'f':
            v_list = []
            L = len(splited_line)
            try:
                for i in range(1, L):
                    splited_face_data = splited_line[i].split('/')
                    v_list.append(int(splited_face_data[0]) - 1 )
                faces_list.append(v_list)
            except ValueError:
                v_list = []
                for i in range(1, L-1):
                    v_list.append(int(splited_line[i]) - 1 )
                faces_list.append(v_list)
        if splited_line[0] == 'vt':
            split_u = splited_line[1].split('\n')
            u = float(split_u[0])
            split_v = splited_line[2].split('\n')
            v = float(split_v[0])
            vertices_list.append([u,v])
        if len(uv_list) > 0:
            self._uv = np.array(uv_list)
    self.make_mesh(vertices_list, faces_list)

# -------------------------------------------------------------------------
#                               Writing
# -------------------------------------------------------------------------

def make_obj_file(self, file_name, overwrite=False):
    """
    Writes the mesh data to an OBJ file.

    This function exports the mesh vertices and faces to an OBJ file, which can be 
    used for visualization or further processing.

    Parameters
    ----------
    file_name : str
        The path to the OBJ file.
    overwrite : bool, optional (default=False)
        Whether to overwrite the file if it already exists.

    Returns
    -------
    file_path : str
        The path to the created OBJ file (without the extension).

    Note
    ----
    This function assumes that the mesh data is valid and manifold.
    The OBJ file format is widely supported for 3D visualization and modeling.

    See Also
    --------
    read_obj_file : Reads mesh data from an OBJ file.
    """
    path = utilities.make_filepath(file_name, 'obj', overwrite)
    obj = open(path, 'w')
    line = ('o {}\n').format(file_name)
    obj.write(line)
    faces = self.faces_list()
    for v in range(self.V):
        vi = self.vertices[v]
        line = ('v {} {} {}\n').format(vi[0], vi[1], vi[2])
        obj.write(line)
    for f in range(self.F):
            obj.write('f ')
            N = len(faces[f])
            for v in range(N - 1):
                vf = str(faces[f][v] + 1)
                #obj.write(unicode(vf + '//' + ' '))
                obj.write(vf + ' ')
            vf = str(faces[f][N - 1] + 1)
            #obj.write(unicode(vf + '//' + '\n'))
            obj.write(vf + '\n')
    obj.close()
    return path.split('.')[0]

# -------------------------------------------------------------------------
#                                 Ordering
# -------------------------------------------------------------------------

def vertex_ring_ordered_halfedges(self):
    """
    Retrieves the half-edges in the vertex ring around each vertex, ordered by traversal.

    This function returns the indices of half-edges that form the vertex ring around 
    each vertex, ordered in a consistent traversal direction.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    ordered_halfedges : numpy array
        The indices of half-edges in the vertex ring, ordered by traversal.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The half-edges are ordered to facilitate consistent traversal around each vertex.

    See Also
    --------
    vertex_ring_vertices_iterators : Iterates over vertices in the vertex ring.
    """
    H = np.copy(self.halfedges)
    i = np.argsort(H[:,0])
    v = H[i,0]
    index = np.arange(H.shape[0])
    _, j = np.unique(v, True)
    v = np.delete(v,j)
    index = np.delete(index,j)
    while v.shape[0] > 0:
        _, j = np.unique(v, True)
        i[index[j]] = H[H[i[index[j] - 1],3],4]
        v = np.delete(v,j)
        index = np.delete(index,j)
    return i

def face_ordered_halfedges(self):
    """
    Retrieves the half-edges in each face, ordered by traversal.

    This function returns the indices of half-edges that form each face, ordered 
    in a consistent traversal direction.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    ordered_halfedges : numpy array
        The indices of half-edges in each face, ordered by traversal.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The half-edges are ordered to facilitate consistent traversal around each face.

    See Also
    --------
    face_vertices_iterators : Iterates over vertices in each face.
    """
    H = np.copy(self.halfedges)
    i = np.argsort(H[:,1])
    i = i[np.where(H[i,1] >= 0)]
    f = H[i,1]
    index = np.arange(i.shape[0])
    _, j = np.unique(f, True)
    f = np.delete(f,j)
    index = np.delete(index, j)
    while f.shape[0] > 0:
        _, j = np.unique(f, True)
        i[index[j]] = H[i[index[j] - 1],2]
        f = np.delete(f, j)
        index = np.delete(index, j)
    return i

# -------------------------------------------------------------------------
#                              Iterators
# -------------------------------------------------------------------------

def vertex_ring_vertices_iterators(self, sort=False, order=False, return_lengths=False):
    """
    Iterates over the vertices in the vertex ring around each vertex.

    This function provides iterators for the vertices in the vertex ring around 
    each vertex, with options for sorting and ordering.

    Parameters
    ----------
    sort : bool, optional (default=False)
        Whether to sort the vertices by vertex index.
    order : bool, optional (default=False)
        Whether to order the vertices by traversal direction.
    return_lengths : bool, optional (default=False)
        Whether to return the lengths of the vertex rings.

    Returns
    -------
    v : numpy array
        The indices of the central vertices.
    vj : numpy array
        The indices of the adjacent vertices in the vertex ring.
    lengths : numpy array, optional
        The lengths of the vertex rings (returned if `return_lengths=True`).

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The vertex ring vertices are useful for local mesh analysis and processing.

    See Also
    --------
    vertex_ring_faces_iterators : Iterates over faces in the vertex ring.

    vertex_ring_edges_iterators : Iterates over edges in the vertex ring.
    """
    H = self.halfedges
    v  = H[:,0]
    vj = H[H[:,4],0]
    if order:
        i  = self.vertex_ring_ordered_halfedges()
        v  = v[i]
        vj = vj[i]
    elif sort:
        i  = np.argsort(v)
        v  = v[i]
        vj = vj[i]
    if return_lengths:
        i  = np.ones(vj.shape[0], dtype=int)
        lj = utilities.sum_repeated(i,v)
        return v, vj, lj
    else:
        return v, vj

def vertex_ring_faces_iterators(self, sort=False, order=False):
    """
    Iterates over the faces in the vertex ring around each vertex.

    This function provides iterators for the faces in the vertex ring around 
    each vertex, with options for sorting and ordering.

    Parameters
    ----------
    sort : bool, optional (default=False)
        Whether to sort the faces by face index.
    order : bool, optional (default=False)
        Whether to order the faces by traversal direction.

    Returns
    -------
    v : numpy array
        The indices of the central vertices.
    fj : numpy array
        The indices of the adjacent faces in the vertex ring.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The vertex ring faces are useful for local mesh analysis and processing.

    See Also
    --------
    vertex_ring_vertices_iterators : Iterates over vertices in the vertex ring.

    vertex_ring_edges_iterators : Iterates over edges in the vertex ring.
    """
    H = self.halfedges
    if order:
        i  = self.vertex_ring_ordered_halfedges()
        v  = H[i,0]
        fj = H[i,1]
    else:
        i  = np.where(H[:,1] >= 0)[0]
        v  = H[i,0]
        fj = H[i,1]
        if sort:
            i  = np.argsort(v)
            v  = v[i]
            fj = fj[i]
    return v, fj

def vertex_ring_edges_iterators(self, sort=False, order=False):
    """
    Iterates over the edges in the vertex ring around each vertex.

    This function provides iterators for the edges in the vertex ring around 
    each vertex, with options for sorting and ordering.

    Parameters
    ----------
    sort : bool, optional (default=False)
        Whether to sort the edges by edge index.
    order : bool, optional (default=False)
        Whether to order the edges by traversal direction.

    Returns
    -------
    v : numpy array
        The indices of the central vertices.
    ej : numpy array
        The indices of the adjacent edges in the vertex ring.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The vertex ring edges are useful for local mesh analysis and processing.

    See Also
    --------
    vertex_ring_vertices_iterators : Iterates over vertices in the vertex ring.

    vertex_ring_faces_iterators : Iterates over faces in the vertex ring.
    """
    H = self.halfedges
    v  = H[:,0]
    ej = H[:,5]
    if order:
        i  = self.vertex_ring_ordered_halfedges()
        v  = v[i]
        ej = ej[i]
    elif sort:
        i  = np.argsort(v)
        v  = v[i]
        ej = ej[i]
    return v, ej

def face_edge_vertices_iterators(self, sort=False, order=False):
    """
    Iterates over the vertices of edges in each face.

    This function provides iterators for the vertices of edges that form each face, 
    with options for sorting and ordering.

    Parameters
    ----------
    sort : bool, optional (default=False)
        Whether to sort the vertices by vertex index.
    order : bool, optional (default=False)
        Whether to order the vertices by traversal direction.

    Returns
    -------
    f : numpy array
        The indices of the faces.
    vi : numpy array
        The indices of the starting vertices of the edges.
    vj : numpy array
        The indices of the ending vertices of the edges.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The face edge vertices are useful for local mesh analysis and processing.

    See Also
    --------
    face_vertices_iterators : Iterates over vertices in each face.

    face_edges_iterators : Iterates over edges in each face.
    """
    H = self.halfedges
    f  = H[:,1]
    vi = H[:,0]
    vj = H[H[:,2],0]
    if order:
        i  = self.face_ordered_halfedges()
        f  = f[i]
        vi = vi[i]
        vj = vj[i]
    else:
        i  = np.where(H[:,1] >= 0)[0]
        f  = f[i]
        vi = vi[i]
        vj = vj[i]
        if sort:
            i  = np.argsort(f)
            vi = vi[i]
            vj = vj[i]
    return f, vi, vj

def face_vertices_iterators(self):
    """
    Iterates over the vertices in each face.

    This function provides iterators for the vertices that form each face.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    f : numpy array
        The indices of the faces.
    vi : numpy array
        The indices of the vertices in each face.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The face vertices are useful for local mesh analysis and processing.

    See Also
    --------
    face_edge_vertices_iterators : Iterates over vertices of edges in each face.

    face_edges_iterators : Iterates over edges in each face.
    """
    H = self.halfedges
    i  = self.face_ordered_halfedges()
    vi = H[i,0]
    f  = H[i,1]
    return f, vi

def face_edges_iterators(self):
    """
    Iterates over the edges in each face.

    This function provides iterators for the edges that form each face.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    f : numpy array
        The indices of the faces.
    ei : numpy array
        The indices of the edges in each face.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The face edges are useful for local mesh analysis and processing.

    See Also
    --------
    face_vertices_iterators : Iterates over vertices in each face.

    face_edge_vertices_iterators : Iterates over vertices of edges in each face.
    """
    H = self.halfedges
    i  = self.face_ordered_halfedges()
    ei = H[i,5]
    f  = H[i,1]
    return f, ei

def edge_vertices_iterators(self, sort=False):
    """
    Iterates over the vertices of each edge.

    This function provides iterators for the vertices that form each edge, 
    with an option for sorting.

    Parameters
    ----------
    sort : bool, optional (default=False)
        Whether to sort the vertices by vertex index.

    Returns
    -------
    e : numpy array
        The indices of the edges.
    vi : numpy array
        The indices of the vertices forming each edge.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The edge vertices are useful for local mesh analysis and processing.

    See Also
    --------
    edge_faces : Retrieves the faces adjacent to each edge.
    """
    H = self.halfedges
    e  = H[:,5]
    vi = H[:,0]
    if sort:
        i  = np.argsort(H[:,5])
        e  = e[i]
        vi = vi[i]
    return e, vi

def vertex_double_ring_vertices_iterators(self):
    """
    Iterates over the vertices in the double ring around each vertex.

    This function provides iterators for the vertices in the double ring around 
    each vertex, which includes the vertex ring and the vertices adjacent to it.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    v : numpy array
        The indices of the central vertices.
    vj : numpy array
        The indices of the vertices in the double ring.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The double ring vertices are useful for local mesh analysis and processing.

    See Also
    --------
    vertex_ring_vertices_iterators : Iterates over vertices in the vertex ring.
    """
    #import time
    #t0 = time.time()
    v, vj = self.vertex_ring_vertices_iterators(sort=True)
    M = coo_matrix((vj, (v, vj)), shape=(self.V, self.V))
    M = M.todense()
    ring = np.copy(M)
    while v.shape[0] > 0:
        vi, j = np.unique(v, True)
        ring[vi] += M[vj[j]]
        v = np.delete(v, j)
        vj = np.delete(vj, j)
    #t4 = time.time()
    #print(t4-t0)
    return ring.nonzero()

# -------------------------------------------------------------------------
#                               Ring lists
# -------------------------------------------------------------------------

def vertex_ring_vertices_list(self):
    """
    Returns a list of vertex rings for each vertex.

    This function constructs a list where each entry corresponds to a vertex and 
    contains the indices of vertices in its vertex ring.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    ring_list : list of lists
        A list where each entry is a list of vertex indices in the vertex ring.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The vertex ring vertices are useful for local mesh analysis and processing.

    See Also
    --------
    vertex_ring_vertices_iterators : Iterates over vertices in the vertex ring.
    """
    ring_list = [[] for i in range(self.V)]
    v, vj = self.vertex_ring_vertices_iterators(order=True)
    for i in range(len(v)):
        ring_list[v[i]].append(vj[i])
    return ring_list

def vertex_double_ring_vertices_list(self):
    """
    Returns a list of double rings for each vertex.

    This function constructs a list where each entry corresponds to a vertex and 
    contains the indices of vertices in its double ring.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    ring_list : list of lists
        A list where each entry is a list of vertex indices in the double ring.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The double ring vertices are useful for local mesh analysis and processing.

    See Also
    --------
    vertex_double_ring_vertices_iterators : Iterates over vertices in the double ring.
    
    vertex_ring_vertices_list : Returns a list of vertex rings for each vertex.
    """
    ring_list = [[] for i in range(self.V)]
    v, vj = self.vertex_double_ring_vertices_iterators()
    for i in range(len(v)):
        ring_list[v[i]].append(vj[i])
    return ring_list

def vertex_ring_edges_list(self):
    """
    Returns a list of edges in the vertex ring for each vertex.

    This function constructs a list where each entry corresponds to a vertex and 
    contains the indices of edges in its vertex ring.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    ring_list : list of lists
        A list where each entry is a list of edge indices in the vertex ring.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The vertex ring edges are useful for local mesh analysis and processing.

    See Also
    --------
    vertex_ring_edges_iterators : Iterates over edges in the vertex ring.
    
    vertex_ring_vertices_list : Returns a list of vertex rings for each vertex.
    """
    ring_list = [[] for i in range(self.V)]
    v, ej = self.vertex_ring_edges_iterators(order=True)
    for i in range(len(v)):
        ring_list[v[i]].append(ej[i])
    return ring_list

def vertex_ring_faces_list(self):
    """
    Returns a list of faces in the vertex ring for each vertex.

    This function constructs a list where each entry corresponds to a vertex and 
    contains the indices of faces in its vertex ring.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    ring_list : list of lists
        A list where each entry is a list of face indices in the vertex ring.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The vertex ring faces are useful for local mesh analysis and processing.

    See Also
    --------
    vertex_ring_faces_iterators : Iterates over faces in the vertex ring.
    
    vertex_ring_vertices_list : Returns a list of vertex rings for each vertex.
    """
    ring_list = [[] for i in range(self.V)]
    v, fj = self.vertex_ring_faces_iterators(order=True)
    for i in range(len(v)):
        ring_list[v[i]].append(fj[i])
    return ring_list

#--------------------------------------------------------------------------
#                                 Faces
#--------------------------------------------------------------------------

def face_lengths(self):
    """
    Computes the number of edges (length) for each face in the mesh.

    This function returns the number of edges for each face, which is equivalent 
    to the number of vertices forming the face.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    face_lengths : numpy array
        The number of edges for each face.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The face lengths are useful for mesh analysis and processing.

    See Also
    --------
    faces_list : Retrieves the list of faces.
    
    face_vertices_iterators : Iterates over vertices of each face.
    """
    H = self.halfedges
    f = H[H[:,1] >= 0,1]
    f = f[np.argsort(f)]
    i = np.ones((f.shape), 'i')
    lengths = utilities.sum_repeated(i, f)
    return lengths

def cell_arrays(self):
    """
    Converts the mesh faces to a cell array format.

    This function returns the faces in a cell array format, where each face is 
    represented by its vertices and a length indicator.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    cells : numpy array
        The cell array representation of the mesh faces.
    cell_types : numpy array
        The types of cells (e.g., triangles, quads).

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The cell arrays are useful for visualization and further processing.

    See Also
    --------
    faces_list : Retrieves the list of faces.

    face_lengths : Computes the number of edges for each face.
    """
    H = self.halfedges
    i  = self.face_ordered_halfedges()
    vi = H[i,0]
    f  = H[i,1]
    i = np.ones((f.shape[0]), 'i')
    j = np.arange(f.shape[0])
    _, k = np.unique(f, True)
    lengths = utilities.sum_repeated(i, f)
    index = j[k]
    cells = np.insert(vi, index, lengths)
    cell_types = lengths - 3
    cell_types[np.where(cell_types[:] > 2)[0]] = 2
    return cells, cell_types

def faces_list(self):
    """
    Retrieves the list of faces in the mesh.

    This function returns the faces of the mesh, where each face is represented 
    by a list of vertex indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    faces : list of lists
        The list of faces, where each face is a list of vertex indices.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The faces list is useful for mesh analysis and processing.

    See Also
    --------
    face_lengths : Computes the number of edges for each face.

    face_vertices_iterators : Iterates over vertices of each face.
    """
    faces_list = [[] for i in range(self.F)]
    fi, vj = self.face_vertices_iterators()
    for i in range(len(fi)):
        faces_list[fi[i]].append(vj[i])
    return faces_list

def face_triangles(self):
    """
    Decomposes each face into triangles.

    This function splits each face of the mesh into triangles, which is useful for 
    rendering or further processing. The resulting triangles are returned along 
    with their corresponding face indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    triangles : numpy array
        An array of triangles, where each triangle is represented by three vertex indices.
    face_indices : numpy array
        The indices of the original faces that each triangle belongs to.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The decomposition is useful for rendering or further processing.

    See Also
    --------
    faces_list : Retrieves the list of faces.

    face_lengths : Computes the number of edges for each face.
    """
    H = np.copy(self.halfedges)
    h = np.argsort(H[:,1])
    h = h[np.where(H[h,1] >= 0)]
    f = H[h,1]
    f_i, j = np.unique(f, True)
    one = np.arange(j.shape[0])
    f = np.delete(f, j)
    f = np.delete(f, j-one)
    f = np.delete(f, j-2*one)
    T = np.column_stack((H[j,0], H[H[j,2],0], H[H[H[j,2],2],0]))
    nex = H[H[H[j,2],2],2]
    face_index = f_i
    offset = 0
    while len(f) > 0:
        f_i, j = np.unique(f, True)
        T_i = np.column_stack((T[offset+f_i,-1], H[nex[f_i],0], T[f_i,0]))
        f = np.delete(f, j)
        nex = H[nex,2]
        T = np.vstack((T, T_i))
        face_index = np.hstack((face_index, f_i))
        offset += len(f_i)
    return T, face_index

# -------------------------------------------------------------------------
#                                  Edges
# -------------------------------------------------------------------------

def edge_vertices(self):
    """
    Retrieves the vertices of each edge in the mesh.

    This function returns the indices of the two vertices that form each edge.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    v1 : numpy array
        The indices of the first vertex of each edge.
    v2 : numpy array
        The indices of the second vertex of each edge.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The edge vertices are useful for local mesh analysis and processing.

    See Also
    --------
    edge_faces : Retrieves the faces adjacent to each edge.

    edge_lengths : Computes the lengths of all edges.
    """
    H  = self.halfedges
    v  = H[np.argsort(H[:,5]),0]
    v1 = v[0::2]
    v2 = v[1::2]
    return v1, v2

def edge_faces(self):
    """
    Retrieves the faces adjacent to each edge in the mesh.

    This function returns the indices of the two faces that share each edge.
    For boundary edges, one of the face indices will be -1.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    f1 : numpy array
        The indices of the first face adjacent to each edge.
    f2 : numpy array
        The indices of the second face adjacent to each edge.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The edge faces are useful for local mesh analysis and processing.

    See Also
    --------
    edge_vertices : Retrieves the vertices of each edge.

    are_boundary_edges : Checks if edges are on the boundary.
    """
    H  = self.halfedges
    f  = H[np.argsort(H[:,5]),1]
    f1 = f[0::2]
    f2 = f[1::2]
    return f1, f2

def vertices_edge_map(self):
    """
    Creates a sparse matrix mapping vertex pairs to edge indices.

    This function constructs a sparse matrix where each entry corresponds to a 
    pair of vertices and contains the index of the edge connecting them.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    edge_map : scipy.sparse.coo_matrix
        A sparse matrix mapping vertex pairs to edge indices.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The edge map is useful for efficient edge lookups and mesh analysis.

    See Also
    --------
    vertices_edge_faces_maps : Creates maps for vertices, edges, and faces.

    edge_vertices : Retrieves the vertices of each edge.
    """
    H  = self.halfedges
    v1 = H[:,0]
    v2 = H[H[:,4],0]
    e  = H[:,5]
    edge_map = coo_matrix((e, (v1,v2)), shape=(self.V, self.V))
    edge_map = edge_map.tocsc()
    return edge_map

def vertices_edge_faces_maps(self):
    """
    Creates sparse matrices mapping vertex pairs to edge and face indices.

    This function constructs sparse matrices where each entry corresponds to a 
    pair of vertices and contains the indices of the edge and adjacent faces.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    f1Map : scipy.sparse.coo_matrix
        A sparse matrix mapping vertex pairs to the first adjacent face.
    f2Map : scipy.sparse.coo_matrix
        A sparse matrix mapping vertex pairs to the second adjacent face.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The maps are useful for efficient lookups and mesh analysis.

    See Also
    --------
    vertices_edge_map : Creates a map for vertices and edges.

    edge_faces : Retrieves the faces adjacent to each edge.
    """
    H  = self.halfedges
    v1 = H[:,0]
    v2 = H[H[:,4],0]
    f1 = H[:,1]
    f2 = H[H[:,4],1]
    f1Map = coo_matrix((f1, (v1,v2)), shape=(self.V, self.V))
    f2Map = coo_matrix((f2, (v1,v2)), shape=(self.V, self.V))
    f1Map = f1Map.tocsc()
    f2Map = f2Map.tocsc()
    return f1Map, f2Map

def edge_halfedge(self, edge_index):
    """
    Retrieves a half-edge associated with a specified edge.

    This function returns one of the half-edges that corresponds to the given edge index.
    Each edge in the half-edge mesh is represented by a pair of half-edges.

    Parameters
    ----------
    edge_index : int
        The index of the edge.

    Returns
    -------
    halfedge_index : int
        The index of a half-edge associated with the specified edge.

    Note
    ----
    This function assumes that the specified edge exists in the mesh.
    The returned half-edge can be used for further mesh traversal and manipulation.

    See Also
    --------
    vertex_halfedge : Retrieves a half-edge originating from a specified vertex.
    """
    H = self.halfedges
    e = np.argsort(H[:,5])
    h1 = e[0::2]
    h2 = e[1::2]
    return h1, h2

# -------------------------------------------------------------------------
#                                 Boundary
# -------------------------------------------------------------------------

def boundary_vertices(self):
    """
    Retrieves the indices of boundary vertices in the mesh.

    This function identifies and returns the vertices that lie on the boundary of the mesh.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    boundary_vertices : numpy array
        The indices of vertices on the boundary.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary vertices are useful for mesh analysis and processing.

    See Also
    --------
    boundary_edges : Retrieves the indices of boundary edges.
    
    boundary_faces : Retrieves the indices of boundary faces.
    """
    H = self.halfedges
    b = np.where(H[:,1] == -1)[0]
    v = H[b,0]
    return v

def inner_halfedges(self):
    """
    Retrieves the indices of inner half-edges in the mesh.

    This function identifies and returns the half-edges that are not on the boundary.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    inner_halfedges : numpy array
        The indices of inner half-edges.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The inner half-edges are useful for mesh analysis and processing.

    See Also
    --------
    boundary_halfedges : Retrieves the indices of boundary half-edges.
    """
    H = self.halfedges
    h = np.where(H[:,1] != -1)[0]
    return h

def boundary_halfedges(self):
    """
    Retrieves the indices of boundary half-edges in the mesh.

    This function identifies and returns the half-edges that lie on the boundary.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    boundary_halfedges : numpy array
        The indices of boundary half-edges.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary half-edges are useful for mesh analysis and processing.

    See Also
    --------
    inner_halfedges : Retrieves the indices of inner half-edges.

    boundary_edges : Retrieves the indices of boundary edges.
    """
    H = self.halfedges
    b = np.where(H[:,1] == -1)[0]
    return b

def boundary_faces(self):
    """
    Retrieves the indices of boundary faces in the mesh.

    This function identifies and returns the faces that lie on the boundary of the mesh.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    boundary_faces : numpy array
        The indices of faces on the boundary.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary faces are useful for mesh analysis and processing.

    See Also
    --------
    boundary_vertices : Retrieves the indices of boundary vertices.

    boundary_edges : Retrieves the indices of boundary edges.
    """
    H = self.halfedges
    b = self.boundary_halfedges()
    e = H[H[b,4],1]
    return e

def inner_vertices(self):
    """
    Retrieves the indices of inner vertices in the mesh.

    This function identifies and returns the vertices that are not on the boundary of the mesh.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    inner_vertices : numpy array
        The indices of vertices that are not on the boundary.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The inner vertices are useful for mesh analysis and processing.

    See Also
    --------
    boundary_vertices : Retrieves the indices of boundary vertices.
    """
    b = self.boundary_vertices()
    v = np.arange(self.V)
    mask = np.invert(np.in1d(v, b))
    v = v[mask]
    return v

def boundary_curves(self, corner_split=False):
    """
    Retrieves the boundary curves of the mesh.

    This function identifies and returns the boundary curves of the mesh, 
    optionally splitting at corners.

    Parameters
    ----------
    corner_split : bool, optional (default=False)
        Whether to split boundary curves at corners.

    Returns
    -------
    boundary_curves : list of numpy arrays
        A list of boundary curves, where each curve is an array of vertex indices.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary curves are useful for mesh analysis and visualization.

    See Also
    --------
    boundary_vertices : Retrieves the indices of boundary vertices.

    boundary_edges : Retrieves the indices of boundary edges.
    """
    H = self.halfedges
    boundaries = []
    boundary_halfedges = []
    for h in range(H.shape[0]):
        if H[h,1] == -1 and h not in boundary_halfedges:
            boundary = []
            h_he = h
            boundary_halfedges.append(h_he)
            boundary.append(H[h_he,0])
            h_he = H[h_he,2]
            while h_he != h:
                boundary_halfedges.append(h_he)
                boundary.append(H[h_he,0])
                h_he = H[h_he,2]
            boundaries.append(np.array(boundary))
    if corner_split:
        corner_boundaries = []
        corners = self.mesh_corners()
        for boundary in boundaries:
            indices = np.arange(len(boundary))
            c = indices[np.in1d(boundary, corners)]
            boundary = np.split(boundary, c)
            for i in range(len(boundary) - 1):
                a = boundary[i]
                boundary[i] = np.insert(a, a.shape ,boundary[i+1][0])
            if len(boundary) > 1:
                boundary[0] = np.hstack((boundary[-1], boundary[0]))
                del boundary[-1]
            corner_boundaries.extend(boundary)
        boundaries = corner_boundaries
    return boundaries

def boundary_curves_halfedges(self, corner_split=False):
    """
    Retrieves the boundary curves of the mesh in terms of half-edges.

    This function identifies and returns the boundary curves of the mesh, 
    optionally splitting at corners, represented by half-edges.

    Parameters
    ----------
    corner_split : bool, optional (default=False)
        Whether to split boundary curves at corners.

    Returns
    -------
    boundary_curves : list of numpy arrays
        A list of boundary curves, where each curve is an array of half-edge indices.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary curves are useful for mesh analysis and processing.

    See Also
    --------
    boundary_curves : Retrieves the boundary curves in terms of vertices.

    boundary_halfedges : Retrieves the indices of boundary half-edges.
    """
    H = self.halfedges
    boundaries = []
    visited = []
    for h in range(H.shape[0]):
        if H[h,1] == -1 and h not in visited:
            boundary = []
            h_he = h
            boundary.append(h_he)
            h_he = H[h_he,2]
            while h_he != h:
                boundary.append(h_he)
                h_he = H[h_he,2]
            boundaries.append(np.array(boundary))
            visited.extend(boundary)
    if corner_split:
        corner_boundaries = []
        corners = self.mesh_corners()
        for boundary in boundaries:
            indices = np.arange(len(boundary))
            c = indices[np.in1d(H[boundary,0], corners)]
            boundary = np.split(boundary, c)
            if len(boundary) > 1:
                boundary[0] = np.hstack((boundary[-1], boundary[0]))
                del boundary[-1]
            corner_boundaries.extend(boundary)
        boundaries = corner_boundaries
    return boundaries

def boundary_polylines(self):
    """
    Converts the boundary curves of the mesh to polyline objects.

    This function creates polyline objects from the boundary curves of the mesh.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    polylines : list of Polyline objects
        A list of polyline objects representing the boundary curves.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The polylines are useful for visualization and further processing.

    See Also
    --------
    boundary_curves : Retrieves the boundary curves in terms of vertices.

    Polyline : A class representing a polyline object.
    """
    polylines = []
    curves = self.boundary_curves(corner_split=False)
    for curve in curves:
        polyline = Polyline(self.vertices[curve,:], closed=True)
        polyline.corner_tolerance = self.corner_tolerance
        polylines.append(polyline)
    return polylines

def mesh_corners(self):
    """
    Identifies the corner vertices of the mesh.

    This function detects and returns the vertices that are considered corners 
    based on the angle between adjacent edges.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    corners : numpy array
        The indices of corner vertices.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The corner vertices are useful for mesh analysis and processing.

    See Also
    --------
    boundary_curves : Retrieves the boundary curves in terms of vertices.

    boundary_curves_halfedges : Retrieves the boundary curves in terms of half-edges.
    """
    H = self.halfedges
    b = np.where(H[:,1] == -1)[0]
    v0 = H[b,0]
    vp = H[H[b,3],0]
    vn = H[H[b,2],0]
    Vp = self.vertices[v0,:] - self.vertices[vp,:]
    Vn = self.vertices[vn,:] - self.vertices[v0,:]
    Vp = Vp / np.linalg.norm(Vp, axis=1, keepdims=True)
    Vn = Vn / np.linalg.norm(Vn, axis=1, keepdims=True)
    C = np.einsum('ij,ij->i', Vp, Vn)
    corners = v0[np.where(C[:] < self.corner_tolerance)[0]]
    return corners

def double_boundary_vertices(self):
    """
    Identifies vertices that are on the boundary of two faces.

    This function detects and returns vertices that are shared by two boundary faces.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    vertices : numpy array
        The indices of vertices that are on the boundary of two faces.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The double boundary vertices are useful for mesh analysis and processing.

    See Also
    --------
    boundary_vertices : Retrieves the indices of boundary vertices.

    boundary_faces : Retrieves the indices of boundary faces.
    """
    bf = self.are_boundary_faces()
    v, fj = self.vertex_ring_faces_iterators(sort=True)
    bf = bf[fj]
    bf = utilities.sum_repeated(bf, v)
    bf = np.where(bf > 0)[0]
    return bf

def boundary_edges(self):
    """
    Retrieves the indices of boundary edges in the mesh.

    This function identifies and returns the edges that lie on the boundary of the mesh.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    boundary_edges : numpy array
        The indices of edges on the boundary.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary edges are useful for mesh analysis and processing.

    See Also
    --------
    boundary_vertices : Retrieves the indices of boundary vertices.

    boundary_faces : Retrieves the indices of boundary faces.
    """
    H = self.halfedges
    ind = np.where(H[:,1] == -1)
    e = H[ind,5]
    e = np.unique(e)
    return e

# -------------------------------------------------------------------------
#                              Global queries
# -------------------------------------------------------------------------

def are_boundary_edges(self):
    """
    Checks if edges are on the boundary of the mesh.

    This function returns a boolean array indicating whether each edge is on the boundary.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    is_boundary : numpy array (bool)
        A boolean array where True indicates that the edge is on the boundary.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary edges are useful for mesh analysis and processing.

    See Also
    --------
    boundary_edges : Retrieves the indices of boundary edges.

    are_boundary_faces : Checks if faces are on the boundary.
    """
    H = self.halfedges
    B = H[np.argsort(H[:,5])]
    B = B[:,1] == -1
    bound = np.logical_or(B[0::2], B[1::2])
    return bound

def are_boundary_faces(self):
    """
    Checks if faces are on the boundary of the mesh.

    This function returns a boolean array indicating whether each face is on the boundary.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    is_boundary : numpy array (bool)
        A boolean array where True indicates that the face is on the boundary.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary faces are useful for mesh analysis and processing.

    See Also
    --------
    boundary_faces : Retrieves the indices of boundary faces.

    are_boundary_edges : Checks if edges are on the boundary.
    """
    H = self.halfedges
    f = np.where(H[:,1] != -1)[0]
    B = H[H[f,4],1] == -1
    i = np.argsort(H[f,1])
    bound = utilities.sum_repeated(B[i], H[f,1])
    return bound

# -------------------------------------------------------------------------
#                                Normals
# -------------------------------------------------------------------------

def face_vector_areas(self):
    """
    Computes the vector areas of all faces in the mesh.

    This function calculates the vector area of each face, which is the cross 
    product of two edge vectors scaled by half. The vector area is a vector 
    perpendicular to the face with a magnitude equal to the face area.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    vector_areas : numpy array
        An array of vector areas for each face.

    Note
    ----
    This function assumes that the mesh faces are planar.
    The vector areas are computed based on the current vertex positions.

    See Also
    --------
    face_areas : Computes the scalar areas of all faces.

    face_normals : Computes the normal vectors of all faces.
    """
    f, v1, v2 = self.face_edge_vertices_iterators(order=True)
    V1 = self.vertices[v1,:]
    V2 = self.vertices[v2,:]
    N  = np.cross(V1,V2)
    normals = utilities.sum_repeated(N, f)
    return 0.5 * normals

def face_normals(self):
    """
    Computes the normal vectors for each face in the mesh.

    This function calculates the face normals by taking the cross product of two edges
    of each face and normalizing the resulting vector.

    The main logic involves:

    1. Iterating over each face to extract its vertices.

    2. Computing the cross product of two edges to obtain the normal vector.

    3. Normalizing the normal vector to ensure unit length.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    normals : numpy array
        An array of face normals, where each normal is a unit vector [nx, ny, nz].

    Note
    ----
    This function assumes that the mesh faces are planar and that the half-edge structure is correctly defined.
    The normals are computed based on the current vertex positions.

    See Also
    --------
    vertex_normals : Computes the normal vectors for each vertex.

    edge_normals : Computes the normal vectors for each edge.
    """
    N = self.face_vector_areas()
    N = N / np.linalg.norm(N, axis=1, keepdims=True)
    return N

def vertex_normals(self):
    """
    Computes the normal vectors for each vertex in the mesh.

    This function calculates the vertex normals by averaging the face normals
    of the faces adjacent to each vertex. The resulting normals are normalized
    to ensure they are unit vectors.

    The main logic involves:

    1. Iterating over each vertex to identify its adjacent faces.

    2. Summing the face normals of these adjacent faces.

    3. Normalizing the summed vector to obtain the vertex normal.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    normals : numpy array
        An array of vertex normals, where each normal is a unit vector [nx, ny, nz].

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The vertex normals are computed based on the current face normals and mesh topology.

    See Also
    --------
    face_normals : Computes the normal vectors for each face.

    edge_normals : Computes the normal vectors for each edge.
    """
    N = self.face_vector_areas()
    v, fi = self.vertex_ring_faces_iterators(sort=True)
    N = N[fi,:]
    normals = utilities.sum_repeated(N, v)
    normals = normals / np.linalg.norm(normals, axis=1, keepdims=True)
    return normals

def edge_normals(self):
    """
    Computes the normal vectors for each edge in the mesh.

    This function calculates the edge normals by averaging the face normals
    of the faces adjacent to each edge. If an edge is on the boundary, its
    normal is computed based on the single adjacent face.

    The main logic involves:

    1. Identifying the faces adjacent to each edge.

    2. Averaging the face normals of these adjacent faces.

    3. Normalizing the resulting vector to obtain the edge normal.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    normals : numpy array
        An array of edge normals, where each normal is a unit vector [nx, ny, nz].

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The edge normals are computed based on the current face normals and mesh topology.

    See Also
    --------
    face_normals : Computes the normal vectors for each face.

    vertex_normals : Computes the normal vectors for each vertex.
    """
    N = self.face_normals()
    N = np.insert(N, N.shape[0], 0, axis=0)
    f1, f2 = self.edge_faces()
    normals = N[f1] + N[f2]
    normals = normals / np.linalg.norm(normals, axis=1, keepdims=True)
    return normals

def boundary_normals(self):
    """
    Computes the normal vectors for boundary edges in the mesh.

    This function calculates the boundary normals by considering the adjacent
    faces and edges. For each boundary edge, it computes a normal vector that
    is perpendicular to the edge and lies in the plane of the adjacent face.

    The main logic involves:

    1. Identifying boundary edges and their adjacent faces.

    2. Computing the cross product of edge vectors and face normals.

    3. Normalizing the resulting vectors to obtain unit normals.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    normals : numpy array
        An array of boundary normals, where each normal is a unit vector [nx, ny, nz].

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary normals are computed based on the current mesh topology and geometry.

    See Also
    --------
    boundary_edges : Retrieves the indices of boundary edges.

    edge_normals : Computes the normal vectors for each edge.
    """
    H = self.halfedges
    b = np.where(H[:,1] == -1)[0]
    face_normals = self.face_normals()
    N1 = face_normals[H[H[b,4],1]]
    N2 = face_normals[H[H[H[b,3],4],1]]
    normals = np.zeros((self.V, 3))
    E1 = self.vertices[H[H[b,2],0]] - self.vertices[H[b,0]]
    E2 = self.vertices[H[b,0]] - self.vertices[H[H[b,3],0]]
    N = np.cross(N1,E1) + np.cross(N2,E2)
    N = N / np.linalg.norm(N, axis=1, keepdims=True)
    normals[H[b,0],:] = N
    return normals

def boundary_tangents(self, normalize=True):
    """
    Computes the tangent vectors for boundary edges in the mesh.

    This function calculates the tangent vectors along the boundary edges.
    The tangents are computed as the direction vectors of the boundary edges.

    The main logic involves:

    1. Identifying boundary edges.

    2. Computing the direction vectors of these edges.

    3. Optionally normalizing the tangent vectors.

    Parameters
    ----------
    normalize : bool, optional (default=True)
        Whether to normalize the tangent vectors to unit length.

    Returns
    -------
    tangents : numpy array
        An array of boundary tangents, where each tangent is a vector [tx, ty, tz].

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The boundary tangents are computed based on the current mesh topology and geometry.

    See Also
    --------
    boundary_edges : Retrieves the indices of boundary edges.

    boundary_normals : Computes the normal vectors for boundary edges.
    """
    H = self.halfedges
    b = np.where(H[:,1] == -1)[0]
    V1 = self.vertices[H[H[b,3],0]]
    V2 = self.vertices[H[H[b,2],0]]
    T = (V2 - V1)
    if normalize:
        T = T / np.linalg.norm(T, keepdims=True)
    else:
        T = T/2
    tangents = np.zeros((self.V, 3))
    tangents[H[b,0],:] = T
    return tangents

# -------------------------------------------------------------------------
#                                  Area
# -------------------------------------------------------------------------

def face_areas(self):
    """
    Computes the areas of all faces in the mesh.

    This function calculates the area of each face using the cross product
    of its edge vectors. For triangular faces, the area is half the magnitude
    of the cross product of two edges.

    The main logic involves:

    1. Iterating over each face to extract its vertices.

    2. Computing the cross product of two edge vectors.
    
    3. Calculating the area as half the magnitude of the cross product.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    areas : numpy array
        An array of face areas.

    Note
    ----
    This function assumes that the mesh faces are planar.
    The face areas are computed based on the current vertex positions.

    See Also
    --------
    area : Computes the total surface area of the mesh.

    face_vector_areas : Computes the vector areas of faces.
    """
    N = self.face_vector_areas()
    A = np.linalg.norm(N, axis=1)
    return A

def area(self):
    """
    Computes the total surface area of the mesh.

    This function calculates the total surface area by summing the areas
    of all faces in the mesh.

    The main logic involves:

    1. Calling `face_areas` to compute the area of each face.

    2. Summing the face areas to obtain the total surface area.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    total_area : float
        The total surface area of the mesh.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The total area is computed based on the current face areas.

    See Also
    --------
    face_areas : Computes the areas of individual faces.
    """
    A = self.face_areas()
    A = np.sum(A)
    return A

def vertex_ring_areas(self):
    """
    Computes the average face areas around each vertex.

    This function calculates the average face area for each vertex by
    summing the areas of adjacent faces and dividing by the number of faces.

    The main logic involves:

    1. Iterating over each vertex to identify its adjacent faces.

    2. Summing the areas of these adjacent faces.

    3. Dividing by the number of adjacent faces to obtain the average area.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    ring_areas : numpy array
        An array of average face areas around each vertex.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The vertex ring areas are computed based on the current face areas and mesh topology.

    See Also
    --------
    face_areas : Computes the areas of individual faces.

    vertex_ring_faces_iterators : Iterates over faces adjacent to each vertex.
    """
    L = self.face_lengths()
    A = self.face_areas()
    v, fi = self.vertex_ring_faces_iterators(sort=True)
    ring_area = A[fi]/L[fi]
    ring_area = utilities.sum_repeated(ring_area, v)
    return ring_area

# -------------------------------------------------------------------------
#                               Closeness
# -------------------------------------------------------------------------

def make_kdtree(self):
    """
    Constructs a k-d tree from the mesh vertices for efficient nearest-neighbor queries.

    This function creates a k-d tree data structure using the vertex positions,
    allowing for fast nearest-neighbor searches.

    The main logic involves:

    1. Creating a k-d tree from the vertex positions.

    2. Storing the k-d tree within the mesh instance for later use.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    None
        The k-d tree is stored within the mesh instance.

    Note
    ----
    This function assumes that the mesh vertices are valid and non-empty.
    The k-d tree is used to accelerate nearest-neighbor queries.

    See Also
    --------
    closest_vertices : Finds the closest vertices to a set of query points.
    """
    kdtree = spatial.cKDTree(self.vertices)
    self._kdtree = kdtree

def closest_vertices(self, points, make_tree=False):
    """
    Finds the closest vertices in the mesh to a set of query points.

    This function uses a k-d tree to efficiently find the closest vertices 
    in the mesh to the given points. If the k-d tree does not exist, it will 
    be constructed automatically.

    Parameters
    ----------
    points : numpy array
        An array of query points [x, y, z].
    make_tree : bool, optional (default=False)
        Whether to force the reconstruction of the k-d tree.

    Returns
    -------
    closest : numpy array
        The indices of the closest vertices in the mesh.

    Note
    ----
    This function assumes that the mesh vertices are valid and non-empty.
    The k-d tree is used to accelerate nearest-neighbor queries.

    See Also
    --------
    make_kdtree : Constructs a k-d tree from the mesh vertices.
    """
    if self._kdtree is None:
        self.make_kdtree()
    elif make_tree:
        self.make_kdtree()
    closest = self._kdtree.query(points)[1]
    return closest

# -------------------------------------------------------------------------
#                                Geometry
# -------------------------------------------------------------------------

def edge_mid_points(self):
    """
    Computes the midpoints of all edges in the mesh.

    This function calculates the midpoint of each edge by averaging the positions
    of its two vertices.

    The main logic involves:

    1. Iterating over all edges.

    2. Computing the midpoint as the average of the two vertex positions.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    midpoints : numpy array
        An array of edge midpoints, where each midpoint is a vector [x, y, z].

    Note
    ----
    This function assumes that the mesh edges are valid and non-degenerate.
    The edge midpoints are computed based on the current vertex positions.

    See Also
    --------
    edge_lengths : Computes the lengths of all edges.

    edge_vectors : Computes the direction vectors of all edges.
    """
    v1, v2 = self.edge_vertices()
    M = 0.5*(self.vertices[v1] + self.vertices[v2])
    return M

def edge_vectors(self):
    """
    Computes the direction vectors of all edges in the mesh.

    This function calculates the direction vector of each edge by subtracting
    the position of the starting vertex from the position of the ending vertex.

    The main logic involves:

    1. Iterating over all edges.

    2. Computing the direction vector as the difference between the two vertex positions.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    edge_vectors : numpy array
        An array of edge direction vectors, where each vector is [dx, dy, dz].

    Note
    ----
    This function assumes that the mesh edges are valid and non-degenerate.
    The edge vectors are computed based on the current vertex positions.

    See Also
    --------
    edge_lengths : Computes the lengths of all edges.

    edge_mid_points : Computes the midpoints of all edges.
    """
    v1, v2 = self.edge_vertices()
    Ev = self.vertices[v1] - self.vertices[v2]
    return Ev

def face_barycenters(self):
    """
    Computes the barycenters (geometric centers) of all faces in the mesh.

    This function calculates the barycenter of each face by averaging the positions
    of its vertices.

    The main logic involves:

    1. Iterating over all faces.

    2. Computing the barycenter as the average of the vertex positions.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    barycenters : numpy array
        An array of face barycenters, where each barycenter is a vector [x, y, z].

    Note
    ----
    This function assumes that the mesh faces are valid and non-degenerate.
    The barycenters are computed based on the current vertex positions.

    See Also
    --------
    face_areas : Computes the areas of all faces.

    face_circum_circles : Computes the circumcircles of all faces.
    """
    H = self.halfedges
    H = H[np.where(H[:,1] >= 0)[0],:]
    i = np.argsort(H[:,1])
    f = H[i,1]
    v = H[i,0]
    V = self.vertices[v,:]
    B = utilities.sum_repeated(V,f)
    L = self.face_lengths()
    L = np.column_stack((L,L,L))
    B = B/L
    return B

def edge_lengths(self):
    """
    Computes the lengths of all edges in the mesh.

    This function calculates the length of each edge using the Euclidean distance
    between its two vertices.

    The main logic involves:

    1. Iterating over all edges.

    2. Computing the length as the Euclidean distance between the two vertex positions.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    lengths : numpy array
        An array of edge lengths.

    Note
    ----
    This function assumes that the mesh edges are valid and non-degenerate.
    The edge lengths are computed based on the current vertex positions.

    See Also
    --------
    edge_vectors : Computes the direction vectors of all edges.

    mean_edge_length : Computes the mean length of all edges.
    """
    v1, v2 = self.edge_vertices()
    V1 = self.vertices[v1]
    V2 = self.vertices[v2]
    V = V1 - V2
    L = np.linalg.norm(V, axis=1)
    return L

def mean_edge_length(self):
    """
    Computes the mean length of all edges in the mesh.

    This function calculates the mean edge length by averaging the lengths of all edges.

    The main logic involves:

    1. Calling `edge_lengths` to compute the lengths of all edges.

    2. Computing the mean length by summing all edge lengths and dividing by the number of edges.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    mean_length : float
        The mean length of all edges in the mesh.

    Note
    ----
    This function assumes that the mesh edges are valid and non-degenerate.
    The mean edge length is computed based on the current edge lengths.

    See Also
    --------
    edge_lengths : Computes the lengths of all edges.
    """
    L = self.edge_lengths()
    M = np.sum(L) / L.shape[0]
    return M

def face_planarity(self, scale_invariant=True):
    """
    Computes the planarity of each face in the mesh.

    This function calculates the planarity of each face by measuring the deviation
    from a perfect plane. The planarity can be computed as either absolute or
    scale-invariant.

    The main logic involves:

    1. Iterating over all faces.

    2. Computing the deviation from planarity using the cross product of edge vectors.
    
    3. Normalizing the deviation if scale-invariant planarity is requested.

    Parameters
    ----------
    scale_invariant : bool, optional (default=True)
        Whether to compute scale-invariant planarity.

    Returns
    -------
    planarity : numpy array
        An array of face planarity values.

    Note
    ----
    This function assumes that the mesh faces are valid and non-degenerate.
    The planarity is computed based on the current vertex positions.

    See Also
    --------
    face_areas : Computes the areas of all faces.

    face_normals : Computes the normal vectors of all faces.
    """
    planarity = np.zeros((self.F))
    f, vi = self.face_vertices_iterators()
    i = np.ones((f.shape[0]),dtype=int)
    j = np.arange(f.shape[0])
    _, k = np.unique(f, True)
    L = utilities.sum_repeated(i, f)
    index = j[k]
    quad = np.where(L > 3)[0]
    shift = 0
    while len(quad) > 0:
        P1 = self.vertices[vi[index[quad] + shift]]
        P2 = self.vertices[vi[index[quad] + shift + 1]]
        P3 = self.vertices[vi[index[quad] + shift + 2]]
        P4 = self.vertices[vi[index[quad] + shift + 3]]
        V1 = P3 - P1
        V2 = P4 - P2
        N  = np.cross(V1,V2)
        eps = np.finfo(float).eps
        P12 = P2 - P1
        norm = ((np.linalg.norm(N, axis=1) + eps))
        d = np.einsum('ij,ij->i', P12, N) / norm
        if scale_invariant:
            d1 = np.linalg.norm(V1, axis=1)
            d2 = np.linalg.norm(V1, axis=1)
            p = np.abs(d) / ((d1 + d2)/2)
        else:
            p = np.abs(d)
        planarity[quad] = np.maximum(p, planarity[quad])
        L -= 1
        shift += 1
        quad = np.where(L > 3)[0]
    return planarity

def edge_versors(self):
    """
    Computes the unit direction vectors (versors) of all edges in the mesh.

    This function normalizes the direction vectors of each edge to obtain unit vectors.

    The main logic involves:

    1. Calling `edge_vectors` to compute the direction vectors of all edges.

    2. Normalizing each direction vector to obtain a unit vector.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    versors : numpy array
        An array of edge versors, where each versor is a unit vector [dx, dy, dz].

    Note
    ----
    This function assumes that the mesh edges are valid and non-degenerate.
    The edge versors are computed based on the current vertex positions.

    See Also
    --------
    edge_vectors : Computes the direction vectors of all edges.

    edge_lengths : Computes the lengths of all edges.
    """
    v1, v2 = self.edge_vertices()
    ver = self.vertices[v2] - self.vertices[v1]
    ver = ver / np.linalg.norm(ver, axis=1, keepdims=True)
    return ver

def bounding_box(self):
    """
    Computes the axis-aligned bounding box of the mesh.

    This function calculates the minimum and maximum coordinates of the mesh
    along each axis (x, y, z).

    The main logic involves:

    1. Finding the minimum and maximum x, y, and z coordinates of all vertices.

    2. Returning the bounding box as a list of ranges for each axis.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    bbox : list of tuples
        The bounding box, represented as [(xmin, xmax), (ymin, ymax), (zmin, zmax)].

    Note
    ----
    This function assumes that the mesh vertices are valid and non-empty.
    The bounding box is computed based on the current vertex positions.

    See Also
    --------
    mesh_center : Computes the geometric center of the mesh.
    """
    Xmin = np.min(self.vertices[:,0])
    Xmax = np.max(self.vertices[:,0])
    Ymin = np.min(self.vertices[:,1])
    Ymax = np.max(self.vertices[:,1])
    Zmin = np.min(self.vertices[:,2])
    Zmax = np.max(self.vertices[:,2])
    return ([Xmin, Xmax], [Ymin, Ymax], [Zmin, Zmax])

def mesh_center(self):
    """
    Computes the geometric center of the mesh.

    This function calculates the center of the mesh by averaging the coordinates
    of its bounding box.

    The main logic involves:

    1. Calling `bounding_box` to compute the axis-aligned bounding box.

    2. Computing the center as the midpoint of each axis range.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    center : numpy array
        The geometric center of the mesh, represented as [x, y, z].

    Note
    ----
    This function assumes that the mesh vertices are valid and non-empty.
    The center is computed based on the current bounding box.

    See Also
    --------
    bounding_box : Computes the axis-aligned bounding box of the mesh.
    """
    B = self.bounding_box()
    x = (B[0][0] + B[0][1])/2
    y = (B[1][0] + B[1][1])/2
    z = (B[2][0] + B[2][1])/2
    return np.array([x,y,z])

def face_circum_circles(self):
    """
    Computes the circumcircles of all faces in the mesh.

    This function calculates the circumcircle of each face by finding the circle
    that passes through its vertices. For triangular faces, the circumcircle is
    defined by the three vertices.

    The main logic involves:

    1. Iterating over all faces.
    
    2. Computing the circumcircle using the vertices of each face.

    3. Returning the center and radius of each circumcircle.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    circumcircles : list of tuples
        A list of circumcircle information for each face, where each entry is 
        (center_x, center_y, center_z, radius).

    Note
    ----
    This function assumes that the mesh faces are valid and non-degenerate.
    The circumcircles are computed based on the current vertex positions.

    See Also
    --------
    face_barycenters : Computes the barycenters of all faces.

    face_areas : Computes the areas of all faces.
    """
    f, vi = self.face_vertices_iterators()
    _, j = np.unique(f,True)
    p1 = self.vertices[vi[j]]
    p2 = self.vertices[vi[j+1]]
    p3 = self.vertices[vi[j+2]]
    circles = circle_three_points(p1, p2, p3)
    return circles

# -------------------------------------------------------------------------
#                                 Topology
# -------------------------------------------------------------------------

def flip_normals(self):
    """
    Flips the normals of the mesh by reversing the half-edge orientation.

    This function inverts the direction of all face normals by swapping the 
    orientation of the half-edges. This effectively reverses the winding order 
    of the faces.

    The main logic involves:

    1. Iterating over all half-edges.

    2. Swapping the origin and destination vertices of each half-edge.

    3. Updating the face normals accordingly.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    None
        The mesh normals are updated in place.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The normals are flipped based on the current half-edge structure.

    See Also
    --------
    face_normals : Computes the normal vectors of all faces.

    vertex_normals : Computes the normal vectors of all vertices.
    """
    H = self.halfedges
    H[:,0] = H[H[:,2],0]
    H[:,[2,3]] = H[:,[3,2]]

def orient_faces(self, vertices_list, faces_list):
    """
    Reorients the faces of the mesh to ensure consistent winding order.

    This function processes the input vertices and faces to ensure that all faces 
    have a consistent orientation. This is necessary for creating a valid half-edge 
    structure.

    The main logic involves:

    1. Iterating over all faces and checking their orientation.

    2. Reversing the vertex order of any face with inconsistent orientation.

    3. Returning the reoriented face list.

    Parameters
    ----------
    vertices_list : list of lists
        The list of vertex coordinates [x, y, z].
    faces_list : list of lists
        The list of faces, where each face is defined by vertex indices.

    Returns
    -------
    oriented_faces : list of lists
        The reoriented face list with consistent winding order.

    Note
    ----
    This function assumes that the input mesh is manifold and orientable.
    The reoriented faces are necessary for constructing a valid half-edge structure.

    See Also
    --------
    make_mesh : Constructs the half-edge mesh from vertices and faces.
    """
    F = len(faces_list)
    V = len(vertices_list)
    fmap = -np.ones((V,V), dtype=int)
    inconsistent = np.zeros((V,V), dtype=int)
    flipped = np.zeros(F, dtype=bool)
    oriented = np.zeros(F, dtype=bool)
    oriented_faces = copy.deepcopy(faces_list)
    for f in range(F):
        face = faces_list[f]
        for j in range(len(face)):
            v0 = face[j-1]
            v1 = face[j]
            if fmap[v0,v1] == -1:
                fmap[v0,v1] = f
            else:
                fmap[v1,v0] = f
                inconsistent[v0,v1] = True
                inconsistent[v1,v0] = True
    ring = [0]
    oriented[0] = True
    i = 1
    while len(ring) > 0:
        next_ring = []
        for f in ring:
            face = faces_list[f]
            for j in range(len(face)):
                flip = False
                v0 = face[j-1]
                v1 = face[j]
                if fmap[v0,v1] == f:
                    v2 = v1
                    v3 = v0
                else:
                    v2 = v0
                    v3 = v1
                if inconsistent[v2,v3] and not flipped[f]:
                    flip = True
                if not inconsistent[v2,v3] and flipped[f]:
                    flip = True
                fi = fmap[v2,v3]
                if fi != -1 and not oriented[fi]:
                    if fi not in next_ring:
                        next_ring.append(fi)
                    if flip:
                        oriented_faces[fi].reverse()
                        flipped[fi] = True
                    i += 1
                    oriented[fi] = True
                    if i == F:
                        return oriented_faces
        ring = next_ring
        if len(ring) == 0:
            try:
                ring = [np.where(oriented == False)[0][0]]
            except:
                return

def vertex_ring_expansion(self, v_index, callback=None, depth=None):
    """
    Expands the vertex ring around a specified vertex up to a given depth.

    This function performs a breadth-first traversal of the mesh starting from 
    the specified vertex and collects all vertices within the given depth.

    The main logic involves:

    1. Initializing the traversal from the specified vertex.

    2. Iteratively expanding the vertex ring by following adjacent vertices.

    3. Optionally applying a callback function to each vertex during traversal.

    Parameters
    ----------
    v_index : int
        The index of the starting vertex.
    callback : function, optional (default=None)
        An optional callback function to apply to each vertex.
    depth : int, optional (default=None)
        The maximum depth of the traversal. If None, expands to the entire mesh.

    Returns
    -------
    expanded_vertices : numpy array
        The indices of all vertices within the specified depth.

    Note
    ----
    This function assumes that the mesh is manifold and connected.
    The vertex ring expansion is based on the current mesh topology.

    See Also
    --------
    vertex_ring_vertices_iterators : Iterates over vertices in the vertex ring.
    """
    vi, vj = self.vertex_ring_vertices_iterators()
    mring = np.full(self.V, False)
    sring = np.full(self.V, False)
    search = np.array([v_index], dtype='i')
    mring[v_index] = True
    if depth is None:
        depth = self.V
    for i in range(depth):
        sring[:] = False
        for v in search:
            ring = vj[vi == v]
            ring = ring[np.invert(mring[ring])]
            mring[ring] = True
            sring[ring] = True
            if callable(callback):
                callback(v, ring)
        search = np.where(sring)[0]
        if np.all(mring):
            return np.where(mring)[0]
    return np.where(mring)[0]

def make_simply_connected(self):
    """
    Modifies the mesh to make it simply connected by removing boundary loops.

    This function cuts the mesh along boundary edges to eliminate multiple 
    boundary loops, resulting in a simply connected mesh.

    The main logic involves:

    1. Identifying all boundary loops in the mesh.

    2. Cutting the mesh along boundary edges to merge boundary loops.

    3. Reconstructing the mesh topology after cutting.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    None
        The mesh topology is updated in place.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The resulting mesh will have a single boundary loop.

    See Also
    --------
    boundary_curves : Retrieves the boundary loops of the mesh.

    cut : Cuts the mesh along a specified boundary edge.
    """
    curves = self.boundary_curves(corner_split=False)
    while len(curves) > 1:
        curve = curves[0]
        i = 0
        v = curve[i]
        while len(self.vertex_ring_vertices(v)) < 3:
            i += 1
            v = curve[i]
        self.cut(v)
        curves = self.boundary_curves(corner_split=False)

# -------------------------------------------------------------------------
#                                 Curves
# -------------------------------------------------------------------------

def mesh_curves(self):
    """
    Extracts the curves from the mesh, including boundary and internal curves.

    This function identifies and extracts curves from the mesh by following 
    edges and vertices. The curves can be used for analysis or visualization.

    The main logic involves:

    1. Identifying boundary and internal edges.

    2. Following edges to form continuous curves.

    3. Returning the extracted curves as lists of vertex indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    curves : list of lists
        A list of curves, where each curve is represented by a list of vertex indices.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The curves are extracted based on the current mesh topology.

    See Also
    --------
    boundary_curves : Retrieves the boundary loops of the mesh.

    mesh_polylines : Converts curves to polyline objects.
    """
    _,_, valence = self.vertex_ring_vertices_iterators(return_lengths=True)
    boundary_vertices = self.boundary_vertices()
    H = self.halfedges
    boundaries = self.boundary_curves_halfedges(True)
    done = []
    curves = []
    for boundary in boundaries:
        family = []
        for h in boundary:
            if H[h,0] not in done:
                curve = [H[h,0]]
                if valence[H[h,0]] <= 3:
                    turn = 1
                else:
                    turn = 2
                for i in range(turn):
                    h = H[H[h,4],2]
                vertex = H[H[h,4],0]
                stop = False
                exclude = False
                if vertex in boundary_vertices:
                    stop = True
                    exclude = True
                while not stop:
                    curve.append(vertex)
                    if vertex in boundary_vertices:
                        stop = True
                        done.append(vertex)
                    if valence[vertex] <= 4:
                        turn = 1
                    else:
                        turn = 2
                    for i in range(turn):
                        h = H[H[H[h,2],4],2]
                    vertex = H[H[h,4],0]
                if not exclude:
                    family.append(curve)
        if len(family) > 0:
            curves.append(family)
    curves.append(self.boundary_curves(True))
    return curves

def mesh_polylines(self):
    """
    Converts the mesh curves to polyline objects for visualization.

    This function takes the extracted curves from the mesh and converts them 
    into polyline objects, which can be used for visualization or further processing.

    The main logic involves:

    1. Calling `mesh_curves` to extract the curves from the mesh.
    
    2. Converting each curve to a polyline object using vertex positions.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    polylines : list of Polyline objects
        A list of polyline objects representing the mesh curves.

    Note
    ----
    This function assumes that the mesh curves are valid and non-degenerate.
    The polylines are created based on the current vertex positions.

    See Also
    --------
    mesh_curves : Extracts the curves from the mesh.

    Polyline : A class representing a polyline object.
    """
    curves = self.mesh_curves()
    polylines = []
    for family in curves:
        poly_family = []
        for curve in family:
            poly_family.append(Polyline(self.vertices[curve,:]))
        polylines.append(poly_family)
    return polylines

# -------------------------------------------------------------------------
#                             Transformations
# -------------------------------------------------------------------------

def move(self, displacement_vector):
    """
    Translates the mesh by a specified displacement vector.

    This function moves the entire mesh by adding the displacement vector to 
    each vertex position.

    The main logic involves:

    1. Adding the displacement vector to each vertex coordinate.

    2. Updating the vertex positions in place.

    Parameters
    ----------
    displacement_vector : list or numpy array
        The displacement vector [dx, dy, dz] to apply to the mesh.

    Returns
    -------
    None
        The mesh vertices are updated in place.

    Note
    ----
    This function assumes that the displacement vector is valid.
    The mesh is translated based on the current vertex positions.

    See Also
    --------
    scale : Scales the mesh by a specified factor.
    """
    self.vertices[:,[0,1,2]] += np.array(displacement_vector)[[0,1,2]]

def scale(self, factor, center=[0, 0, 0]):
    """
    Scales the mesh by a specified factor relative to a center point.

    This function scales the mesh by multiplying each vertex position by the 
    specified factor, relative to a given center point.

    The main logic involves:

    1. Translating the mesh to the origin by subtracting the center point.

    2. Scaling the mesh by the specified factor.

    3. Translating the mesh back to the original center point.

    Parameters
    ----------
    factor : float
        The scaling factor to apply to the mesh.
    center : list or numpy array, optional (default=[0, 0, 0])
        The center point relative to which the mesh is scaled.

    Returns
    -------
    None
        The mesh vertices are updated in place.

    Note
    ----
    This function assumes that the scaling factor is non-zero.
    The mesh is scaled based on the current vertex positions.

    See Also
    --------
    move : Translates the mesh by a specified displacement vector.
    """
    self.vertices[:,:] *= factor
    self.vertices[:,0] -= center[0]
    self.vertices[:,1] -= center[1]
    self.vertices[:,2] -= center[2]

# -------------------------------------------------------------------------
#                     Discrete Differential Geometry
# -------------------------------------------------------------------------

def vertex_ring_parametrization(self):
    """
    Computes a parametrization of the vertex ring around each vertex.

    This function assigns a parametric coordinate (U, V) to each vertex in the 
    vertex ring, based on the angular position around the central vertex.

    The main logic involves:

    1. Iterating over each vertex and its adjacent vertices (vertex ring).

    2. Computing the angular position of each adjacent vertex.

    3. Assigning parametric coordinates (U, V) based on the angular position.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    v : numpy array
        The indices of the central vertices.
    vj : numpy array
        The indices of the adjacent vertices in the vertex ring.
    U : numpy array
        The U-coordinates of the parametrization.
    V : numpy array
        The V-coordinates of the parametrization.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The parametrization is based on the current mesh topology and geometry.

    See Also
    --------
    vertex_ring_vertices_iterators : Iterates over vertices in the vertex ring.
    """
    v, vj, l = self.vertex_ring_vertices_iterators(sort=True,
                                                return_lengths=True)
    index = np.arange(v.shape[0])
    step = np.zeros(v.shape[0])
    _, unique = np.unique(v, return_index=True)
    vertices = np.delete(v, unique)
    index = np.delete(index, unique)
    value = 0
    while len(vertices) > 0:
        value += 1
        _, unique = np.unique(vertices, return_index=True)
        step[index[unique]] = value
        vertices = np.delete(vertices, unique)
        index = np.delete(index, unique)
    phi = 2*np.pi*step / l[v]
    U = np.sin(phi)
    V = np.cos(phi)
    return v, vj, U, V

def vertex_local_frame(self):
    """
    Computes a local coordinate frame (origin, x, y, z) for each vertex.

    This function assigns a local frame to each vertex, where:
    - The origin is the vertex position.
    - The z-axis is aligned with the vertex normal.
    - The x and y axes are orthogonal to the z-axis and each other.

    The main logic involves:

    1. Computing the vertex normals.

    2. Finding orthogonal vectors to the vertex normals.

    3. Constructing the local frames.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    frame : Frame object
        A Frame object containing the local frames for each vertex.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The local frames are computed based on the current vertex positions and normals.

    See Also
    --------
    vertex_normals : Computes the normal vectors of all vertices.

    Frame : A class representing a local coordinate frame.
    """
    o = self.vertices
    z = self.vertex_normals()
    x = utilities.orthogonal_vectors(z)
    y = np.cross(z,x)
    frame = Frame(o, x, y, z)
    return frame

def edge_angle_vectors(self):
    """
    Computes the angle vectors for each edge in the mesh.

    This function calculates the angle vectors by considering the dihedral angle 
    between adjacent faces. The angle vector is proportional to the angle between 
    the face normals.

    The main logic involves:

    1. Iterating over each edge and its adjacent faces.
    
    2. Computing the cross product of the face normals.

    3. Normalizing the resulting vector to obtain the angle vector.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    angle_vectors : numpy array
        An array of angle vectors for each edge.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The angle vectors are computed based on the current face normals.

    See Also
    --------
    edge_sine_vectors : Computes the sine vectors for each edge.

    face_normals : Computes the normal vectors of all faces.
    """
    sin = self.edge_sine_vectors()
    beta = np.arcsin(np.linalg.norm(sin, axis=1))
    beta = np.array([beta,beta,beta]).T * utilities.normalize(sin)
    return beta

def edge_angles(self):
    """
    Computes the dihedral angles for each edge in the mesh.

    This function calculates the dihedral angle between adjacent faces for each edge.
    The dihedral angle is the angle between the face normals.

    The main logic involves:

    1. Iterating over each edge and its adjacent faces.
    
    2. Computing the dot product of the face normals.

    3. Calculating the dihedral angle using the arccosine function.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    angles : numpy array
        An array of dihedral angles for each edge.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The dihedral angles are computed based on the current face normals.

    See Also
    --------
    edge_angle_vectors : Computes the angle vectors for each edge.

    face_normals : Computes the normal vectors of all faces.
    """
    sin = self.edge_sine_vectors()
    beta = np.arcsin(np.linalg.norm(sin, axis=1))
    return beta

def edge_sine_vectors(self):
    """
    Computes the sine vectors for each edge in the mesh.

    This function calculates the sine vectors by considering the sine of the 
    dihedral angle between adjacent faces. The sine vector is proportional to 
    the cross product of the face normals.

    The main logic involves:

    1. Iterating over each edge and its adjacent faces.

    2. Computing the cross product of the face normals.

    3. Normalizing the resulting vector to obtain the sine vector.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    sine_vectors : numpy array
        An array of sine vectors for each edge.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The sine vectors are computed based on the current face normals.

    See Also
    --------
    edge_angle_vectors : Computes the angle vectors for each edge.

    face_normals : Computes the normal vectors of all faces.
    """
    v, ej = self.vertex_ring_edges_iterators(sort=True)
    normals = self.face_normals()
    v1, v2 = self.edge_vertices()
    f1, f2 = self.edge_faces()
    bf1 = np.where(f1 == -1)[0]
    bf2 = np.where(f2 == -1)[0]
    f1[bf1] = f2[bf1]
    f2[bf2] = f1[bf2]
    F1 = normals[f1,:]
    F2 = normals[f2,:]
    sin = np.cross(F1, F2)
    return sin

def extended_shape_operator(self, area_normalization=False, use_sine=False):
    """
    Computes the extended shape operator for each vertex in the mesh.

    This function calculates the extended shape operator, which is a matrix 
    representing the curvature properties of the mesh at each vertex.

    The main logic involves:

    1. Iterating over each vertex and its adjacent edges.

    2. Computing the curvature contributions from each edge.

    3. Summing the contributions to form the extended shape operator.

    Parameters
    ----------
    area_normalization : bool, optional (default=False)
        Whether to normalize the shape operator by the vertex area.
    use_sine : bool, optional (default=False)
        Whether to use sine vectors instead of angle vectors.

    Returns
    -------
    shape_operator : numpy array
        An array of extended shape operators for each vertex.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The shape operator is computed based on the current mesh geometry and topology.

    See Also
    --------
    principal_curvatures : Computes the principal curvatures from the shape operator.
    
    edge_angle_vectors : Computes the angle vectors for each edge.
    
    edge_sine_vectors : Computes the sine vectors for each edge.
    """
    v, ej = self.vertex_ring_edges_iterators(sort=True)
    v1, v2 = self.edge_vertices()
    if use_sine:
        B = self.edge_sine_vectors()
    else:
        B = self.edge_angle_vectors()
    V12 = self.vertices[v1[ej]] - self.vertices[v2[ej]]
    B12 = B[ej] / 2
    W = np.einsum('ij,ik -> ijk', B12, V12)
    W = utilities.sum_repeated(W, v)
    if area_normalization:
        A = self.vertex_ring_areas()
        A = np.array([[A,A,A],[A,A,A],[A,A,A]]).T
        W = W/A
    return W

def principal_curvatures(self, area_normalization=False, use_sine=False):
    """
    Computes the principal curvatures and directions for each vertex in the mesh.

    This function calculates the principal curvatures and their corresponding 
    directions by analyzing the extended shape operator at each vertex.

    The main logic involves:

    1. Computing the extended shape operator for each vertex.

    2. Diagonalizing the shape operator to obtain the principal curvatures and directions.

    Parameters
    ----------
    area_normalization : bool, optional (default=False)
        Whether to normalize the shape operator by the vertex area.
    use_sine : bool, optional (default=False)
        Whether to use sine vectors instead of angle vectors.

    Returns
    -------
    k1 : numpy array
        The first principal curvature at each vertex.
    k2 : numpy array
        The second principal curvature at each vertex.
    D1 : numpy array
        The direction of the first principal curvature.
    D2 : numpy array
        The direction of the second principal curvature.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The principal curvatures are computed based on the current mesh geometry and topology.

    See Also
    --------
    extended_shape_operator : Computes the extended shape operator for each vertex.
    
    gaussian_curvature : Computes the Gaussian curvature from principal curvatures.
    
    mean_curvature : Computes the mean curvature from principal curvatures.
    """
    W = self.extended_shape_operator(area_normalization, use_sine)
    try:
        eig = np.linalg.eigh(W)
        srt = np.argsort(np.abs(eig[0]), axis=1)
        i = np.arange(self.V)
        i1 = srt[:,1]
        i2 = srt[:,2]
        k1 = eig[0][i,i1]
        k2 = eig[0][i,i2]
        V1 = eig[1][i,:,i1]
        V2 = eig[1][i,:,i2]
        N = self.vertex_normals()
        D1 = utilities.normalize(np.cross(V1,N))
        D2 = utilities.normalize(np.cross(V2,N))
    except:
        V = self.V
        return (np.ones(V), np.ones(V), np.ones((V,3)), np.ones((V,3)))
    return (k1, k2, D1, D2)

def curvature_ratios(self):
    """
    Computes the curvature ratios for each vertex in the mesh.

    This function calculates the ratio of principal curvatures at each vertex, 
    which provides insight into the local shape characteristics (e.g., whether 
    the surface is more cylindrical or saddle-like).

    The main logic involves:

    1. Computing the principal curvatures using `principal_curvatures`.

    2. Calculating the ratio of the absolute values of the principal curvatures.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    ratios : numpy array
        The curvature ratios at each vertex.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The curvature ratios are computed based on the principal curvatures.

    See Also
    --------
    principal_curvatures : Computes the principal curvatures for each vertex.
    """
    K = self.principal_curvatures(True)
    k1 = K[0]
    k2 = K[1]
    R = np.maximum(np.abs(k2/(k1 + 1e-10)),
                            np.abs(k1/(k2 + 1e-10))) * np.sign(k1*k2)
    return R

def gaussian_curvature(self):
    """
    Computes the Gaussian curvature for each vertex in the mesh.

    This function calculates the Gaussian curvature as the product of the 
    principal curvatures at each vertex.

    The main logic involves:

    1. Computing the principal curvatures using `principal_curvatures`.

    2. Calculating the product of the principal curvatures.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    gaussian_curv : numpy array
        The Gaussian curvature at each vertex.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The Gaussian curvature is computed based on the principal curvatures.

    See Also
    --------
    principal_curvatures : Computes the principal curvatures for each vertex.
   
    mean_curvature : Computes the mean curvature for each vertex.
    """
    K = self.principal_curvatures(True, use_sine=True)
    return K[0]*K[1]

def mean_curvature(self):
    """
    Computes the mean curvature for each vertex in the mesh.

    This function calculates the mean curvature as the average of the 
    principal curvatures at each vertex.

    The main logic involves:

    1. Computing the principal curvatures using `principal_curvatures`.
    
    2. Calculating the average of the principal curvatures.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    mean_curv : numpy array
        The mean curvature at each vertex.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The mean curvature is computed based on the principal curvatures.

    See Also
    --------
    principal_curvatures : Computes the principal curvatures for each vertex.
    
    gaussian_curvature : Computes the Gaussian curvature for each vertex.
    """
    K = self.principal_curvatures(True)
    H =  0.5*(K[0] + K[1])
    return H

def edge_cotangents_weigths(self):
    """
    Computes the cotangent weights for each edge in the mesh.

    This function calculates the cotangent weights based on the angles 
    between adjacent edges. Cotangent weights are commonly used in 
    discrete differential geometry for various computations.

    The main logic involves:

    1. Iterating over each edge and its adjacent faces.

    2. Computing the cotangent of the angles between adjacent edges.

    3. Summing the cotangent weights for each edge.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    cotangent_weights : numpy array
        The cotangent weights for each edge.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The cotangent weights are computed based on the current mesh topology and geometry.

    See Also
    --------
    mean_curvature_normal : Computes the mean curvature normal using cotangent weights.
    """
    H = self.halfedges
    V0 = self.vertices[H[:,0]]
    V1 = self.vertices[H[H[:,2],0]]
    V2 = self.vertices[H[H[:,3],0]]
    E1 = (V0 - V2) / np.linalg.norm(V0 - V2, axis=1, keepdims=True)
    E2 = (V1 - V2) / np.linalg.norm(V1 - V2, axis=1, keepdims=True)
    cos = np.einsum('ij,ij->i', E1, E2)
    sin = np.linalg.norm(np.cross(E1, E2, axis=1), axis=1)
    cotan = cos / sin
    h = np.argsort(H[:,5])
    cotan = cotan[h]
    e = H[h,5]
    cotan = utilities.sum_repeated(cotan, e)
    return cotan / 2

def mean_curvature_normal(self):
    """
    Computes the mean curvature normal for each vertex in the mesh.

    This function calculates the mean curvature normal using cotangent weights 
    and vertex positions. The mean curvature normal is a vector field that 
    represents the direction and magnitude of the mean curvature.

    The main logic involves:

    1. Computing the cotangent weights using `edge_cotangents_weights`.

    2. Summing the weighted edge vectors to obtain the mean curvature normal.

    3. Normalizing the resulting vectors.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    mean_curvature_normals : numpy array
        The mean curvature normal vectors for each vertex.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The mean curvature normal is computed based on the current mesh geometry and topology.

    See Also
    --------
    edge_cotangents_weights : Computes the cotangent weights for each edge.

    mean_curvature : Computes the mean curvature for each vertex.
    """
    b = self.boundary_vertices()
    v, e = self.vertex_ring_edges_iterators(sort=True)
    v, vj = self.vertex_ring_vertices_iterators(sort=True)
    cotan = self.edge_cotangents_weigths()
    cotan = np.array([cotan, cotan, cotan]).T
    M = cotan[e] * (self.vertices[vj] - self.vertices[v])
    M = utilities.sum_repeated(M, v)
    M[b,:] = np.array([0,0,0])
    return M

# -------------------------------------------------------------------------
#                                  Cutting
# -------------------------------------------------------------------------

def cut(self, vertex_index):
    """
    Cuts the mesh along a boundary edge starting from a specified vertex.

    This function modifies the mesh topology by splitting a boundary edge and
    creating new vertices and half-edges. The cut is performed along the boundary
    until another boundary vertex is reached.

    The main logic involves:

    1. Identifying the boundary edge starting from the specified vertex.

    2. Iteratively splitting the boundary edge and updating the half-edge structure.
    
    3. Creating new vertices and half-edges to maintain mesh consistency.

    Parameters
    ----------
    vertex_index : int
        The index of the vertex where the cut starts.

    Returns
    -------
    Hcouple : numpy array
        An array of half-edge pairs representing the cut edges.

    Note
    ----
    This function assumes that the specified vertex is on the boundary of the mesh.
    The mesh must be manifold and orientable for the cut operation to be valid.

    See Also
    --------
    make_simply_connected : Simplifies the mesh topology by making it simply connected.
    """
    boundary = self.boundary_vertices()
    hc1 = []
    hc2 = []
    if not np.in1d(vertex_index, boundary)[0]:
        return False
    V = self.V
    E = self.E
    W = 2*E
    H = self.halfedges
    v, vj, l = self.vertex_ring_vertices_iterators(return_lengths=True)
    v0 = vertex_index
    vvv = [v0]
    if l[v0] < 3:
        return False
    h = np.where(H[:,0] == v0)[0]
    h = h[np.where(H[h,1] == -1)[0]]
    h0 = copy.copy(h)
    H[h0,0] = V
    n_v1 = np.copy(self.vertices[v0])
    for i in range(l[v0]//2):
        h = H[H[h,4],2]
        H[h,0] = V
    hc1.append(h); hc2.append(H[h,4])
    v1 = H[H[h,4],0][0]
    n_h1 = np.array([[V+1, -1, h0, W+2, h, E]], 'i')
    n_h2 = np.array([[v0, -1, W+3, H[h0,3], H[h,4], H[h,5]]], 'i')
    H[H[h,4],4] = W+1
    H[h,0] = V
    H[h,4] = W
    H[h,5] = E
    H[H[h0,3],2] = W+1
    H[h0,3] = W
    H = np.vstack((H, n_h1, n_h2))
    self._vertices = np.vstack((self.vertices, n_v1))
    W += 2
    E += 1
    V += 1
    while not np.in1d(v1, boundary)[0]:
        v0 = v1
        h = H[h,2]
        H[h,0] = V
        n_v1 = np.copy(self.vertices[v0])
        for i in range(int(l[v0]//2 - 1)):
            h = H[H[h,4],2]
            H[h,0] = V
        hc1.append(h); hc2.append(H[h,4])
        v1 = H[H[h,4],0][0]
        vvv.append(v0)
        n_h1 = np.array([[V+1, -1, W-2, W+2, h, E]], 'i')
        n_h2 = np.array([[v0, -1, W+3, W-1, H[h,4], H[h,5]]], 'i')
        H[H[h,4],4] = W+1
        H[h,0] = V
        H[h,4] = W
        H[h,5] = E
        H = np.vstack((H, n_h1, n_h2))
        self._vertices = np.vstack((self.vertices, n_v1))
        W += 2
        E += 1
        V += 1
    W -= 2
    n_v1 = np.copy(self.vertices[v0])
    while H[h,1] != -1:
        h = H[H[h,2],4]
        H[H[h,4],0] = V
    H[W+1,2] = copy.copy(H[h,2])
    H[H[h,2],3] = W+1
    H[h,2] = W
    H[W,3] = h
    self._halfedges = H
    n_v1 = np.copy(self.vertices[v1])
    self._vertices = np.vstack((self.vertices, n_v1))
    self.topology_update()
    Hcouple = np.array(hc2)
    '''
    from geometrylab import vtkplot
    pl_m = vtkplot.Edges(self)
    pl_p = vtkplot.Points(self.vertices[vvv])
    vtkplot.view([pl_m,pl_p])
    vtkplot.check(self)
    '''
    return Hcouple

# -------------------------------------------------------------------------
#                          Edit global connectivity
# -------------------------------------------------------------------------

def is_triangular_mesh(self):
    """
    Checks if the mesh is a triangular mesh.

    This function verifies whether all faces in the mesh are triangles.

    The main logic involves:

    1. Iterating over all faces.

    2. Checking the number of vertices per face.

    3. Returning True if all faces are triangles, otherwise False.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    is_triangular : bool
        Whether the mesh is a triangular mesh.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The check is based on the current mesh topology.

    See Also
    --------
    loop : Applies the Loop subdivision algorithm (requires a triangular mesh).
    """
    l = self.face_lengths()
    if len(np.where(l != 3)[0] > 0):
        return False
    else:
        return True

def loop(self, steps=1):
    """
    Applies the Loop subdivision algorithm to the mesh.

    This function refines the mesh by subdividing each triangle into four smaller triangles.
    The new vertices are positioned based on the Loop subdivision rules.

    The main logic involves:

    1. Identifying the vertices and edges of each triangle.

    2. Computing the new vertex positions using the Loop subdivision weights.

    3. Creating new faces and updating the mesh topology.

    Parameters
    ----------
    steps : int, optional (default=1)
        The number of subdivision steps to apply.

    Returns
    -------
    None
        The mesh is updated within the class instance.

    Note
    ----
    This function assumes that the mesh is triangular and manifold.
    The Loop subdivision algorithm is only applicable to triangular meshes.

    See Also
    --------
    catmull_clark : Applies the Catmull-Clark subdivision algorithm.
    """
    if not self.is_triangular_mesh():
        return False
    def _loop(self):
        V = self.V
        H = self.halfedges
        _, h1 = np.unique(H[:,1], True)
        h1 = np.delete(h1, np.where(H[h1,1] == -1))
        h2 = H[h1,2]
        h3 = H[h1,3]
        F0 = np.array((H[h1,5]+V, H[h2,5]+V, H[h3,5]+V)).T
        F1 = np.array((H[h1,0], H[h1,5]+V, H[H[h1,3],5]+V)).T
        F2 = np.array((H[h2,0], H[h2,5]+V, H[H[h2,3],5]+V)).T
        F3 = np.array((H[h3,0], H[h3,5]+V, H[H[h3,3],5]+V)).T
        new_faces = np.vstack((F0, F1, F2, F3)).tolist()
        v, vj, l = self.vertex_ring_vertices_iterators(sort=True,
                                                    return_lengths=True)
        c = 1./l * (5./8 - (3./8 + 1./4*np.cos(2*np.pi*l**(-1.)))**2)
        d = 1 - l*c
        c = np.array([c,c,c]).T
        d = np.array([d,d,d]).T
        ring = utilities.sum_repeated(self.vertices[vj], v)
        V0 = c*ring + d*self.vertices
        _, e = np.unique(H[:,5], True)
        v1 = self.vertices[H[e,0]]
        v2 = self.vertices[H[H[e,4],0]]
        v3 = self.vertices[H[H[e,3],0]]
        v4 = self.vertices[H[H[H[e,4],3],0]]
        be = self.boundary_edges()
        v3[be] = v1[be]
        v4[be] = v2[be]
        V1 = 3./8*v1 + 3./8*v2 + 1./8*v3 + 1./8*v4
        bh = np.where(H[:,1] == -1)[0]
        v0 = self.vertices[H[bh,0]]
        v5 = self.vertices[H[H[bh,3],0]]
        v6 = self.vertices[H[H[bh,2],0]]
        V0[H[bh,0]] = 1./8*v6 + 1./8*v5 + 3./4*v0
        new_vertices = np.vstack((V0, V1))
        self.make_mesh(new_vertices, new_faces)
    for i in range(steps):
        self._loop()

def catmull_clark(self, steps=1):
    """
    Applies the Catmull-Clark subdivision algorithm to the mesh.

    This function refines the mesh by subdividing each face into smaller faces 
    using the Catmull-Clark rules. It works for both triangular and quadrilateral meshes.

    The main logic involves:

    1. Iterating over each face and its vertices.

    2. Computing new vertex positions using the Catmull-Clark weights.

    3. Creating new faces and updating the mesh topology.

    Parameters
    ----------
    steps : int, optional (default=1)
        The number of subdivision steps to apply.

    Returns
    -------
    None
        The mesh is updated within the class instance.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The Catmull-Clark algorithm is applicable to both triangular and quadrilateral meshes.

    See Also
    --------
    loop : Applies the Loop subdivision algorithm (for triangular meshes).
    """
    def _catmull_clark(self):
        V = self.V
        E = self.E
        H = self.halfedges
        h = self.inner_halfedges()
        v1 = H[h,0]
        v2 = H[h,5] + V
        v3 = H[h,1] + V + E
        v4 = H[H[h,3],5] + V
        new_faces = (np.array([v1, v2, v3, v4]).T).tolist()
        v  = H[:,0]
        i  = np.argsort(v)
        v  = v[i]
        v1 = H[H[:,4],0][i]
        v2 = H[H[H[:,2],4],0][i]
        s = np.ones(len(v))
        l = utilities.sum_repeated(s,v)
        d = 1. / (4*l**2)
        c = 3. / (2*l**2)
        e = 1 - 7./(4*l)
        c = np.array([c,c,c]).T
        d = np.array([d,d,d]).T
        e = np.array([e,e,e]).T
        v1 = utilities.sum_repeated(self.vertices[v1], v)
        v2 = utilities.sum_repeated(self.vertices[v2], v)
        old_vertices = c*v1 + d*v2 + e*self.vertices
        _, e = np.unique(H[:,5], True)
        v1 = self.vertices[H[e,0]]
        v2 = self.vertices[H[H[e,4],0]]
        v3 = self.vertices[H[H[e,3],0]]
        v4 = self.vertices[H[H[H[e,4],3],0]]
        v5 = self.vertices[H[H[H[e,2],2],0]]
        v6 = self.vertices[H[H[H[H[e,4],2],2],0]]
        be = self.boundary_edges()
        v3[be] = v5[be] = v1[be]
        v4[be] = v6[be] = v2[be]
        mid_points = 3./8*(v1 + v2) + 1./16*(v3 + v4 + v5 + v6)
        bh = np.where(H[:,1] == -1)[0]
        v0 = self.vertices[H[bh,0]]
        v5 = self.vertices[H[H[bh,3],0]]
        v6 = self.vertices[H[H[bh,2],0]]
        old_vertices[H[bh,0]] = 1./8*v6 + 1./8*v5 + 3./4*v0
        barycenters = self.face_barycenters()
        new_vertices = np.vstack((old_vertices, mid_points, barycenters))
        self.make_mesh(new_vertices, new_faces)
    for i in range(steps):
        self._catmull_clark()


def dual_mesh(self, make_boundary=True):
    """
    Constructs the dual mesh of the current mesh.

    This function creates a new mesh where each face in the original mesh corresponds
    to a vertex in the dual mesh, and each vertex in the original mesh corresponds
    to a face in the dual mesh.

    The main logic involves:

    1. Computing the barycenters of the original faces to define the dual vertices.
    
    2. Creating new faces in the dual mesh based on the original vertex connectivity.
    
    3. Handling boundary cases if `make_boundary` is True.

    Parameters
    ----------
    make_boundary : bool, optional (default=True)
        Whether to create boundary faces in the dual mesh.

    Returns
    -------
    dual_mesh : Mesh
        The dual mesh instance.

    Note
    ----
    This function assumes that the input mesh is manifold and orientable.
    The dual mesh is constructed based on the current mesh topology and geometry.

    See Also
    --------
    exploded_mesh : Creates an exploded version of the mesh.
    """
    H = self.halfedges
    HD = np.copy(H)
    b = np.where(H[:,1] == -1)[0]
    B = b.shape[0]
    fb = np.arange(B) + self.F
    hb1 = np.arange(B) + 2*self.E
    hb2 = np.arange(B) + B + 2*self.E
    eb = np.arange(B) + self.E
    HD[b,1] = np.copy(fb)
    HD[b,3] = np.copy(hb1)
    HD[H[b,3],2] = np.copy(hb2)
    Hb1 = np.zeros((B,6), dtype=int)
    Hb2 = np.zeros((B,6), dtype=int)
    Hb1[:,0] = -1
    Hb2[:,0] = H[b,0]
    Hb1[:,1] = fb
    Hb2[:,1] = HD[H[b,3],1]
    Hb1[:,2] = b
    Hb2[:,2] = HD[H[b,3],3]#HD[H[b,2],2]##
    Hb1[:,3] = HD[b,2]
    Hb2[:,3] = H[b,3]
    Hb1[:,4] = hb2
    Hb2[:,4] = hb1
    Hb1[:,5] = eb
    Hb2[:,5] = eb
    HD = np.vstack((HD, Hb1, Hb2))
    HR = np.copy(HD)
    HD[:,0] = HR[HR[:,4],1]
    HD[:,1] = HR[:,0]
    HD[:,2] = HR[HR[:,3],4]
    HD[:,3] = HR[HR[:,4],2]
    dual_vertices = self.face_barycenters()
    face_normals = self.face_normals()
    edge_vec = self.vertices[H[b,0]] - self.vertices[H[H[b,2],0]]
    normals = np.cross(edge_vec, face_normals[H[H[b,4],1]])
    new_vertices = self.edge_mid_points()[H[b,5]] + normals/2
    dual_vertices = np.vstack((dual_vertices, new_vertices))
    dual_mesh = Mesh()
    dual_mesh._halfedges = HD
    dual_mesh._vertices = dual_vertices
    dual_mesh.topology_update()
    return dual_mesh

def delete_faces(self, faces):
    """
    Deletes specified faces from the mesh.

    This function removes the specified faces and updates the mesh topology.
    It also removes any half-edges and vertices that are no longer connected.

    The main logic involves:

    1. Identifying the half-edges and vertices associated with the faces to be deleted.
    
    2. Removing the faces and updating the half-edge structure.
    
    3. Cleaning up any unconnected vertices and edges.

    Parameters
    ----------
    faces : int or list of int
        The indices of the faces to be deleted.

    Returns
    -------
    None
        The mesh is updated within the class instance.

    Note
    ----
    This function assumes that the specified faces exist in the mesh.
    The mesh topology is updated to maintain consistency after deletion.

    See Also
    --------
    delete_unconnected_vertices : Removes unconnected vertices from the mesh.
    """
    if type(faces) is int:
        faces = [faces]
    H = self.halfedges
    self._open_trash()
    self._cull_faces(faces)
    hf = np.arange(H.shape[0])[np.in1d(H[:,1], faces)]
    bh = hf[H[H[hf,4],1] == -1]
    self._cull_halfedges(bh)
    self._cull_halfedges(H[bh,4])
    self._cull_edges(H[bh,5])
    dh = hf[np.in1d(H[hf,4], hf)]
    self._cull_halfedges(dh)
    self._cull_halfedges(H[dh,4])
    self._cull_edges(H[dh,5])
    H[hf,1] = -1
    self._clean()
    self.delete_unconnected_vertices()

def delete_unconnected_vertices(self):
    """
    Removes unconnected vertices from the mesh.

    This function identifies and deletes vertices that are not part of any face.
    It updates the mesh topology accordingly.

    The main logic involves:

    1. Identifying vertices that are not referenced by any half-edge.

    2. Removing these vertices from the mesh.

    3. Updating the vertex indices in the half-edge structure.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    None
        The mesh is updated within the class instance.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The mesh topology is updated to maintain consistency after deletion.

    See Also
    --------
    delete_faces : Deletes specified faces from the mesh.
    """
    H = self.halfedges
    v = np.arange(len(self.vertices))
    cull = np.invert(np.in1d(v, H[:,0]))
    self._open_trash()
    self._cull_vertices(v[cull])
    self._clean()

def exploded_mesh(self):
    """
    Creates an exploded version of the mesh.

    This function generates a new mesh where each face is separated from the others,
    effectively "exploding" the mesh into individual faces.

    The main logic involves:

    1. Duplicating the vertices for each face.

    2. Creating new faces using the duplicated vertices.

    3. Constructing the exploded mesh topology.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    exploded_mesh : Mesh
        The exploded mesh instance.

    Note
    ----
    This function assumes that the input mesh is manifold and orientable.
    The exploded mesh is constructed based on the current mesh topology and geometry.

    See Also
    --------
    dual_mesh : Constructs the dual mesh of the current mesh.
    """
    f, v = self.face_vertices_iterators()
    vertices = self.vertices[v]
    k = np.arange(len(v))
    faces_list = [[] for i in range(self.F)]
    for i in range(len(f)):
        faces_list[f[i]].append(k[i])
    exp_mesh = Mesh()
    exp_mesh.make_mesh(vertices, faces_list)
    return exp_mesh

# -------------------------------------------------------------------------
#                            Edit edge connectivity
# -------------------------------------------------------------------------

def delete_edge(self, edge_index):
    """
    Deletes a specified edge from the mesh.

    This function removes the specified edge and updates the mesh topology.
    It also removes any half-edges and vertices that are no longer connected.

    The main logic involves:

    1. Identifying the half-edges associated with the edge to be deleted.

    2. Removing the edge and updating the half-edge structure.

    3. Cleaning up any unconnected vertices and edges.

    Parameters
    ----------
    edge_index : int
        The index of the edge to be deleted.

    Returns
    -------
    None
        The mesh is updated within the class instance.

    Note
    ----
    This function assumes that the specified edge exists in the mesh.
    The mesh topology is updated to maintain consistency after deletion.

    See Also
    --------
    delete_faces : Deletes specified faces from the mesh.

    delete_unconnected_vertices : Removes unconnected vertices from the mesh.
    """
    h = self.edge_halfedge(edge_index)
    self._open_trash()
    self._delete_halfedge(h)
    self._clean()

def flip_edge(self, edge_index):
    """
    Flips the specified edge in the mesh.

    This function reverses the direction of the specified edge by swapping its 
    adjacent faces. This operation is useful for optimizing mesh quality.

    The main logic involves:

    1. Identifying the half-edges associated with the edge.

    2. Swapping the adjacent faces of the edge.

    3. Updating the half-edge structure.

    Parameters
    ----------
    edge_index : int
        The index of the edge to be flipped.

    Returns
    -------
    success : bool
        Whether the edge flip was successful.

    Note
    ----
    This function assumes that the specified edge exists in the mesh.
    The edge flip is only performed if it results in a valid mesh topology.

    See Also
    --------
    split_edge : Splits an edge into two edges.

    collapse_edge : Collapses an edge into a single vertex.
    """
    h = self.edge_halfedge(edge_index)
    if not self.is_halfedge_bounding_tri_faces(h):
        return False
    self._flip_halfedge(h)
    self.topology_update()

def split_edge(self, edge_index):
    """
    Splits the specified edge into two edges by inserting a new vertex.

    This function subdivides the specified edge by adding a new vertex at its midpoint.
    The new vertex is connected to the original vertices of the edge.

    The main logic involves:

    1. Identifying the half-edges associated with the edge.

    2. Inserting a new vertex at the midpoint of the edge.

    3. Creating new half-edges and updating the mesh topology.

    Parameters
    ----------
    edge_index : int
        The index of the edge to be split.

    Returns
    -------
    success : bool
        Whether the edge split was successful.

    Note
    ----
    This function assumes that the specified edge exists in the mesh.
    The edge split is only performed if it results in a valid mesh topology.

    See Also
    --------
    flip_edge : Flips an edge by swapping its adjacent faces.

    collapse_edge : Collapses an edge into a single vertex.
    """
    h = self.edge_halfedge(edge_index)
    if not self.is_halfedge_bounding_tri_faces(h):
        return False
    self._expand_arrays()
    self._split_halfedge(h)
    self._cull_arrays()

def collapse_edge(self, edge_index):
    """
    Collapses the specified edge into a single vertex.

    This function merges the two vertices of the specified edge into a single vertex.
    The resulting vertex replaces the original vertices in the mesh topology.

    The main logic involves:

    1. Identifying the half-edges associated with the edge.

    2. Merging the two vertices of the edge into a single vertex.

    3. Updating the half-edge structure to maintain mesh consistency.

    Parameters
    ----------
    edge_index : int
        The index of the edge to be collapsed.

    Returns
    -------
    success : bool
        Whether the edge collapse was successful.

    Note
    ----
    This function assumes that the specified edge exists in the mesh.
    The edge collapse is only performed if it results in a valid mesh topology.

    See Also
    --------
    split_edge : Splits an edge into two edges.

    flip_edge : Flips an edge by swapping its adjacent faces.
    """
    h = self.edge_halfedge(edge_index)
    if not self.is_halfedge_bounding_tri_faces(h):
        return False
    self._open_trash()
    self._collapse_halfedge(h)
    self._clean()

# -------------------------------------------------------------------------
#                                   Remesh
# -------------------------------------------------------------------------

def equalize_valences(self):
    """
    Equalizes the valences of vertices in a triangular mesh.

    This function attempts to balance the number of edges connected to each vertex
    by performing local mesh operations such as edge flips, splits, and collapses.
    The goal is to achieve a more uniform vertex valence distribution.

    The main logic involves:

    1. Iterating over all edges and checking their valences.

    2. Performing edge flips, splits, or collapses to balance valences.

    3. Repeating the process until the desired valence distribution is achieved.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    success : bool
        Whether the valence equalization was successful.

    Note
    ----
    This function assumes that the mesh is a triangular mesh.
    The mesh topology is updated to achieve a more uniform valence distribution.

    See Also
    --------
    split_edges : Splits edges to achieve a desired maximum length.

    collapse_edges : Collapses edges to achieve a desired minimum length.
    """
    if not self.is_triangular_mesh():
        return False
    H = self.halfedges
    _, he = np.unique(H[:,5], return_index=True)
    t = np.repeat(6, self.V)
    t[self.boundary_vertices()] = 4
    t[self.mesh_corners()] = 3
    _, _, l = self.vertex_ring_vertices_iterators(True, False, True)
    for h in he:
        a = H[h,0]
        b = H[H[h,4],0]
        c = H[H[h,3],0]
        d = H[H[H[h,4],3],0]
        deviation_0  = (l[a] - t[H[h,0]])**2
        deviation_0 += (l[b] - t[H[H[h,4],0]])**2
        deviation_0 += (l[c] - t[H[H[h,3],0]])**2
        deviation_0 += (l[d] - t[H[H[H[h,4],3],0]])**2
        deviation_1  = (l[a] - t[H[h,0]] - 1)**2
        deviation_1 += (l[b] - t[H[H[h,4],0]] - 1)**2
        deviation_1 += (l[c] - t[H[H[h,3],0]] + 1)**2
        deviation_1 += (l[d] - t[H[H[H[h,4],3],0]] + 1)**2
        if deviation_1 < deviation_0:
            if self._flip_halfedge(h):
                l[a] -= 1
                l[b] -= 1
                l[c] += 1
                l[d] += 1
    self.topology_update()
    return True

def split_edges(self, max_length):
    """
    Splits edges in the mesh that exceed a specified maximum length.

    This function iterates over all edges and splits those that are longer than 
    the specified maximum length. New vertices are inserted at the midpoints of 
    these edges, and the mesh topology is updated accordingly.

    The main logic involves:

    1. Iterating over all edges and checking their lengths.

    2. Splitting edges that exceed the maximum length.

    3. Updating the mesh topology with new vertices and edges.

    Parameters
    ----------
    max_length : float
        The maximum allowed length for edges.

    Returns
    -------
    success : bool
        Whether the edge splitting was successful.

    Note
    ----
    This function assumes that the mesh is a triangular mesh.
    The mesh topology is updated to ensure all edges are within the specified length.

    See Also
    --------
    collapse_edges : Collapses edges that are shorter than a specified minimum length.
    
    equalize_valences : Balances vertex valences in the mesh.
    """
    if not self.is_triangular_mesh():
        return False
    H = self.halfedges
    _, he = np.unique(H[:,5], return_index=True)
    self._expand_arrays(len(he))
    for h in he:
        if self.halfedge_length(h) > max_length:
            self._split_halfedge(h)
    self._cull_arrays()
    return True

def collapse_edges(self, min_length):
    """
    Collapses edges in the mesh that are shorter than a specified minimum length.

    This function iterates over all edges and collapses those that are shorter than 
    the specified minimum length. The two vertices of each collapsed edge are merged 
    into a single vertex, and the mesh topology is updated accordingly.

    The main logic involves:

    1. Iterating over all edges and checking their lengths.

    2. Collapsing edges that are shorter than the minimum length.

    3. Updating the mesh topology with merged vertices.

    Parameters
    ----------
    min_length : float
        The minimum allowed length for edges.

    Returns
    -------
    success : bool
        Whether the edge collapse was successful.

    Note
    ----
    This function assumes that the mesh is a triangular mesh.
    The mesh topology is updated to ensure all edges are longer than the specified length.

    See Also
    --------
    split_edges : Splits edges that exceed a specified maximum length.

    equalize_valences : Balances vertex valences in the mesh.
    """
    if not self.is_triangular_mesh():
        return False
    H = self.halfedges
    _, he = np.unique(H[:,5], return_index=True)
    self._open_trash()
    for h in he:
        if self._is_culled_halfedge(h):
            h = H[h,4]
        if self.halfedge_length(h) < min_length:
            self._collapse_halfedge(h)
    self._clean()
    return True

# -------------------------------------------------------------------------
#                             Local connectivity
# -------------------------------------------------------------------------

def edge_halfedge(self, edge_index):
    """
    Retrieves a half-edge associated with a specified edge.

    This function returns one of the half-edges that corresponds to the given edge index.
    Each edge in the half-edge mesh is represented by a pair of half-edges.

    Parameters
    ----------
    edge_index : int
        The index of the edge.

    Returns
    -------
    halfedge_index : int
        The index of a half-edge associated with the specified edge.

    Note
    ----
    This function assumes that the edge index is valid and exists in the mesh.
    The returned half-edge is one of the two half-edges that form the edge.

    See Also
    --------
    vertex_halfedge : Retrieves a half-edge associated with a specified vertex.
    """
    H = self.halfedges
    h = np.where(H[:,5] == edge_index)[0][0]
    return h

def vertex_halfedge(self, vertex_index):
    """
    Retrieves a half-edge originating from a specified vertex.

    This function returns one of the half-edges that starts at the specified vertex.
    The returned half-edge can be used to traverse the mesh topology.

    The main logic involves:

    1. Searching for a half-edge that originates from the specified vertex.

    2. Returning the index of the found half-edge.

    Parameters
    ----------
    vertex_index : int
        The index of the vertex.

    Returns
    -------
    halfedge_index : int
        The index of a half-edge originating from the specified vertex.

    Note
    ----
    This function assumes that the specified vertex exists in the mesh.
    The returned half-edge can be used for further mesh traversal and manipulation.

    See Also
    --------
    halfedge_ring : Retrieves the half-edge ring around a specified half-edge.
    
    vertex_ring_vertices : Retrieves the vertices in the vertex ring around a specified vertex.
    """
    H = self.halfedges
    v = np.where(H[:,0] == vertex_index)[0][0]
    return v

def halfedge_ring(self, halfedge_index):
    """
    Retrieves the half-edge ring around a specified half-edge.

    This function returns the sequence of half-edges that form a loop around the 
    specified half-edge. The half-edge ring can be used to traverse the mesh topology.

    The main logic involves:
    
    1. Starting from the specified half-edge.

    2. Iteratively following the next half-edge in the loop until the starting half-edge is reached.
    
    3. Returning the sequence of half-edges in the loop.

    Parameters
    ----------
    halfedge_index : int
        The index of the starting half-edge.

    Returns
    -------
    ring : list of int
        The indices of half-edges in the half-edge ring.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The half-edge ring can be used for further mesh traversal and manipulation.

    See Also
    --------
    vertex_halfedge : Retrieves a half-edge originating from a specified vertex.
    
    halfedge_ring_vertices : Retrieves the vertices in the half-edge ring.
    """
    H = self.halfedges
    h0 = halfedge_index
    ring = [h0]
    h = H[H[h0,3],4]
    while h != h0:
        ring.append(h)
        h = H[H[h,3],4]
    return ring

def vertex_ring_vertices(self, vertex_index):
    """
    Retrieves the vertices in the vertex ring around a specified vertex.

    This function returns the sequence of vertices that are directly connected to 
    the specified vertex. The vertex ring can be used to traverse the mesh topology.

    The main logic involves:

    1. Starting from a half-edge originating from the specified vertex.
    
    2. Iteratively following the next half-edge in the loop to collect adjacent vertices.
    
    3. Returning the sequence of vertices in the vertex ring.

    Parameters
    ----------
    vertex_index : int
        The index of the vertex.

    Returns
    -------
    ring_vertices : list of int
        The indices of vertices in the vertex ring.

    Note
    ----
    This function assumes that the specified vertex exists in the mesh.
    The vertex ring can be used for further mesh traversal and manipulation.

    See Also
    --------
    vertex_halfedge : Retrieves a half-edge originating from a specified vertex.
    
    halfedge_ring : Retrieves the half-edge ring around a specified half-edge.
    """
    h = self.vertex_halfedge(vertex_index)
    ring = self.halfedge_ring_vertices(h)
    return ring

def vertex_multiple_ring_vertices(self, vertex_index, depth=1):
    """
    Retrieves the vertices within a specified depth from a given vertex.

    This function performs a breadth-first traversal starting from the specified 
    vertex and collects all vertices within the specified depth. The result 
    includes the starting vertex and all vertices reachable within the given depth.

    The main logic involves:
    
    1. Initializing a queue with the starting vertex.
    
    2. Iteratively expanding the vertex ring up to the specified depth.
   
    3. Collecting all visited vertices.

    Parameters
    ----------
    vertex_index : int
        The index of the starting vertex.
    depth : int, optional (default=1)
        The maximum depth of traversal.

    Returns
    -------
    ring_vertices : numpy array
        The indices of vertices within the specified depth.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The traversal is based on the current mesh topology.

    See Also
    --------
    vertex_ring_vertices : Retrieves the vertices in the immediate vertex ring.
    
    ertex_ring_expansion : Expands the vertex ring with optional callbacks.
    """
    vi, vj = self.vertex_ring_vertices_iterators()
    ring = np.array([], dtype='i')
    search = np.array([vertex_index], dtype='i')
    for i in range(int(depth)):
        vring = np.array([], dtype='i')
        for v in search:
            vring = np.hstack((vj[vi == v], vring))
        vring = np.unique(vring)
        vring = vring[np.invert(np.in1d(vring, ring))]
        search = vring
        ring = np.hstack((ring, vring))
        if len(ring) == self.V:
            return ring
    return np.unique(ring)

def halfedge_ring_vertices(self, halfedge_index):
    """
    Retrieves the vertices in the half-edge ring around a specified half-edge.

    This function returns the sequence of vertices that form the half-edge ring 
    starting from the specified half-edge. The vertices are returned in the order 
    they appear around the ring.

    Parameters
    ----------
    halfedge_index : int
        The index of the starting half-edge.

    Returns
    -------
    vertices : list of int
        The indices of vertices in the half-edge ring.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The half-edge ring vertices are useful for local mesh analysis and processing.

    See Also
    --------
    halfedge_ring : Retrieves the half-edge ring around a specified half-edge.
    
    vertex_ring_vertices_iterators : Iterates over vertices in the vertex ring.
    """
    H = self.halfedges
    ring = self.halfedge_ring(halfedge_index)
    vertices = H[H[ring,2],0]
    return vertices

def halfedge_ring_faces(self, halfedge_index):
    """
    Retrieves the faces adjacent to the half-edge ring of a specified half-edge.

    This function returns the sequence of faces that are adjacent to the half-edge 
    ring starting from the specified half-edge. The faces are returned in the order 
    they appear around the ring.

    The main logic involves:
    
    1. Starting from the specified half-edge.
    
    2. Iteratively following the next half-edge in the ring to collect adjacent faces.
    
    3. Returning the sequence of faces.

    Parameters
    ----------
    halfedge_index : int
        The index of the starting half-edge.

    Returns
    -------
    ring_faces : list of int
        The indices of faces adjacent to the half-edge ring.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The faces are returned in the order they appear around the ring.

    See Also
    --------
    halfedge_ring : Retrieves the half-edge ring around a specified half-edge.
    
    face_vertices_iterators : Iterates over vertices of all faces.
    """
    H = self.halfedges
    ring = self.halfedge_ring(halfedge_index)
    faces = H[H[ring,2],1]
    return faces

def halfedge_face_vertices(self, halfedge_index):
    """
    Retrieves the vertices of the face associated with a specified half-edge.

    This function returns the sequence of vertices that form the face associated 
    with the specified half-edge. The vertices are returned in the order they 
    appear around the face.

    The main logic involves:

    1. Starting from the specified half-edge.

    2. Iteratively following the next half-edge in the face loop to collect vertices.
    
    3. Returning the sequence of vertices.

    Parameters
    ----------
    halfedge_index : int
        The index of the half-edge.

    Returns
    -------
    face_vertices : list of int
        The indices of vertices forming the face.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The vertices are returned in the order they appear around the face.

    See Also
    --------
    halfedge_ring_vertices : Retrieves the vertices in the half-edge ring.
    
    face_vertices_iterators : Iterates over vertices of all faces.
    """
    H = self.halfedges
    ring = self.halfedge_face_ring(halfedge_index)
    vertices = H[ring,0]
    return vertices

# -------------------------------------------------------------------------
#                             Local queries
# -------------------------------------------------------------------------

def is_boundary_halfedge_ring(self, ring):
    """
    Checks if a half-edge ring is on the boundary of the mesh.

    This function verifies whether the specified half-edge ring is part of the 
    mesh boundary. It returns True if the ring is on the boundary, False otherwise.

    The main logic involves:

    1. Iterating over the half-edge ring.

    2. Checking if any half-edge in the ring is a boundary half-edge.

    Parameters
    ----------
    ring : list of int
        The indices of half-edges in the ring.

    Returns
    -------
    is_boundary : bool
        Whether the half-edge ring is on the boundary.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The check is based on the current mesh topology and half-edge definitions.

    See Also
    --------
    is_halfedge_bounding_tri_faces : Checks if a half-edge bounds triangular faces.
    """
    H = self.halfedges
    for h in ring:
        if H[h,1] == -1:
            v0 = H[h,0]
            v1 = H[H[h,2],0]
            v2 = H[H[h,3],0]
            E1 = (self.vertices[v1] - self.vertices[v0])
            E2 = (self.vertices[v0] - self.vertices[v2])
            E1 = E1 / (E1[0]**2 + E1[1]**2 + E1[2]**2 + 1e-10)**0.5
            E2 = E2 / (E2[0]**2 + E2[1]**2 + E2[2]**2 + 1e-10)**0.5
            dot = E1[0]*E2[0] + E1[1]*E2[1] + E1[2]*E2[2]
            if dot < self.corner_tolerance:
                return 2
            else:
                return 1
    return 0

def is_halfedge_bounding_tri_faces(self, halfedge_index):
    """
    Checks if a specified half-edge bounds triangular faces.

    This function verifies whether the specified half-edge is part of a triangular 
    face or a boundary edge. It returns True if the half-edge bounds triangular faces, 
    False otherwise.

    The main logic involves:

    1. Checking the face associated with the half-edge.

    2. Verifying if the face is triangular or a boundary face.

    Parameters
    ----------
    halfedge_index : int
        The index of the half-edge.

    Returns
    -------
    is_bounding : bool
        Whether the half-edge bounds triangular faces.

    Note
    ----
    This function assumes that the mesh is manifold and orientable.
    The check is based on the current mesh topology and face definitions.

    See Also
    --------
    is_boundary_halfedge_ring : Checks if a half-edge ring is on the boundary.
    """
    H = self.halfedges
    h = halfedge_index
    for i in range(2):
        counter = 1
        h0 = h
        h = H[h,2]
        while h != h0:
            h = H[h,2]
            counter += 1
            if counter > 3:
                return False
        h = H[halfedge_index,4]
    return True

# -------------------------------------------------------------------------
#                             Local geometry
# -------------------------------------------------------------------------

def halfedge_length(self, halfedge_index):
    """
    Computes the length of a specified half-edge.

    This function calculates the Euclidean distance between the origin and destination 
    vertices of the specified half-edge.

    The main logic involves:

    1. Retrieving the origin and destination vertices of the half-edge.
    
    2. Computing the Euclidean distance between the two vertices.

    Parameters
    ----------
    halfedge_index : int
        The index of the half-edge.

    Returns
    -------
    length : float
        The length of the specified half-edge.

    Note
    ----
    This function assumes that the specified half-edge exists in the mesh.
    The length is computed based on the current vertex positions.

    See Also
    --------
    edge_lengths : Computes the lengths of all edges in the mesh.
    """
    H = self.halfedges
    h = halfedge_index
    V1 = self.vertices[H[h,0]]
    V2 = self.vertices[H[H[h,4],0]]
    E = V1 - V2
    return (E[0]**2 + E[1]**2 + E[2]**2)**0.5

# -------------------------------------------------------------------------
#                         Edit half-edge connectivity
# -------------------------------------------------------------------------