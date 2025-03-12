"""
forked and built on geometrylab/geometry/meshpy.py

meshpy.py --> quadrings.py --> gridshell_new.py(GeolabGUI)--> gui_basic.py -->
guidedprojection_net.py + opt_din.py --> read_file
"""

__author__ = 'Hui Wang'


#------------------------------------------------------------------------------
import numpy as np
#------------------------------------------------------------------------------
from geometrylab.geometry.meshpy import Mesh

from archgeolab.archgeometry.orient import  orient_rings

from archgeolab.archgeometry.conicSection import interpolate_sphere
    
from archgeolab.archgeometry.curves import mesh_polylines,make_polyline_from_endpoints,\
        make_multiple_polylines_from_endpoints,\
        get_isoline_between_2bdry,get_diagonal_polyline_from_2points
        
from archgeolab.archgeometry.make_mesh import make_quad_mesh_pieces,get_barycenters_mesh
#------------------------------------------------------------------------------

# -------------------------------------------------------------------------
#                        Supplementary functions
# -------------------------------------------------------------------------

def regular_rectangle_patch(self):
    """
    Creates a patch matrix for a regular rectangular mesh.

    This function constructs a matrix representing the connectivity of a regular 
    rectangular mesh, where each entry corresponds to a vertex index.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    None
        The patch matrix is stored within the class instance.

    Note
    ----
    This function assumes that the mesh is a regular rectangular patch.
    The patch matrix is useful for organizing and accessing mesh vertices.

    See Also
    --------
    regular_rotational_patch : Creates a patch matrix for a rotational mesh.
    """
    H = self.halfedges
    b = np.where(H[:,1]==-1)[0]
    c = np.where(H[H[H[H[b,4],3],4],1]==-1)[0]
    corner = H[H[b[c],2],0]
    e0 = b[c[0]] # choose the first corner-edge
    def row(e): # vertices' index --- up--down
        row = []
        while H[e,0] not in corner:
            row.append(e)
            e = H[e,3]
        row.append(e)
        return row

    def column(e): # vertices' index --- left--right
        col = H[e,0]
        ecol = [e]
        e = H[H[H[e,4],2],2]
        while H[e,4] not in b:
            col = np.r_[col,H[e,0]]
            ecol.append(e)
            e = H[H[H[e,4],2],2]
        col = np.r_[col,H[e,0]]
        ecol.append(e)
        return col,ecol
    
    left = row(e0)
    v1,e1 = column(e0)
    v0 = H[H[e1,4],0]
    for e in left:  
        vi,_ = column(e)
        v0 = np.c_[v0,vi]
        
    v,v1,v2,v3,v4 = self.ver_star_matrix[0,:]
    if np.where(v0==v)[1] != np.where(v0==v1)[1]:
        self._patch_matrix = v0.T
    else:
        "v1,v,v3 from up to down; v2,v,v4 from left to right"
        self._patch_matrix = v0 
        
def regular_rotational_patch(self):
    """
    Creates a patch matrix for a regular rotational mesh.

    This function constructs a matrix representing the connectivity of a regular 
    rotational mesh, where each entry corresponds to a vertex index.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    None
        The patch matrix is stored within the class instance.

    Note
    ----
    This function assumes that the mesh is a regular rotational patch.
    The patch matrix is useful for organizing and accessing mesh vertices.

    See Also
    --------
    regular_rectangle_patch : Creates a patch matrix for a rectangular mesh.
    """
    H = self.halfedges
    eb1 = np.where(H[:,1]==-1)[0][0]
    erow1 = [eb1] # index for edges
    eb = eb1
    while H[H[eb,2],0] != H[eb1,0]:
        eb = H[eb,2]
        erow1.append(eb)
    erow1.append(H[eb,2]) # last should be equal to the first
    vm = np.array([],dtype=int)
    for e in erow1:
        eci = [e] # index for edge
        while H[H[e,4],1]!=-1:
            e = H[H[H[e,4],2],2]
            eci.append(e)
        vm = np.r_[vm, H[eci,0]]
    vM = (vm.reshape(len(erow1),len(eci))).T
    self._rot_patch_matrix = vM 


def get_patch_and_rot_matrix_ind(self):
    """
    Retrieves the indices of vertices in the patch and rotational matrices.

    This function identifies the vertices in the patch and rotational matrices 
    and returns their indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    indices : list
        A list containing the indices of vertices in the patch and rotational matrices.

    Note
    ----
    This function assumes that the mesh is a regular patch or rotational mesh.
    The indices are useful for accessing specific vertices in the mesh.

    only work for patch/annulus srf
    from mesh rrvertex, return [v1-v-v3] in column direction of indM
                                [v2-v-v4] in row direction of indM
    indM_inner = [[0,1,2,3,4,5],
                    [6,7,8,9,10,11]
                    [......]
                    []]
    [v,v1,v2,v3,v4] ~ indM_inn
    [v[indv],v1[indv],v2[indv],v3[indv],v4[indv]] ~ indM_inner
    indv: indices of M_inn in rrv4f4
    [ivf1,ivf2,ivf3,ivf4]: quadfaces indices of M_inn_inn in rrv4f4


    See Also
    --------
    regular_rectangle_patch : Creates a patch matrix for a rectangular mesh.
    regular_rotational_patch : Creates a patch matrix for a rotational mesh.
    """
    v,v1,v2,v3,v4 = self.rrv4f4

    is_patch_or_rotation = True
    "shared by project self.mesh, pleated, aotu"
    if len(self.corner)!=0:
        "[0,1,2,3,4,5]"
        M = self.patch_matrix
        M_inn = M[1:-1,1:-1]
    else:
        "[0,1,2,3,4,5; 0]"
        M = self.rot_patch_matrix ##M[:,-1]==M[:,0]
        M_inn = M[1:-1,:]
        is_patch_or_rotation = False

    self._is_patch_or_rot = is_patch_or_rotation
    
    i,j = np.where(M==v[0])
    "note i,j may equal to [1 1] [ 0 80]"
    i,j = i[0], j[0]
    if M[i-1,j]==v1[0] or M[i-1,j]==v3[0]:
        "choose the vertical of M being [v1-v-v3], then horizon is [v2-v-v4]"
        pass
    elif M[i-1,j]==v2[0] or M[i-1,j]==v4[0]:
        M_inn = M_inn.T
        M = M.T
    else:
        print("there is a bug of the corresponding between rrvstar and indM")

    self._vMatrix = M

    indv = []
    for iv in M_inn.flatten():
        j = np.where(v==iv)[0][0]
        indv.append(j)
    indv = np.array(indv)
    #Vc = V[v[indv]] ##note will not update later optimization
    #V1,V2,V3,V4 = V[v1[indv]],V[v2[indv]],V[v3[indv]],V[v4[indv]]
    
    indvM = indv.reshape(M_inn.shape) ##~M_inn
    ivf1 = indvM[:-1,:-1].flatten()
    ivf2 = indvM[1:,:-1].flatten()
    ivf3 = indvM[1:,1:].flatten()
    ivf4 = indvM[:-1,1:].flatten() 
    vf1,vf2,vf3,vf4 = v[ivf1],v[ivf2],v[ivf3],v[ivf4]

    indM = np.arange(M.shape[0]*M.shape[1]).reshape(-1,M.shape[1])
    indM_inn = np.arange(M_inn.shape[0]*M_inn.shape[1]).reshape(-1,M_inn.shape[1])

    f1 = M[:-1,:-1].flatten()
    f2 = M[1:,:-1].flatten()
    f3 = M[1:,1:].flatten()
    f4 = M[:-1,1:].flatten()  
    self._v3D_f1234 = [f1,f2,f3,f4]
    
    if1 = indM_inn[:-1,:-1].flatten()
    if2 = indM_inn[1:,:-1].flatten()
    if3 = indM_inn[1:,1:].flatten()
    if4 = indM_inn[:-1,1:].flatten()  
    self._vind_f1234 = [len(ivf1),[vf1,vf2,vf3,vf4],[if1,if2,if3,if4]]##~rrv4f4~binormal

    num_row,num_col = M.shape
    Vnum_2D = num_row * num_col
    self._Vnum_2D = Vnum_2D

    if len(self.corner)!=0:
        arr = np.arange(Vnum_2D)
        M2D = arr.reshape(num_row,-1)
    else:
        arr = np.arange(Vnum_2D-num_row)
        Minn = arr.reshape(num_row,-1)
        col = Vnum_2D-num_row + np.arange(num_row)
        M2D = np.insert(Minn,Minn.shape[1],col,axis=1)
    self._vMatrix2d = M2D
    f1 = M2D[:-1,:-1].flatten()
    f2 = M2D[1:,:-1].flatten()
    f3 = M2D[1:,1:].flatten()
    f4 = M2D[:-1,1:].flatten() 
    self._v2D_f1234 = [f1,f2,f3,f4]
    
    return [[v,v1,v2,v3,v4], [v[indv],v1[indv],v2[indv],v3[indv],v4[indv]],\
        [M, M_inn], [indM,indM_inn], indv]

def get_rregular_split_list_for_polysegmnet(self, diag=False, is_poly=False):
    """
    Splits the regular polyline into segments.

    This function divides the regular polyline into segments based on the 
    specified criteria (diagonal or non-diagonal).

    Parameters
    ----------
    diag : bool, optional (default=False)
        Whether to split based on diagonal edges.
    is_poly : bool, optional (default=False)
        Whether to return the result as a polyline.

    Returns
    -------
    segments : list
        A list of segments, where each segment is a list of vertex indices.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The segments are useful for further processing and analysis.

    filter/divide the rregular-polyline into 2 class-lists
    This is the moset safe one, suitible for any mesh with even singularity

    See Also
    --------
    nonsingular : Identifies non-singular vertices in the mesh.
    """
    if diag:
        v0,v1,v2,v3,v4 = self.rr_star_corner
    else:
        v0,v1,v2,v3,v4 = self.rrv4f4
    
    def _nonmultiple(v1,v3):
        vl = np.r_[v1,v0]
        vr = np.r_[v0,v3]

        multi = []
        for i in range(len(vl)):
            j = np.where(vl==vl[i])[0]
            for jj in j:
                if jj>i and vr[jj]==vr[i]:
                    multi.append(jj)
            
            j = np.where(vr==vl[i])[0]
            for jj in j:
                if jj>i and vl[jj]==vr[i]:
                    multi.append(jj)
                    
        left = np.delete(np.arange(len(vl)),np.array(multi))            
        return vl[left], vr[left]

    v13l,v13r = _nonmultiple(v1,v3)
    v24l,v24r = _nonmultiple(v2,v4)

    bdry = self.boundary_vertices()
    _,_,lj = self.vertex_ring_vertices_iterators(return_lengths=True)
    singular = np.where(lj>4)[0]
    nonregular = np.unique(np.r_[bdry,singular])
    ##print(singular,nonregular)
    
    def _onelist(vl,vr):
        """
        left ~ vl ~ vr
        vil ~ vir ~ ck1 ~ ck2
        elist ~ sign
        seglist=[[],[],[]...] ~ signlist=[[],[],[]...]
        """
        vvl,vvr = vl,vr
        left = np.arange(len(vvl)) ##num of segments
        seglist = []
        signlist = []
        while len(left)!=0:
            elist = [left[0]] ## i-segment
            sign = [1]
            left = np.delete(left,0)
            vil,vir = vl[0],vr[0]
            vl = np.delete(vl,0)
            vr = np.delete(vr,0)
            vends = np.unique(np.r_[vl,vr])
            ck1 = vil in nonregular or vil not in vends
            ck2 = vir in nonregular or vir not in vends
            if ck1 and ck2:
                seglist.append(elist)
                signlist.append(sign)
            save_vlist = [vil,vir]
            ##print(save_vlist)
            while not (ck1 and ck2):    
                js1 = np.where(vl==vil)[0] ##start in leftstart
                js2 = np.where(vr==vil)[0] ##start in leftend
                je1 = np.where(vl==vir)[0] ##end in leftstart
                je2 = np.where(vr==vir)[0] ##end in leftend   
                if len(np.r_[js1,js2])!=0:
                    "expand in left-direction"
                    if len(js1)!=0:
                        j = js1[0]
                        vil = vr[j]
                        sign.insert(0,-1) ##=[-1,...]
                        elist.insert(0,left[j])
                    elif len(js2)!=0:
                        j = js2[0]
                        vil = vl[j]
                        sign.insert(0,1) ##[+1,...]
                        elist.insert(0,left[j])
                    # if len(sign)!=len(elist):
                    #     print(j,left[j])
                    left = np.delete(left,j)
                    vl = np.delete(vl,j)
                    vr = np.delete(vr,j)
                    vends = np.unique(np.r_[vl,vr])
                    ck1 = vil in nonregular or vil not in vends
                    ck3 = vil in save_vlist or vil == vir
                    if ck1 and ck2 or ck3 or len(left)==0:
                        #print('break1')
                        break
                    ##save_vlist.append(vil)
                    #print('l',vil)

                if len(np.r_[je1,je2])!=0:
                    "expand in right-direction"
                    je1 = np.where(vl==vir)[0] ##end in leftstart
                    je2 = np.where(vr==vir)[0] ##end in leftend
                    if len(je1)!=0:
                        j = je1[0]
                        vir = vr[j]
                        sign.append(1) ##[...,1]
                        elist.append(left[j])
                    elif len(je2)!=0:
                        j = je2[0]
                        vir = vl[j]
                        sign.append(-1) ##[...,-1]
                        elist.append(left[j])
                    # if len(sign)!=len(elist):
                    #     print(j,left[j])
                    left = np.delete(left,j)
                    vl = np.delete(vl,j)
                    vr = np.delete(vr,j)
                    vends = np.unique(np.r_[vl,vr])
                    ck2 = vir in nonregular or vir not in vends
                    ck3 = vir in save_vlist or vil == vir
                    if ck1 and ck2 or ck3 or len(left)==0:
                        #print('break2')
                        break
                    ##save_vlist.append(vir)
                    #print('r',vir)
                    
            #print(save_vlist)
            if ck1 and ck2 or ck3 or len(left)==0:
                seglist.append(elist)
                signlist.append(sign)
                
        #print('---',len(seglist))
        #print('---',len(signlist))

        alllist=[]
        for i in range(len(seglist)):
            elist = seglist[i]
            sign = signlist[i]
            ilist=[]
            for j in range(len(elist)):
                n,s = elist[j],sign[j]
                if s==-1:
                    ilist.append(vvr[n])
                    ilist.append(vvl[n])
                else:
                    ilist.append(vvl[n])
                    ilist.append(vvr[n])  
            jlist = ilist[::2]
            jlist.append(ilist[-1])
            alllist.append(jlist)
            
        ##print('===',len(alllist),alllist)
        
        seg_list = []
        for vlist in alllist:
            "break at singular vertex"
            if len(vlist)>2:
                ilist = np.array(vlist[1:-1])
                k = np.intersect1d(ilist,singular)
                if len(k)!=0:
                    div = []
                    for i in k:
                        j = np.where(ilist==i)[0][0]
                        div.append(j+1)
                    m = 0
                    for n in div:
                        if len(vlist[m:n+1])!=0:
                            seg_list.append(vlist[m:n+1])
                        m += n
                    if len(vlist[n:])!=0:
                        seg_list.append(vlist[n:])
                    
                else:
                    seg_list.append(vlist)
            else:
                seg_list.append(vlist)
                
        ##print('---',len(seg_list),seg_list) ##right: show all loop for one-poly
        vllist,vrlist,vstart,vend = [],[],[],[]
        innlist,iinnlist, = [],[]
        inn_in_vlist_arr,iinn_in_vlist_arr = [],[]
        k = 0
        for vlist in seg_list:
            vstart.append(vlist[0])  ##note: for rotation-onelist, vstart doesn't include endpts
            vend.append(vlist[-1])
            vllist.append(vlist[:-1])
            vrlist.append(vlist[1:])

            if len(vlist)>=4:
                innlist.append(vlist[1:-1]) ## 2-dim
                inn_in_vlist_arr.append(k)
                inn = []
                for v in vlist[1:-1]:
                    if v in v0:
                        inn.append(np.where(v0==v)[0][0])
                if len(inn) !=0:
                    iinnlist.append(inn) ## 2-dim
                    iinn_in_vlist_arr.append(k)
            k += 1

        vstart,vend = np.array(vstart), np.array(vend)     
        if False: ##ordered arranged polys; need to choose True or False
            """ new added after AAG-project for AoTu-project
            above can't make sure the eachpoly inorder,i.e. everysecond maybe wrong
                reorder the list: output the polylines in order 
                from left to right of the srf
                but only work for non-loop-vlist, 
                for loops without bdry, skip, todo in the furture.
            """
            order = self.get_index_in_polyline(vstart,patch_or_annulus=True) ##only work for patch+annulus; others need todo
            vstart, vend = vstart[order], vend[order]
            seg,l,r,inn,iinn = [],[],[],[],[]
            for i in order:
                seg.append(seg_list[i])
                l.append(vllist[i])
                r.append(vrlist[i])
                if i in inn_in_vlist_arr:
                    j = np.where(inn_in_vlist_arr==i)[0][0]
                    inn.append(innlist[j])
                if i in iinn_in_vlist_arr:
                    j = np.where(iinn_in_vlist_arr==i)[0][0]
                    iinn.append(iinnlist[j])
            seg_list,vllist,vrlist,innlist,iinnlist = seg,l,r,inn,iinn 
        return [seg_list,[vllist,vrlist],[innlist,iinnlist],[vstart,vend]]            
        #-----------------------------------------------------------------
        
    if is_poly:
        V = self.vertices
        
        if 0:
            ##mesh_polylines have bug, not work.
            seglists1 = _onelist(v13l,v13r)[0]
            pl1 = mesh_polylines(V, seglists1)
            seglists2 = _onelist(v24l,v24r)[0]
            pl2 = mesh_polylines(V, seglists2)
        
        
        vllist,vrlist = _onelist(v13l,v13r)[1]
        vl = np.array(sum(vllist,[]),dtype=int)
        vr = np.array(sum(vrlist,[]),dtype=int)
        
        pl1 = make_polyline_from_endpoints(V[vl],V[vr])
        
        V1ll = 2*V[vl]/3+V[vr]/3
        V1lr = V[vl]/3+2*V[vr]/3

        vllist,vrlist = _onelist(v24l,v24r)[1]
        vl = np.array(sum(vllist,[]),dtype=int)
        vr = np.array(sum(vrlist,[]),dtype=int)
        
        pl2 = make_polyline_from_endpoints(V[vl],V[vr])
        V2ll = 2*V[vl]/3+V[vr]/3
        V2lr = V[vl]/3+2*V[vr]/3
        return pl1,pl2,V1ll,V1lr,V2ll,V2lr

    return [_onelist(v13l,v13r), _onelist(v24l,v24r)]

def nonsingular(self):
    """
    Identifies non-singular vertices in the mesh.

    This function detects vertices with valence 4 (non-singular) and stores 
    their indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    None
        The non-singular vertices are stored within the class instance.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    Non-singular vertices are useful for mesh analysis and processing.
    nonsingular(=regular) vertices v in increased order

    See Also
    --------
    nonsingular_star_matrix : Creates a star matrix for non-singular vertices.
    """
    self._vi,self._vj,lj = self.vertex_ring_vertices_iterators(sort=True,return_lengths=True)
    order = np.where(lj==4)[0]
    self._ver_regular = order
    self._num_regular = len(order)
    self._corner = np.where(lj==2)[0]
    self._inner = np.setdiff1d(np.arange(self.V),self.boundary_vertices())

def nonsingular_star_matrix(self):
    """
    Creates a star matrix for non-singular vertices.

    This function constructs a matrix representing the connectivity of 
    non-singular vertices, where each row corresponds to a vertex and its 
    neighboring vertices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    star_matrix : numpy array
        The star matrix for non-singular vertices.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The star matrix is useful for organizing and accessing mesh vertices.

    See Also
    --------
    nonsingular : Identifies non-singular vertices in the mesh.
    """
    order = self.ver_regular
    ring = [[] for i in range(len(order))]
    for i in range(len(order)):
        v = order[i]
        ring[i]= self.ringlist[v]
    star = np.c_[order.T,ring]
    return star

def quadfaces(self):
    """
    Identifies quad faces in the mesh.

    This function detects quad faces and stores their indices and connectivity.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    None
        The quad faces are stored within the class instance.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    Quad faces are useful for mesh analysis and processing.

    See Also
    --------
    regular_vertex_regular_quad : Identifies regular vertices and quad faces.
    """
    f, v1, v2 = self.face_edge_vertices_iterators(order=True)
    f4,vi = [],[]
    for i in range(self.F):
        ind = np.where(f==i)[0]
        if len(ind)==4:
            f4.extend([i,i,i,i])
            vi.extend(v1[ind])
            #vj.extend(v2[ind])
    self._num_quadface = len(f4) // 4
    #v1,v2,v3,v4 = vi[::4],vi[1::4],vi[2::4],vi[3::4]
    self._quadface = np.array(vi,dtype=int)
    self._quadface_order = np.unique(f4)

def regular_vertex_regular_quad(self,delete_multi=True):
    """
    This function processes half-edge data structures to identify and organize 
    regular vertices and regular quadrilateral faces in a mesh. It returns 
    structured information about the regular vertices, regular quads, and their 
    connectivity.

    The function's main logic involves:

    1. Iterating over each vertex to identify oriented quad faces using the half-edge structure.

    2. Handling boundary cases where quad faces may be incomplete or partially defined.

    3. Optionally deleting multiple faces that overlap, ensuring a unique representation of quad faces.

    4. Returning detailed information about regular vertices, regular quads, and their connectivity.

    Parameters
    ----------
    halfedges : numpy int array
        The half-edge data structure representing the mesh. 
        Each row corresponds to a half-edge with columns representing: 
        [origin vertex, twin half-edge, next half-edge, prev half-edge, face].
    ver_star_matrix : numpy int array
        A matrix where each row corresponds to a vertex and contains the indices of its neighboring vertices.
    delete_multi : bool, optional (default True)
        whether to deletes multiple overlapping faces to ensure unique quad faces or not.

    Returns
    -------
    num_rrf : int
        The number of regular quad faces.
    rr_quadface : numpy int array
        The regular quad faces, each defined by four vertices.
    rr_quadface_order :  numpy int array
        The order of the quad faces.
    num_rrv : int
        The number of regular vertices.
    rr_star : numpy int array
        The star matrix for regular vertices.
    rr_4quad_vers : list
        A list of all oriented quad faces, including boundary faces.


    Note:
    --------
    This function assumes the input mesh is manifold and orientable. It also assumes that the half-edge structure is correctly defined, with each half-edge pointing to its twin, next, and previous half-edges.
        

    See Also: 
    --------
    quadfaces()

    
    Examples
    --------
    ```python
    from quadrings import regular_vertex_regular_quad
    num_rrf, rr_quadface, rr_quadface_order, num_rrv, rr_star, rr_4quad_vers = regular_vertex_regular_quad(halfedges, ver_star_matrix)
    ```
    """
    H = self.halfedges
    ##starM = np.array(orient_rings(self))
    starM = self.ver_star_matrix
    num = len(starM)
    f4 = []
    for i in range(num):
        "multiple oriented quad faces"
        v,v1,v2,v3,v4 = starM[i,:]
        ei = np.where(H[:,0]==v)[0]
        ej = H[H[H[H[ei,2],2],2],2]
        e1 = ei[np.where(H[H[ei,4],0]==v1)[0]]
        e2 = ei[np.where(H[H[ei,4],0]==v2)[0]]
        e3 = ei[np.where(H[H[ei,4],0]==v3)[0]]
        e4 = ei[np.where(H[H[ei,4],0]==v4)[0]]
        if any(list(ej-ei)): # whose neighbor include not quad face
            if H[e1,1]==-1 and H[H[e2,4],1]==-1:
                f4.append([v2,v,v1,-1])
                f4.append([H[H[H[e2,2],2],0][0],v3,v,v2])
                f4.append([v3,H[H[H[e3,2],2],0][0],v4,v])
                f4.append([v,v4,H[H[H[e4,2],2],0][0],v1])
            elif H[e2,1]==-1 and H[H[e3,4],1]==-1:
                f4.append([v2,v,v1,H[H[H[e1,2],2],0][0]])
                f4.append([-1,v3,v,v2])
                f4.append([v3,H[H[H[e3,2],2],0][0],v4,v])
                f4.append([v,v4,H[H[H[e4,2],2],0][0],v1])
            elif H[e3,1]==-1 and H[H[e4,4],1]==-1:
                f4.append([v2,v,v1,H[H[H[e1,2],2],0][0]])
                f4.append([H[H[H[e2,2],2],0][0],v3,v,v2])
                f4.append([v3,-1,v4,v])
                f4.append([v,v4,H[H[H[e4,2],2],0][0],v1])
            elif H[e4,1]==-1 and H[H[e1,4],1]==-1:
                f4.append([v2,v,v1,H[H[H[e1,2],2],0][0]])
                f4.append([H[H[H[e2,2],2],0][0],v3,v,v2])
                f4.append([v3,H[H[H[e3,2],2],0][0],v4,v])
                f4.append([v,v4,-1,v1])
        else:
            if H[H[H[H[e1,2],2],2],0]==v2:
                "one quad face [v2,v,v1,x]"
                f4.append([v2,v,v1,H[H[H[e1,2],2],0][0]])
            if H[H[H[H[e2,2],2],2],0]==v3:
                "one quad face [x,v3,v,v2]"
                f4.append([H[H[H[e2,2],2],0][0],v3,v,v2])
            if H[H[H[H[e3,2],2],2],0]==v4:
                "one quad face [v3,x,v4,v]"
                f4.append([v3,H[H[H[e3,2],2],0][0],v4,v])
            if H[H[H[H[e4,2],2],2],0]==v1:
                "one quad face [v,v4,x,v1]"
                f4.append([v,v4,H[H[H[e4,2],2],0][0],v1])


    farr = np.unique(f4,axis=0)
    a,b = np.where(farr==-1)
    farr = np.delete(farr,a,axis=0)
    forder=np.array([],dtype=int)
    for f in farr:
        e1=np.where(H[:,0]==f[0])[0]
        e2=np.where(H[H[:,4],0]==f[1])[0]
        e = np.intersect1d(e1,e2)
        forder = np.r_[forder, H[e,1]]

    f4list = np.array(f4)
    if delete_multi: # delete multiple-left faces
        #forder, ind = np.unique(forder,return_index=True) # changed order
        ind=[]
        multi=[]
        for i in range(len(forder)):
            f = forder[i]
            if f not in forder[ind]:
                ind.append(i)
            else:
                j = np.where(forder[ind]==f)[0][0]
                k = forder[ind][j]
                l = np.setdiff1d(np.where(forder==k)[0], np.array([i]))[0]
                multi.append(list(farr[l]))
        forder = forder[ind]
        farr = farr[ind]
        for f in multi:
            index=np.array([],dtype=int)
            e1,e2,e3,e4 = f
            a,b,c = [e4,e1,e2,e3],[e3,e4,e1,e2],[e2,e3,e4,e1]
            if a in f4:
                ind,_ = np.where(f4list==a)
                index = np.r_[index,ind]
            if b in f4:
                ind,_ = np.where(f4list==b)
                index = np.r_[index,ind]
            if c in f4:
                ind,_ = np.where(f4list==c)
                index = np.r_[index,ind]
            f4list[index]=np.array(f)

    self._num_rrf = len(farr)
    self._rr_quadface = farr
    self._rr_quadface_order = forder
    #print(len(farr),len(f4list),num,len(index))
    self._num_rrv = num # same with num_regular
    #order = np.setdiff1d(np.arange(num),index//4)
    self._rr_star = starM # same with starM
    self._rr_4quad_vers = f4 #rr_4quad_vers
    #return f4, farr, starM

def index_of_rrquad_face_vertex_with_rrv(self):
    """
    Retrieves the indices of vertices in regular quad faces that are also regular vertices.

    This function identifies the vertices in regular quad faces that are also 
    regular vertices and returns their indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    indices : numpy array
        The indices of vertices in regular quad faces that are also regular vertices.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The indices are useful for accessing specific vertices in the mesh.
    self.rr_quadface() may include vertices that are not in rr_star; this one works

    See Also
    --------
    regular_vertex_regular_quad : Identifies regular vertices and quad faces.
    """
    v1,v2,v3,v4  = self.rr_quadface.T
    rrv = self.ver_rrv4f4
    ind = []
    for i in range(len(v1)):
        if v1[i] in rrv and v2[i] in rrv and v3[i] in rrv and v4[i] in rrv:
            ind.append(i)
    self._ind_rr_quadface_with_rrv = np.array(ind)
    
def get_rr_quadface_boundaryquad_index(self, vb_quad=True):
    """
    Retrieves the indices of boundary quad faces in the regular quad face list.

    This function identifies the boundary quad faces and returns their indices 
    in the regular quad face list.

    Parameters
    ----------
    vb_quad : bool, optional (default=True)
        Whether to include quad faces that have vertices on the boundary.

    Returns
    -------
    inner_indices : numpy array
        The indices of inner quad faces.
    boundary_indices : numpy array
        The indices of boundary quad faces.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The indices are useful for accessing specific quad faces in the mesh.

    See Also
    --------
    regular_vertex_regular_quad : Identifies regular vertices and quad faces.
    """
    v1,v2,v3,v4  = self.rr_quadface.T
    forder = self.rr_quadface_order
    fb = self.boundary_faces()
    _,out,_ = np.intersect1d(forder,fb, return_indices=True)
    inn = np.setdiff1d(np.arange(len(v1)), out)
    if vb_quad:
        "also including quad'vertex belong to boundary"
        boundary = self.boundary_vertices()
        vb = []
        for i in inn:
            if v1[i] in boundary or v2[i] in boundary:
                vb.append(i)
            if v3[i] in boundary or v4[i] in boundary:
                vb.append(i)
        inn = np.setdiff1d(inn,np.array(vb))
        out = np.r_[out,np.array(vb,dtype=int)]
    return inn, out

def get_rr_vs_4face_corner(self):
    """
    Retrieves the corner vertices of regular quad faces.

    This function identifies the corner vertices of regular quad faces and 
    returns their indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    corner_vertices : numpy array
        The indices of corner vertices of regular quad faces.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The corner vertices are useful for mesh analysis and processing.

    self.rr_star_corner (should be same with get_vs_diagonal_v)
    vertex star' 4faces' 4 corner vertex [a,b,c,d]
        a   1    d
        2   v    4
        b   3    c

    See Also
    --------
    regular_vertex_regular_quad : Identifies regular vertices and quad faces.
    """
    H = self.halfedges
    v,v1,v2,v3,v4 = self.rrv4f4
    va,vb,vc,vd = [],[],[],[]
    for i in range(len(v)):
        e1=np.intersect1d(np.where(H[:,0]==v[i])[0],np.where(H[H[:,4],0]==v1[i])[0])
        va.append(H[H[H[e1,2],2],0])
        e2=np.intersect1d(np.where(H[:,0]==v[i])[0],np.where(H[H[:,4],0]==v2[i])[0])
        vb.append(H[H[H[e2,2],2],0])            
        e3=np.intersect1d(np.where(H[:,0]==v[i])[0],np.where(H[H[:,4],0]==v3[i])[0])
        vc.append(H[H[H[e3,2],2],0])
        e4=np.intersect1d(np.where(H[:,0]==v[i])[0],np.where(H[H[:,4],0]==v4[i])[0])
        vd.append(H[H[H[e4,2],2],0])      
    self._rr_star_corner = np.c_[v,va,vb,vc,vd].T    

        
def get_rr_vs_4face_centers(self):
    """
    Computes the centers of regular quad faces in the mesh.

    This function calculates the geometric centers (barycenters) of regular quad 
    faces in the mesh. The centers are computed as the average of the vertex 
    positions of each quad face.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    face_centers : numpy array
        The coordinates of the centers of regular quad faces.

    Note
    ----
    This function assumes that the mesh is a regular mesh with quad faces.
    The face centers are useful for mesh analysis, processing, and visualization.

    Example
    -------
    >>> mesh = Mesh()  # Assuming Mesh is a class that contains the mesh data
    >>> face_centers = mesh.get_rr_vs_4face_centers()
    >>> print(face_centers)
    Array of face center coordinates

    See Also
    --------
    get_rr_vs_bounary : Separates regular vertices into inner and boundary vertices.
    get_rr_quadface_boundaryquad_index : Retrieves indices of boundary quad faces.
    """
    V = self.vertices
    f4 = self.rr_4quad_vers
    f41,f42,f43,f44 = f4[::4],f4[1::4],f4[2::4],f4[3::4]
    def _face_center(f4i):
        fic = (V[f4i[::1]]+V[f4i[1::4]]+V[f4i[2::4]]+V[f4i[3::4]]) / 4.0
        return fic
    f1c = _face_center(f41)
    f2c = _face_center(f42)
    f3c = _face_center(f43)
    f4c = _face_center(f44)
    return f1c,f2c,f3c,f4c

def get_rr_vs_bounary(self):
    """
    Separates regular vertices into inner and boundary vertices.

    This function identifies which regular vertices are on the boundary and 
    which are in the interior of the mesh.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    inner_indices : numpy array
        Indices of regular vertices that are in the interior.
    boundary_indices : numpy array
        Indices of regular vertices that are on the boundary.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The separation of vertices is useful for mesh analysis and processing.

    See Also
    --------
    boundary_vertices : Retrieves the indices of boundary vertices.
    """
    fb = self.boundary_faces()
    vsb = np.array([],dtype=int)
    for f in fb:
        vsb = np.r_[vsb,np.array(self.faces_list()[f])]
    vsb = np.unique(vsb)
    rrv = self.ver_rrv4f4
    inn = np.setdiff1d(rrv,vsb)
    rrb = np.intersect1d(rrv,vsb)
    idi,idb = [],[]
    for i in inn:
        idi.append(np.where(rrv==i)[0][0])
    for j in rrb:
        idb.append(np.where(rrv==j)[0][0])
    return np.array(idi), np.array(idb)

def index_of_4quad_face_order_at_regular_vs(self):
    """
    Retrieves the indices of quad faces adjacent to regular vertices.

    This function identifies the quad faces adjacent to each regular vertex 
    and returns their indices in the quad face list.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    indices : numpy array
        Indices of quad faces adjacent to each regular vertex.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The indices are useful for accessing specific quad faces in the mesh.

    most strong regular case: regular vertex & regular quads; 
    star = self.rr_star
    star = star[self.ind_rr_star_v4f4]
    if1,if2,if3,if4 = self.ind_rr_quadface_order.T

    See Also
    --------
    regular_vertex_regular_quad : Identifies regular vertices and quad faces.
    """
    H = self.halfedges
    star = self.rr_star#[self.ind_rr_star_v4f4] ##NOTE: NEED TO CHECK IF [self.ind_rr_star_v4f4]
    forder = self.rr_quadface_order
    flist = []
    ind = []
    for i in range(len(star)):
        v,v1,v2,v3,v4 = star[i,:]
        e=np.where(H[:,0]==v)[0]
        e1=np.where(H[H[:,4],0]==v1)[0]
        e2=np.where(H[H[:,4],0]==v2)[0]
        e3=np.where(H[H[:,4],0]==v3)[0]
        e4=np.where(H[H[:,4],0]==v4)[0]
        i1 = np.intersect1d(e,e1)
        i2 = np.intersect1d(e,e2)
        i3 = np.intersect1d(e,e3)
        i4 = np.intersect1d(e,e4)
        f1,f2,f3,f4 = H[i1,1],H[i2,1],H[i3,1],H[i4,1]
        if f1!=-1 and f2!=-1 and f3!=-1 and f4!=-1:
            "for v whose faces are not 4, de-select it"
            if1 = np.where(forder==f1)[0]
            if2 = np.where(forder==f2)[0]
            if3 = np.where(forder==f3)[0]
            if4 = np.where(forder==f4)[0]
            if len(if1)!=0 and len(if2)!=0 and len(if3)!=0 and len(if4)!=0:
                "for face who does belong to forder"
                ind.append(i)
                flist.append([if1[0],if2[0],if3[0],if4[0]])
    self._ind_rr_star_v4f4 = np.array(ind) # ind_rr_star_v4f4
    self._ind_rr_quadface_order = np.array(flist) 


def orient(self, S0, A, B, C, D, Nv4):
    """
    Orients the mesh based on the given parameters.

    This function adjusts the orientation of the mesh by modifying the 
    vertex normals and related parameters.

    Parameters
    ----------
    S0 : numpy array
        The center point of the mesh.
    A : numpy array
        The first parameter for orientation.
    B : numpy array
        The second parameter for orientation.
    C : numpy array
        The third parameter for orientation.
    D : numpy array
        The fourth parameter for orientation.
    Nv4 : numpy array
        The normal vectors of the mesh.

    Returns
    -------
    B : numpy array
        The modified second parameter.
    C : numpy array
        The modified third parameter.
    D : numpy array
        The modified fourth parameter.
    Nv4 : numpy array
        The modified normal vectors.
    x_orient : numpy array
        The orientation factor.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The orientation is useful for mesh analysis and processing.

    See Also
    --------
    new_vertex_normals : Computes the vertex normals of the mesh.
    """
    orientrn = self.new_vertex_normals()
    ##print(Nv4.shape, orientrn.shape)

    ind1 = np.where(np.einsum('ij,ij->i',orientrn,Nv4) < 0)[0]
    if len(ind1)!=0:
        Nv4[ind1] = -Nv4[ind1]
        "new_c = v[ind1]+r[ind1]*n[ind1], --> (b,c,d)[ind1] = -2*a*(v+rn)[ind1]"
        x,y,z = S0.T
        bb = -4*A*x-B
        cc = -4*A*y-C
        dd = -4*A*z-D
        B[ind1] = bb[ind1]
        C[ind1] = cc[ind1]
        D[ind1] = dd[ind1]
        
    x_orient = np.sqrt(np.abs(np.einsum('ij,ij->i',Nv4,orientrn)))        
    return B,C,D,Nv4,x_orient

def new_vertex_normals(self, is_orient_to_spherecenter=True, is_rrv4f4=True):
    """
    Computes the vertex normals of the mesh.

    This function calculates the vertex normals and optionally orients them 
    towards the sphere center.

    Parameters
    ----------
    is_orient_to_spherecenter : bool, optional (default=True)
        Whether to orient the normals towards the sphere center.
    is_rrv4f4 : bool, optional (default=True)
        Whether to use the regular vertices and quad faces.

    Returns
    -------
    normals : numpy array
        The vertex normals of the mesh.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The vertex normals are useful for mesh analysis and processing.
    It makes vertex_normals orient non-convex direction

    See Also
    --------
    vertex_normals : Computes the vertex normals of the mesh.
    """
    v0,v1,v2,v3,v4 = self.rrv4f4 #(self.rr_star).T
    V = self.vertices
    VN = np.copy(self.vertex_normals())
    
    orientrn = self.get_v4_orient_unit_normal(rregular=is_rrv4f4)[1]
    id0 = np.where(np.einsum('ij,ij->i', orientrn, VN[v0]) < 0)[0]
    if len(id0)!=0:
        orientrn[id0] = -orientrn[id0]
    VN[v0] = orientrn
        
    if is_orient_to_spherecenter:
        S0,S1,S2,S3,S4 = V[v0],V[v1],V[v2],V[v3],V[v4]
        _,radius,coeff,Nv4 = interpolate_sphere(S0,S1,S2,S3,S4)#Nv4 is from vertex to center
        id1 = np.where(np.einsum('ij,ij->i', VN[v0], Nv4) < 0)[0] 
        if len(id1)>int(len(v0)/2-1):
            orientn = -VN
        else:
            orientn = VN
    else:
        Ns = self.vertex_normals()[v0]
        uE1 = (V[v1]-V[v0]) / np.linalg.norm(V[v1]-V[v0],axis=1)[:,None]
        uE2 = (V[v2]-V[v0]) / np.linalg.norm(V[v2]-V[v0],axis=1)[:,None]
        uE3 = (V[v3]-V[v0]) / np.linalg.norm(V[v3]-V[v0],axis=1)[:,None]
        uE4 = (V[v4]-V[v0]) / np.linalg.norm(V[v4]-V[v0],axis=1)[:,None]
        N1, N2 = uE1 + uE3, uE2 + uE4
        id1 = np.where(np.einsum('ij,ij->i', Ns, N1) < 0)[0] # non-convex direction
        id2 = np.where(np.einsum('ij,ij->i', Ns, N2) < 0)[0] # non-convex direction
        number = max(len(id1),len(id2))
        if number > len(v0)-len(id1)-len(id2):
            orientn = self.vertex_normals()
        else:
            orientn = -self.vertex_normals()
            
    if is_rrv4f4:
        return orientn[v0]
    else:
        return orientn

def boundary_vertex_3neibs(self):
    """
    Retrieves the neighbors of boundary vertices with valence 3.

    This function identifies the neighbors of boundary vertices that have 
    valence 3 and returns their indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    boundary_vertices : numpy array
        Indices of boundary vertices with valence 3.
    left_neighbors : numpy array
        Indices of the left neighbors.
    right_neighbors : numpy array
        Indices of the right neighbors.
    inner_neighbors : numpy array
        Indices of the inner neighbors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The neighbors are useful for mesh analysis and processing.

    bdry(include corner) --> return [vl,vr,vinn];
    vl,vr in bdry, vinn in inner
    bdry_3neibs[0]:= 2nd_boundary_vertices, including corners
    bdry_3neibs[1]:= left of bdry_vertices
    bdry_3neibs[2]:= right of bdry_vertices

    See Also
    --------
    boundary_vertices : Retrieves the indices of boundary vertices.
    """
    H = self.halfedges
    ibdry = self.boundary_vertices()
    vl,vr,vinn = [],[],[]
    for v in ibdry:
        ie = np.where(H[:,0]==v)[0]
        if len(ie)==3:
            "bdry vertex of valence 3"
            inn = np.intersect1d(np.where(H[ie,1]!=-1)[0],np.where(H[H[ie,4],1]!=-1)[0])
            j = ie[inn][0]
            vinn.append(H[H[j,4],0])
            vl.append(H[H[j,3],0])
            vr.append(H[H[H[H[j,4],2],4],0])
        elif len(ie)==2:
            "corner with valence 2"          
            if H[ie[0],1]==-1:
                j = ie[0][0]
            else:
                j = ie[1][0]
            vl.append(H[H[j,4],0])
            vr.append(H[H[j,3],0])
            vinn.append([H[H[H[j,4],3],0]])
    self._ver_bdry_3neibs = [ibdry,np.array(vinn),np.array(vl),np.array(vr)]
    
def vertex_valence3_neib(self, boundary_vertex=True, corner=False):
    """
    Retrieves the neighbors of vertices with valence 3.

    This function identifies the neighbors of vertices that have valence 3 
    and returns their indices.

    Parameters
    ----------
    boundary_vertex : bool, optional (default=True)
        Whether to consider boundary vertices.
    corner : bool, optional (default=False)
        Whether to consider corner vertices.

    Returns
    -------
    vertices : numpy array
        Indices of vertices with valence 3.
    neighbors : numpy array
        Indices of the neighbors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The neighbors are useful for mesh analysis and processing.

    See Also
    --------
    boundary_vertex_3neibs : Retrieves neighbors of boundary vertices with valence 3.
    """
    v,vj,lj = self.vertex_ring_vertices_iterators(sort=True,return_lengths=True)
    ##o2 = np.where(lj==2)[0]
    o3 = np.where(lj==3)[0]
    ##o4 = np.where(lj==4)[0]
    boundary = self.boundary_vertices()
    ##self.bvalen2 = np.intersect1d(o2, boundary)

    if boundary_vertex:
        bvalen3 = np.intersect1d(o3, boundary)
        if len(bvalen3)==0:
            pass
        else:
            H = self.halfedges
            vl,vr = [],[]
            for v in bvalen3:
                "in order"
                ie = np.where(H[:,0]==v)[0]
                i = np.intersect1d(np.where(H[ie,1]!=-1)[0],np.where(H[H[ie,4],1]!=-1)[0])
                vl.append(H[H[ie[i],3],0][0])
                vr.append(H[H[H[H[ie[i],4],2],4],0][0])                    
            self._ver_bod_valence3_neib = [bvalen3,np.array(vl),np.array(vr)]
    else:
        invalen3 = np.setdiff1d(o3,boundary)
        if len(invalen3)==0:
            pass
        neib = []
        for v in invalen3:
            neib.append(self.ringlist[v])
        self._ver_inn_valence3_neib = [invalen3, np.array(neib)]
    
def vertex_corner_valence3_neib(self):
    """
    Retrieves the neighbors of corner vertices with valence 3.

    This function identifies the neighbors of corner vertices that have 
    valence 3 and returns their indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    vertices : numpy array
        Indices of corner vertices with valence 3.
    left_neighbors : numpy array
        Indices of the left neighbors.
    right_neighbors : numpy array
        Indices of the right neighbors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The neighbors are useful for mesh analysis and processing.

    See Also
    --------
    vertex_valence3_neib : Retrieves neighbors of vertices with valence 3.
    """
    v,vl,vr = self.ver_bod_valence3_neib
    c = self.corner
    ic = []
    for i in range(len(v)):
        if vl[i] in c or vr[i] in c:
            ic.append(i)
    ic = np.array(ic)
    return v[ic],vl[ic],vr[ic]   
    
def vertex_valence5_neib(self):
    """
    Retrieves the neighbors of vertices with valence 5.

    This function identifies the neighbors of vertices that have valence 5 
    and returns their indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    vertices : numpy array
        Indices of vertices with valence 5.
    neighbors : numpy array
        Indices of the neighbors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The neighbors are useful for mesh analysis and processing.

    See Also
    --------
    vertex_valence3_neib : Retrieves neighbors of vertices with valence 3.
    """
    v,vj,lj = self.vertex_ring_vertices_iterators(sort=True,return_lengths=True)
    o5 = np.where(lj==5)[0]
    boundary = self.boundary_vertices()
    inv5 = np.setdiff1d(o5, boundary)   
    if len(inv5)==0:
        pass
    else:
        neib = []
        for v in inv5:
            neib.append(self.ringlist[v])
        self._ver_inn_valence5_neib = [inv5, np.array(neib)]
    
def vertex_valence6_neib(self):
    """
    Retrieves the neighbors of vertices with valence 6.

    This function identifies the neighbors of vertices that have valence 6 
    and returns their indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    vertices : numpy array
        Indices of vertices with valence 6.
    neighbors : numpy array
        Indices of the neighbors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The neighbors are useful for mesh analysis and processing.

    See Also
    --------
    vertex_valence3_neib : Retrieves neighbors of vertices with valence 3.
    """
    v,vj,lj = self.vertex_ring_vertices_iterators(sort=True,return_lengths=True)
    o6 = np.where(lj==6)[0]
    boundary = self.boundary_vertices()
    inv6 = np.setdiff1d(o6, boundary)   
    if len(inv6)==0:
        pass
    else:
        neib = []
        for v in inv6:
            neib.append(self.ringlist[v])
        self._ver_inn_valence6_neib = [inv6, np.array(neib)]     
        
def vertex_valence4_neib(self, corner=True):
    """
    Retrieves the neighbors of vertices with valence 4.

    This function identifies the neighbors of vertices that have valence 4 
    and returns their indices.

    Parameters
    ----------
    corner : bool, optional (default=True)
        Whether to consider corner vertices.

    Returns
    -------
    vertices : numpy array
        Indices of vertices with valence 4.
    neighbors : numpy array
        Indices of the neighbors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The neighbors are useful for mesh analysis and processing.

    v, neib=[v1,v2,v3,v4,va,vb,vc,vd]

    See Also
    --------
    vertex_valence3_neib : Retrieves neighbors of vertices with valence 3.
    """
    _,_,lj = self.vertex_ring_vertices_iterators(sort=True,return_lengths=True)
    o4 = np.where(lj==4)[0]
    boundary = self.boundary_vertices()
    inv4 = np.setdiff1d(o4, boundary)   
    if len(inv4)==0:
        pass
    else:
        H = self.halfedges
        vc4 = []
        neib = np.array([],dtype=int)
        for v in inv4:
            if len(np.intersect1d(self.ringlist[v],boundary))==2:
                vc4.append(v)
                ie = np.where(H[:,0]==v)[0]
                abcd = H[H[H[ie,2],2],0]
                neib = np.r_[neib,abcd,np.array(self.ringlist[v])]
        neib = neib.reshape(-1,8)
        self._ver_corner_valence4_neib = [np.array(vc4), neib]      
    
def vertex_corner_neib(self):
    """
    Retrieves the neighbors of corner vertices.

    This function identifies the neighbors of corner vertices and returns 
    their indices.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    left_neighbors : numpy array
        Indices of the left neighbors.
    center_neighbors : numpy array
        Indices of the center neighbors.
    right_neighbors : numpy array
        Indices of the right neighbors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The neighbors are useful for mesh analysis and processing.

    See Also
    --------
    vertex_valence3_neib : Retrieves neighbors of vertices with valence 3.
    """
    H = self.halfedges
    corner = self.corner
    el = []
    for i in range(len(corner)):
        c = corner[i]
        e = np.intersect1d(np.where(H[:,0]==c)[0],np.where(H[:,1]==-1)[0])[0]
        el.append(e)
    va,vb = H[H[el,2],0],H[H[H[el,2],2],0]
    v1,v2 = H[H[el,3],0],H[H[H[el,3],3],0]
    vl,vc,vr = np.r_[corner,corner],np.r_[va,v1],np.r_[vb,v2]
    self._ver_corner_neib = [vl,vc,vr]
    

def get_a_boundary_L_strip(self, direction):
    """
    Retrieves an L-shaped strip of boundary vertices.

    This function identifies an L-shaped strip of boundary vertices based on 
    the specified direction.

    Parameters
    ----------
    direction : bool
        The direction of the L-strip (True for one direction, False for the other).

    Returns
    -------
    strip : numpy array
        The indices of vertices forming the L-strip.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The L-strip is useful for mesh analysis and processing.

    AG-net: only rectangular-patch shape

    See Also
    --------
    get_cylinder_annulus_mesh_diagonal_oriented_vertices : Retrieves diagonal-oriented vertices.
    """
    H = self.halfedges
    v,va,vb,vc,vd = self.rr_star_corner
    corner = self.corner
    vii = vb if direction else va
    i = np.intersect1d(vii,corner)[0]

    def _get_v13(i):
        "should have at least one"
        e = np.intersect1d(np.where(H[:,1]==-1)[0],np.where(H[:,0]==i)[0])[0]
        e0 = e
        ib1,ib2 = [],[]
        while H[H[e,2],0] not in corner:
            ib1.append(e)
            e = H[e,2]
        ib1.append(e)
        
        #---------
        # e = H[e,2]
        # while H[H[e,2],0] not in corner:
        #     e = H[e,2]
        #     ib2.append(e)  
            
        # v1,v3 = H[ib1,0], H[H[H[ib1,4],3],0]
        # v1 = np.r_[v1,H[H[H[H[ib2,4],3],3],0]]
        # v3 = np.r_[v3,H[H[ib2,4],0]]

        #----------------
        ib1 = ib1[::-1]
        v1,v3 = H[ib1,0], H[H[H[ib1,4],3],0]
        while H[H[e0,3],0] not in corner:
            e0 = H[e0,3]
            ib2.append(e0)
        ib2.append(H[e0,3])
        
        v1 = np.r_[v1,H[H[ib2[1:],4],0]]
        v3 = np.r_[v3,H[H[H[H[ib2[1:],4],3],3],0]]  
        return np.c_[v1,v3]
    "return L-shape boundary quads (Lx1): v and only 1diag v"
    return _get_v13(i)

def get_cylinder_annulus_mesh_diagonal_oriented_vertices(self, direction):
    """
    Retrieves diagonal-oriented vertices for a cylinder or annulus mesh.

    This function identifies vertices oriented along the diagonal direction 
    for a cylinder or annulus mesh.

    Parameters
    ----------
    direction : bool
        The direction of the diagonal (True for one direction, False for the other).

    Returns
    -------
    vertices : numpy array
        The indices of diagonal-oriented vertices.

    Note
    ----
    This function assumes that the mesh is a cylinder or annulus mesh.
    The diagonal-oriented vertices are useful for mesh analysis and processing.

    AG-net: only work for cylinder-annulus-shape

    return loop boundary quads (Lx1): v and only 1diag v
    

    See Also
    --------
    get_a_boundary_L_strip : Retrieves an L-shaped strip of boundary vertices.
    """
    H = self.halfedges
    vb,_ = self.get_i_boundary_vertex_indices(0)
    vnext = []
    for v in vb:
        if direction:
            e = np.intersect1d(np.where(H[:,1]==-1)[0],np.where(H[:,0]==v)[0])[0]
            vdiag = H[H[H[e,4],3],0]
        else:
            e = np.intersect1d(np.where(H[:,1]==-1)[0],np.where(H[:,0]==v)[0])[0]
            e = H[e,3]
            vdiag = H[H[H[H[e,4],2],2],0]
        vnext.append(vdiag)
    return np.c_[vb, vnext]
    

## using this first for choosing mesh polys
def get_both_isopolyline(self,diagpoly=False,is_one_or_another=False,
                            is_poly=False,only_inner=False,interval=1,
                            is_demultiple=False):
    """
    Retrieves two families of isopolyline curves from the mesh.

    This function identifies and returns two families of isopolyline curves 
    based on the specified parameters.

    Parameters
    ----------
    diagpoly : bool, optional (default=False)
        Whether to use diagonal polyline.
    is_one_or_another : bool, optional (default=False)
        Whether to select one family or the other.
    is_poly : bool, optional (default=False)
        Whether to return the result as polylines.
    only_inner : bool, optional (default=False)
        Whether to consider only inner vertices.
    interval : int, optional (default=1)
        The interval between selected vertices.
    is_demultiple : bool, optional (default=False)
        Whether to remove multiple edges.

    Returns
    -------
    polylines : list
        A list of isopolyline curves.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The isopolyline curves are useful for mesh analysis and processing.

    AG-net:works for any geometry; but the polys are not orderly arranged

    menubar_basicplot: show_isoline1,show_isoline2; no multiple edges

    See Also
    --------
    get_isoline_vertex_list : Retrieves isoline vertices.
    """
    _,_,lj = self.vertex_ring_vertices_iterators(return_lengths=True)
    o56 = np.where(lj>4)[0]
    o4 = np.where(lj==4)[0]
    boundary = self.boundary_vertices()
    inv4 = np.intersect1d(o4, boundary)  ## if !=0, bdry has V4
    if len(o56)==0 and len(inv4)==0: ## no singular inside, no v4 in bdry
        "regular-patch-shape or rotational-shape"
        ipllist = []
        vl = vr = np.array([],dtype=int)
        if diagpoly:
            if len(self.corner)!=0:
                "surface shape: rectangular-patch"
                vf1 = self.get_a_boundary_L_strip(direction=True)
                vf2 = self.get_a_boundary_L_strip(direction=False)
            else:
                "surface shape: cylinder-annulus"
                vf1 = self.get_cylinder_annulus_mesh_diagonal_oriented_vertices(True)
                vf2 = self.get_cylinder_annulus_mesh_diagonal_oriented_vertices(False)
            allplv = vf1 if is_one_or_another else vf2
            if interval ==1:
                for i, k in enumerate(allplv):
                    if i%interval==0:
                        iv,_ = get_diagonal_polyline_from_2points(self,k,is_poly=False)   
                        if len(iv)!=0:
                            ipllist.append(iv)
                            vl = np.r_[vl,iv[:-1]]
                            vr = np.r_[vr,iv[1:]]
            else:
                ##AAG_2_4x4patch_16x16_add1loopbdry +3
                # for i in range(len(allplv)+3): 
                #     if i ==0 or i==1 or i==2:
                #         continue
                #     else:
                #         k = allplv[i-3]
                for i in range(len(allplv)+1): 
                    if i ==0 :
                        continue
                    else:
                        k = allplv[i-1]
                        if i%interval==0:
                            iv,_ = get_diagonal_polyline_from_2points(self,k,is_poly=False)   
                            if len(iv)!=0:
                                ipllist.append(iv)
                                vl = np.r_[vl,iv[:-1]]
                                vr = np.r_[vr,iv[1:]]
        else:
            if len(self.corner)!=0:
                M = self.patch_matrix
            else:
                M = self.rot_patch_matrix  
            if not is_one_or_another:
                M = M.T
            if only_inner:
                #M = M[1:-1,1:-1]
                M = M[1:-1,:] ## to except the bdry strips
            allplv = M.tolist() 
            for i, iv in enumerate(allplv):
                if i%interval==0:
                    ipllist.append(iv)
                    if len(iv)!=0:
                        vl = np.r_[vl,iv[:-1]]
                        vr = np.r_[vr,iv[1:]]
    # elif len(o56)==0 and len(inv4)!=0:
    #     "patch"
    #     vb1,_ = self.get_i_boundary_vertex_indices(0) # i=0,1,2,3
    #     try:
    #         vb2,_ = self.get_i_boundary_vertex_indices(1)# i=0,1,2,3 # to check if 1 or 2
    #     except:
    #         vb2,_ = self.get_i_boundary_vertex_indices(1)
    #     vb = vb1 if is_one_or_another else vb2
    #     allplv = get_isoline_between_2bdry(self,vb)   
    #     #ipllist = allplv
        
    #     if only_inner:
    #         if allplv[0][0] in self.corner:
    #             allplv.pop(0)
    #         if allplv[-1][0] in self.corner:
    #             allplv.pop()
    elif len(o56)==0 and len(self.corner)!=0:
        "schwarzh_02_diag_unitscale_AAG_AAG"        
        ipl1,ipl2 = self.get_2families_polyline_from_1closed_bdry(diag=diagpoly,
                                                        interval=interval,
                                                        inner=False) ## need to choose True or False
        allplv = ipl1[1] if is_one_or_another else ipl2[1]
        
        ipllist = []
        vl = vr = np.array([],dtype=int)
        for i, iv in enumerate(allplv):
            if i%interval==0:
                ipllist.append(iv)
                if len(iv)!=0:
                    vl = np.r_[vl,iv[:-1]]
                    vr = np.r_[vr,iv[1:]]

    if is_poly:
        Vl,Vr = self.vertices[vl], self.vertices[vr]
        return make_polyline_from_endpoints(Vl,Vr)
    else:
        if is_demultiple:
            "for ZIG-bdry mesh, if others need to check again"
            from huilab.huimesh.curves import remove_multiple_repetitive_lists
            ipllist = remove_multiple_repetitive_lists(ipllist)
        return ipllist

def get_isoline_vertex_list(self, interval, another_direction=True):
    """
    Retrieves the vertices forming isoline curves.

    This function identifies and returns the vertices forming isoline curves 
    based on the specified interval and direction.

    Parameters
    ----------
    interval : int
        The interval between selected vertices.
    another_direction : bool, optional (default=True)
        Whether to consider another direction.

    Returns
    -------
    all_v0 : list
        All vertices forming the isoline curves.
    all_vs_v0 : list
        All neighboring vertices forming the isoline curves.
    select_v0 : list
        Selected vertices forming the isoline curves.
    select_vs_v0 : list
        Selected neighboring vertices forming the isoline curves.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The isoline vertices are useful for mesh analysis and processing.

    along one bdry for a patch-shape; two bdry for a star-shape

    See Also
    --------
    get_both_isopolyline : Retrieves two families of isopolyline curves.
    """
    v,v1,v2,v3,v4 = self.rrv4f4
    if another_direction:
        vl,vr = v2,v4
    else:
        vl,vr = v1,v3
    if False:
        "patch-shape"
        M = self.patch_matrix
        if another_direction:
            M = M.T
        vb = M[-1,1:-1]
    else:
        "more general, one boundary v"
        vb1,_ = self.get_i_boundary_vertex_indices(0) # i=0,1,2,3
        vb2,_ = self.get_i_boundary_vertex_indices(1) # i=0,1,2,3
        if len(np.intersect1d(np.r_[vl,vr],vb1))!=0:
            vb = vb1
        elif len(np.intersect1d(np.r_[vl,vr],vb2))!=0:
            vb = vb2
        else:
            "need to check and rewrite"
            if len(self.corner)!=0:
            #try:
                "patch-shape"
                M = self.patch_matrix
            else:
            #except:
                M = self.rot_patch_matrix
            
            if another_direction:
                M = M.T
            vb = M[-1,:]
        #vb = vb[1:-1]##[::interval] # Hui comment, may have problem later.
        if len(np.intersect1d(vr,vb))!=0:
            "make sure vl in vb, vr not in vb"
            vl,vr = vr,vl    
            
    alls_v0, alls_vs_v0 = [],[]
    select_v0,select_vs_v0 = [],[]
    allplv = get_isoline_between_2bdry(self,vb)
    for k in range(len(vb)):
        iv = allplv[k]
        iv0,jv0 = [],[]
        for i in iv:
            if i in v:
                j = np.where(v==i)[0][0]
                iv0.append(i)
                jv0.append(j)  
        if len(jv0)>1:#!=0:
            alls_v0.append(iv0)
            alls_vs_v0.append(jv0)
            if k % interval ==0:
                select_v0.append(iv0)
                select_vs_v0.append(jv0)
    return alls_v0,alls_vs_v0,select_v0,select_vs_v0  
###--------------------------------------------------------

def get_index_in_polyline(self, vstartlist, patch_or_annulus=True):
    """
    Retrieves the indices of vertices in a polyline.

    This function identifies and returns the indices of vertices in a polyline 
    based on the specified start vertices and mesh type.

    Parameters
    ----------
    vstartlist : list
        The list of start vertices.
    patch_or_annulus : bool, optional (default=True)
        Whether the mesh is a patch or an annulus.

    Returns
    -------
    indices : numpy array
        The indices of vertices in the polyline.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The indices are useful for accessing specific vertices in the mesh.

    reorder for a continuous but messy order vertex-list;

    vstartlist in boundary_curves

    See Also
    --------
    get_polyline_from_an_edge : Retrieves a polyline from an edge.
    """
    if patch_or_annulus: ##only work for these two cases. should be compatible with self.indM
        "also works for inner loops without mesh-bdry-pts in vstartlist"
        if len(self.corner)!=0:
            "patch shape: [0,1,2,3,4,5]"
            M = self.patch_matrix.T ##need check if .T
            #vlist = np.r_[M[0,:],M[-1,:]]
        else:
            "rotational shape: [0,1,2,3,4,5; 0]"
            M = self.rot_patch_matrix ##M[:,-1]==M[:,0]
            # vlist = np.r_[M[0,:-1],M[-1,:-1]]
            M = M[:,:-1]
        for vlist in M:
            if len(np.intersect1d(vlist, vstartlist))>1:
                "vstartlist include mesh-boundary-pts"
                order = []
                for v in vlist:
                    if v in vstartlist:
                        order.append(np.where(vstartlist==v)[0][0])
                break
            
        for vlist in M.T:
            if len(np.intersect1d(vlist, vstartlist))>1:
                "vstartlist include mesh-boundary-pts"
                order = []
                for v in vlist:
                    if v in vstartlist:
                        order.append(np.where(vstartlist==v)[0][0])
                break
        return np.array(order)    
    ### below works bad. no more use now.
    elif False:
        ### problem: edge happens at inner bdry, vlist is the whole bdry of the patch
        H = self.halfedges
        i = 0
        edge = None
        while edge is None:
            ei = np.where(H[:,0]==vstartlist[i])[0]
            for e in ei:
                if H[H[e,4],0] in vstartlist:
                    edge = e
                    break
            i += 1

        vlist,_ = self.get_polyline_from_an_edge(edge,is_bdry=True)##note, Vlist is full-poly, right!
    elif False:
        ### same problem as the above
        vlist = self.boundary_vertices()
    elif False:
        ### note, vlist=[1,2,3,4, 4,5,6,7,8, 8,9,10,11, ..., 1]
        ### all these three ways have the prelimarary that vstartlist in bdry
        brdies = self.boundary_curves(True)
        vlist = np.array([],dtype=int)
        for brr in brdies:
            vlist = np.r_[vlist,brr] ##has multiple vertices
            
    #print(vlist,vstartlist,len(np.intersect1d(vlist, vstartlist)))

def get_polylines_from_edges(self, es, is_poly=True):
    """
    Retrieves polylines from a list of edges.

    This function identifies and returns polylines based on the specified edges.

    Parameters
    ----------
    es : list
        The list of edge indices.
    is_poly : bool, optional (default=True)
        Whether to return the result as polylines.

    Returns
    -------
    ivs : numpy array
        The indices of vertices forming the polylines.
    V : numpy array
        The coordinates of vertices forming the polylines.
    polys : list
        The list of polylines (if is_poly=True).

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The polylines are useful for mesh analysis and processing.

    See Also
    --------
    get_polyline_from_an_edge : Retrieves a polyline from an edge.
    """
    V = self.vertices
    ivs = np.array([],dtype=int)
    ns = []
    for e in es:
        iv,VV = self.get_polyline_from_an_edge(e,is_halfedge=False)
        ivs = np.r_[ivs, iv]
        ns.append(len(VV)-1)
    ns = np.array(ns) 
    
    if is_poly:
        polys = make_multiple_polylines_from_endpoints(V[ivs], ns)
        return ivs, V[ivs], polys
    return ivs, V[ivs]

def get_polyline_from_an_edge(self, e, is_bdry=False, is_halfedge=True, is_poly=False):
    """
    Retrieves a polyline from a specified edge.

    This function identifies and returns a polyline based on the specified edge.

    Parameters
    ----------
    e : int
        The index of the edge.
    is_bdry : bool, optional (default=False)
        Whether the edge is on the boundary.
    is_halfedge : bool, optional (default=True)
        Whether to consider half-edges.
    is_poly : bool, optional (default=False)
        Whether to return the result as a polyline.

    Returns
    -------
    iv : numpy array
        The indices of vertices forming the polyline.
    VV : numpy array
        The coordinates of vertices forming the polyline.
    poly : Polyline object
        The polyline object (if is_poly=True).

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The polyline is useful for mesh analysis and processing.

    See Also
    --------
    get_polylines_from_edges : Retrieves polylines from a list of edges.
    """
    H = self.halfedges
    vl = vr = np.array([],dtype=int)
    if H[e,1]==-1:
        if is_bdry:
            "loop case"
            while H[e,0] not in vl:
                vl = np.r_[vl,H[e,0]]
                e = H[e,2]
            iv = vl
        else:
            e0 = e
            while H[H[e,4],0] not in self.corner:
                vr = np.r_[vr, H[H[e,4],0]]
                e = H[e,2]
            vr = np.r_[vr, H[H[e,4],0]]
            e = e0
            while H[H[e,3],0] not in self.corner:
                vl = np.r_[vl, H[H[e,3],0]]
                e = H[e,3]
            vl = np.r_[vl, H[H[e,3],0]]
            iv = np.r_[vl[::-1],H[e0,0],vr]
            #print('H[e,1]==-1',iv)
        VV = self.vertices[iv]
    elif H[H[e,4],1]==-1:
        if is_bdry:
            "loop case"
            e = H[e,4]
            while H[e,0] not in vl:
                vl = np.r_[vl,H[e,0]]
                e = H[e,2]
            iv = vl
        else:
            e0 = H[e,4]
            e = e0
            while H[H[e,4],0] not in self.corner:
                vr = np.r_[vr, H[H[e,4],0]]
                e = H[e,2]
            vr = np.r_[vr, H[H[e,4],0]]
            e = e0
            while H[H[e,3],0] not in self.corner:
                vl = np.r_[vl, H[H[e,3],0]]
                e = H[e,3]
            vl = np.r_[vl, H[H[e,3],0]]
            iv = np.r_[vl[::-1],H[e0,0],vr]
            #print('H[H[e,4],1]==-1',iv)
        VV = self.vertices[iv]
    else:
        if is_halfedge:
            il,ir = e,H[e,4]
        else:
            ##print(np.where(H[:,5]==e),H[e,0],H[H[e,4],0])
            il,ir = np.where(H[:,5]==e)[0] # should be two
        vl = np.r_[vl,H[il,0]]
        while H[H[H[il,2],4],1]!=-1 and H[H[H[H[H[il,2],4],2],4],1]!=-1:
            il = H[H[H[il,2],4],2]
            if H[il,0] in vl:
                break
            vl = np.r_[vl, H[il,0]]
        vl = np.r_[vl,H[H[il,4],0]]
        while H[H[H[ir,2],4],1]!=-1 and H[H[H[H[H[ir,2],4],2],4],1]!=-1:
            ir = H[H[H[ir,2],4],2]
            if  H[H[ir,4],0] in vr:
                break
            vr = np.r_[vr, H[H[ir,4],0]]     
        "Vl = self.vertices[vl[::-1]]; Vr = self.vertices[vr]"
        iv = np.r_[vl[::-1],vr]
        VV = self.vertices[iv]
        if is_poly:
            poly = make_polyline_from_endpoints(VV[:-1,:],VV[1:,:])
            return iv,VV,poly
        
        if False:#has problem
            """above shows a problem for the chebyshev-sphere mesh, 
                seem has extra1pt from bdry (in fact the first one), 
                remove such point
            """
            j = []
            for i, v in enumerate(iv):
                if False:
                    "work, [v0,vwrong,v1,..]"
                    if i<len(iv)-1:
                        e = None
                        ei = np.where(H[:,0]==v)[0]
                        ej = np.where(H[H[:,4],0]==iv[i+1])[0]
                        e = np.intersect1d(ei,ej)[0]
                        if e is None:
                            j.append(i)
                            #print('1',j)
                        if H[H[H[e,4],3],0] in iv:
                            j.append(i+1)
                            #print('2',j)
                        if H[H[e,3],0] in iv:
                            j.append(i+1)
                            # print('3',j)
                else:
                    if i<len(iv)-1 and i>0:
                        #print(iv[i-1],iv[i],iv[i+1])
                        e1 = np.where(H[:,0]==iv[i-1])[0]
                        e2 = np.where(H[:,0]==iv[i])[0]
                        e3 = np.where(H[:,0]==iv[i+1])[0]
                        #print(e1,e2,e3)
                        e12,e23 = None,None
                        e12 = np.intersect1d(H[e2,4],e1)[0]
                        e23 = np.intersect1d(H[e2,4],e3)[0]
                        #print(e12,e23)
                        if e12 is None or e23 is None:
                            j.append(i)
                            print('123',e12,e23,i)
            j = np.unique(j)
            #print(iv,j)
            if len(j) !=0:
                iv = np.delete(iv,j)
                VV = self.vertices[iv]   
            #print(iv)
    return iv,VV


def get_polylines_fair_index(self, es, vvv, vva, vvb):
    """
    Retrieves indices of fair polylines.

    This function identifies and returns the indices of fair polylines based 
    on the specified edges and vertex lists.

    Parameters
    ----------
    es : list
        The list of edge indices.
    vvv : list
        The list of vertex indices.
    vva : list
        The list of neighboring vertex indices.
    vvb : list
        The list of neighboring vertex indices.

    Returns
    -------
    indices : numpy array
        The indices of fair polylines.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The indices are useful for accessing specific polylines in the mesh.

    vvv=[v,v],vva=[v1,v2],vvb=[v3,v4]

    See Also
    --------
    get_polyline_from_an_edge : Retrieves a polyline from an edge.
    """
    #v,v1,v2,v3,v4 = self.ver_regular_star.T
    ivs,_,_ = self.get_polylines_from_edges(es)
    idvl = []
    for iv in ivs:
        iia = np.where(vva==iv)[0]
        if len(iia)!=0:
            for i in iia:
                if vvv[i] in ivs and vvb[i] in ivs:
                    idvl.append(i)
    return np.array(idvl,dtype=int)    
    
def get_quadface_diagonal(self):
    """
    Computes the diagonal vectors of quad faces.

    This function calculates the diagonal vectors of quad faces and returns 
    their lengths and unit vectors.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    ld1 : numpy array
        The lengths of the first diagonal.
    ld2 : numpy array
        The lengths of the second diagonal.
    ud1 : numpy array
        The unit vectors of the first diagonal.
    ud2 : numpy array
        The unit vectors of the second diagonal.

    Note
    ----
    This function assumes that the mesh is a regular mesh with quad faces.
    The diagonal vectors are useful for mesh analysis and processing.

    for isometry_checkerboard

    See Also
    --------
    get_quad_diagonal : Computes the diagonal vectors for a given set of vertices.
    """
    vi = self.quadface
    v1,v2,v3,v4 = vi[::4],vi[1::4],vi[2::4],vi[3::4]
    ld1 = np.linalg.norm(self.vertices[v1]-self.vertices[v3],axis=1)
    ld2 = np.linalg.norm(self.vertices[v2]-self.vertices[v4],axis=1)
    ud1 = (self.vertices[v1]-self.vertices[v3]) / ld1[:,None]
    ud2 = (self.vertices[v2]-self.vertices[v4]) / ld2[:,None]
    return ld1,ld2,ud1,ud2

def get_quad_diagonal(self, V, plot=False):
    """
    Computes the diagonal vectors of quad faces for a given set of vertices.

    This function calculates the diagonal vectors of quad faces and returns 
    their lengths and unit vectors. Optionally, it can plot the diagonals.

    Parameters
    ----------
    V : numpy array
        The vertex positions.
    plot : bool, optional (default=False)
        Whether to plot the diagonal vectors.

    Returns
    -------
    ld1 : numpy array
        The lengths of the first diagonal.
    ld2 : numpy array
        The lengths of the second diagonal.
    ud1 : numpy array
        The unit vectors of the first diagonal.
    ud2 : numpy array
        The unit vectors of the second diagonal.
    anchor : numpy array
        The anchor points for the diagonals.

    Note
    ----
    This function assumes that the mesh is a regular mesh with quad faces.
    The diagonal vectors are useful for mesh analysis and processing.

    isogonal_ck_based: get diagonal edge_length / unit_vector

    See Also
    --------
    get_quadface_diagonal : Computes the diagonal vectors of quad faces.
    """
    v1,v2,v3,v4 = self.rr_quadface.T # in odrder
    ld1 = np.linalg.norm(V[v3]-V[v1],axis=1)
    ld2 = np.linalg.norm(V[v4]-V[v2],axis=1)
    ud1 = (V[v3]-V[v1]) / ld1[:,None]
    ud2 = (V[v4]-V[v2]) / ld2[:,None]
    anchor = (V[v1]+V[v2]+V[v3]+V[v4]) / 4
    if plot:
        a,b = ld1/4.0, ld2/4.0
        Vr,Vl = anchor+ud1*a[:,None], anchor-ud1*a[:,None]
        pl1 = make_polyline_from_endpoints(Vl,Vr)
        Vr,Vl = anchor+ud2*b[:,None], anchor-ud2*b[:,None] 
        pl2 = make_polyline_from_endpoints(Vl,Vr)
        return pl1,pl2          
    return ld1,ld2,ud1,ud2,anchor

def get_quad_midpoint_cross_vectors(self, plot=False, scale=1/3):
    """
    Computes the midpoint cross vectors of quad faces.

    This function calculates the midpoint cross vectors of quad faces and 
    returns their lengths and unit vectors. Optionally, it can plot the vectors.

    Parameters
    ----------
    plot : bool, optional (default=False)
        Whether to plot the midpoint cross vectors.
    scale : float, optional (default=1/3)
        The scaling factor for the vectors.

    Returns
    -------
    lp : numpy array
        The lengths of the first cross vector.
    lq : numpy array
        The lengths of the second cross vector.
    m13 : numpy array
        The unit vectors of the first cross vector.
    m24 : numpy array
        The unit vectors of the second cross vector.
    anchor : numpy array
        The anchor points for the cross vectors.

    Note
    ----
    This function assumes that the mesh is a regular mesh with quad faces.
    The midpoint cross vectors are useful for mesh analysis and processing.

    isogonal_face_based: get quadface midpoint edge_length / unit_vector

    See Also
    --------
    get_quad_diagonal : Computes the diagonal vectors of quad faces.
    """
    v1,v2,v3,v4 = self.rr_quadface.T # in odrder
    e1 = self.vertices[v3]+self.vertices[v4]-self.vertices[v1]-self.vertices[v2]
    e2 = self.vertices[v1]+self.vertices[v4]-self.vertices[v2]-self.vertices[v3]
    ld1 = np.linalg.norm(e1,axis=1)
    ld2 = np.linalg.norm(e2,axis=1)
    ud1 = e1 / ld1[:,None]
    ud2 = e2 / ld2[:,None]
    anchor = (self.vertices[v1]+self.vertices[v2]+self.vertices[v3]+self.vertices[v4]) / 4
    if plot:
        a,b = ld1*scale, ld2*scale
        Vr,Vl = anchor+ud1*a[:,None], anchor-ud1*a[:,None]
        pl1 = make_polyline_from_endpoints(Vl,Vr)
        Vr,Vl = anchor+ud2*b[:,None], anchor-ud2*b[:,None] 
        pl2 = make_polyline_from_endpoints(Vl,Vr)
        return pl1,pl2            
    return ld1,ld2,ud1,ud2,anchor

def get_quad_parallelogram(self, is_theta_equal_theta0=False, is_plot=False):
    """
    Computes the properties of parallelogram faces in the mesh.

    This function calculates the properties of parallelogram faces, including 
    lengths, angles, and similarity metrics. Optionally, it can plot the results.

    Parameters
    ----------
    is_theta_equal_theta0 : bool, optional (default=False)
        Whether to assume equal angles.
    is_plot : bool, optional (default=False)
        Whether to plot the results.

    Returns
    -------
    properties : numpy array
        The computed properties of the parallelogram faces.

    Note
    ----
    This function assumes that the mesh is a regular mesh with parallelogram faces.
    The properties are useful for mesh analysis and processing.

    I-net with similar parallelograms

    See Also
    --------
    get_quad_parallelogram_similarity : Computes the similarity of parallelogram faces.
    """
    la,lb,d1,d2,_ = self.get_quad_diagonal(self.vertices)
    lp,lq,m13,m24,anchor = self.get_quad_midpoint_cross_vectors()
    La,Lb,Lp,Lq = la**2,lb**2,lp**2,lq**2
    c1,c2 = Lb/La, Lq/Lp
    lamda,mu = np.mean(np.sqrt(c1)),np.mean(np.sqrt(c2))
    delt1 = np.sqrt(mu-1)
    delt2 = lamda*mu+lamda-mu+1
    delt3 = -lamda*mu+lamda+mu+1
    A = np.mean(c1)
    delt4 = np.sqrt(np.abs(A-np.sqrt(2)+1))
    delt5 = np.sqrt(np.abs(np.sqrt(2)+1-A))
    
    if is_plot:
        v1,v2,_,_ = self.rr_quadface.T
        an1 = (self.vertices[v1]+self.vertices[v2])*0.5
        return [an1,d1,d2, anchor, m13,m24] 
        
    if is_theta_equal_theta0:
        return np.r_[La,Lb,A,delt4,delt5]
    else:
        return np.r_[La,Lb,Lp,Lq,c1,c2,lamda,mu,delt1,delt2,delt3]
    
def get_quad_parallelogram_similarity(self):
    """
    Computes the similarity of parallelogram faces in the mesh.

    This function calculates the similarity metrics of parallelogram faces, 
    including lambda and mu values.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    lambda : numpy array
        The lambda values for the parallelogram faces.
    mu : numpy array
        The mu values for the parallelogram faces.

    Note
    ----
    This function assumes that the mesh is a regular mesh with parallelogram faces.
    The similarity metrics are useful for mesh analysis and processing.

    See Also
    --------
    get_quad_parallelogram : Computes the properties of parallelogram faces.
    """
    V = self.vertices
    v1,v2,v3,v4 = self.rr_quadface.T # in odrder
    la = np.linalg.norm(V[v1]-V[v3],axis=1)
    lb = np.linalg.norm(V[v2]-V[v4],axis=1)
    lp = np.linalg.norm(V[v3]+V[v4]-V[v1]-V[v2],axis=1)
    lq = np.linalg.norm(V[v1]+V[v4]-V[v2]-V[v3],axis=1)
    lamda,mu = lb/la, lq/lp
    return lamda, mu
    

def get_quad_midline(self):
    """
    Computes the midlines of quad faces in the mesh.

    This function calculates the midlines of quad faces and returns the 
    corresponding polylines.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    pl1 : Polyline object
        The first midline polyline.
    pl2 : Polyline object
        The second midline polyline.

    Note
    ----
    This function assumes that the mesh is a regular mesh with quad faces.
    The midlines are useful for mesh analysis and processing.

    See Also
    --------
    get_quad_diagonal : Computes the diagonal vectors of quad faces.
    """
    V = self.vertices
    v1,v2,v3,v4 = self.rr_quadface.T
    V12 = (V[v1]+V[v2])/2.0
    V23 = (V[v2]+V[v3])/2.0
    V34 = (V[v3]+V[v4])/2.0
    V41 = (V[v4]+V[v1])/2.0
    pl1 = make_polyline_from_endpoints(V12,V34)
    pl2 = make_polyline_from_endpoints(V23,V41)
    return pl1,pl2   

def quadface_neighbour_star(self, boundary=True):
    """
    Computes the neighbor star matrix for quad faces.

    This function calculates the neighbor star matrix for quad faces, 
    representing the connectivity of adjacent faces.

    Parameters
    ----------
    boundary : bool, optional (default=True)
        Whether to include boundary faces.

    Returns
    -------
    tristar : numpy array
        The neighbor star matrix for quad faces.

    Note
    ----
    This function assumes that the mesh is a regular mesh with quad faces.
    The neighbor star matrix is useful for mesh analysis and processing.

    rr_quadface_tristar: [f,fi,fj]; get the star metrix of forde/rr_quadface

    f in forder, its neighbour f1234 also in forder; 

    get the corresponding order of forder

    See Also
    --------
    get_quadface_diagonal : Computes the diagonal vectors of quad faces.
    """
    H = self.halfedges
    forder = self.rr_quadface_order
    farr = self.rr_quadface
    numf = self.num_rrf
    id00,id01,idl,idr,idu,idd=[],[],[],[],[],[]
    if boundary:
        iddd,iddl,iddr = [],[],[]
    if0,if1,if2,if3,if4 = [],[],[],[],[]
    for i in range(numf):
        first = farr[i]
        e1=np.where(H[:,0]==first[0])[0]
        e2=np.where(H[H[:,4],0]==first[1])[0]
        e = np.intersect1d(e1,e2)
        if H[e,1]==-1:
                e = H[e,4]
        e1 = H[H[e,3],4] # left
        e2 = H[e,4] # down
        e3 = H[H[e,2],4] #right
        e4 = H[H[H[e,2],2],4] #up
        if H[e1,1]!=-1 and H[e3,1]!=-1 and H[e1,1] in forder and H[e3,1] in forder:
            id00.append(i)
            idl.append(np.where(forder==H[e1,1])[0][0])
            idr.append(np.where(forder==H[e3,1])[0][0])
        if H[e2,1]!=-1 and H[e4,1]!=-1 and H[e2,1] in forder and H[e4,1] in forder:
            id01.append(i)
            idd.append(np.where(forder==H[e2,1])[0][0])
            idu.append(np.where(forder==H[e4,1])[0][0])
            
        if H[e1,1]!=-1 and H[e3,1]!=-1 and H[e1,1] in forder and H[e3,1] in forder:
            if H[e2,1]!=-1 and H[e4,1]!=-1 and H[e2,1] in forder and H[e4,1] in forder:
                if0.append(i)
                if1.append(np.where(forder==H[e1,1])[0][0])
                if3.append(np.where(forder==H[e3,1])[0][0])
                if2.append(np.where(forder==H[e2,1])[0][0])
                if4.append(np.where(forder==H[e4,1])[0][0])  
            
        if boundary:
            "edge e has four positon: up,down,left,right"
            fd,fu = H[H[e,4],1], H[H[H[H[e,2],2],4],1]
            fl,fr = H[H[H[e,3],4],1], H[H[H[e,2],4],1]
            if fd==-1 and (fl!=-1 or fr!=-1):
                try:
                    a = np.where(forder==fl)[0][0]
                    b = np.where(forder==fr)[0][0]
                    iddd.append(i)
                    iddl.append(a)
                    iddr.append(b)
                except:
                    pass
            if fu==-1 and (fl!=-1 or fr!=-1):
                try:
                    a = np.where(forder==fl)[0][0]
                    b = np.where(forder==fr)[0][0]
                    iddd.append(i)
                    iddl.append(b)
                    iddr.append(a)
                except:
                    pass
            if fl==-1 and (fd!=-1 or fu!=-1):
                try:
                    a = np.where(forder==fd)[0][0]
                    b = np.where(forder==fu)[0][0]
                    iddd.append(i)
                    iddl.append(b)
                    iddr.append(a)
                except:
                    pass
            if fr==-1 and (fd!=-1 or fu!=-1):
                try:
                    a = np.where(forder==fd)[0][0]
                    b = np.where(forder==fu)[0][0]
                    iddd.append(i)
                    iddl.append(a)
                    iddr.append(b)
                except:
                    pass     
    if boundary:
        ###print(iddd,iddl,iddr)
        tristar = np.vstack((np.r_[id00,id01,iddd],np.r_[idl,idu,iddl],np.r_[idr,idd,iddr])).T
    else:
        tristar = np.vstack((np.r_[id00,id01],np.r_[idl,idu],np.r_[idr,idd])).T
    self._rr_quadface_tristar = tristar
    self._rr_quadface_4neib = np.c_[if0,if1,if2,if3,if4]

# -------------------------------------------------------------------------
#                       all faces / vertices
# -------------------------------------------------------------------------
def get_checker_select_vertex(self, first=[0]):
    """
    Retrieves the checkerboard selection of vertices.

    This function identifies and returns the indices of vertices in a 
    checkerboard pattern, starting from the specified vertex.

    Parameters
    ----------
    first : list, optional (default=[0])
        The starting vertex index.

    Returns
    -------
    indices : list
        The indices of vertices in the checkerboard pattern.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The checkerboard selection is useful for mesh analysis and processing.

    self.checker_vertex: all vertex...blue ~ 0 ; red ~ 1; left ~ -1

    See Also
    --------
    get_checkerboard_black : Retrieves the black vertices in a checkerboard pattern.
    get_checkerboard_white : Retrieves the white vertices in a checkerboard pattern.
    """
    vall,vallij=self.vertex_ring_vertices_iterators(sort=True, order=False)
    #vall,vallij = self._vi,self._vj
    arr = -np.ones(self.V)
    left = np.arange(self.V)
    i0 = first       ## f0
    arr[i0] = 0      ## blue
    ind = np.where(vall==i0)[0]
    ij = vallij[ind] ## fneib
    left = np.setdiff1d(left,np.r_[i0,ij])
    a = 0
    while len(left)!=0:
        a = a % 2
        vneib=[]
        for i in ij:
            if i !=-1:
                arr[i]=1-a
                ind = np.where(vall==i)[0]
                vneib.extend(list(vallij[ind]))
        left = np.setdiff1d(left,np.r_[ij,vneib])
        ij=[]
        for v in vneib:
            if v!=-1 and arr[v]==-1:
                arr[v]=a
                ij.append(v)
        a += 1
    blue = np.where(arr==0)[0]
    red = np.where(arr==1)[0]
    self._checker_vertex = [blue,red]
    
# -------------------------------------------------------------------------
#                        Diagonal mesh
# -------------------------------------------------------------------------
def make_quad_mesh_from_indices(self, V, v1, v2, v3, v4):
    """
    Creates a quad mesh from the given vertex indices.

    This function constructs a quad mesh using the specified vertex indices and 
    their corresponding positions.

    Parameters
    ----------
    V : numpy array
        The vertex positions.
    v1, v2, v3, v4 : numpy array
        The indices of vertices forming the quad faces.

    Returns
    -------
    quad_mesh : Mesh object
        The constructed quad mesh.

    Note
    ----
    This function assumes that the input indices form valid quad faces.
    The resulting quad mesh is useful for further analysis and processing.

    See Also
    --------
    get_diagonal_mesh : Creates a diagonal mesh from the current mesh.
    """
    vrr = np.unique(np.r_[v1,v2,v3,v4])
    numf = len(v1)
    vlist = V[vrr]
    farr = np.array([],dtype=int)
    for i in range(numf):
        i1 = np.where(vrr==v1[i])[0]
        i2 = np.where(vrr==v2[i])[0]
        i3 = np.where(vrr==v3[i])[0]
        i4 = np.where(vrr==v4[i])[0]
        farr = np.r_[farr,i1,i2,i3,i4]
    flist = farr.reshape(-1,4).tolist()
    ck = Mesh()
    ck.make_mesh(vlist,flist)
    return ck

def get_diagonal_mesh(self, sort=True, blue=False, whole=False):
    """
    Creates a diagonal mesh from the current mesh.

    This function constructs a diagonal mesh by connecting the vertices of the 
    current mesh in a specific pattern. The resulting mesh can be sorted and 
    colored based on the specified parameters.

    Parameters
    ----------
    sort : bool, optional (default=True)
        Whether to sort the vertices.
    blue : bool, optional (default=False)
        Whether to color the vertices blue.
    whole : bool, optional (default=False)
        Whether to include the whole mesh.

    Returns
    -------
    diagonal_mesh : Mesh object
        The constructed diagonal mesh.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The diagonal mesh is useful for mesh analysis and processing.

    See Also
    --------
    make_quad_mesh_from_indices : Creates a quad mesh from vertex indices.
    """
    V = self.vertices
    if whole:
        num = self.V
        bary = self.face_barycenters()
        dV = np.vstack((V,bary))
        H = self.halfedges
        i1 = np.where(H[:,1] >= 0)[0]
        i2 = np.where(H[H[:,4],1] >= 0)[0]
        i = np.intersect1d(i1,i2)
        e = np.array([i[0]])
        for j in i[1:]:
            if H[H[j,4],1] not in H[e,1]:
                e = np.r_[e,j]
        v1, v3 = H[e,0], H[H[e,4],0]
        v2, v4 = num+H[H[e,4],1], num+H[e,1]
        dallv = np.unique(np.r_[v1,v2,v3,v4])
        vlist = dV[dallv]
        iv1 = [np.argwhere(dallv == item)[0][0] for item in v1]
        iv2 = [np.argwhere(dallv == item)[0][0] for item in v2]
        iv3 = [np.argwhere(dallv == item)[0][0] for item in v3]
        iv4 = [np.argwhere(dallv == item)[0][0] for item in v4]
        flist = (np.array([iv1,iv2,iv3,iv4]).T).tolist()
        dmesh = Mesh()
        dmesh.make_mesh(vlist,flist)
        return dmesh
    else:
        v,v1,v2,v3,v4 = self.ver_regular_star.T
        if sort:
            def _ind_red_regular_vertex(blue):
                vblue,vred = self.checker_vertex # depends on all vertex-checker
                order = self.ver_regular
                idb,idr = [],[]
                for i in range(len(order)):
                    v = order[i]
                    if v in vblue:
                        idb.append(i)
                    elif v in vred:
                        idr.append(i)
                if blue:
                    return np.array(idb)
                return np.array(idr)
            order = _ind_red_regular_vertex(blue)
        else:
            order = np.arange(self.num_regular)[::2]
        # V1,V2,V3,V4 = V[v1[order]],V[v2[order]],V[v3[order]],V[v4[order]]
        # dmesh = make_quad_mesh_pieces(V1,V2,V3,V4)
        v1,v2,v3,v4 = v1[order],v2[order],v3[order],v4[order]
        dmesh = self.make_quad_mesh_from_indices(V,v1,v2,v3,v4)
        return dmesh

def get_midline_mesh(self):
    """
    Creates a midline mesh from the current mesh.

    This function constructs a midline mesh by connecting the midpoints of the 
    edges of the current mesh. The resulting mesh represents the midlines of the 
    original mesh.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    midline_mesh : Mesh object
        The constructed midline mesh.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The midline mesh is useful for mesh analysis and processing.

    (remesh of principal mesh by midline-ortho + planarity; the quadface should be planar)

    See Also
    --------
    get_diagonal_mesh : Creates a diagonal mesh from the current mesh.
    """
    V = self.vertices
    vi = self.quadface
    v1,v2,v3,v4 = vi[::4],vi[1::4],vi[2::4],vi[3::4]
    bary = (V[v1]+V[v2]+V[v3]+V[v4])/4.0
    arr = self.V + np.arange(len(v1))
    vlist = np.vstack((V,bary))
    
    edges,edge12,edge23,edge34,edge41 = [],[],[],[],[]
    H = self.halfedges
    for i in range(len(v1)):
        e1i = np.where(H[:,0]== v1[i])[0]
        e1j = np.where(H[H[:,4],0]== v2[i])[0]
        e = np.intersect1d(e1i,e1j)[0]
        e12, e23 = H[e,5],H[H[e,2],5]
        e34, e41 = H[H[H[e,2],2],5], H[H[e,3],5]
        edge12.append(e12)
        edge23.append(e23)
        edge34.append(e34)
        edge41.append(e41)
        edges.append(e12)
        edges.append(e23)
        edges.append(e34)
        edges.append(e41)
    
    edges = np.unique(edges)
    for j in edges:
        jl,jr = np.where(H[:,5]==j)[0]
        Vm = (V[H[jl,0]] + V[H[jr,0]])/2
        vlist = np.vstack((vlist,Vm))
    
    i12,i23,i34,i41 = [],[],[],[]
    for i in range(len(v1)):
        i12.append(np.where(edges==edge12[i])[0])
        i23.append(np.where(edges==edge23[i])[0])
        i34.append(np.where(edges==edge34[i])[0])
        i41.append(np.where(edges==edge41[i])[0])
    m1 = np.array(i12)+ self.V + len(v1)
    m2 = np.array(i23)+ self.V + len(v1)
    m3 = np.array(i34)+ self.V + len(v1)
    m4 = np.array(i41)+ self.V + len(v1)
    f1 = np.tile(arr,4)
    f2 = np.r_[m4,m1,m2,m3]
    f3 = np.r_[v1,v2,v3,v4]
    f4 = np.r_[m1,m2,m3,m4]
    flist = (np.c_[f1,f2,f3,f4]).tolist()
    mmesh = Mesh()
    mmesh.make_mesh(vlist,flist)
    return mmesh

def get_checkerboard_black(self, is_rr=False, is_ck_n=False, is_GI=False):
    """
    Retrieves the black vertices in a checkerboard pattern.

    This function identifies and returns the black vertices in a checkerboard 
    pattern, optionally considering regular vertices and normal vectors.

    Parameters
    ----------
    is_rr : bool, optional (default=False)
        Whether to consider regular vertices.
    is_ck_n : bool, optional (default=False)
        Whether to return normal vectors.
    is_GI : bool, optional (default=False)
        Whether to compute the gradient of the normal vectors.

    Returns
    -------
    black_vertices : numpy array
        The indices of black vertices in the checkerboard pattern.
    normals : numpy array, optional
        The normal vectors of the black vertices (if is_ck_n=True).
    gradient : numpy array, optional
        The gradient of the normal vectors (if is_GI=True).

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The black vertices are useful for mesh analysis and processing.

    See Also
    --------
    get_checker_select_vertex : Retrieves the checkerboard selection of vertices.
    get_checkerboard_white : Retrieves the white vertices in a checkerboard pattern.
    """
    if is_rr:
        "oriented order of the vertices of quadfaces"
        order = self.rr_quadface_order
        v1,v2,v3,v4 = self.rr_quadface.T # in odrder
    else:
        vi = self.quadface
        v1,v2,v3,v4 = vi[::4],vi[1::4],vi[2::4],vi[3::4]
        order = self.quadface_order
    
    V = self.vertices
    
    if is_ck_n:
        an = (V[v1]+V[v2]+V[v3]+V[v4])/4.0
        fn = np.cross(V[v1]-V[v3], V[v2]-V[v4])
        fn = fn / np.linalg.norm(fn, axis=1)[:,None]
        return an, fn
    if is_GI:
        "below works the same as self.get_face_normal_GI"
        all_fn = self.face_normals()
        
        fn = np.cross(V[v1]-V[v3], V[v2]-V[v4])
        fn = fn / np.linalg.norm(fn, axis=1)[:,None]

        order_fn = all_fn[order]
        dot = np.einsum('ij,ij->i',order_fn, fn)
        idd = np.where(dot<0)[0]
        fn[idd] = -fn[idd]
        
        all_fn[order] = fn
        
        GI = get_barycenters_mesh(self,verlist=all_fn)  
        return GI
    
    P1 = (V[v1] + V[v2]) * 0.5
    P2 = (V[v2] + V[v3]) * 0.5
    P3 = (V[v3] + V[v4]) * 0.5
    P4 = (V[v4] + V[v1]) * 0.5
    ck = make_quad_mesh_pieces(P1,P2,P3,P4)
    return ck

def get_checkerboard_white(self):
    """
    Retrieves the white vertices in a checkerboard pattern.

    This function identifies and returns the white vertices in a checkerboard 
    pattern.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    white_vertices : numpy array
        The indices of white vertices in the checkerboard pattern.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The white vertices are useful for mesh analysis and processing.

    See Also
    --------
    get_checker_select_vertex : Retrieves the checkerboard selection of vertices.
    get_checkerboard_black : Retrieves the black vertices in a checkerboard pattern.
    """
    V = self.vertices
    v,v1,v2,v3,v4 = self.ver_regular_star.T
    P1 = (V[v1] + V[v]) * 0.5
    P2 = (V[v2] + V[v]) * 0.5
    P3 = (V[v3] + V[v]) * 0.5
    P4 = (V[v4] + V[v]) * 0.5
    ck = make_quad_mesh_pieces(P1,P2,P3,P4)
    return ck
    
#--------------------------------------------------------------------------
#                 Plot polylines: isolines + diagonals
#-------------------------------------------------------------------------- 

def get_quad_mesh_1family_isoline(self, diagnet=False, direction=True, edge=False):
    """
    Retrieves one family of isolines from the quad mesh.

    This function identifies and returns one family of isolines from the quad mesh, 
    optionally considering diagonal edges and edge data.

    Parameters
    ----------
    diagnet : bool, optional (default=False)
        Whether to consider diagonal edges.
    direction : bool, optional (default=True)
        The direction of the isolines.
    edge : bool, optional (default=False)
        Whether to return edge data.

    Returns
    -------
    isolines : Polyline object
        The isolines in the specified direction.
    edge_data : numpy array, optional
        The edge data (if edge=True).

    Note
    ----
    This function assumes that the mesh is a quad mesh.
    The isolines are useful for mesh analysis and visualization.

    See Also
    --------
    get_quad_mesh_2family_isolines : Retrieves two families of isolines.
    """    
    V = self.vertices
    v1,v2 = self.get_1family_oriented_polyline(diagnet,poly2=direction)
    pl = make_polyline_from_endpoints(V[v1],V[v2])
    if edge:
        "edge_data = np.zeros(self.E)"
        H = self.halfedges
        e,ib,eb = [],[],[]
        for i in range(len(v1)):
            a,b = v1[i],v2[i]
            j = np.where(H[:,0]==a)[0]
            k = np.where(H[H[:,4],0]==b)[0]
            m = np.intersect1d(j,k)[0]
            e.append(H[m,5])
            if H[m,1]==-1 or H[H[m,4],1]==-1:
                ib.append(i)
                eb.append(H[m,5])
        return pl, np.array(e),[np.array(ib),np.array(eb)]
    else:
        an,vec = V[v1], V[v2]-V[v1]
        return pl,an,vec 

def get_1family_oriented_polyline(self, diagnet=False, poly2=True):
    """
    Retrieves one family of oriented polylines from the mesh.

    This function identifies and returns one family of oriented polylines, 
    optionally considering diagonal edges and polyline direction.

    Parameters
    ----------
    diagnet : bool, optional (default=False)
        Whether to consider diagonal edges.
    poly2 : bool, optional (default=True)
        The direction of the polylines.

    Returns
    -------
    polyline : numpy array
        The oriented polyline.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The oriented polylines are useful for mesh analysis and processing.

    (!still have problem for the demultiple and oriendted quad faces,bad for thickness)

    See Also
    --------
    get_2families_polyline_from_1closed_bdry : Retrieves two families of polylines.
    """
    if diagnet:
        v,v1,v2,v3,v4 = self.rr_star_corner# in diagonal direction
    else:
        v,v1,v2,v3,v4 = self.ver_star_matrix.T
    num = len(v)
    if poly2:
        vl,vr = v2,v4
    else:
        vl,vr = v1,v3
    va,vb = vl,v
    for i in range(num):
        if v[i] not in np.r_[vl,v] or vr[i] not in np.r_[vl,v]:
            va = np.r_[va,v[i]]
            vb = np.r_[vb,vr[i]]
        else:
            i1 = np.where(va==v[i])[0]
            i2 = np.where(vb==vr[i])[0]
            if len(np.intersect1d(i1,i2))==0:
                va = np.r_[va,v[i]]
                vb = np.r_[vb,vr[i]]
    if demultiple:
        "remove multiple edges:"
        ind_del = []
        ck = []
        for i in range(len(va)):
            i1 = np.where(vb==va[i])[0]
            if len(i1) !=0:
                for j in i1:
                    if j not in ck and va[j] == vb[i]:
                        ind_del.append(i)
                        #print(i,va[i],vb[i],va[j],vb[j])
                        ck.append(i)
        if len(ind_del)!=0:
            ind = np.array(ind_del)   
            va = np.delete(va,ind)
            vb = np.delete(vb,ind)  
    return va,vb
    
def get_2families_polyline_from_1closed_bdry(self, diagnet=False, direction=True):
    """
    Retrieves two families of polylines from a closed boundary.

    This function identifies and returns two families of polylines from a closed 
    boundary, optionally considering diagonal edges and direction.

    Parameters
    ----------
    diagnet : bool, optional (default=False)
        Whether to consider diagonal edges.
    direction : bool, optional (default=True)
        The direction of the polylines.

    Returns
    -------
    polyline1 : numpy array
        The first family of polylines.
    polyline2 : numpy array
        The second family of polylines.

    Note
    ----
    This function assumes that the mesh has a closed boundary.
    The polylines are useful for mesh analysis and processing.

    (along one bdry for a patch-shape; two bdry for a star-shape)

    See Also
    --------
    get_1family_oriented_polyline : Retrieves one family of oriented polylines.
    """
    if diag:
        v0,v1,v2,v3,v4 = self.rr_star_corner
    else:
        v0,v1,v2,v3,v4 = self.rrv4f4
    v13,v24 = np.unique(np.r_[v1,v3]),np.unique(np.r_[v2,v4])
    
    "more general, one closed-boundary v"
    vb,_ = self.get_a_closed_boundary()
    
    H = self.halfedges
    "boundary are divided/filter into 2class based on v13 or v24 direction"
    vb1,vb2 = [],[]
    vfmly1,vfmly2 = [],[]
    i1l=i1r=i2l=i2r = np.array([],dtype=int)
    for i, v in enumerate(vb):
        if i%interval==0:
            if v in v13:
                "filter to 1st polyline-bdry-vertices"
                vb1.append(v)
                if v in v1:
                    j = np.where(v1==v)[0]
                elif v in v3:
                    j = np.where(v3==v)[0]
                vx = v0[j]
                if diag:
                    vvx = np.r_[v,vx]
                    vpl,_ = get_diagonal_polyline_from_2points(self,vvx,is_poly=False)
                else:
                    e = np.intersect1d(np.where(H[:,0]==v)[0],np.where(H[H[:,4],0]==vx)[0])[0]
                    vpl,_ = self.get_polyline_from_an_edge(e)
                if inner:
                    vpl = vpl[1:-1]
                vfmly1.append(vpl)
                i1l = np.r_[i1l,vpl[:-1]]
                i1r = np.r_[i1r,vpl[1:]]
            if v in v24:
                "filter to 2nd polyline-bdry-vertices"
                vb2.append(v)
                if v in v2:
                    j = np.where(v2==v)[0]
                elif v in v4:
                    j = np.where(v4==v)[0]
                vx = v0[j]
                if diag:
                    vvx = np.r_[v,vx]
                    vpl,_ = get_diagonal_polyline_from_2points(self,vvx,is_poly=False)
                else:
                    e = np.intersect1d(np.where(H[:,0]==v)[0],np.where(H[H[:,4],0]==vx)[0])[0]
                    vpl,_ = self.get_polyline_from_an_edge(e)
                if inner:
                    vpl = vpl[1:-1]
                vfmly2.append(vpl)
                i2l = np.r_[i2l,vpl[:-1]]
                i2r = np.r_[i2r,vpl[1:]]
    if is_poly:
        V = self.vertices
        pl1 = make_polyline_from_endpoints(V[i1l],V[i1r])
        pl2 = make_polyline_from_endpoints(V[i2l],V[i2r])
        return pl1,pl2
    if is_seg:
        return [i1l,i1r],[i2l,i2r]
    return [vb1,vfmly1],[vb2,vfmly2]

def get_i_boundary_vertices(self, e, by_closed=False, by_corner2=False):
    """
    Retrieves the boundary vertices from a specified edge.

    This function identifies and returns the boundary vertices starting from 
    the specified edge, optionally considering closed boundaries or corner vertices.

    Parameters
    ----------
    e : int
        The index of the starting edge.
    by_closed : bool, optional (default=False)
        Whether to consider closed boundaries.
    by_corner2 : bool, optional (default=False)
        Whether to consider corner vertices.

    Returns
    -------
    boundary_vertices : numpy array
        The indices of boundary vertices.
    boundary_positions : numpy array
        The positions of boundary vertices.

    Note
    ----
    This function assumes that the mesh has a boundary.
    The boundary vertices are useful for mesh analysis and processing.

    See Also
    --------
    get_i_boundary_vertex_indices : Retrieves the indices of boundary vertices.
    """
    H = self.halfedges
    il,ir = np.where(H[:,5]==e)[0] # should be two
    if H[il,1]==-1:
        e = il
    else:
        e = ir
    if by_closed:
        "no corner vertex of valence 2, which forms closed boundry"
        eb = []
        e1=e
        e2 = H[e,2]
        while H[e1,0] not in H[np.array(eb,dtype=int),0]:
            eb.append(e1)
            e1 = H[e1,3]
        #eb.append(e1)
        while H[e2,0] not in H[np.array(eb,dtype=int),0]:
            eb.append(e2)
            e2 = H[e2,2]
        #eb.append(e2)
        vb = H[np.array(eb,dtype=int),0]
        "above should be no duplicate, then below no use"
        #u, ind = np.unique(vb, return_index=True)
        #vb = u[np.argsort(ind)]
        return vb, self.vertices[vb]
    elif by_corner2: 
        "split by corner vertex of valence 2"
        eb = []
        e1=e
        e2 = H[e,2]
        _,_,lj = self.vertex_ring_vertices_iterators(sort=True,return_lengths=True)
        corner = np.where(lj==2)[0]
        while H[e1,0] not in corner:
            eb.append(e1)
            e1 = H[e1,3]
        eb.append(e1)
        while H[e2,0] not in corner:
            eb.append(e2)
            e2 = H[e2,2]
        eb.append(e2)
        vb = H[np.array(eb,dtype=int),0]
        return vb, self.vertices[vb]
    else:
        "all boundary vertices"
        vb = self.boundary_curves()
        for v in vb:
            if H[e,0] in v:
                return v, self.vertices[v]

def get_i_boundary_vertex_indices(self, i):
    """
    Retrieves the indices of boundary vertices for a specified boundary.

    This function identifies and returns the indices of boundary vertices 
    for the specified boundary.

    Parameters
    ----------
    i : int
        The index of the boundary.

    Returns
    -------
    boundary_vertices : numpy array
        The indices of boundary vertices.
    boundary_positions : numpy array
        The positions of boundary vertices.

    Note
    ----
    This function assumes that the mesh has multiple boundaries.
    The boundary vertices are useful for mesh analysis and processing.
    (work for rectangular-patch + cylinder-annulus)

    See Also
    --------
    get_i_boundary_vertices : Retrieves the boundary vertices from a specified edge.
    """
    #i = self.show_i_boundary
    _,_,lj = self.vertex_ring_vertices_iterators(sort=True,return_lengths=True)
    vcorner = np.where(lj==2)[0]
    H = self.halfedges
    ec = []
    if len(vcorner)!=0:
        for v in vcorner:
            e = np.intersect1d(np.where(H[:,0]==v)[0],np.where(H[:,1]==-1)[0])[0]
            ec.append(e)
        ec = H[np.array(ec,dtype=int),5]
        v,B = self.get_i_boundary_vertices(ec[i],by_corner2=True) 
    else:
        "rotational-case"
        i = int(i%2) # i always equal 0 or 1
        bi = self.boundary_curves(corner_split=True)
        for b in bi:
            v = b[0]
            e = np.intersect1d(np.where(H[:,0]==v)[0],np.where(H[:,1]==-1)[0])[0]
            ec.append(e)
        ec = H[np.array(ec,dtype=int),5]     
        v,B = self.get_i_boundary_vertices(ec[i],by_closed=True)         
    return v,B

def get_all_boundary_vertices(self, order=False):
    """
    Retrieves all boundary vertices in the mesh.

    This function identifies and returns all boundary vertices in the mesh, 
    optionally in a specific order.

    Parameters
    ----------
    order : bool, optional (default=False)
        Whether to return the boundary vertices in a specific order.

    Returns
    -------
    boundary_vertices : numpy array
        The indices of boundary vertices.
    boundary_positions : numpy array
        The positions of boundary vertices.

    Note
    ----
    This function assumes that the mesh has a boundary.
    The boundary vertices are useful for mesh analysis and processing.

    See Also
    --------
    boundary_vertices : Retrieves the indices of boundary vertices.
    """
    v = self.boundary_vertices()
    if order:
        # v = np.array([],dtype=int)
        v = self.boundary_curves(corner_split=False)[0]
        # for i in vi:
        #     v = np.r_[v,i]
    B = self.vertices[v]
    return v,B

# -------------------------------------------------------------------------
#                        For ploting
# -------------------------------------------------------------------------
def get_v4_unit_edge(self, rregular=True):
    """
    Computes the unit edge vectors for vertices with valence 4.

    This function calculates the unit edge vectors for vertices with valence 4, 
    optionally considering regular vertices.

    Parameters
    ----------
    rregular : bool, optional (default=True)
        Whether to consider regular vertices.

    Returns
    -------
    v : numpy array
        The indices of vertices.
    l1, l2, l3, l4 : numpy array
        The lengths of the edge vectors.
    e1, e2, e3, e4 : numpy array
        The unit edge vectors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The unit edge vectors are useful for mesh analysis and processing.

    See Also
    --------
    get_v4_unit_tangents : Computes the unit tangent vectors for vertices with valence 4.
    """
    V = self.vertices
    if rregular:
        v,v1,v2,v3,v4 = self.rrv4f4
    else:
        ## v,v1,v2,v3,v4 = self.ver_regular_star.T
        v,v1,v2,v3,v4 = self.ver_star_matrix.T
    E1 = V[v1]-V[v]
    E2 = V[v2]-V[v]
    E3 = V[v3]-V[v]
    E4 = V[v4]-V[v]
    l1 = np.linalg.norm(E1, axis=1)
    l2 = np.linalg.norm(E2, axis=1)
    l3 = np.linalg.norm(E3, axis=1)
    l4 = np.linalg.norm(E4, axis=1)
    e1 = E1 / l1[:,None]
    e2 = E2 / l2[:,None]
    e3 = E3 / l3[:,None]
    e4 = E4 / l4[:,None]
    return v,l1,l2,l3,l4,e1,e2,e3,e4

def get_v4_unit_tangents(self, plot=False, rregular=True):
    """
    Computes the unit tangent vectors for vertices with valence 4.

    This function calculates the unit tangent vectors for vertices with valence 4, 
    optionally plotting the results and considering regular vertices.

    Parameters
    ----------
    plot : bool, optional (default=False)
        Whether to plot the tangent vectors.
    rregular : bool, optional (default=True)
        Whether to consider regular vertices.

    Returns
    -------
    lt1, lt2 : numpy array
        The lengths of the tangent vectors.
    ut1, ut2 : numpy array
        The unit tangent vectors.
    anchor : numpy array
        The anchor points for the tangent vectors.
    angle : numpy array
        The angles between the tangent vectors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The unit tangent vectors are useful for mesh analysis and processing.

    (only for vertices of valence 4, not depends on X)

    See Also
    --------
    get_v4_unit_edge : Computes the unit edge vectors for vertices with valence 4.
    """
    v,l1,l2,l3,l4,e1,e2,e3,e4 = self.get_v4_unit_edge(rregular)
    #v = self.ver_star_matrix[:,0]
    anchor = self.vertices[v]
    t1 = (e1-e3)
    t2 = (e2-e4)
    lt1 = np.linalg.norm(t1,axis=1)
    lt2 = np.linalg.norm(t2,axis=1)
    ut1 = t1 / lt1[:,None]
    ut2 = t2 / lt2[:,None]
    angle = np.arccos(np.einsum('ij,ij->i', ut1, ut2))*180/np.pi
    if plot:
        a,b = (l1+l3)/5.0, (l2+l4)/5.0
        Vl,Vr = anchor+ut1*a[:,None], anchor-ut1*a[:,None]
        pl1 = make_polyline_from_endpoints(Vl,Vr)
        Vl,Vr = anchor+ut2*b[:,None], anchor-ut2*b[:,None] 
        pl2 = make_polyline_from_endpoints(Vl,Vr)
        return pl1,pl2
    return lt1,lt2,ut1,ut2,anchor,angle

def get_v4_diag_unit_edge(self):
    """
    Computes the unit edge vectors for diagonal edges of vertices with valence 4.

    This function calculates the unit edge vectors for diagonal edges of vertices 
    with valence 4.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    v : numpy array
        The indices of vertices.
    l1, l2, l3, l4 : numpy array
        The lengths of the diagonal edge vectors.
    e1, e2, e3, e4 : numpy array
        The unit diagonal edge vectors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The unit diagonal edge vectors are useful for mesh analysis and processing.

    See Also
    --------
    get_v4_unit_edge : Computes the unit edge vectors for vertices with valence 4.
    """
    V = self.vertices
    v,v1,v2,v3,v4 = self.rr_star_corner
    E1 = V[v1]-V[v]
    E2 = V[v2]-V[v]
    E3 = V[v3]-V[v]
    E4 = V[v4]-V[v]
    l1 = np.linalg.norm(E1, axis=1)
    l2 = np.linalg.norm(E2, axis=1)
    l3 = np.linalg.norm(E3, axis=1)
    l4 = np.linalg.norm(E4, axis=1)
    e1 = E1 / l1[:,None]
    e2 = E2 / l2[:,None]
    e3 = E3 / l3[:,None]
    e4 = E4 / l4[:,None]     
    return v,l1,l2,l3,l4,e1,e2,e3,e4

def get_v4_diag_unit_tangents(self, plot=False):
    """
    Computes the unit tangent vectors for diagonal edges of vertices with valence 4.

    This function calculates the unit tangent vectors for diagonal edges of vertices 
    with valence 4, optionally plotting the results.

    Parameters
    ----------
    plot : bool, optional (default=False)
        Whether to plot the tangent vectors.

    Returns
    -------
    lt1, lt2 : numpy array
        The lengths of the tangent vectors.
    ut1, ut2 : numpy array
        The unit tangent vectors.
    anchor : numpy array
        The anchor points for the tangent vectors.
    angle : numpy array
        The angles between the tangent vectors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The unit tangent vectors are useful for mesh analysis and processing.

    (only for valence 4, not depends on X)

    See Also
    --------
    get_v4_unit_tangents : Computes the unit tangent vectors for vertices with valence 4.
    """
    v,l1,l2,l3,l4,e1,e2,e3,e4 = self.get_v4_diag_unit_edge()
    #v,_,_,_,_ = self.rr_star_corner
    anchor = self.vertices[v]
    t1 = (e1-e3)
    t2 = (e2-e4)
    lt1 = np.linalg.norm(t1,axis=1)
    lt2 = np.linalg.norm(t2,axis=1)
    ut1 = t1 / lt1[:,None]
    ut2 = t2 / lt2[:,None]
    angle = np.arccos(np.einsum('ij,ij->i', ut1, ut2))*180/np.pi
    if plot:
        a,b = (l1+l3)/6.0, (l2+l4)/6.0
        Vl,Vr = anchor+ut1*a[:,None], anchor-ut1*a[:,None]
        pl1 = make_polyline_from_endpoints(Vl,Vr)
        Vl,Vr = anchor+ut2*b[:,None], anchor-ut2*b[:,None] 
        pl2 = make_polyline_from_endpoints(Vl,Vr)
        return pl1,pl2           
    return lt1,lt2,ut1,ut2,anchor,angle

def get_v4_unit_normal(self, diag=False, rregular=True):
    """
    Computes the unit normal vectors for vertices with valence 4.

    This function calculates the unit normal vectors for vertices with valence 4, 
    optionally considering diagonal edges and regular vertices.

    Parameters
    ----------
    diag : bool, optional (default=False)
        Whether to consider diagonal edges.
    rregular : bool, optional (default=True)
        Whether to consider regular vertices.

    Returns
    -------
    v : numpy array
        The indices of vertices.
    anchor : numpy array
        The anchor points for the normal vectors.
    unit_normal : numpy array
        The unit normal vectors.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The unit normal vectors are useful for mesh analysis and processing.

    See Also
    --------
    get_v4_unit_tangents : Computes the unit tangent vectors for vertices with valence 4.
    """
    v = self.ver_rrv4f4
    if diag:
        _,_,t1,t2,an,_ = self.get_v4_diag_unit_tangents()
    else:
        _,_,t1,t2,an,_ = self.get_v4_unit_tangents(rregular=rregular)
    n = np.cross(t1,t2)
    un = n / np.linalg.norm(n,axis=1)[:,None]
    return v,an, un

def get_v4_orient_unit_normal(self, diag=False, rregular=True):
    """
    Computes the oriented unit normal vectors for vertices with valence 4.

    This function calculates the oriented unit normal vectors for vertices with 
    valence 4, optionally considering diagonal edges and regular vertices.

    Parameters
    ----------
    diag : bool, optional (default=False)
        Whether to consider diagonal edges.
    rregular : bool, optional (default=True)
        Whether to consider regular vertices.

    Returns
    -------
    anchor : numpy array
        The anchor points for the normal vectors.
    unit_normal : numpy array
        The oriented unit normal vectors.
    angle : numpy array
        The angles between the normal vectors and the reference direction.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The oriented unit normal vectors are useful for mesh analysis and processing.

    (updated for each time, orientn; defined at rrv4f4)
    
    See Also
    --------
    get_v4_unit_normal : Computes the unit normal vectors for vertices with valence 4.
    """
    v,an,vN = self.get_v4_unit_normal(diag,rregular)
    Nv = self.vertex_normals()[v]
    i = np.where(np.einsum('ij,ij->i',Nv,vN)<0)[0]
    vN[i] = -vN[i]
    a = np.sqrt(np.abs(np.einsum('ij,ij->i',vN,Nv))) ##vN*Nv=a^2
    return [an,vN,a]

def get_net_crossing_angle(self, ut1, ut2):
    """
    Computes the crossing angles between two sets of tangent vectors.

    This function calculates the crossing angles between two sets of tangent vectors.

    Parameters
    ----------
    ut1 : numpy array
        The first set of unit tangent vectors.
    ut2 : numpy array
        The second set of unit tangent vectors.

    Returns
    -------
    angle : numpy array
        The crossing angles between the tangent vectors.

    Note
    ----
    This function assumes that the input tangent vectors are unit vectors.
    The crossing angles are useful for mesh analysis and processing.

    (for orthogonal / isogonal case)

    See Also
    --------
    get_v4_unit_tangents : Computes the unit tangent vectors for vertices with valence 4.
    """
    cos1 = np.einsum('ij,ij->i', ut1,ut2)
    A = np.arccos(cos1)*180/np.pi
    print('----- net crossing angles : -------')
    print('max=', '%.2g' % np.max(A))
    print('mean=', '%.2g' % np.mean(A))
    print('min=', '%.2g' % np.min(A))
    print('----------------------------------')
    #return np.max(A),np.min(A),np.mean(A),np.median(A)

# -------------------------------------------------------------------------
#                     Discrete Differential Geometry
# -------------------------------------------------------------------------
def get_curvature(self, mesh, order=None):
    """
    Computes the curvature properties of the mesh.

    This function calculates the curvature properties of the mesh, including 
    principal curvatures, Gaussian curvature, and mean curvature.

    Parameters
    ----------
    mesh : Mesh object
        The input mesh.
    order : numpy array, optional (default=None)
        The order of vertices to consider.

    Returns
    -------
    ratio : list
        The ratio of principal curvatures.
    K : list
        The Gaussian curvature values.
    H : list
        The mean curvature values.
    D1, D2 : numpy array
        The principal curvature directions.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The curvature properties are useful for mesh analysis and processing.

    (from davide' eigen of shape-operator)

    See Also
    --------
    get_curvature_libigl : Computes the curvature properties using libigl.
    """
    k1,k2, D1, D2 = mesh.principal_curvatures(True, use_sine=True)
    K = mesh.gaussian_curvature()
    H = mesh.mean_curvature()
    eps = np.finfo(float).eps
    if order is not None:
        K,H = K[order],H[order]
        D1,D2 = D1[order],D2[order]
        k1,k2 = k1[order],k2[order]
    ratio = [np.min(k1/(k2+eps)),np.mean(k1/(k2+eps)),np.max(k1/(k2+eps))]
    return ratio,[np.min(K),np.mean(K),np.max(K)],[np.min(H),np.mean(H),np.max(H)],D1,D2

def get_curvature_libigl(self, mesh, evalue=False):
    """
    Computes the curvature properties of the mesh using libigl.

    This function calculates the curvature properties of the mesh using libigl, 
    including principal curvatures, Gaussian curvature, and mean curvature.

    Parameters
    ----------
    mesh : Mesh object
        The input mesh.
    evalue : bool, optional (default=False)
        Whether to return the eigenvalues of the shape operator.

    Returns
    -------
    vertices : numpy array
        The vertex positions.
    K : list
        The Gaussian curvature values.
    H : list
        The mean curvature values.
    D1, D2 : numpy array
        The principal curvature directions.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The curvature properties are useful for mesh analysis and processing.

    (via quadric fitting (Panozzo, 2010))

    See Also
    --------
    get_curvature : Computes the curvature properties of the mesh.
    """ 
    import igl
    trv = mesh.vertices
    trf,_ = mesh.face_triangles()
    D1,D2,k1,k2 = igl.principal_curvature(trv,trf)  
    if evalue:
        K,H = k1*k2, (k1+k2)/2
        return trv,[np.min(K),np.mean(K),np.max(K)],[np.min(H),np.mean(H),np.max(H)],D1,D2
    return trv, k2, k1, D2, D1 #rearange in [min,max]

# -------------------------------------------------------------------------
#                         Multinets # Xinye
# -------------------------------------------------------------------------

def Get_Inner_TwinsHalfedge_Index(self):
    """
    Retrieves the indices of inner twin half-edges.

    This function identifies and returns the indices of inner twin half-edges 
    in the mesh.

    goal: find the inner twins halfedges indexes (unrepeatable)
    
    method:

    1. give the sequence number of all halfedges;

    2. delete the boundary twin halfedges;

    3. delete the inner repeated sequence number

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    indices : numpy array
        The indices of inner twin half-edges.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The inner twin half-edges are useful for mesh analysis and processing.


    See Also
    --------
    Get_Diagonals_of_Multinets : Retrieves the diagonals of multinets.
    """
    H=self.halfedges  
    number_of_halfedges = len(H[:,0])
    Boundary_SequenceNumber = []
    Boundary_Twin_SequenceNumber = []
    Inner_Repeated_SequenceNumber =[]
    Bool_SequenceNumber = np.ones(number_of_halfedges)
    SequenceNumber = list(range(0,number_of_halfedges))   
    _,_, valence = self.vertex_ring_vertices_iterators(return_lengths=True)
    
    for i in SequenceNumber :
        if H[i,1] == -1:
            Boundary_SequenceNumber.append(i)
            Bool_SequenceNumber[i] = 0
        if H[i,4] < i:
            Inner_Repeated_SequenceNumber.append(i)
            Bool_SequenceNumber[i] = 0           
        
    for i in Boundary_SequenceNumber :
        n = H[i,4]
        Boundary_Twin_SequenceNumber.append(n)
        Bool_SequenceNumber[n] = 0
        
    for i in SequenceNumber :
        if H[H[H[i,2],4],1]!=-1 and H[H[H[H[i,2],2],4],1]!=-1 and valence[H[i,0]] != 4 :
            Bool_SequenceNumber[i] = 0
        if H[H[H[i,2],4],1]!=-1 and H[H[H[H[i,2],2],4],1]!=-1 and valence[H[H[i,2],0]] != 4 :
            Bool_SequenceNumber[i] = 0
        if H[H[H[i,2],4],1]!=-1 and H[H[H[H[i,2],2],4],1]!=-1 and valence[H[H[H[i,2],2],0]] != 4 :
            Bool_SequenceNumber[i] = 0
        if H[H[H[i,2],4],1]!=-1 and H[H[H[H[i,2],2],4],1]!=-1 and valence[H[H[i,3],0]] != 4 :
            Bool_SequenceNumber[i] = 0
        
    self._Bool_SequenceNumber = Bool_SequenceNumber


def Get_Diagonals_of_Multinets(self):
    """
    Retrieves the diagonals of multinets in the mesh.

    This function identifies and returns the diagonals of multinets in the mesh. 
    It computes the diagonals based on the inner twin half-edges.

    Parameters
    ----------
    None
        The function operates on the mesh data within the class instance.

    Returns
    -------
    pl1 : Polyline object
        The first set of diagonal polylines.
    pl2 : Polyline object
        The second set of diagonal polylines.

    Note
    ----
    This function assumes that the mesh is a regular mesh.
    The diagonals of multinets are useful for mesh analysis and processing.

    See Also
    --------
    Get_Inner_TwinsHalfedge_Index : Retrieves the indices of inner twin half-edges.
    """
    H=self.halfedges  
    number_of_halfedges = len(H[:,0])
    V=self.vertices
    Bool_List = self.Bool_SequenceNumber
    Halfedge_1_NextNext = np.array([[0,0,0]])
    Halfedge_1_Prev = np.array([[0,0,0]])
    Halfedge_2_NextNext = np.array([[0,0,0]])
    Halfedge_2_Prev = np.array([[0,0,0]])
    # a=V[[H[H[H[1,2],2],0]]]
    # print(a)
    for i in range(number_of_halfedges) :
        if Bool_List[i] == 1:
            Halfedge_1_NextNext = np.r_[Halfedge_1_NextNext,V[[H[H[H[i,2],2],0]]]]
            Halfedge_1_Prev = np.r_[Halfedge_1_Prev,V[[H[H[i,3],0]]]]
            Halfedge_2_NextNext = np.r_[Halfedge_2_NextNext,V[[H[H[H[H[i,4],2],2],0]]]]
            Halfedge_2_Prev = np.r_[Halfedge_2_Prev,V[[H[H[H[i,4],3],0]]]]
    Halfedge_1_NextNext = np.delete(Halfedge_1_NextNext,0,axis=0)
    Halfedge_1_Prev = np.delete(Halfedge_1_Prev,0,axis=0)
    Halfedge_2_NextNext = np.delete(Halfedge_2_NextNext,0,axis=0)
    Halfedge_2_Prev = np.delete(Halfedge_2_Prev,0,axis=0)
    pl1 = make_polyline_from_endpoints(Halfedge_1_NextNext,Halfedge_2_NextNext)
    pl2 = make_polyline_from_endpoints(Halfedge_1_Prev,Halfedge_2_Prev)
    return pl1,pl2

