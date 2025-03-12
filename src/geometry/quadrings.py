
import numpy as np


def regular_vertex_regular_quad(halfedges,ver_star_matrix,
                                delete_multi=True):
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
        

    See Also: quadfaces()

    
    Examples
    --------
    ```python
    from quadrings import regular_vertex_regular_quad
    num_rrf, rr_quadface, rr_quadface_order, num_rrv, rr_star, rr_4quad_vers = regular_vertex_regular_quad(halfedges, ver_star_matrix)
    ```
    """

    H, starM = halfedges,ver_star_matrix
    
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

    num_rrf = len(farr)
    rr_quadface = farr
    rr_quadface_order = forder
    num_rrv = num # same with num_regular
    #order = np.setdiff1d(np.arange(num),index//4)
    rr_star = starM # same with starM
    rr_4quad_vers = f4 #rr_4quad_vers

    return num_rrf, rr_quadface, rr_quadface_order, num_rrv, rr_star, rr_4quad_vers