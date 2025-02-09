
#include <queue>
#include <set>
#include <unordered_map>
#include <iostream>

#include "../geometry/halfedge.h"
#include "debug.h"

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it does not want to perform the operation for
    whatever reason (e.g. you don't want to allow the user to erase the last vertex).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementaiton, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/

/*
    Compute the normal of a face, given the coordinates of its vertices
*/
Vec3 normal_of_vecs(std::vector<Vec3> positions) {
    Vec3 n;
    int n_verts = positions.size();
    for (int i = 0; i < n_verts; ++i) {
        n += cross(positions[i], positions[(i + 1) % n_verts]);
    }
    return n.unit();
}

Vec3 barycenter_of_vecs(std::vector<Vec3> positions) {
    Vec3 avg;
    int n_verts = positions.size();
    for (auto p : positions) avg += p;
    return avg * (1.0 / n_verts);
}

/*
    This method should replace the given vertex and all its neighboring
    edges and faces with a single face, returning the new face.
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_vertex(Halfedge_Mesh::VertexRef v) {
    HalfedgeRef hi = v->halfedge();
    // `f` will be the merged face 
    FaceRef f = hi->face();
    auto vhe = v->neighborhood_halfedges();
    int n_hes = vhe.size();
    for (int i = 0; i < n_hes; ++i) {
        HalfedgeRef he_nxt = vhe[i], he_cur = vhe[(i + 1) % n_hes];
        VertexRef v_cur = he_cur->twin()->vertex();
        v_cur->_halfedge = he_cur->next();
        
        HalfedgeRef he = he_cur;
        while (he->next() != he_nxt->twin()) {
            he = he->next();
            he->_face = f;
        };
        he->_next = he_nxt->next();
    }
    f->_halfedge = hi->next();

    for (auto he : vhe) {
        erase(he->edge()), erase(he), erase(he->twin());
        if (he->face() != f)
           erase(he->face());
    }
    erase(v);

    return f;
}

/*
    This method should erase the given edge and return an iterator to the
    merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_edge(Halfedge_Mesh::EdgeRef e) {
    HalfedgeRef he_1 = e->halfedge(), he_2 = he_1->twin();
    // Refuse to remove if `he_1` and `he_2` have `next` relation
    if (he_1 == he_2->next() || he_2 == he_1->next())
        return std::nullopt;
    HalfedgeRef he_1_nxt = he_1->next(), he_2_nxt = he_2->next();
    VertexRef v_1 = he_1->vertex(), v_2 = he_2->vertex();
    // `f_1` will be the merged face
    FaceRef f_1 = he_1->face(), f_2 = he_2->face();
    // Refuse to remove if the two sides of `e` connect to the same face
    if (f_1 == f_2)
        return std::nullopt;

    HalfedgeRef he_1_prev = he_1;
    while (he_1_prev->next() != he_1) he_1_prev = he_1_prev->next();
    HalfedgeRef he_2_prev = he_2;
    while (he_2_prev->next() != he_2) he_2_prev = he_2_prev->next();
    
    he_2_prev->_next = he_1_nxt, he_1_prev->_next = he_2_nxt;
    HalfedgeRef he = he_1_nxt;
    do {
        he->_face = f_1;
        he = he->next();
    } while (he != he_1_nxt);
    v_1->_halfedge = he_2_nxt, v_2->_halfedge = he_1_nxt;
    f_1->_halfedge = he_1_nxt;
    f_1->boundary |= f_2->boundary;
    
    erase(e), erase(he_1), erase(he_2), erase(f_2);
    return f_1;
}

/*
    This method should collapse the given edge and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(Halfedge_Mesh::EdgeRef e) {

    /*
       | he_1p   he_1n |
       |               |
       |      he_1     |     he_1 : v_1 -> v_2
      v_1 ----------- v_2
       |      he_2     |     he_2 : v_2 -> v_1
       |               |
       | he_2n   he_2p |

       v_2 will be erased
       if `he_1` is in a triangle, `he_1n, he_1p, he_1n->_edge` will be erased
       if `he_2` is in a triangle, `he_2n, he_2p, he_2n->_edge` will be erased
    */
    HalfedgeRef he_1 = e->halfedge(), he_2 = he_1->twin(), he_1n = he_1->next(), he_2n = he_2->next();
    VertexRef v_1 = he_1->vertex(), v_2 = he_2->vertex();
    HalfedgeRef he_1p, he_2p;
    {
        HalfedgeRef he = he_1;
        while (he->next() != he_1)
            he = he->next();
        he_1p = he;
    }
    {
        HalfedgeRef he = he_2;
        while (he->next() != he_2)
            he = he->next();
        he_2p = he;
    }

    // Reassign the `vertex` field of halfedges starting from `v_2`
    std::vector v_2_nhe = v_2->neighborhood_halfedges();
    for (auto he : v_2_nhe)
        he->_vertex = v_1;

    HalfedgeRef he_1n_twin = he_1n->_twin, he_1p_twin = he_1p->_twin;
    he_1n->_edge->_halfedge = he_1n_twin, he_1p->_edge->_halfedge = he_1p_twin; // In case `he_1n, he_1p` are is removed
    if (he_1n->next() == he_1p) { //if triangle
        he_1n_twin->_twin = he_1p_twin, he_1p_twin->_twin = he_1n_twin; //set twins
        erase(he_1n), erase(he_1p), erase(he_1n->_edge); //erase inner halfedges
        erase(he_1n->_face);
        he_1n_twin->_edge = he_1p->_edge; //reassign edge
        he_1p->_vertex->_halfedge = he_1n_twin; //reassign halfedge
    } else {
        he_1p->_next = he_1n, he_1p->_face->_halfedge = he_1p;
    }

    //repeat for other side
    HalfedgeRef he_2n_twin = he_2n->_twin, he_2p_twin = he_2p->_twin;
    he_2n->_edge->_halfedge = he_2n_twin, he_2p->_edge->_halfedge = he_2p_twin; // In case `he_2n, he_2p` are is removed
    if (he_2n->next() == he_2p) {
        he_2n_twin->_twin = he_2p_twin, he_2p_twin->_twin = he_2n_twin;
        erase(he_2n), erase(he_2p), erase(he_2n->_edge);
        erase(he_2n->_face);
        he_2n_twin->_edge = he_2p->_edge;
        he_2p->_vertex->_halfedge = he_2n_twin;
    } else {
        he_2p->_next = he_2n, he_2p->_face->_halfedge = he_2p;
    }

    v_1->_halfedge = he_2p_twin;
    v_1->pos = (v_1->pos + v_2->pos) * 0.5;

    erase(he_1), erase(he_2), erase(v_2), erase(he_1->_edge);

    // If both sides of `he_1n_twin->_edge` (`he_2n_twin->_edge`) are boundary faces, the edges can be removed
    if (he_1n_twin->_face->boundary && he_1n_twin->_twin->_face->boundary)
        erase_edge(he_1n_twin->_edge);
    if (he_2n_twin->_face->boundary && he_2n_twin->_twin->_face->boundary)
        erase_edge(he_2n_twin->_edge);
    return v_1;
}

/*
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(Halfedge_Mesh::EdgeRef e) {
    //return std::nullopt;
    //collect edges and halfedges
    std::vector<HalfedgeRef> h;
    std::vector<EdgeRef> eH;
    std::vector<HalfedgeRef> kill;
    std::vector<FaceRef> killF;

    HalfedgeRef fX = e -> halfedge() -> next();
    bool killfirst = false;
    if(fX->next()->next()->next()==fX)
    {
        killfirst = true;
        kill.push_back(fX->next());
        killF.push_back(fX -> face());
    }
    HalfedgeRef hcurrent = fX;
    h.push_back(fX);
    eH.push_back(fX->edge());
    hcurrent = hcurrent -> twin() -> next();
    while (hcurrent != fX)
    {
        h.push_back(hcurrent);
        eH.push_back(hcurrent->edge());
        hcurrent = hcurrent -> twin() -> next();
    }

    std::vector<HalfedgeRef> r;
    std::vector<EdgeRef> eR;


    bool killsend = false;
    fX = e -> halfedge() -> twin() -> next();
    if(fX->next()->next()->next()==fX)
    {
        kill.push_back(fX->next());
        killF.push_back(fX->face());
        killsend = true;
    }

    hcurrent = fX;

    r.push_back(fX);
    eR.push_back(fX->edge());
    hcurrent = hcurrent -> twin() -> next();
    while (hcurrent != r.front())
    {
        r.push_back(hcurrent);
        eR.push_back(hcurrent->edge());
        hcurrent = hcurrent -> twin() -> next();
    }

    VertexRef v0 = e->halfedge()->vertex();
    VertexRef v1 = e->halfedge()->twin()->vertex();

    //reassign


    for (unsigned long i = 0; i < h.size(); i++)
    {
        h.at(i) -> vertex() = v0;
        eH.at(i) -> halfedge() = h.at(i); // this shouldnt matter
    }

    for (unsigned long i = 0; i < r.size(); i++)
    {
        r.at(i) -> vertex() = v0; //this should already be set
        eR.at(i) -> halfedge() = r.at(i); //this should already be set
    }
    
    h.at(0) -> face() -> halfedge() = h.at(0); //ensure face isnt part of deleted halfedge
    r.at(0) -> face() -> halfedge() = r.at(0);
    h.at(0) -> next() -> twin() -> face() -> halfedge() = h.at(0) -> next() -> twin() -> next();
    r.at(0) -> next() -> twin() -> face() -> halfedge() = r.at(0) -> next() -> twin() -> next();
    v0->halfedge() = h.at(0);
    
    if(killfirst)
    {
    h.at(0) -> next() = r.at(r.size()-2) -> next(); //this is the next of second to last element; the element (not the next) will be deleted
    h.at(0) -> next() -> vertex() -> halfedge() = r.at(r.size()-2) -> next();
    h.at(0) -> face() = r.at(r.size()-2) -> next() -> face(); //assign it to new face bc that face will be deleated
    h.at(h.size()-3) -> twin() -> next() = r.at(0);
    }
    else
    {
        h.at(0) -> face() -> halfedge() = h.at(0); //ensure face isnt part of deleted halfedge
        r.at(r.size()-2) -> twin() -> next() = h.at(0);
    }

    if(killsend)
    {
    r.at(0) -> next() = h.at(h.size()-2) -> next(); //this is the next of second to last element; the element (not the next) will be deleted
    r.at(0) -> next() -> vertex() -> halfedge() = h.at(h.size()-2) -> next();
    r.at(0) -> face() = h.at(h.size()-2) -> next() -> face(); //assign it to new face bc that face will be deleated
    r.at(r.size()-3) -> twin() -> next() = h.at(0);
    }
    else
    {
        h.at(h.size()-2) -> twin() -> next() = r.at(0);
    }

    v0 -> pos = .5* (v0 -> pos + v1 -> pos);
    //now delete
    erase(v1);
    erase(e);
    erase(e -> halfedge());
    erase(e -> halfedge() -> twin());

    for (unsigned long i = 0; i < kill.size(); i++)
    {
        erase(kill.at(i) -> twin());
        erase(kill.at(i) -> edge());
        erase(kill.at(i));
        erase(killF.at(i));
    }
    


    return v0;
}
*/

/*
    This method should collapse the given face and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(Halfedge_Mesh::FaceRef f) {

    (void)f;
    return std::nullopt;
}

/*
    This method should flip the given edge and return an iterator to the
    flipped edge.
*/
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(Halfedge_Mesh::EdgeRef e) {

    if(e->on_boundary()) return std::nullopt;
    std::vector<HalfedgeRef> h;

    std::vector<VertexRef> v;
    std::vector<EdgeRef> eR;
    int l1edges = 0;
    HalfedgeRef hcurrent = e->halfedge();
    //collect all of the half edges, vertices, edges on one face
    h.push_back(hcurrent);
    v.push_back(hcurrent->vertex());
    eR.push_back(e);
    hcurrent = hcurrent->next();
    l1edges +=1;
    while (h.front() != hcurrent)
    {
        h.push_back(hcurrent);
        v.push_back(hcurrent->vertex());
        eR.push_back(hcurrent->edge());
        hcurrent = hcurrent->next();
        l1edges++;
    }
    
    //and then on the other face; hcurrent = h on edge
    hcurrent = hcurrent->twin();
    HalfedgeRef htw = hcurrent;
    h.push_back(hcurrent);
    if(std::find(v.begin(), v.end(),hcurrent->vertex())==v.end()) //if the vertex is not currently in the vector; should happen twice
        v.push_back(hcurrent->vertex());
    if(std::find(eR.begin(), eR.end(),hcurrent->edge())==eR.end()) //if the edge is not currently in the vector; should happen once
        eR.push_back(hcurrent->edge());
    hcurrent = hcurrent->next();
 
     while (htw != hcurrent)
    {
        h.push_back(hcurrent);
        if(std::find(v.begin(), v.end(),hcurrent->vertex())==v.end()) //if the vertex is not currently in the vector; should happen twice
            v.push_back(hcurrent->vertex());
        if(std::find(eR.begin(), eR.end(),hcurrent->edge())==eR.end()) //if the edge is not currently in the vector; should happen once
            eR.push_back(hcurrent->edge());
        hcurrent = hcurrent->next();

    } 

    
    
    //collect faces
    FaceRef f0 = e->halfedge()->face();
    FaceRef f1 = e->halfedge()->twin()->face();
    
    //reassign
    h.at(0) -> next() = h.at(2);
    h.at(0) -> vertex() = v.at(l1edges);
    h.at(0) -> twin() = h.at(l1edges);
    h.at(0) -> edge() = e; //unchanged
    h.at(0) -> face() = f0; //unchanged

    h.at(l1edges) -> next() = h.at(l1edges + 2) ;
    h.at(l1edges) -> vertex() = v.at(2);
    h.at(l1edges) -> twin() = h.at(0);
    h.at(l1edges) -> edge() = e;
    h.at(l1edges) -> face() = f1;

    h.at(l1edges-1) -> next() = h.at(l1edges+1);
    h.at(l1edges-1) -> face() = f0; //this should stay the same


    h.at(l1edges+1) -> next() = h.at(0);
    h.at(l1edges+1) -> face() = f0;

    h.back() -> next() = h.at(1);
    h.back() -> face() = f1; //this should be unchanged

    h.at(1) -> next() = h.at(l1edges);
    h.at(1) -> face() = f1;
    
    v.at(0) -> halfedge() = h.at(l1edges+1);
    v.at(1) -> halfedge() = h.at(1);

    e -> halfedge() = h.at(0);
    f0 -> halfedge() = h.at(0);
    f1 -> halfedge() = h.at(l1edges);



 
    
    

    //TODO: Reassignment
    return e;
}

/*
    This method should split the given edge and return an iterator to the
    newly inserted vertex. The halfedge of this vertex should point along
    the edge that was split, rather than the new edges.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(Halfedge_Mesh::EdgeRef e) {

    if (e->on_boundary())
    {
        HalfedgeRef hcurrent = e->halfedge();
            //collect
        std::vector<HalfedgeRef> h;
        std::vector<VertexRef> v;
        std::vector<EdgeRef> eR;

            //collect all of the half edges, vertices, edges on one face
        h.push_back(hcurrent);
        v.push_back(hcurrent->vertex());
        eR.push_back(e);
        hcurrent = hcurrent->next();

        while (h.front() != hcurrent)
        {
            h.push_back(hcurrent);
            v.push_back(hcurrent->vertex());
            hcurrent = hcurrent->next();

        }
    }

    HalfedgeRef hcurrent = e->halfedge();
    // Check both sides are on triangles
    if (hcurrent->next()->next()->next() != hcurrent)
        return std::nullopt;
    if (hcurrent->twin()->next()->next()->next() != hcurrent->twin())
        return std::nullopt;

    
    //collect
    std::vector<HalfedgeRef> h;
    std::vector<VertexRef> v;
    std::vector<EdgeRef> eR;

    //collect all of the half edges, vertices, edges on one face
    h.push_back(hcurrent);
    v.push_back(hcurrent->vertex());
    eR.push_back(e);
    hcurrent = hcurrent->next();

    while (h.front() != hcurrent)
    {
        h.push_back(hcurrent);
        v.push_back(hcurrent->vertex());
        hcurrent = hcurrent->next();

    }
    
    //and then on the other face; hcurrent = h on edge
    hcurrent = hcurrent->twin();
    HalfedgeRef htw = hcurrent;
    h.push_back(hcurrent);
    if(std::find(v.begin(), v.end(),hcurrent->vertex())==v.end()) //if the vertex is not currently in the vector; should happen twice
        v.push_back(hcurrent->vertex());

    hcurrent = hcurrent->next();
 
     while (htw != hcurrent)
    {
        h.push_back(hcurrent);
        if(std::find(v.begin(), v.end(),hcurrent->vertex())==v.end()) //if the vertex is not currently in the vector; should happen twice
            v.push_back(hcurrent->vertex());
        hcurrent = hcurrent->next();

    } 

    
    //collect faces
    FaceRef f0 = e->halfedge()->face();
    FaceRef f1 = e->halfedge()->twin()->face();


    //create new vertex
    v.push_back(new_vertex());

    //create new half edges
    for (int i = 0; i < 6; i++)
    {
        h.push_back(new_halfedge());
    }
    
    //create new edges
    for (int i = 0; i < 3; i++)
    {
        eR.push_back(new_edge());
    }
    
    // create new faces
    FaceRef f2 = new_face();
    FaceRef f3 = new_face();

    //H should have 12 elements;

 
    //now assign like hell
    h.at(0) -> next() = h.at(3);
    h.at(0) -> vertex() = v.at(0); // unchanged
    h.at(0) -> twin() = h.at(11);
    h.at(0) -> edge() = eR.at(0); //this is e
    h.at(0) -> face() = f0; //unchanged

    h.at(1) -> next() = h.at(6);
    h.at(1) -> face() = f2;

    h.at(2) -> next() = h.at(0);
    h.at(2) -> face() = f0; //unchanged

    h.at(3) -> next() = h.at(2);
    h.at(3) -> vertex() = v.at(4);
    h.at(3) -> twin() = h.at(6);
    h.at(3) -> face() = f0;
    h.at(3) -> edge() = eR.at(2);

    h.at(4) -> next() = h.at(10);
    h.at(4) -> face() = f1; //unchanged

    h.at(5) -> next() = h.at(8);
    h.at(5) -> face() = f3; 

    h.at(6) -> set_neighbors(h.at(7),h.at(3), v.at(2), eR.at(2),f2);
    h.at(7) -> set_neighbors(h.at(1), h.at(8), v.at(4), eR.at(1), f2);
    h.at(8) -> set_neighbors(h.at(9),h.at(7),v.at(1),eR.at(1), f3);
    h.at(9) -> set_neighbors(h.at(5),h.at(10),v.at(4),eR.at(3),f3);
    h.at(10) -> set_neighbors(h.at(11),h.at(9),v.at(3),eR.at(3),f1);
    h.at(11) -> set_neighbors(h.at(4),h.at(0),v.at(4),eR.at(0),f1);

    //set vertex
    v.at(0) -> halfedge() = h.at(4);
    v.at(1) -> halfedge() = h.at(1);
    v.at(2) -> halfedge() = h.at(2);
    v.at(3) -> halfedge() = h.at(5);
    v.at(4) -> halfedge() = h.at(3);
    v.at(4) -> pos = .5* (v.at(0) -> pos + v.at(1) -> pos);

    eR.at(0) -> halfedge() = h.at(0);
    eR.at(1) -> halfedge() = h.at(7);
    eR.at(2) -> halfedge() = h.at(3);
    eR.at(3) -> halfedge() = h.at(9);
    
    f0 -> halfedge() = h.at(0);
    f1 -> halfedge() = h.at(4);
    f3 -> halfedge() = h.at(5);
    f2 -> halfedge() = h.at(1);

    return v.at(4);
}

/* Note on the beveling process:

    Each of the bevel_vertex, bevel_edge, and bevel_face functions do not represent
    a full bevel operation. Instead, they should update the _connectivity_ of
    the mesh, _not_ the positions of newly created vertices. In fact, you should set
    the positions of new vertices to be exactly the same as wherever they "started from."

    When you click on a mesh element while in bevel mode, one of those three functions
    is called. But, because you may then adjust the distance/offset of the newly
    beveled face, we need another method of updating the positions of the new vertices.

    This is where bevel_vertex_positions, bevel_edge_positions, and
    bevel_face_positions come in: these functions are called repeatedly as you
    move your mouse, the position of which determins the normal and tangent offset
    parameters. These functions are also passed an array of the original vertex
    positions: for  bevel_vertex, it has one element, the original vertex position,
    for bevel_edge,  two for the two vertices, and for bevel_face, it has the original
    position of each vertex in halfedge order. You should use these positions, as well
    as the normal and tangent offset fields to assign positions to the new vertices.

    Finally, note that the normal and tangent offsets are not relative values - you
    should compute a particular new position from them, not a delta to apply.
*/

/*
    This method should replace the vertex v with a face, corresponding to
    a bevel operation. It should return the new face.  NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_vertex_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(Halfedge_Mesh::VertexRef v) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)v;
    return std::nullopt;
}

/*
    This method should replace the edge e with a face, corresponding to a
    bevel operation. It should return the new face. NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_edge_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(Halfedge_Mesh::EdgeRef e) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)e;
    return std::nullopt;
}

/*
    This method should replace the face f with an additional, inset face
    (and ring of faces around it), corresponding to a bevel operation. It
    should return the new face.  NOTE: This method is responsible for updating
    the *connectivity* of the mesh only---it does not need to update the vertex
    positions. These positions will be updated in
    Halfedge_Mesh::bevel_face_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_face(Halfedge_Mesh::FaceRef f) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."
    assert(!f->boundary); // f must not be a boundary face

    /* 
       Construct `hes_1, hes_2, hes_3, hes_4, hes_n, vs, nfs, res, ies`
       * `he_1` are all the halfedges of the face `f`
       * For each `i`, `vs[i]` is the copy of `hes_1[i]->vertex()`
       * For each `i`, `nfs[i]` is the new `ring` face associated with `hes_1[i]`
       * For each `i`, `hes_1[i]` is on the original face `f`
       * For each `i`, `hes_n[i]` is on the inset face
       * For each `i`, `res[i]` is the edge associated with `hes_2[i]`
       * For each `i`, `ies[i]` is the edge associated with `hes_3[i]`
                  he_n
              <-----------
             |    he_3    |
             | he_4  he_2 |
             |    he_1    |
              ----------->
    */
    FaceRef inset_face = f;
    std::vector<HalfedgeRef> hes_1, hes_2, hes_3, hes_4, hes_n;
    std::vector<EdgeRef> res, ies;
    std::vector<VertexRef> vs;
    std::vector<FaceRef> nfs;
    HalfedgeRef he_start = f->halfedge();
    he_start->vertex();
    // `he` iterates over all half edges of `f`
    HalfedgeRef he = he_start;
    do {
        hes_1.push_back(he);
        hes_2.push_back(new_halfedge());
        hes_3.push_back(new_halfedge());
        hes_4.push_back(new_halfedge());
        hes_n.push_back(new_halfedge());
        res.push_back(new_edge());
        ies.push_back(new_edge());
        vs.push_back(new_vertex());
        nfs.push_back(new_face());
        he = he->next();
    } while (he != he_start);

    // Wrap around
    int n_verts = vs.size();
    hes_1.push_back(hes_1[0]), hes_2.push_back(hes_2[0]);
    hes_3.push_back(hes_3[0]), hes_4.push_back(hes_4[0]);
    hes_n.push_back(hes_n[0]);
    res.push_back(res[0]), ies.push_back(ies[0]);
    vs.push_back(vs[0]), nfs.push_back(nfs[0]);

    // Update mesh
    for (int i = 0; i < n_verts; ++i) {
        /*
                      he_n
               v4 <----------- v3
                 |    he_3    |
                 | he_4  he_2 |
                 |    he_1    |
               v1 -----------> v2
        */
        HalfedgeRef he_1 = hes_1[i], he_2 = hes_2[i], he_3 = hes_3[i], he_4 = hes_4[i], he_n = hes_n[i];
        VertexRef v1 = he_1->vertex(), v2 = he_1->next()->vertex(), v3 = vs[i + 1], v4 = vs[i];
        EdgeRef e_2 = res[i + 1], e_3 = ies[i], e_4 = res[i];
        FaceRef nf = nfs[i];

        // Halfedge
        he_1->_next = he_2, he_1->_face = nf;
        he_2->set_neighbors(he_3, hes_4[i + 1], v2, e_2, nf);
        he_3->set_neighbors(he_4, he_n, v3, e_3, nf);
        HalfedgeRef he_4_twin = (i == 0) ? hes_2[n_verts - 1] : hes_2[i - 1];
        he_4->set_neighbors(he_1, he_4_twin, v4, e_4, nf);
        he_n->set_neighbors(hes_n[i + 1], he_3, v4, e_3, inset_face);
        // Vertex
        v2->_halfedge = he_2;
        v3->pos = v2->pos, v3->_halfedge = he_3;
        // Edge
        e_2->_halfedge = he_2, e_3->_halfedge = he_3;
        v4->pos = v1->pos;
        // Face
        nf->_halfedge = he_1;
    }

    inset_face->_halfedge = hes_n[0];
    return inset_face;
}

/*
    Compute new vertex positions for the vertices of the beveled vertex.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the original vertex position and its associated outgoing edge
    to compute a new vertex position along the outgoing edge.
*/
void Halfedge_Mesh::bevel_vertex_positions(const std::vector<Vec3>& start_positions,
                                           Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled edge.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the orig array) to compute an offset vertex position.

    Note that there is a 1-to-1 correspondence between halfedges in
    newHalfedges and vertex positions
    in orig.  So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vector3D pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_edge_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled face.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the start_positions array) to compute an offset vertex
    position.

    Note that there is a 1-to-1 correspondence between halfedges in
    new_halfedges and vertex positions
    in orig. So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vec3 pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_face_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset,
                                         float normal_offset) {

    if(flip_orientation) normal_offset = -normal_offset;
    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while (h != face->halfedge());

    int n_verts = new_halfedges.size();
    //Vec3 start_norm = normal_of_vecs(start_positions);
    Vec3  start_center = barycenter_of_vecs(start_positions);
    for (int i = 0; i < n_verts; ++i) {
        VertexRef v = new_halfedges[i]->vertex();
        Vec3 pos = start_positions[i];
        pos -= normal_offset * normal_of_vecs(start_positions);
        pos += tangent_offset * (start_positions[i] - start_center);
        v->pos = pos;
    }
}

void Halfedge_Mesh::triangulate_face(FaceRef f) {
    // Do not triangulate virtual boundary face
    if (f->boundary)
        return;
    HalfedgeRef hi = f->halfedge();
    // Do not triangulate triangular faces
    if (hi->next()->next()->next() == hi)
        return;
    
    HalfedgeRef he;
    /*
      We use the "fan" triangulation method

      Suppose `f` is a `n`-gon, we will have `n - 3` line segments
      starting from the `base` vertex that divides `f` into `n - 2` triangles

    vs[i] ---------- vs[i + 1]
         \ hes_f[i] /
          \        /
           \      / hes_radial[i + 1]
            \    /  es_radial[i + 1]
             \  /
              \/
             base
    
      `hes_twin[i]` is the twin of `hes_radial[i]`
    */
    VertexRef base = hi->vertex();
    std::vector<HalfedgeRef> hes_radial, hes_twin, hes_f;
    std::vector<VertexRef> vs;
    std::vector<EdgeRef> es_radial;
    std::vector<FaceRef> fs_n;
    he = hi;
    hes_radial.push_back(he), hes_twin.push_back(he->twin());
    vs.push_back(he->next()->vertex());
    es_radial.push_back(he->edge());
    fs_n.push_back(f);
    do {
        if (he != hi) {
            hes_f.push_back(he);
            if (he != hi->next()) {
                hes_radial.push_back(new_halfedge());
                hes_twin.push_back(new_halfedge());
                vs.push_back(he->vertex());
                es_radial.push_back(new_edge());
                fs_n.push_back(new_face());
            }
        }
        he = he->next();
    } while (he->next() != hi);
    hes_radial.push_back(he->twin()), hes_twin.push_back(he);
    vs.push_back(he->vertex());
    es_radial.push_back(he->edge());

    // Update mesh
    /*
            e_2
    v_2 ------------ v_1
       \    he_2    /
        \ he_3     /
         \   he_1 /
      e_3 \      / e_1
           \    /
            \  /
             \/
            base
    `nf` is the `FaceRef` of the above triangle
    `he_3` points from `base` to `v_2`, `he_2` points from `v_2` to `v_1`
    */
    int n_trgs = fs_n.size();
    for (int i = 0; i < n_trgs; ++i) {
        HalfedgeRef he_1 = hes_twin[i + 1], he_2 = hes_f[i], he_3 = hes_radial[i];
        VertexRef v_1 = vs[i + 1];
        EdgeRef e_1 = es_radial[i + 1], e_3 = es_radial[i];
        FaceRef nf = fs_n[i];

        // Halfedge
        he_1->set_neighbors(he_3, hes_radial[i + 1], v_1, e_1, nf);
        he_2->_next = he_1, he_2->_face = nf;
        he_3->set_neighbors(he_2, hes_twin[i], base, e_3, nf);
        // Edge
        e_1->_halfedge = he_1, e_3->_halfedge = he_3;
        // Face
        nf->_halfedge = he_2;
    }
}

/*
    Splits all non-triangular faces into triangles.
*/
void Halfedge_Mesh::triangulate() {
    std::vector<FaceRef> faces_copy;
    for (FaceRef f = faces.begin(); f != faces.end(); ++f)
        faces_copy.push_back(f);
    for (FaceRef f : faces_copy)
        triangulate_face(f);
}

/* Note on the quad subdivision process:

        Unlike the local mesh operations (like bevel or edge flip), we will perform
        subdivision by splitting *all* faces into quads "simultaneously."  Rather
        than operating directly on the halfedge data structure (which as you've
        seen is quite difficult to maintain!) we are going to do something a bit nicer:
           1. Create a raw list of vertex positions and faces (rather than a full-
              blown halfedge mesh).
           2. Build a new halfedge mesh from these lists, replacing the old one.
        Sometimes rebuilding a data structure from scratch is simpler (and even
        more efficient) than incrementally modifying the existing one.  These steps are
        detailed below.

  Step I: Compute the vertex positions for the subdivided mesh.
        Here we're going to do something a little bit strange: since we will
        have one vertex in the subdivided mesh for each vertex, edge, and face in
        the original mesh, we can nicely store the new vertex *positions* as
        attributes on vertices, edges, and faces of the original mesh. These positions
        can then be conveniently copied into the new, subdivided mesh.
        This is what you will implement in linear_subdivide_positions() and
        catmullclark_subdivide_positions().

  Steps II-IV are provided (see Halfedge_Mesh::subdivide()), but are still detailed
  here:

  Step II: Assign a unique index (starting at 0) to each vertex, edge, and
        face in the original mesh. These indices will be the indices of the
        vertices in the new (subdivided mesh).  They do not have to be assigned
        in any particular order, so long as no index is shared by more than one
        mesh element, and the total number of indices is equal to V+E+F, i.e.,
        the total number of vertices plus edges plus faces in the original mesh.
        Basically we just need a one-to-one mapping between original mesh elements
        and subdivided mesh vertices.

  Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
        the element indices defined above. In other words, each new quad should be
        of the form (i,j,k,l), where i,j,k and l are four of the indices stored on
        our original mesh elements.  Note that it is essential to get the orientation
        right here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces
        should circulate in the same direction as old faces (think about the right-hand
        rule).

  Step IV: Pass the list of vertices and quads to a routine that clears
        the internal data for this halfedge mesh, and builds new halfedge data from
        scratch, using the two lists.
*/

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    simple linear interpolation, e.g., the edge midpoints and face
    centroids.
*/
void Halfedge_Mesh::linear_subdivide_positions() {



    // For each vertex, assign Vertex::new_pos to
    // its original position, Vertex::pos.

    for (VertexRef v = vertices_begin(); v != vertices_end(); v++)
    {
        v -> new_pos = v -> pos;
    }
    

    // For each edge, assign the midpoint of the two original
    // positions to Edge::new_pos.

    for (EdgeRef e = edges_begin(); e != edges_end(); e++)
    {
        e -> new_pos = e -> center();
    }
    
    for (FaceRef f = faces_begin(); f != faces_end() ; f++)
    {
        f -> new_pos = f -> center();
    }
    

    // For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::new_pos. Note
    // that in general, NOT all faces will be triangles!
}

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    the Catmull-Clark rules for subdivision.

    Note: this will only be called on meshes without boundary
*/
void Halfedge_Mesh::catmullclark_subdivide_positions() {

    // The implementation for this routine should be
    // a lot like Halfedge_Mesh:linear_subdivide_positions:(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules. (These rules are outlined in the Developer Manual.)

    // Faces
    for (FaceRef f = faces_begin(); f != faces_end() ; f++)
    {
        f -> new_pos = f -> center();
    }


    // Edges
    for (EdgeRef e = edges_begin(); e != edges_end(); e++)
    {
        e -> new_pos = .5*(e -> center() 
        + e->halfedge() -> face() -> center() * .5 
        + e->halfedge() -> twin() -> face() -> center() * .5);
    }
    

    // Vertices
    for (VertexRef v = vertices_begin(); v != vertices_end(); v++)
    {
        //get face positions
        Vec3 Q(0);
        Vec3 R;
        int total = 0;
        HalfedgeRef hcurrent = v -> halfedge();
        HalfedgeRef first = hcurrent;
        Q += hcurrent -> face() -> center();
        R += hcurrent -> edge() -> center();
        total++;
        hcurrent = hcurrent -> twin() -> next();
        
        while (hcurrent != first)
        {
            Q += hcurrent -> face() -> center();
            R += hcurrent -> edge() -> center();
            total++;
            hcurrent = hcurrent -> twin() -> next();
        }
        
        Q = Q/total;
        R = R/total;

        v -> new_pos = (Q + 2*R + (total-3)*v->pos)/total;

    }
}

/*
        This routine should increase the number of triangles in the mesh
        using Loop subdivision. Note: this is will only be called on triangle meshes.
*/
void Halfedge_Mesh::loop_subdivide() {

    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::new_pos.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh. Use Vertex::is_new for this.
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::new_pos.
    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::is_new. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    // -> Now flip any new edge that connects an old and new vertex.
    // -> Finally, copy the new vertex positions into final Vertex::pos.

    // Each vertex and edge of the original surface can be associated with a
    // vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions
    // using the connectivity of the original (coarse) mesh; navigating this mesh
    // will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We
    // will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute updated positions for all the vertices in the original mesh, using
    // the Loop subdivision rule.

    // Next, compute the updated vertex positions associated with edges.

    // Next, we're going to split every edge in the mesh, in any order. For
    // future reference, we're also going to store some information about which
    // subdivided edges come from splitting an edge in the original mesh, and
    // which edges are new.
    // In this loop, we only want to iterate over edges of the original
    // mesh---otherwise, we'll end up splitting edges that we just split (and
    // the loop will never end!)

    // Finally, flip any new edge that connects an old and new vertex.

    // Copy the updated vertex positions to the subdivided mesh.
}

/*
    Isotropic remeshing. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if this is not a triangle mesh)
*/
bool Halfedge_Mesh::isotropic_remesh() {

    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}

const float invertible_threshold = 1e-6;

/* Helper type for quadric simplification */
struct Edge_Record {
    Edge_Record() {
    }
    Edge_Record(std::unordered_map<Halfedge_Mesh::VertexRef, Mat4>& vertex_quadrics,
                Halfedge_Mesh::EdgeRef e)
        : edge(e) {

        // Compute the combined quadric from the edge endpoints.
        // -> Build the 3x3 linear system whose solution minimizes the quadric error
        //    associated with these two endpoints.
        // -> Use this system to solve for the optimal position, and store it in
        //    Edge_Record::optimal.
        // -> Also store the cost associated with collapsing this edge in
        //    Edge_Record::cost.
        Halfedge_Mesh::HalfedgeRef he = e->halfedge();
        Halfedge_Mesh::VertexRef v_1 = he->vertex(), v_2 = he->twin()->vertex();
        Mat4 edge_quadric = vertex_quadrics[v_1] + vertex_quadrics[v_2], A = edge_quadric;
        Vec3 b(A[3][0], A[3][1], A[3][2]);
        A[3][0] = A[3][1] = A[3][2] = A[0][3] = A[1][3] = A[2][3] = 0, A[3][3] = 1;
        if (A.det() > invertible_threshold * pow(e->length(), 3.0)) {
            // `A` is invertible
            optimal = -1 * (A.inverse() * b);
            Vec4 optimal_(optimal, 1);
            cost = dot(optimal_, edge_quadric * optimal_);
        } else {
            // `A` is (approximately) singular
            float cost_1 = dot(Vec4(v_1->pos, 1.0), edge_quadric * Vec4(v_1->pos, 1.0));
            float cost_2 = dot(Vec4(v_2->pos, 1.0), edge_quadric * Vec4(v_2->pos, 1.0));
            Vec3 mid = (v_1->pos + v_2->pos) * 0.5;
            float cost_mid = dot(Vec4(mid, 1.0), edge_quadric * Vec4(mid, 1.0));
            // Find `t` such that `(1 - t) v_1 + t v_2` is best
            // Suppose `cost(t) = a t^2 + b t + c`, then
            // `c = cost_1, a / 4 + b / 2 + c = cost_mid, a + b + c = cost_2`
            float a = 2 * (cost_2 - 2 * cost_mid + cost_1), b = cost_2 - cost_1 - a, c = cost_1;
            float t = (b == 0) ? 0.5 : - a / (2 * b);
            if (t < 0) t = 0.0;
            if (t > 1) t = 1.0;
            optimal = (1 - t) * v_1->pos + t * v_2->pos;
            cost = a * t * t + b * t + c;
        }
    }
    Halfedge_Mesh::EdgeRef edge;
    Vec3 optimal;
    float cost;
};

/* Comparison operator for Edge_Records so std::set will properly order them */
bool operator<(const Edge_Record& r1, const Edge_Record& r2) {
    if(r1.cost != r2.cost) {
        return r1.cost < r2.cost;
    }
    Halfedge_Mesh::EdgeRef e1 = r1.edge;
    Halfedge_Mesh::EdgeRef e2 = r2.edge;
    return &*e1 < &*e2;
}

/** Helper type for quadric simplification
 *
 * A PQueue is a minimum-priority queue that
 * allows elements to be both inserted and removed from the
 * queue.  Together, one can easily change the priority of
 * an item by removing it, and re-inserting the same item
 * but with a different priority.  A priority queue, for
 * those who don't remember or haven't seen it before, is a
 * data structure that always keeps track of the item with
 * the smallest priority or "score," even as new elements
 * are inserted and removed.  Priority queues are often an
 * essential component of greedy algorithms, where one wants
 * to iteratively operate on the current "best" element.
 *
 * PQueue is templated on the type T of the object
 * being queued.  For this reason, T must define a comparison
 * operator of the form
 *
 *    bool operator<( const T& t1, const T& t2 )
 *
 * which returns true if and only if t1 is considered to have a
 * lower priority than t2.
 *
 * Basic use of a PQueue might look
 * something like this:
 *
 *    // initialize an empty queue
 *    PQueue<myItemType> queue;
 *
 *    // add some items (which we assume have been created
 *    // elsewhere, each of which has its priority stored as
 *    // some kind of internal member variable)
 *    queue.insert( item1 );
 *    queue.insert( item2 );
 *    queue.insert( item3 );
 *
 *    // get the highest priority item currently in the queue
 *    myItemType highestPriorityItem = queue.top();
 *
 *    // remove the highest priority item, automatically
 *    // promoting the next-highest priority item to the top
 *    queue.pop();
 *
 *    myItemType nextHighestPriorityItem = queue.top();
 *
 *    // Etc.
 *
 *    // We can also remove an item, making sure it is no
 *    // longer in the queue (note that this item may already
 *    // have been removed, if it was the 1st or 2nd-highest
 *    // priority item!)
 *    queue.remove( item2 );
 *
 */
template<class T> struct PQueue {
    void insert(const T& item) {
        queue.insert(item);
    }
    void remove(const T& item) {
        if(queue.find(item) != queue.end()) {
            queue.erase(item);
        }
    }
    const T& top(void) const {
        return *(queue.begin());
    }
    void pop(void) {
        queue.erase(queue.begin());
    }
    size_t size() {
        return queue.size();
    }

    std::set<T> queue;
};

Mat4 face_quadric(Halfedge_Mesh::FaceRef f) {
    assert (!f->is_boundary());
    assert (f->halfedge()->next()->next()->next() == f->halfedge());
    Vec3 norm = f->normal();
    Vec4 norm4(norm, -dot(norm, f->halfedge()->vertex()->pos));
    return outer(norm4, norm4);
}

Mat4 vertex_quadric(Halfedge_Mesh::VertexRef v, std::unordered_map<Halfedge_Mesh::FaceRef, Mat4>& face_quadrics) {
    Halfedge_Mesh::HalfedgeRef hi = v->halfedge(), he = hi;
    Mat4 ret(Mat4::Zero);
    do {
        Halfedge_Mesh::FaceRef f = he->face();
        if (!f->is_boundary())
            ret += face_quadrics[f];
        he = he->twin()->next();
    } while (he != hi);
    return ret;
}

const int simplification_factor = 4;

/*
    An edge `e` is collapsable iff collapsing it will result in a good mesh
*/
bool edge_collapsable(Halfedge_Mesh::EdgeRef e) {
    // `he_1` and `he_2` are the two halfedges of `e`
    // The vertex of `he_1` is `v_1` and the vertex of `he_2` is `v_2`
    Halfedge_Mesh::HalfedgeRef he_1 = e->halfedge(), he_2 = he_1->twin();
    Halfedge_Mesh::VertexRef v_1 = he_1->vertex(), v_2 = he_2->vertex();
    // If the two vertices of `e` are identical
    if (he_1->vertex() == he_1->next()->vertex())
        return false;
    // If `he_1` or `he_2` are within some `2-gon`
    if (he_1->next()->next() == he_1 || he_2->next()->next() == he_2)
        return false;
    // If there are two faces sharing two edges, `he_1->edge()` and `he_1->next()->edge()`
    if (he_1->next()->twin()->next() == he_1->twin())
        return false;
    // If there are two faces sharing two edges, `he_2->edge()` and `he_2->next()->edge()`
    if (he_2->next()->twin()->next() == he_2->twin())
        return false;
    std::unordered_map<Halfedge_Mesh::VertexRef, Halfedge_Mesh::HalfedgeRef> m_1 = v_1->neighborhood_map(), m_2 = v_2->neighborhood_map();

    for (auto vhe_3 : m_1) {
        Halfedge_Mesh::VertexRef v_3 = vhe_3.first;
        // If both `v_1` and `v_2` connect to the same vertex `v_3`
        if (m_2.find(v_3) != m_2.end()) {
            Halfedge_Mesh::HalfedgeRef he_13 = vhe_3.second, he_23 = m_2[v_3];
            // If `v_1, v_2, v_3` does not form a triangle, then collapsing
            // `e` will result in `v_1 v_3` and `v_2 v_3` to be incident to more than two faces
            bool v_123 = (he_13->twin()->next() == he_1) && (he_1->next() == he_23);
            bool v_321 = (he_23->twin()->next() == he_2) && (he_2->next() == he_13);
            if (!v_123 && !v_321)
                return false;
            // If the next halfedge of `v_1 v_3` and `v_2 v_3` are on the same edge,
            // then collapsing `e` will result in the edge `(v_1 v_2) v_3` to
            // be on two faces that are the same around the vertex `v_3`
            if (he_13->next()->twin()->next() == he_23->twin() || he_23->next()->twin()->next() == he_13->twin())
                return false;
            // If edge `v_1 v_3` and `v_2 v_3` are simultaneously on two identical faces,
            // then collapsing `e` will result in the two sides of the edge `(v_1 v_2) v_3` to
            // be on the same face
            if (he_13->twin()->face() == he_23->face() && he_23->twin()->face() == he_13->face())
                return false;
        }
    }
    return true;
}

/*
    Mesh simplification. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if you can't simplify the mesh any
    further without destroying it.)
*/
bool Halfedge_Mesh::simplify() {

    std::unordered_map<VertexRef, Mat4> vertex_quadrics;
    std::unordered_map<FaceRef, Mat4> face_quadrics;
    std::unordered_map<EdgeRef, Edge_Record> edge_records;
    PQueue<Edge_Record> edge_queue;

    for (Face f : faces)
        // If there are non-triangular faces, refuse to simplify
        if (!f.is_boundary() && f.halfedge()->next()->next()->next() != f.halfedge())
            return false;

    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in face_quadrics
    for (auto f = faces.begin(); f != faces.end(); ++f)
        if (!f->is_boundary())
            face_quadrics[f] = face_quadric(f);
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics
    for (auto v = vertices.begin(); v != vertices.end(); ++v)
        vertex_quadrics[v] = vertex_quadric(v, face_quadrics);
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above PQueue<Edge_Record> for this.
    for (auto e = edges.begin(); e != edges.end(); ++e) {
        Edge_Record er(vertex_quadrics, e);
        edge_records[e] = er;
        edge_queue.insert(er);
    }
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

    size_t target = faces.size() - (face_quadrics.size() - face_quadrics.size() / simplification_factor);
    bool collapsed = false;
    while (faces.size() > target && edge_queue.size()) {
        Edge_Record top = edge_queue.top();
        edge_queue.pop();
        if (!edge_collapsable(top.edge))
            continue;
        // Erase the two vertices of `top.edge` from `vertex_quadrics`
        VertexRef v1 = top.edge->halfedge()->vertex(), v2 = top.edge->halfedge()->twin()->vertex();
        Mat4 new_quardric = vertex_quadrics[v1] + vertex_quadrics[v2];
        std::vector<VertexRef> two_vertices({v1, v2});
        for (VertexRef v : two_vertices) {
            HalfedgeRef hi = v->halfedge(), he = hi;
            vertex_quadrics.erase(v);
            do {
                if (edge_records.find(he->edge()) != edge_records.end()) {
                    edge_queue.remove(edge_records[he->edge()]);
                    edge_records.erase(he->edge());
                }
                he = he->twin()->next();
            } while (he != hi);
        }
        auto v_collapse_option = collapse_edge_erase(top.edge);
        collapsed = true;
        assert (v_collapse_option.has_value());
        VertexRef v_collapse = v_collapse_option.value();
        {
            HalfedgeRef hi = v_collapse->halfedge(), he = hi;
            vertex_quadrics[v_collapse] = new_quardric;
            do {
                EdgeRef er = he->edge();
                if (edge_records.find(er) == edge_records.end()) {
                    Edge_Record e_rec(vertex_quadrics, er);
                    edge_records[er] = e_rec;
                    edge_queue.insert(e_rec);
                }
                he = he->twin()->next();
            } while (he != hi);
        }
    }

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.
    
    return collapsed;
}
