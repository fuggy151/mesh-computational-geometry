# Mesh Computational Geometry

This work comes from practical assignments at École Centrale de Lyon, course *Meshes and Computational Geometry*.

**Authors:** Danylo Mazurak, Anna Samsonenko  
**Professor:** Raphaëlle Chaine

---

## How to run the program

This assignment was completed using the gKit framework, studied during the practical sessions with Jean-Claude Iehl, for the visualization of results.

Before compiling, the following libraries may need to be installed:
```bash
sudo apt install libglew-dev libglew2.2
sudo apt install libsdl2-image-2.0-0 libsdl2-image-dev
sudo apt install libglfw3 libglfw3-dev
```

The program is compiled using the provided `Makefile`. From the project directory, run:
```bash
make
./main
```

At startup, the program prompts the user to choose between four modes:
```
1 - create meshes
2 - load mesh and open viewer (TP2-3)
3 - split operations (TP4)
4 - Delaunay and Voronoi
```

**Choosing `1`** creates and validates the three elementary meshes (tetrahedron, pyramid, and bounding box) and all test meshes, runs the full suite of unit tests (split, flip, orientation predicates, naive insertion, Delaunay), and saves results as `.off` files in the `data/` directory.

**Choosing `2`** loads a previously saved (mesh `data/queen.off` or any other OFF mesh can be used instead) and opens an interactive viewer with heat diffusion and curvature modes inherited from TP2–3.

**Choosing `3`** opens an interactive viewer with a mesh where mesh operations can be applied step by step using the keyboard.

**Choosing `4`** opens the Delaunay viewer for interactive point insertion and Voronoi diagram visualization.

---

### Controls — Mode `2` viewer

Inside the viewer, three display modes can be switched using the keyboard:

- `1` — plain mesh rendering (blue fill, black edges)
- `2` — heat diffusion simulation with color visualization
- `3` — mean curvature estimation with HSV color visualization

---

### Controls — Mode `3` viewer

- `T` — triangle split (split the face by inserting its centroid)
- `E` — edge split (split edge of the face by inserting its midpoint)
- `F` — edge flip
- `R` — reset to the original mesh

---

### Controls — Mode `4` viewer

- `I` — naive insertion of the point
- `D` — Delaunay insertion of the point (incremental Lawson)
- `P` — perform one Lawson flip step on a non-Delaunay edge
- `R` — reset the triangulation
- `N` — toggle the Voronoi diagram overlay

---

## Data structure

The `MeshSewn` data structure stores a triangular mesh as a list of `Vertex` records and a list of `Face` records. Each `Face` holds three vertex indices `v[0..2]` and three neighbor face indices `neigh[0..2]`, where `neigh[i]` is the face adjacent along edge $(v_i, v_{i+1 \bmod 3})$, or $-1$ on a boundary. Each `Vertex` stores a 3D position and a reference `oneFace` to any incident face, sufficient for local neighborhood traversal. The sewing relation is mutual: if `neigh[ei]` of face $f_i$ points to $f_j$, then $f_j$ must have a neighbor entry pointing back to $f_i$ across the same edge. All operations in this assignment preserve this invariant, which is verified after each modification by `check_sewing_verbose`.

## Implementation (TP2–3)

### Cotangent Laplacian

The Laplacian of a scalar function $u$ defined on the mesh vertices is computed using the cotangent formula:

$$
(Lu)_i = \frac{1}{2A_i} \sum_j (\cot \alpha_{ij} + \cot \beta_{ij})(u_j - u_i)
$$

where $A_i$ is estimated as one third of the total area of the triangles incident to vertex $i$. The implementation iterates over all faces once: for each triangle, the cotangent at each vertex is computed and its contribution is accumulated symmetrically to both endpoints of the opposite edge, so each edge is handled in a single pass without requiring explicit neighbor traversal.

---

### Heat diffusion

Heat diffusion on the surface is simulated using the explicit Euler scheme:

$$
u(t + \delta t) = u(t) + \delta t \cdot \Delta u
$$

A single heat source is placed at vertex $0$, whose temperature is held constant at $1.0$ throughout the simulation by resetting it after every integration step. All other vertices are initialized to $0.0$. At each frame, one time step is performed using the cotangent Laplacian. The simulation evolves until a stationary state is reached, where the temperature distribution no longer changes.

![Heat diffusion early stage](images/1.png)
*Heat diffusion at an early stage. The heat source at vertex 0 is visible as a small red region. The surrounding area is still cold.*

![Heat diffusion late stage](images/2.png)
*Heat diffusion at a later stage. The heat has spread further across the surface,following the intrinsic geometry of the mesh independently of any parametrization.*

#### Color normalization

For visualization, the temperature field is mapped to RGB colors: hot regions appear red and cold regions appear blue. Initially, a linear normalization was used:

$$
t_{\text{norm}} = \frac{u_i - u_{\min}}{u_{\max} - u_{\min}}
$$

However, since most vertices have very small temperature values compared to the source, the linear mapping produced a nearly uniform blue mesh with only the source visible. To improve the visual result, a power normalization was adopted:

$$
t_{\text{norm}} = \left(\frac{u_i - u_{\min}}{u_{\max} - u_{\min}}\right)^{0.3}
$$

This compresses the lower end of the range, making the gradient visible across the entire surface. The simulation data itself is not affected — only the color mapping changes.

![Heat diffusion linear normalization](images/4.png)
*Heat diffusion with linear normalization. Most of the mesh appears uniformly blue, with little visible gradient away from the source.*

---

### Mean curvature

The mean curvature at each vertex is estimated by applying the Laplacian operator to each coordinate function. For a vertex $s$ with position $(x, y, z)$, the Laplacian vector $\Delta s = (\Delta x, \Delta y, \Delta z)$ satisfies $\Delta s = 2Hn$, where $H$ is the mean curvature and $n$ is the surface normal. The scalar curvature is thus:

$$
H = \frac{\|\Delta s\|}{2}
$$

In practice, $\Delta x$, $\Delta y$, and $\Delta z$ are each computed by a separate call to `compute_laplacian`, passing the corresponding coordinate of each vertex as the scalar function.

The curvature field is computed once when mode `3` is activated and does not change over time. For visualization, curvature values are normalized and encoded in HSV coordinates. The hue ranges from blue (low curvature, flat regions) to red (high curvature, sharp features). Saturation and value are fixed at $1.0$.

To handle outliers in the curvature distribution (range up to 320), a percentile clamp is applied: only values between the 5th and 95th percentiles are used for normalization. The displayed range for the queen mesh is $[0.88, 69.92]$.

![Mean curvature front](images/5.png)
![Mean curvature side](images/6.png)
*Mean curvature visualization of the queen mesh. Flat regions such as the forehead and cheeks appear blue, while sharp features such as the ears, lips, and hair details appear red.*

---


## General mesh operations (TP4)

### Triangle split

The `split_triangle(fi, pos)` operation inserts a new vertex at position `pos` inside face `fi`, replacing it with three new triangles that share the inserted vertex. The original face $(a, b, c)$ is reused as $(a, b, m)$, and two new faces $(b, c, m)$ and $(c, a, m)$ are appended, where $m$ is the index of the new vertex. The three inner edges connecting $m$ to the original vertices are sewn together. The three outer edges, previously connected to `fi`, are redistributed to the appropriate new faces by scanning the neighbor lists of the formerly adjacent faces and updating any back-pointer that still pointed to `fi`.

After the operation, the total face count increases by 2 and the vertex count by 1. The sewing check confirms that all neighbor relations remain mutual.

<table>
  <tr>
    <td><img src="images/split_triangle0.png" width="100%"/></td>
    <td><img src="images/split_triangle1.png" width="100%"/></td>
  </tr>
  <tr>
    <td><img src="images/split_triangle2.png" width="100%"/></td>
    <td><img src="images/split_triangle3.png" width="100%"/></td>
  </tr>
</table>

*Triangle split.*

---

### Edge split

The `split_edge(fi, ei, pos)` operation inserts a new vertex $m$ at position `pos` on edge `ei` of face `fi`, splitting both incident faces into two. Two cases are handled.

When the edge is a boundary (no adjacent face), only the face `fi` = $(a, b, c)$ is split into $(a, m, c)$ and $(m, b, c)$; the shared inner edge $m$–$c$ is sewn, and the former outer neighbors are redistributed.

When the edge is interior, both faces `fi` = $(a, b, c)$ and `fj` = $(b, a, d)$ are involved. They are rebuilt as $(a, m, c)$, $(m, b, c)$, $(a, d, m)$, and $(b, m, d)$. Four inner edges are sewn: $a$–$m$, $m$–$c$, $d$–$m$, and $m$–$b$. The four outer neighbors are then restored to their respective new faces by updating back-pointers in formerly adjacent faces.

After the operation, the face count increases by 2 and the vertex count by 1, regardless of whether the edge was interior or on the boundary.

<table>
  <tr>
    <td><img src="images/split_edge0.png" width="100%"/></td>
    <td><img src="images/split_edge2.png" width="100%"/></td>
  </tr>
  <tr>
    <td><img src="images/split_edge1.png" width="100%"/></td>
    <td><img src="images/split_edge3.png" width="100%"/></td>
  </tr>
</table>

*Edge split. A shared edge of the pyramid is replaced by two edges meeting at the inserted midpoint vertex. Screenshots from different angles.*

---

### Edge flip

The `flip_edge(fi, ei)` operation replaces the diagonal shared by two adjacent triangular faces with the other diagonal of their common quadrilateral. Given face `fi` = $(a, b, c)$ and its neighbor `fj` across edge `ei` with opposite vertex $d$, the two faces are rebuilt as $(a, d, c)$ and $(b, c, d)$. All four outer neighbors (two per original face) are preserved by reassigning them to the appropriate reconstructed face. The new shared edge $c$–$d$ (which becomes `neigh[1]` in both new faces) is then sewn.

No back-pointer correction of outer neighbors is needed here as a separate step, because the outer neighbors of `fi` and `fj` themselves do not change which face they point to — only the geometry of those faces changes. However, two of the four outer neighbors formerly pointing to one face now belong to the other, so their indices are corrected via the `fix_neigh` lambda.

**To avoid breaking the mesh, edge flips are applied only to the base edges of the pyramid, not to its side edges, and this restriction is used only in this visualization.**

<table>
  <tr>
    <td><img src="images/flip_edge0.png" width="100%"/></td>
    <td><img src="images/flip_edge1.png" width="100%"/></td>
  </tr>
</table>

*Edge flip.*

---

## Operations for 2D triangulations 

### Orientation test

The `orientation(A, B, C)` test evaluates the signed area of the triangle formed by three 2D points using the $2 \times 2$ determinant of the edge vectors:

$$
\mathrm{orient}(A, B, C) = (B_x - A_x)(C_y - A_y) - (B_y - A_y)(C_x - A_x)
$$

A positive result indicates counter-clockwise orientation, a negative result indicates clockwise orientation, and zero indicates collinearity. The predicate is used as the building block for all subsequent 2D geometric tests.

---

### "In triangle" test

The `in_triangle(A, B, C, P)` test checks whether a point $P$ lies inside, on the boundary of, or outside a triangle $(A, B, C)$. The idea is straightforward: if $P$ is inside the triangle, it should be on the same side of every edge. We check this by calling `orientation` three times, once for each edge:

$$
d_0 = \mathrm{orient}(A, B, P), \quad
d_1 = \mathrm{orient}(B, C, P), \quad
d_2 = \mathrm{orient}(C, A, P)
$$

If all three results have the same sign, then $P$ is on the same side of all three edges, which means it is inside the triangle. If any $d_i = 0$, the point lies exactly on that edge, so the function returns $0.0$ (on boundary). If all signs match and none is zero, the function returns $1.0$ (strictly inside). If the signs do not all match, $P$ crossed to the wrong side of at least one edge and is outside: the function returns $-1.0$.

---

### Naive insertion of a point in a 2D triangulation

The `insert_point(P)` function inserts a 2D point into the triangulation without any quality guarantee. Insertion proceeds by scanning all faces linearly in `find_triangle` to locate the face containing $P$. If $P$ falls strictly inside a face, `split_triangle` is called. If $P$ lies exactly on an edge, `split_edge` is called instead. The boundary-edge index returned by `find_triangle` determines which case applies. If no face contains $P$, the point is outside the bounding box and insertion is rejected.

<table>
  <tr>
    <td><img src="images/naive_insertion0.png" width="100%"/></td>
    <td><img src="images/naive_insertion1.png" width="100%"/></td>
    <td><img src="images/naive_insertion2.png" width="100%"/></td>
  </tr>
</table>

*First image: Initial mesh without inserted points. Second: Naive insertion of 2 points. Non-Delaunay faces are highlighted in red. Third image: Naive insertion of 25 points.*

---

## Delaunay triangulation

### InCircle predicate

The `in_circle(A, B, C, D)` predicate determines whether a point $D$ lies inside the circumcircle of a counter-clockwise triangle $(A, B, C)$. This is evaluated using the standard $3 \times 3$ determinant after translating $D$ to the origin:

$$
\mathrm{InCircle}(A,B,C,D) = \begin{vmatrix}
A_x - D_x & A_y - D_y & (A_x-D_x)^2 + (A_y-D_y)^2 \\
B_x - D_x & B_y - D_y & (B_x-D_x)^2 + (B_y-D_y)^2 \\
C_x - D_x & C_y - D_y & (C_x-D_x)^2 + (C_y-D_y)^2
\end{vmatrix}
$$

A positive result means $D$ is inside the circumcircle; zero means $D$ is on the circle; a negative result means $D$ is outside.

---

### Local Delaunay check

The `is_locally_delaunay(fi, ei)` function checks whether edge `ei` of face `fi` satisfies the Delaunay condition. Boundary edges (with no adjacent face) are trivially Delaunay and always return `true`. For interior edges, the function retrieves the opposite vertex $D$ of the neighboring face and evaluates $\mathrm{InCircle}(A, B, C, D)$ where $A$, $B$ are the shared edge vertices and $C$ is the opposite vertex in `fi`. The edge is locally Delaunay if and only if $\mathrm{InCircle} \leq 0$, i.e. $D$ does not lie strictly inside the circumcircle of `fi`.

---

### Lawson's algorithm — `make_delaunay`

The `make_delaunay` function converts an arbitrary valid triangulation into a Delaunay triangulation by repeatedly scanning all edges and flipping any edge that violates the local Delaunay condition, until no such edge remains. Each flip strictly increases the minimum angle across the triangulation, so the algorithm always terminates. The outer loop iterates until a full pass produces no flip, at which point the triangulation is globally Delaunay.

<table>
  <tr>
    <td><img src="images/delaunay_vs_naive0.png" width="100%"/></td>
    <td><img src="images/delaunay_vs_naive1.png" width="100%"/></td>
  </tr>
</table>

*Left, naive insertion of 10 points, non-Delaunay edges in red, and right, Delaunay triangulation of non-Delaunay triangles.*

Each step of the Lawson algorithm:

<table>
  <tr>
    <td><img src="images/delaunay_vs_naive0.png" width="100%"/></td>
    <td><img src="images/Lawson1.png" width="100%"/></td>
    <td><img src="images/Lawson2.png" width="100%"/></td>
    <td><img src="images/Lawson3.png" width="100%"/></td>
  </tr>
  <tr>
    <td><img src="images/Lawson4.png" width="100%"/></td>
    <td><img src="images/Lawson5.png" width="100%"/></td>
    <td><img src="images/Lawson6.png" width="100%"/></td>
    <td><img src="images/Lawson7.png" width="100%"/></td>
  </tr>
  <tr>
    <td><img src="images/Lawson8.png" width="100%"/></td>
    <td><img src="images/Lawson9.png" width="100%"/></td>
    <td><img src="images/Lawson10.png" width="100%"/></td>
    <td><img src="images/Lawson11.png" width="100%"/></td>
  </tr>
  <tr>
    <td><img src="images/Lawson12.png" width="100%"/></td>
    <td><img src="images/Lawson13.png" width="100%"/></td>
    <td><img src="images/Lawson14.png" width="100%"/></td>
    <td><img src="images/delaunay_vs_naive1.png" width="100%"/></td>
  </tr>
</table>

---

### Incremental Delaunay insertion — `insert_point_delaunay`

The `insert_point_delaunay(P)` function maintains the Delaunay property incrementally after each point insertion. The point is first inserted naively using `insert_point`. Then, rather than running Lawson over the entire mesh, only edges not incident to the newly inserted vertex $m$ are examined, since edges incident to $m$ are always locally Delaunay by construction. Flips are applied one at a time; after each flip the scan restarts from the beginning of the face list until no non-Delaunay edge remains among the faces adjacent to $m$. This restricts the algorithmic work to the neighborhood of the insertion point.

<table>
  <tr>
    <td><img src="images/delaunay_insert0.png" width="100%"/></td>
    <td><img src="images/delaunay_insert1.png" width="100%"/></td>
    <td><img src="images/delaunay_insert2.png" width="100%"/></td>
  </tr>
  <tr>
    <td><img src="images/delaunay_insert3.png" width="100%"/></td>
    <td><img src="images/delaunay_insert4.png" width="100%"/></td>
    <td><img src="images/delaunay_insert5.png" width="100%"/></td>
  </tr>
  <tr>
    <td><img src="images/delaunay_insert6.png" width="100%"/></td>
    <td><img src="images/delaunay_insert7.png" width="100%"/></td>
    <td><img src="images/delaunay_insert8.png" width="100%"/></td>
  </tr>
</table>

*Step-by-step insertion of the 8 Delaunay points. All faces are blue.*


## Bonus: Voronoi diagram

### Circumcenter computation

For each triangular face $(A, B, C)$ the circumcenter is computed analytically. Translating to the local frame $\mathbf{a} = B - A$, $\mathbf{b} = C - A$ and letting $D = 2(\mathbf{a}_x \mathbf{b}_y - \mathbf{a}_y \mathbf{b}_x)$, the circumcenter offset from $A$ is:

$$
u_x = \frac{\mathbf{b}_y \|\mathbf{a}\|^2 - \mathbf{a}_y \|\mathbf{b}\|^2}{D}, \qquad
u_y = \frac{\mathbf{a}_x \|\mathbf{b}\|^2 - \mathbf{b}_x \|\mathbf{a}\|^2}{D}
$$

so the circumcenter is $(A_x + u_x,\; A_y + u_y,\; 0)$. When $|D| < 10^{-10}$ (degenerate triangle), the centroid is returned as a fallback.

---

### Building Voronoi edges

The `build_voronoi_edges` function iterates over all interior edges of the triangulation. For each sewn pair of faces $(f_i, f_j)$ with $j > i$ (to avoid duplicates), the circumcenters $c_i$ and $c_j$ are computed and a line segment $(c_i, c_j)$ is added to the output. By duality, this segment is the Voronoi edge separating the Voronoi cells of the two vertices shared by the Delaunay edge between $f_i$ and $f_j$.

---

### Rendering

Voronoi edges are uploaded to a dedicated VBO (`vbo_voronoi`) and rendered as `GL_LINES` with a line width of 2 pixels in a golden-yellow color (`vec3(1.0, 0.85, 0.0)`), overlaid on top of the Delaunay triangulation mesh. The Voronoi VAO is separate from the mesh VAO so neither rendering pass interferes with the other.

<table>
  <tr>
    <td><img src="images/voronoi1.png" width="100%"/></td>
    <td><img src="images/voronoi2.png" width="100%"/></td>
  </tr>
  <tr>
    <td><img src="images/voronoi5.png" width="100%"/></td>
    <td><img src="images/voronoi25.png" width="100%"/></td>
  </tr>
</table>

*Voronoi diagram (yellow) overlaid on the Delaunay triangulation (blue) for 1, 2, 5 and 25 random points.*


## Visualization with gKit

### TP2–3 viewer

Visualization was implemented using the gKit framework through an interactive viewer class derived from `AppCamera`.

A first approach using the gKit `Mesh` class was initially implemented. However, this required calling `mesh.release()` and `create_buffers()` at every color update, which caused repeated GPU memory allocation and deallocation — up to 12 times per second during heat diffusion. This also produced unwanted log messages from the gKit internals. A manual OpenGL approach was therefore adopted instead.

Two separate vertex buffer objects are maintained:

- `vbo_position` — stores vertex positions, allocated once with `GL_STATIC_DRAW`   since geometry never changes.
- `vbo_color` — stores per-vertex RGBA colors, allocated with `GL_DYNAMIC_DRAW` and updated at runtime using `glBufferSubData`, which uploads only the color data without touching the position buffer.

Both buffers are bound in a single VAO. The shader `mesh_color.glsl` uses a uniform integer flag `useVertexColor` to switch between two rendering modes: when set to `0`, a constant `baseColor` is used (mode `1`, plain blue mesh); when set to `1`, per-vertex colors from the color buffer are used (modes `2` and `3`). Each frame is rendered in two passes: a filled polygon pass for the surface, followed by an edge pass using `GL_POLYGON_OFFSET_LINE` to avoid z-fighting between the fill and the edges.

For the heat diffusion mode, the color buffer is updated every five frames. For the curvature mode, it is computed once and never updated again. In both cases, only the color buffer is reuploaded — the position buffer remains untouched throughout the entire session.

### TP4 — `DelaunayViewer`

The `DelaunayViewer` class, derived from `AppCamera`, manages interactive insertion and rendering in mode `4`. The same manual OpenGL approach used in TP2–3 is applied here: a position VBO and a color VBO are maintained separately, and `glBufferData` is called after every insertion or flip to reflect the updated mesh geometry and coloring.

After each insertion, `highlight_non_delaunay` scans all faces and colors any face red that has at least one non-Delaunay edge, leaving Delaunay faces blue. This gives immediate visual feedback on the state of the triangulation. Pressing `P` performs exactly one Lawson flip on the first non-Delaunay edge found, so the convergence of the algorithm can be observed step by step.

The rendering pipeline uses the same two-pass approach as in TP2–3: a filled polygon pass for face colors, followed by a wireframe pass with `GL_POLYGON_OFFSET_LINE` to avoid z-fighting on the edges.