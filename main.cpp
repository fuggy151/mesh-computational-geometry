#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <algorithm>
#include <ctime>
#include "wavefront.h"
#include "draw.h"
#include "app_camera.h"
#include "program.h"
#include "uniforms.h"
#include "window.h"

struct Vec3 { double x=0, y=0, z=0; };
struct Tri { std::array<int,3> v; };

struct Face {
    std::array<int,3> v;
    std::array<int,3> neigh;
    Face() : v{0,0,0}, neigh{-1,-1,-1} {}
    Face(int a,int b,int c) : v{a,b,c}, neigh{-1,-1,-1} {}
};

struct Vertex {
    Vec3 p;
    int oneFace = -1;
};

struct MeshSimple {
    std::vector<Vec3> vertices;
    std::vector<Tri> tris;
};

double dot(const Vec3& a, const Vec3& b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vec3 operator-(const Vec3& a, const Vec3& b)
{
    return {a.x-b.x, a.y-b.y, a.z-b.z};
}

Vec3 cross(const Vec3& a, const Vec3& b)
{
    return {
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

double norm(const Vec3& a)
{
    return std::sqrt(dot(a,a));
}

struct MeshSewn {
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    void sew(int fA, int eA, int fB, int eB) {
        faces[fA].neigh[eA] = fB;
        faces[fB].neigh[eB] = fA;
    }

    void update_oneFace_per_vertex() {
        for (auto& v : vertices) v.oneFace = -1;
        for (int fi = 0; fi < (int)faces.size(); ++fi) {
            for (int k = 0; k < 3; ++k) {
                int vid = faces[fi].v[k];
                if (vertices[vid].oneFace == -1) vertices[vid].oneFace = fi;
            }
        }
    }

    static int find_edge_index(const Face& f, int a, int b) {
        for (int i = 0; i < 3; ++i) {
            int u = f.v[i], v = f.v[(i+1)%3];
            if ((u==a && v==b) || (u==b && v==a)) return i;
        }
        return -1;
    }

    bool check_sewing() const {
        for (int fi = 0; fi < (int)faces.size(); ++fi) {
            for (int ei = 0; ei < 3; ++ei) {
                int fj = faces[fi].neigh[ei];
                if (fj == -1) continue;
                int a = faces[fi].v[ei];
                int b = faces[fi].v[(ei+1)%3];
                int ej = find_edge_index(faces[fj], a, b);
                if (ej == -1) return false;
                if (faces[fj].neigh[ej] != fi) return false;
            }
        }
        return true;
    }

    int split_triangle(int fi, Vec3 pos)
    {
        Face old = faces[fi];
        int a = old.v[0], b = old.v[1], c = old.v[2];
        int n0 = old.neigh[0], n1 = old.neigh[1], n2 = old.neigh[2];

        int m = (int)vertices.size();
        vertices.push_back(Vertex{pos});

        faces[fi] = Face(a, b, m);
        int f1 = (int)faces.size(); faces.push_back(Face(b, c, m));
        int f2 = (int)faces.size(); faces.push_back(Face(c, a, m));

        sew(fi, 1, f1, 2);
        sew(f1, 1, f2, 2);
        sew(f2, 1, fi, 2);

        faces[fi].neigh[0] = n0;
        faces[f1].neigh[0] = n1;
        faces[f2].neigh[0] = n2;

        auto fix_neigh = [&](int fj, int old_fi, int new_fi) {
            if (fj == -1) return;
            for (int e = 0; e < 3; ++e)
                if (faces[fj].neigh[e] == old_fi)
                    faces[fj].neigh[e] = new_fi;
        };
        fix_neigh(n1, fi, f1);
        fix_neigh(n2, fi, f2);

        update_oneFace_per_vertex();
        return m;
    }

    int split_edge(int fi, int ei, Vec3 pos)
    {
        int a = faces[fi].v[ei];
        int b = faces[fi].v[(ei+1)%3];
        int c = faces[fi].v[(ei+2)%3];

        int fj = faces[fi].neigh[ei];
        int fi_bc = faces[fi].neigh[(ei+1)%3];
        int fi_ca = faces[fi].neigh[(ei+2)%3];

        int m = (int)vertices.size();
        vertices.push_back(Vertex{pos});

        if (fj == -1) {
            // fi=(a,b,c) → fi=(a,m,c) + f1=(m,b,c)
            int f1 = (int)faces.size(); faces.push_back(Face(m, b, c));
            faces[fi] = Face(a, m, c);

            sew(fi, 1, f1, 2);      // m-c shared

            faces[fi].neigh[2] = fi_ca;
            faces[f1].neigh[1] = fi_bc;
            // edge a-m and m-b remain open (-1)

            auto fix = [&](int fk, int old_f, int new_f) {
                if (fk == -1) return;
                for (int e = 0; e < 3; ++e)
                    if (faces[fk].neigh[e] == old_f)
                        faces[fk].neigh[e] = new_f;
            };
            fix(fi_bc, fi, f1);

            update_oneFace_per_vertex();
            return m;
        }

        int ej = find_edge_index(faces[fj], a, b);
        int d = faces[fj].v[(ej+2)%3];
        int fj_bd = faces[fj].neigh[(ej+1)%3];
        int fj_da = faces[fj].neigh[(ej+2)%3];

        // New faces
        int f1 = (int)faces.size(); faces.push_back(Face(m, b, c));
        int f2 = (int)faces.size(); faces.push_back(Face(b, m, d));

        // Rebuild original faces
        faces[fi] = Face(a, m, c);
        faces[fj] = Face(a, d, m);

        // Internal sews
        sew(fi, 0, fj, 2);   // a-m ↔ m-a
        sew(fi, 1, f1, 2);   // m-c ↔ c-m
        sew(fj, 1, f2, 1);   // d-m ↔ m-d
        sew(f1, 0, f2, 0);   // m-b ↔ b-m

        // External neighbours
        faces[fi].neigh[2] = fi_ca;
        faces[f1].neigh[1] = fi_bc;
        faces[fj].neigh[0] = fj_da;
        faces[f2].neigh[2] = fj_bd;

        auto fix = [&](int fk, int old_f, int new_f) {
            if (fk == -1) return;
            for (int e = 0; e < 3; ++e)
                if (faces[fk].neigh[e] == old_f)
                    faces[fk].neigh[e] = new_f;
        };
        fix(fi_bc, fi, f1);
        fix(fj_bd, fj, f2);

        update_oneFace_per_vertex();
        return m;
    }

    void flip_edge(int fi, int ei)
    {
        int fj = faces[fi].neigh[ei];
        if (fj == -1) return;

        // Get vertices of face fi
        int a = faces[fi].v[ei];
        int b = faces[fi].v[(ei+1)%3];
        int c = faces[fi].v[(ei+2)%3];

        // Find edge (a,b) in neighboring face fj
        int ej = find_edge_index(faces[fj], a, b);
        int d = faces[fj].v[(ej+2)%3];

        // Store all neighbors before modification
        int fi_ac = faces[fi].neigh[(ei+2)%3];
        int fi_bc = faces[fi].neigh[(ei+1)%3];
        int fj_da = faces[fj].neigh[(ej+2)%3];
        int fj_bd = faces[fj].neigh[(ej+1)%3];

        faces[fi] = Face(a, d, c);
        faces[fj] = Face(b, c, d);

        faces[fi].neigh[0] = fj_da;
        faces[fi].neigh[1] = fj;
        faces[fi].neigh[2] = fi_ac;

        faces[fj].neigh[0] = fi_bc;
        faces[fj].neigh[1] = fi;
        faces[fj].neigh[2] = fj_bd;

        sew(fi, 1, fj, 1);

        auto fix_neigh = [&](int fk, int old_f, int new_f) {
            if (fk == -1) return;
            for (int e = 0; e < 3; ++e) {
                if (faces[fk].neigh[e] == old_f) {
                    faces[fk].neigh[e] = new_f;
                }
            }
        };

        fix_neigh(fi_bc, fi, fj);
        fix_neigh(fj_da, fj, fi);
        for(auto& f : faces)
            f.neigh = {-1, -1, -1};

        build_adjacency();
        update_oneFace_per_vertex();
    }

    double orientation(const Vec3& A, const Vec3& B, const Vec3& C) const
    {
        return (B.x - A.x) * (C.y - A.y)
            - (B.y - A.y) * (C.x - A.x);
    }

    double in_triangle(const Vec3& A, const Vec3& B, const Vec3& C, const Vec3& P) const
    {
        double d0 = orientation(A, B, P);
        double d1 = orientation(B, C, P);
        double d2 = orientation(C, A, P);

        bool all_pos = (d0 >= 0) && (d1 >= 0) && (d2 >= 0);
        bool all_neg = (d0 <= 0) && (d1 <= 0) && (d2 <= 0);

        if (!all_pos && !all_neg) return -1.0;  // outside

        if (d0 == 0 || d1 == 0 || d2 == 0) return 0.0;

        return 1.0;
    }

    void make_bounding_box(const std::vector<Vec3>& pts, double margin = 1.0)
    {
        double xmin = pts[0].x, xmax = pts[0].x;
        double ymin = pts[0].y, ymax = pts[0].y;
        for (const auto& p : pts) {
            xmin = std::min(xmin, p.x); xmax = std::max(xmax, p.x);
            ymin = std::min(ymin, p.y); ymax = std::max(ymax, p.y);
        }

        // Add margin
        xmin -= margin; ymin -= margin;
        xmax += margin; ymax += margin;

        vertices.clear();
        faces.clear();
        vertices.push_back(Vertex{{xmin, ymin, 0}});  // 0
        vertices.push_back(Vertex{{xmax, ymin, 0}});  // 1
        vertices.push_back(Vertex{{xmax, ymax, 0}});  // 2
        vertices.push_back(Vertex{{xmin, ymax, 0}});  // 3

        faces.push_back(Face(0, 1, 2));
        faces.push_back(Face(0, 2, 3));

        sew(0, 2, 1, 0);

        update_oneFace_per_vertex();
    }

    int find_triangle(const Vec3& P, int& boundary_edge) const
    {
        boundary_edge = -1;
        for (int fi = 0; fi < (int)faces.size(); ++fi) {
            const Face& f = faces[fi];
            const Vec3& A = vertices[f.v[0]].p;
            const Vec3& B = vertices[f.v[1]].p;
            const Vec3& C = vertices[f.v[2]].p;

            double d = in_triangle(A, B, C, P);
            if (d < 0) continue;  // outside this face

            if (d == 0) {
                // On boundary — find which edge
                double d0 = orientation(A, B, P);
                double d1 = orientation(B, C, P);
                double d2 = orientation(C, A, P);
                if (d0 == 0) boundary_edge = 0;
                else if (d1 == 0) boundary_edge = 1;
                else if (d2 == 0) boundary_edge = 2;
            }
            return fi;
        }
        return -1;
    }

    int insert_point(const Vec3& P)
    {
        int boundary_edge = -1;
        int fi = find_triangle(P, boundary_edge);

        if (fi == -1) {
            std::cout << "[insert_point] Point outside bounding box!\n";
            return -1;
        }

        if (boundary_edge == -1) {
            return split_triangle(fi, P);
        } else {
            return split_edge(fi, boundary_edge, P);
        }
    }

    // InCircle test (2D)
    double in_circle(const Vec3& A, const Vec3& B, const Vec3& C, const Vec3& D) const
    {
        double ax = A.x - D.x, ay = A.y - D.y;
        double bx = B.x - D.x, by = B.y - D.y;
        double cx = C.x - D.x, cy = C.y - D.y;

        return ax * (by * (cx*cx + cy*cy) - cy * (bx*bx + by*by))
            - ay * (bx * (cx*cx + cy*cy) - cx * (bx*bx + by*by))
            + (ax*ax + ay*ay) * (bx*cy - by*cx);
    }

    // Check if edge (fi, ei) is locally Delaunay
    bool is_locally_delaunay(int fi, int ei) const
    {
        int fj = faces[fi].neigh[ei];
        if (fj == -1) return true;  

        const Vec3& A = vertices[faces[fi].v[ei]].p;
        const Vec3& B = vertices[faces[fi].v[(ei+1)%3]].p;
        const Vec3& C = vertices[faces[fi].v[(ei+2)%3]].p;

        int ej = find_edge_index(faces[fj], faces[fi].v[ei], faces[fi].v[(ei+1)%3]);
        const Vec3& D = vertices[faces[fj].v[(ej+2)%3]].p;

        return in_circle(A, B, C, D) <= 0;
    }

    void make_delaunay()
    {
        bool flipped = true;
        while (flipped)
        {
            flipped = false;
            for (int fi = 0; fi < (int)faces.size(); ++fi)
            {
                for (int ei = 0; ei < 3; ++ei)
                {
                    if (!is_locally_delaunay(fi, ei))
                    {
                        flip_edge(fi, ei);
                        flipped = true;
                    }
                }
            }
        }
    }

    int insert_point_delaunay(const Vec3& P)
    {
        // naive insertion
        int new_vi = insert_point(P);
        if (new_vi == -1) return -1;

        // Collect all faces incident to new_vi
        bool flipped = true;
        while (flipped)
        {
            flipped = false;
            for (int fi = 0; fi < (int)faces.size(); ++fi)
            {
                const Face& f = faces[fi];
                // Check if this face contains new_vi
                bool has_new = (f.v[0] == new_vi ||
                                f.v[1] == new_vi ||
                                f.v[2] == new_vi);
                if (!has_new) continue;

                for (int ei = 0; ei < 3; ++ei)
                {
                    // Skip edges that touch new_vi
                    if (f.v[ei] == new_vi || f.v[(ei+1)%3] == new_vi) continue;

                    if (!is_locally_delaunay(fi, ei))
                    {
                        flip_edge(fi, ei);
                        flipped = true;
                        break;
                    }
                }
                if (flipped) break;
            }
        }

        return new_vi;
    }

    Vec3 circumcenter(int fi) const
    {
        const Vec3& A = vertices[faces[fi].v[0]].p;
        const Vec3& B = vertices[faces[fi].v[1]].p;
        const Vec3& C = vertices[faces[fi].v[2]].p;

        double ax = B.x - A.x, ay = B.y - A.y;
        double bx = C.x - A.x, by = C.y - A.y;

        double D = 2.0 * (ax * by - ay * bx);
        if (std::abs(D) < 1e-10)
            return { (A.x + B.x + C.x) / 3.0,
                    (A.y + B.y + C.y) / 3.0, 0.0 };

        double ux = (by * (ax*ax + ay*ay) - ay * (bx*bx + by*by)) / D;
        double uy = (ax * (bx*bx + by*by) - bx * (ax*ax + ay*ay)) / D;

        return { A.x + ux, A.y + uy, 0.0 };
    }

std::vector<std::pair<Vec3, Vec3>> build_voronoi_edges() const
{
    std::vector<std::pair<Vec3, Vec3>> edges;

    for (int fi = 0; fi < (int)faces.size(); ++fi)
    {
        for (int ei = 0; ei < 3; ++ei)
        {
            int fj = faces[fi].neigh[ei];
            if (fj == -1) continue;
            if (fj <= fi) continue;

            Vec3 ci = circumcenter(fi);
            Vec3 cj = circumcenter(fj);
            edges.push_back({ci, cj});
        }
    }

    return edges;
}

    bool save_off(const std::string& path) const {
        std::ofstream out(path);
        if (!out) return false;
        out << "OFF\n";
        out << vertices.size() << " " << faces.size() << " 0\n";
        for (const auto& v : vertices) out << v.p.x << " " << v.p.y << " " << v.p.z << "\n";
        for (const auto& f : faces) out << "3 " << f.v[0] << " " << f.v[1] << " " << f.v[2] << "\n";
        return true;
    }

    void print_neighbors(const std::string& name) const {
        std::cout << "\n== " << name << " ==\n";
        for (int fi = 0; fi < (int)faces.size(); ++fi) {
            const auto& f = faces[fi];
            std::cout << "F" << fi
                      << " v=(" << f.v[0] << "," << f.v[1] << "," << f.v[2] << ")"
                      << " neigh=(" << f.neigh[0] << "," << f.neigh[1] << "," << f.neigh[2] << ")\n";
        }
    }

    bool check_sewing_verbose(const std::string& name) const {
        bool ok = true;
        for (int fi = 0; fi < (int)faces.size(); ++fi) {
            for (int ei = 0; ei < 3; ++ei) {
                int fj = faces[fi].neigh[ei];
                if (fj == -1) continue;
                int a = faces[fi].v[ei];
                int b = faces[fi].v[(ei + 1) % 3];
                int ej = find_edge_index(faces[fj], a, b);
                if (ej == -1 || faces[fj].neigh[ej] != fi) {
                    std::cout << "[FAIL " << name << "] F" << fi
                              << " edge(" << a << "," << b
                              << ") -> neigh F" << fj
                              << " does not point back\n";
                    ok = false;
                }
            }
        }
        if (ok) std::cout << "[OK " << name << "] mutual sewing\n";
        return ok;
    }

    void build_adjacency() {
        std::map<std::pair<int,int>, std::pair<int,int>> edgeMap;
        for (int fi = 0; fi < (int)faces.size(); ++fi) {
            for (int ei = 0; ei < 3; ++ei) {
                int a = faces[fi].v[ei];
                int b = faces[fi].v[(ei+1)%3];
                std::pair<int,int> key = (a < b) ? std::make_pair(a,b)
                                                 : std::make_pair(b,a);
                auto it = edgeMap.find(key);
                if (it == edgeMap.end()) {
                    edgeMap[key] = std::make_pair(fi, ei);
                } else {
                    int fj = it->second.first;
                    int ej = it->second.second;
                    sew(fi, ei, fj, ej);
                }
            }
        }
    }

    bool load_off(const std::string& path) {
        std::ifstream in(path);
        if (!in) return false;

        std::string header;
        in >> header;
        if (header != "OFF") return false;

        int nV, nF, nE;
        in >> nV >> nF >> nE;

        vertices.clear();
        faces.clear();
        vertices.resize(nV);
        faces.resize(nF);

        for (int i = 0; i < nV; ++i)
            in >> vertices[i].p.x >> vertices[i].p.y >> vertices[i].p.z;

        for (int i = 0; i < nF; ++i) {
            int k, a, b, c;
            in >> k >> a >> b >> c;
            faces[i] = Face(a, b, c);
        }

        build_adjacency();
        update_oneFace_per_vertex();
        return true;
    }

    double triangle_area(int fi) const
    {
        const Face& f = faces[fi];
        const Vec3& a = vertices[f.v[0]].p;
        const Vec3& b = vertices[f.v[1]].p;
        const Vec3& c = vertices[f.v[2]].p;
        return 0.5 * norm(cross(b - a, c - a));
    }

    double vertex_area(int vi) const
    {
        double area = 0.0;
        for(int fi = 0; fi < (int)faces.size(); ++fi)
        {
            const Face& f = faces[fi];
            if(f.v[0]==vi || f.v[1]==vi || f.v[2]==vi)
                area += triangle_area(fi) / 3.0;
        }
        return area;
    }

    std::vector<double> compute_laplacian(const std::vector<double>& u) const
    {
        const int n = vertices.size();
        std::vector<double> L(n, 0.0);
        std::vector<double> A(n, 0.0);

        for(int fi = 0; fi < (int)faces.size(); ++fi)
        {
            const Face& f = faces[fi];
            const Vec3& p0 = vertices[f.v[0]].p;
            const Vec3& p1 = vertices[f.v[1]].p;
            const Vec3& p2 = vertices[f.v[2]].p;
            double area = 0.5 * norm(cross(p1 - p0, p2 - p0));
            A[f.v[0]] += area / 3.0;
            A[f.v[1]] += area / 3.0;
            A[f.v[2]] += area / 3.0;
        }

        for(int fi = 0; fi < (int)faces.size(); ++fi)
        {
            const Face& f = faces[fi];
            int i0 = f.v[0], i1 = f.v[1], i2 = f.v[2];
            const Vec3& p0 = vertices[i0].p;
            const Vec3& p1 = vertices[i1].p;
            const Vec3& p2 = vertices[i2].p;

            double cot0 = dot(p1-p0, p2-p0) / norm(cross(p1-p0, p2-p0));
            double cot1 = dot(p2-p1, p0-p1) / norm(cross(p2-p1, p0-p1));
            double cot2 = dot(p0-p2, p1-p2) / norm(cross(p0-p2, p1-p2));

            L[i1] += cot0 * (u[i2] - u[i1]);
            L[i2] += cot0 * (u[i1] - u[i2]);
            L[i2] += cot1 * (u[i0] - u[i2]);
            L[i0] += cot1 * (u[i2] - u[i0]);
            L[i0] += cot2 * (u[i1] - u[i0]);
            L[i1] += cot2 * (u[i0] - u[i1]);
        }

        for(int i = 0; i < n; ++i)
            if(A[i] > 1e-12)
                L[i] /= (2.0 * A[i]);

        return L;
    }

    // Ds = (Lx, Ly, Lz) = 2Hn => H = ||Ds|| / 2
    std::vector<double> compute_mean_curvature() const
    {
        const int n = vertices.size();

        std::vector<double> ux(n), uy(n), uz(n);
        for(int i = 0; i < n; ++i)
        {
            ux[i] = vertices[i].p.x;
            uy[i] = vertices[i].p.y;
            uz[i] = vertices[i].p.z;
        }

        std::vector<double> Lx = compute_laplacian(ux);
        std::vector<double> Ly = compute_laplacian(uy);
        std::vector<double> Lz = compute_laplacian(uz);

        std::vector<double> H(n);
        for(int i = 0; i < n; ++i)
            H[i] = std::sqrt(Lx[i]*Lx[i] + Ly[i]*Ly[i] + Lz[i]*Lz[i]) / 2.0;

        return H;
    }
};

MeshSimple make_tetra_simple() {
    MeshSimple m;
    m.vertices = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
    m.tris = { Tri{{0,2,1}}, Tri{{0,1,3}}, Tri{{1,2,3}}, Tri{{2,0,3}} };
    return m;
}

MeshSewn make_tetra_sewn() {
    MeshSewn m;
    m.vertices = { Vertex{{0,0,0}}, Vertex{{1,0,0}}, Vertex{{0,1,0}}, Vertex{{0,0,1}} };
    m.faces = { Face(0,2,1), Face(0,1,3), Face(1,2,3), Face(2,0,3) };
    m.sew(0,2, 1,0);
    m.sew(0,1, 2,0);
    m.sew(0,0, 3,0);
    m.sew(1,1, 2,2);
    m.sew(1,2, 3,1);
    m.sew(2,1, 3,2);
    m.update_oneFace_per_vertex();
    return m;
}

MeshSewn make_pyramid_sewn() {
    MeshSewn m;
    m.vertices = {
        Vertex{{-1,-1,0}}, Vertex{{ 1,-1,0}}, Vertex{{ 1, 1,0}}, Vertex{{-1, 1,0}},
        Vertex{{ 0, 0,1}}
    };
    m.faces = {
        Face(0,1,2), Face(0,2,3),
        Face(0,1,4), Face(1,2,4), Face(2,3,4), Face(3,0,4)
    };
    m.sew(0,2, 1,0);
    m.sew(0,0, 2,0);
    m.sew(0,1, 3,0);
    m.sew(1,1, 4,0);
    m.sew(1,2, 5,0);
    m.sew(2,1, 3,2);
    m.sew(3,1, 4,2);
    m.sew(4,1, 5,2);
    m.sew(5,1, 2,2);
    m.update_oneFace_per_vertex();
    return m;
}

MeshSewn make_bbox_infinite_sewn() {
    MeshSewn m;
    m.vertices = {
        Vertex{{0,0,0}}, Vertex{{1,0,0}}, Vertex{{1,1,0}}, Vertex{{0,1,0}},
        Vertex{{0.5,0.5,-5.0}}
    };
    m.faces = {
        Face(0,1,2), Face(0,2,3),
        Face(0,1,4), Face(1,2,4), Face(2,3,4), Face(3,0,4)
    };
    m.sew(0,2, 1,0);
    m.sew(0,0, 2,0);
    m.sew(0,1, 3,0);
    m.sew(1,1, 4,0);
    m.sew(1,2, 5,0);
    m.sew(2,1, 3,2);
    m.sew(3,1, 4,2);
    m.sew(4,1, 5,2);
    m.sew(5,1, 2,2);
    m.update_oneFace_per_vertex();
    return m;
}

// Viewer
enum DisplayMode {
    MODE_MESH = 0,
    MODE_HEAT = 1,
    MODE_CURVATURE = 2
};

class Viewer : public AppCamera
{
public:
    Viewer(MeshSewn& m)
        : AppCamera(1600, 900),
        ref(m),
        mode(MODE_MESH),
        vao(0), vbo_position(0), vbo_color(0), vertex_count(0),
        dt(5e-8),
        heat_initialized(false),
        curvature_initialized(false)

    {}

    int init()
    {
        camera().lookat(Point(0,0,0), 2.f);
        glEnable(GL_DEPTH_TEST);
        glClearColor(0.2f, 0.2f, 0.2f, 1.f);
        glDisable(GL_CULL_FACE);

        program = read_program("shaders/mesh_color.glsl");
        program_print_errors(program);

        init_gl_buffers();

        return 0;
    }

    int render()
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // keyboard controls
        static bool key1_pressed = false;
        static bool key2_pressed = false;
        static bool key3_pressed = false;

        if(key_state('1'))
        {
            if(!key1_pressed)
            {
                mode = MODE_MESH;
                reset_colors_to_blue();
                std::cout << "Mode: MESH\n";
                key1_pressed = true;
            }
        }
        else key1_pressed = false;

        if(key_state('2'))
        {
            if(!key2_pressed)
            {
                mode = MODE_HEAT;
                heat_initialized = false;
                std::cout << "Mode: HEAT\n";
                key2_pressed = true;
            }
        }
        else key2_pressed = false;

        if(key_state('3'))
        {
            if(!key3_pressed)
            {
                mode = MODE_CURVATURE;
                curvature_initialized = false;
                std::cout << "Mode: CURVATURE\n";
                key3_pressed = true;
            }
        }
        else key3_pressed = false;

        // simulation / curvature 
        if(mode == MODE_HEAT)
        {
            if(!heat_initialized)
                init_heat();

            static int frame = 0;
            frame++;
            heat_step();

            if(frame % 5 == 0)
                rebuild_mesh_with_temperature();
        }
        else if(mode == MODE_CURVATURE)
        {
            if(!curvature_initialized)
            {
                rebuild_mesh_with_curvature();
                curvature_initialized = true;
            }
        }
            
        static bool keyT_pressed = false;
        static bool keyE_pressed = false;
        static bool keyR_pressed = false;
        static bool keyF_pressed = false;

        if(key_state('t') || key_state('T'))
        {
            if(!keyT_pressed)
            {
                if(split_face_index < (int)ref.faces.size())
                {
                    const Face& f = ref.faces[split_face_index];
                    const Vec3& p0 = ref.vertices[f.v[0]].p;
                    const Vec3& p1 = ref.vertices[f.v[1]].p;
                    const Vec3& p2 = ref.vertices[f.v[2]].p;
                    Vec3 centroid = {
                        (p0.x + p1.x + p2.x) / 3.0,
                        (p0.y + p1.y + p2.y) / 3.0,
                        (p0.z + p1.z + p2.z) / 3.0
                    };
                    ref.split_triangle(split_face_index, centroid);
                    rebuild_position_buffer();
                    std::cout << "[T] split_triangle face " << split_face_index
                            << ": " << ref.faces.size() << " faces, "
                            << ref.vertices.size() << " vertices\n";
                    split_face_index++;
                }
                else
                    std::cout << "[T] всі грані вже розділені, натисни R для reset\n";

                keyT_pressed = true;
            }
        }
        else keyT_pressed = false;

        if(key_state('e') || key_state('E'))
        {
            if(!keyE_pressed)
            {
                if(split_face_index < (int)ref.faces.size())
                {
                    int fi = split_face_index, ei = 0;
                    int a = ref.faces[fi].v[ei];
                    int b = ref.faces[fi].v[(ei+1)%3];
                    Vec3 pa = ref.vertices[a].p;
                    Vec3 pb = ref.vertices[b].p;
                    Vec3 mid = {
                        (pa.x + pb.x) / 2.0,
                        (pa.y + pb.y) / 2.0,
                        (pa.z + pb.z) / 2.0
                    };
                    ref.split_edge(fi, ei, mid);
                    rebuild_position_buffer();
                    std::cout << "[E] split_edge face " << split_face_index
                            << ": " << ref.faces.size() << " faces, "
                            << ref.vertices.size() << " vertices\n";
                    split_face_index++;
                }
                else
                    std::cout << "[E] всі грані вже розділені, натисни R для reset\n";

                keyE_pressed = true;
            }
        }
        else keyE_pressed = false;

        if(key_state('r') || key_state('R'))
        {
            if(!keyR_pressed)
            {
                ref = make_pyramid_sewn();
                split_face_index = 0;
                flip_mode = false;
                flip_done = false;
                rebuild_position_buffer();
                std::cout << "[R] reset: "
                        << ref.faces.size() << " faces, "
                        << ref.vertices.size() << " vertices\n";
                keyR_pressed = true;
            }
        }
        else keyR_pressed = false;

        if (key_state('f') || key_state('F'))
        {
            if (!keyF_pressed)
            {
                int target_fi = -1, target_ei = -1;
                for (int fi = 0; fi < (int)ref.faces.size() && target_fi == -1; ++fi)
                {
                    const Face& f = ref.faces[fi];
                    for (int ei = 0; ei < 3 && target_fi == -1; ++ei)
                    {
                        int fj = f.neigh[ei];
                        if (fj == -1) continue;

                        auto base_vertex = [&](int v){ return ref.vertices[v].p.z == 0.0; };
                        bool fi_base = base_vertex(f.v[0]) && base_vertex(f.v[1]) && base_vertex(f.v[2]);
                        const Face& nf = ref.faces[fj];
                        bool fj_base = base_vertex(nf.v[0]) && base_vertex(nf.v[1]) && base_vertex(nf.v[2]);
                        if (fi_base && fj_base)
                        {
                            target_fi = fi;
                            target_ei = ei;
                        }
                    }
                }

                if (target_fi != -1)
                {
                    ref.flip_edge(target_fi, target_ei);
                    rebuild_position_buffer();
                    std::cout << "[F] Flipped base diagonal (face " << target_fi
                            << ", edge " << target_ei << ")\n";
                }
                else
                    std::cout << "[F] No base edge to flip.\n";

                keyF_pressed = true;
            }
        }
        else keyF_pressed = false;

        // transforms 
        Transform view = camera().view();
        Transform proj = camera().projection(window_width(), window_height(), 45);
        Transform mvp = proj * view * Identity();

        // shader setup
        glUseProgram(program);
        program_uniform(program, "mvpMatrix", mvp);
        program_uniform(program, "useVertexColor", (mode != MODE_MESH) ? 1 : 0);
        program_uniform(program, "baseColor", vec3(0.1f, 0.3f, 1.0f));

        // fill pass
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, vertex_count);
        glBindVertexArray(0);

        // wireframe pass
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_POLYGON_OFFSET_LINE);
        glPolygonOffset(-1.f, -1.f);
        program_uniform(program, "useVertexColor", 0);
        program_uniform(program, "baseColor", vec3(0.0f, 0.0f, 0.0f));
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, vertex_count);
        glBindVertexArray(0);
        glDisable(GL_POLYGON_OFFSET_LINE);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        return 1;
    }

    int quit()
    {
        glDeleteBuffers(1, &vbo_position);
        glDeleteBuffers(1, &vbo_color);
        glDeleteVertexArrays(1, &vao);
        release_program(program);
        return 0;
    }

private:

    // GPU buffers
    void init_gl_buffers()
    {
        vertex_count = (int)ref.faces.size() * 3;

        std::vector<float> positions;
        positions.reserve(vertex_count * 3);
        for(const auto& f : ref.faces)
            for(int i = 0; i < 3; i++)
            {
                const auto& p = ref.vertices[f.v[i]].p;
                positions.push_back((float)p.x);
                positions.push_back((float)p.y);
                positions.push_back((float)p.z);
            }

        color_buffer.assign(vertex_count * 4, 0.0f);
        for(int i = 0; i < vertex_count; i++)
        {
            color_buffer[i*4 + 0] = 0.1f;
            color_buffer[i*4 + 1] = 0.3f;
            color_buffer[i*4 + 2] = 1.0f;
            color_buffer[i*4 + 3] = 1.0f;
        }

        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        // VBO of positions — static
        glGenBuffers(1, &vbo_position);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_position);
        glBufferData(GL_ARRAY_BUFFER,
                     positions.size() * sizeof(float),
                     positions.data(),
                     GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);

        // VBO of colors — dynamic
        glGenBuffers(1, &vbo_color);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
        glBufferData(GL_ARRAY_BUFFER,
                     color_buffer.size() * sizeof(float),
                     color_buffer.data(),
                     GL_DYNAMIC_DRAW);
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(3);

        glBindVertexArray(0);
    }

    void update_color_buffer()
    {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
        glBufferSubData(GL_ARRAY_BUFFER, 0,
                        color_buffer.size() * sizeof(float),
                        color_buffer.data());
    }

    int split_face_index = 0;
    void rebuild_position_buffer()
    {
        vertex_count = (int)ref.faces.size() * 3;

        std::vector<float> positions;
        positions.reserve(vertex_count * 3);
        for(const auto& f : ref.faces)
            for(int i = 0; i < 3; i++)
            {
                const auto& p = ref.vertices[f.v[i]].p;
                positions.push_back((float)p.x);
                positions.push_back((float)p.y);
                positions.push_back((float)p.z);
            }

        color_buffer.assign(vertex_count * 4, 0.0f);
        for(int i = 0; i < vertex_count; i++)
        {
            color_buffer[i*4 + 0] = 0.1f;
            color_buffer[i*4 + 1] = 0.3f;
            color_buffer[i*4 + 2] = 1.0f;
            color_buffer[i*4 + 3] = 1.0f;
        }

        glBindBuffer(GL_ARRAY_BUFFER, vbo_position);
        glBufferData(GL_ARRAY_BUFFER,
                    positions.size() * sizeof(float),
                    positions.data(),
                    GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
        glBufferData(GL_ARRAY_BUFFER,
                    color_buffer.size() * sizeof(float),
                    color_buffer.data(),
                    GL_DYNAMIC_DRAW);
    }

    void reset_colors_to_blue()
    {
        for(int i = 0; i < vertex_count; i++)
        {
            color_buffer[i*4 + 0] = 0.1f;
            color_buffer[i*4 + 1] = 0.3f;
            color_buffer[i*4 + 2] = 1.0f;
            color_buffer[i*4 + 3] = 1.0f;
        }
        update_color_buffer();
    }

    // Heat diffusion
    void init_heat()
    {
        int n = ref.vertices.size();
        temperature.assign(n, 0.0);
        vertex_colors.assign(n, Color(0, 0, 1));
        temperature[0] = 1.0;
        heat_initialized = true;
    }

    void heat_step()
    {
        std::vector<double> L = ref.compute_laplacian(temperature);
        for(int i = 0; i < (int)temperature.size(); ++i)
            temperature[i] += dt * L[i];

        temperature[0] = 1.0;
    }

    void rebuild_mesh_with_temperature()
    {
        double minv = *std::min_element(temperature.begin(), temperature.end());
        double maxv = *std::max_element(temperature.begin(), temperature.end());

        for(int i = 0; i < (int)temperature.size(); i++)
        {
            double normt = (maxv - minv > 1e-12)
                // ? (temperature[i] - minv) / (maxv - minv)
                ? std::pow((temperature[i] - minv) / (maxv - minv), 0.3)
                : 0.0;
            vertex_colors[i] = Color((float)normt,
                                      0.2f * (1.0f - (float)normt),
                                      1.0f - (float)normt);
        }

        int idx = 0;
        for(const auto& f : ref.faces)
            for(int i = 0; i < 3; i++)
            {
                int vid = f.v[i];
                color_buffer[idx*4 + 0] = vertex_colors[vid].r;
                color_buffer[idx*4 + 1] = vertex_colors[vid].g;
                color_buffer[idx*4 + 2] = vertex_colors[vid].b;
                color_buffer[idx*4 + 3] = 1.0f;
                idx++;
            }

        update_color_buffer();
    }

    // Mean curvature
    Color hsv_to_rgb(float h, float s, float v)
    {
        float r=0, g=0, b=0;
        int i = (int)(h * 6.0f);
        float f = h * 6.0f - i;
        float p = v * (1.0f - s);
        float q = v * (1.0f - f * s);
        float t = v * (1.0f - (1.0f - f) * s);
        switch(i % 6)
        {
            case 0: r=v; g=t; b=p; break;
            case 1: r=q; g=v; b=p; break;
            case 2: r=p; g=v; b=t; break;
            case 3: r=p; g=q; b=v; break;
            case 4: r=t; g=p; b=v; break;
            case 5: r=v; g=p; b=q; break;
        }
        return Color(r, g, b);
    }

    void rebuild_mesh_with_curvature()
    {
        std::cout << "Computing mean curvature...\n";
        std::vector<double> H = ref.compute_mean_curvature();

        std::vector<double> sorted_H = H;
        std::sort(sorted_H.begin(), sorted_H.end());
        double clamp_min = sorted_H[(int)(sorted_H.size() * 0.05)];
        double clamp_max = sorted_H[(int)(sorted_H.size() * 0.95)];
        std::cout << "Curvature range: [" << clamp_min << ", " << clamp_max << "]\n";

        int idx = 0;
        for(const auto& f : ref.faces)
            for(int i = 0; i < 3; i++)
            {
                int vid = f.v[i];
                double normH = (clamp_max - clamp_min > 1e-12)
                    ? std::clamp((H[vid] - clamp_min) / (clamp_max - clamp_min), 0.0, 1.0)
                    : 0.0;

                float hue = 0.66f * (1.0f - (float)normH);
                Color c = hsv_to_rgb(hue, 1.0f, 1.0f);

                color_buffer[idx*4 + 0] = c.r;
                color_buffer[idx*4 + 1] = c.g;
                color_buffer[idx*4 + 2] = c.b;
                color_buffer[idx*4 + 3] = 1.0f;
                idx++;
            }

        update_color_buffer();
        std::cout << "Curvature mesh built.\n";
    }

    void highlight_non_delaunay_edges()
    {
        // Reset all to blue
        for (int i = 0; i < vertex_count; i++) {
            color_buffer[i*4+0] = 0.1f;
            color_buffer[i*4+1] = 0.3f;
            color_buffer[i*4+2] = 1.0f;
            color_buffer[i*4+3] = 1.0f;
        }

        int idx = 0;
        for (int fi = 0; fi < (int)ref.faces.size(); ++fi) {
            bool bad = false;
            for (int ei = 0; ei < 3; ++ei)
                if (!ref.is_locally_delaunay(fi, ei)) bad = true;

            float r = bad ? 1.0f : 0.1f;
            float g = bad ? 0.0f : 0.3f;
            float b = bad ? 0.0f : 1.0f;

            for (int i = 0; i < 3; i++) {
                color_buffer[idx*4+0] = r;
                color_buffer[idx*4+1] = g;
                color_buffer[idx*4+2] = b;
                color_buffer[idx*4+3] = 1.0f;
                idx++;
            }
        }
        update_color_buffer();
    }

    // Private members
    MeshSewn& ref;
    GLuint program;
    DisplayMode mode;

    // GPU buffers
    GLuint vao;
    GLuint vbo_position;
    GLuint vbo_color;
    int vertex_count;
    std::vector<float> color_buffer;

    // heat
    std::vector<double> temperature;
    std::vector<Color> vertex_colors;
    double dt;
    bool heat_initialized;

    // curvature
    bool curvature_initialized;

    // flip mode
    bool flip_mode = false;
    bool flip_done = false;
    std::vector<Vec3> flip_pts = {
        {0.2, 0.2, 0}, {0.7, 0.1, 0}, {0.5, 0.7, 0},
        {0.3, 0.5, 0}, {0.6, 0.4, 0}, {0.1, 0.8, 0},
        {0.8, 0.6, 0}, {0.4, 0.3, 0}
    };
};

class DelaunayViewer : public AppCamera
{
public:
    DelaunayViewer()
        : AppCamera(1600, 900),
        vao(0), vbo_position(0), vbo_color(0),
        vertex_count(0), naive_index(0), delaunay_index(0), delaunay_done(false),
        vbo_voronoi(0), vao_voronoi(0), voronoi_vertex_count(0), show_voronoi(false)
    {
        srand(static_cast<unsigned>(time(nullptr)));
        for (int i = 0; i < 25; i++)
            pts.push_back({
                0.1 + 0.8 * (rand() / (double)RAND_MAX),
                0.1 + 0.8 * (rand() / (double)RAND_MAX),
                0.0
            });
    }

    int init()
    {
        camera().lookat(Point(0.5f, 0.5f, 0.f), 1.5f);
        glEnable(GL_DEPTH_TEST);
        glClearColor(0.2f, 0.2f, 0.2f, 1.f);
        glDisable(GL_CULL_FACE);

        program = read_program("shaders/mesh_color.glsl");
        program_print_errors(program);

        reset_mesh();
        init_gl_buffers();

        // Voronoi VBO
        glGenBuffers(1, &vbo_voronoi);
        glGenVertexArrays(1, &vao_voronoi);
        glBindVertexArray(vao_voronoi);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_voronoi);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
        return 0;
    }

    int render()
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        static bool keyI_pressed = false;
        if (key_state('i') || key_state('I'))
        {
            if (!keyI_pressed)
            {
                if (point_index < (int)pts.size())
                {
                    mesh.insert_point(pts[point_index]);
                    std::cout << "[I] Naive insert point " << point_index
                            << " | faces: " << mesh.faces.size() << "\n";
                    point_index++;
                    delaunay_done = false;
                    rebuild_buffers();
                    highlight_non_delaunay();
                }
                else
                    std::cout << "[I] All points inserted. Press R to reset.\n";
                
                if (show_voronoi) rebuild_voronoi_buffer();
                keyI_pressed = true;
            }
        }
        else keyI_pressed = false;
    
        static bool keyD_pressed = false;
        if (key_state('d') || key_state('D'))
        {
            if (!keyD_pressed)
            {
                if (point_index < (int)pts.size())
                {
                    mesh.insert_point_delaunay(pts[point_index]);
                    std::cout << "[D] Delaunay insert point " << point_index
                            << " | faces: " << mesh.faces.size() << "\n";
                    point_index++;
                    delaunay_done = false;
                    rebuild_buffers();
                    highlight_non_delaunay();
                }
                else
                    std::cout << "[D] All points inserted. Press R to reset.\n";
                
                if (show_voronoi) rebuild_voronoi_buffer();
                keyD_pressed = true;
            }
        }
        else keyD_pressed = false;

        static bool keySpace_pressed = false;
        if (key_state(' '))
        {
            if (!keySpace_pressed)
            {
                if (delaunay_done)
                {
                    std::cout << "[Space] Already Delaunay.\n";
                }
                else
                {
                    bool found = false;
                    for (int fi = 0; fi < (int)mesh.faces.size() && !found; ++fi)
                        for (int ei = 0; ei < 3 && !found; ++ei)
                            if (!mesh.is_locally_delaunay(fi, ei))
                            {
                                mesh.flip_edge(fi, ei);
                                std::cout << "[Space] Flipped edge (face " << fi
                                          << ", edge " << ei << ")\n";
                                found = true;
                            }

                    if (!found)
                    {
                        delaunay_done = true;
                        std::cout << "[Space] Delaunay complete!\n";
                    }

                    rebuild_buffers();
                    highlight_non_delaunay();
                }
                keySpace_pressed = true;
            }
        }
        else keySpace_pressed = false;

        // R — reset
        static bool keyR_pressed = false;
        if (key_state('r') || key_state('R'))
        {
            if (!keyR_pressed)
            {
                reset_mesh();
                rebuild_buffers();
                std::cout << "[R] Reset.\n";
                keyR_pressed = true;
            }
        }
        else keyR_pressed = false;

        static bool keyV_pressed = false;
        if (key_state('n') || key_state('N'))
        {
            if (!keyV_pressed)
            {
                show_voronoi = !show_voronoi;
                if (show_voronoi) rebuild_voronoi_buffer();
                std::cout << "[V] Voronoi " << (show_voronoi ? "ON" : "OFF") << "\n";
                keyV_pressed = true;
            }
        }
        else keyV_pressed = false;

        // Render
        Transform view = camera().view();
        Transform proj = camera().projection(window_width(), window_height(), 45);
        Transform mvp = proj * view * Identity();

        glUseProgram(program);
        program_uniform(program, "mvpMatrix", mvp);
        program_uniform(program, "useVertexColor", 1);
        program_uniform(program, "baseColor", vec3(0.1f, 0.3f, 1.0f));

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, vertex_count);
        glBindVertexArray(0);

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_POLYGON_OFFSET_LINE);
        glPolygonOffset(-1.f, -1.f);
        program_uniform(program, "useVertexColor", 0);
        program_uniform(program, "baseColor", vec3(0.0f, 0.0f, 0.0f));
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, vertex_count);
        glBindVertexArray(0);
        glDisable(GL_POLYGON_OFFSET_LINE);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // Voronoi overlay
        if (show_voronoi && voronoi_vertex_count > 0)
        {
            program_uniform(program, "useVertexColor", 0);
            program_uniform(program, "baseColor", vec3(1.0f, 0.85f, 0.0f));

            glBindVertexArray(vao_voronoi);
            glLineWidth(2.0f);
            glDrawArrays(GL_LINES, 0, voronoi_vertex_count);
            glLineWidth(1.0f);
            glBindVertexArray(0);
        }
        return 1;
    }

    int quit()
    {
        glDeleteBuffers(1, &vbo_voronoi);
        glDeleteBuffers(1, &vbo_position);
        glDeleteBuffers(1, &vbo_color);
        glDeleteVertexArrays(1, &vao);
        glDeleteVertexArrays(1, &vao_voronoi);
        release_program(program);
        return 0;
    }

private:
    void reset_mesh()
    {
        mesh.make_bounding_box(pts, 0.5);
        point_index = 0;
        delaunay_done = false;
    }

    void init_gl_buffers()
    {
        vertex_count = (int)mesh.faces.size() * 3;

        std::vector<float> positions;
        for (const auto& f : mesh.faces)
            for (int i = 0; i < 3; i++) {
                const auto& p = mesh.vertices[f.v[i]].p;
                positions.push_back((float)p.x);
                positions.push_back((float)p.y);
                positions.push_back((float)p.z);
            }

        color_buffer.assign(vertex_count * 4, 0.0f);
        for (int i = 0; i < vertex_count; i++) {
            color_buffer[i*4+0] = 0.1f;
            color_buffer[i*4+1] = 0.3f;
            color_buffer[i*4+2] = 1.0f;
            color_buffer[i*4+3] = 1.0f;
        }

        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        glGenBuffers(1, &vbo_position);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_position);
        glBufferData(GL_ARRAY_BUFFER, positions.size()*sizeof(float), positions.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);

        glGenBuffers(1, &vbo_color);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
        glBufferData(GL_ARRAY_BUFFER, color_buffer.size()*sizeof(float), color_buffer.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(3);

        glBindVertexArray(0);
    }

    void rebuild_buffers()
    {
        vertex_count = (int)mesh.faces.size() * 3;

        std::vector<float> positions;
        for (const auto& f : mesh.faces)
            for (int i = 0; i < 3; i++) {
                const auto& p = mesh.vertices[f.v[i]].p;
                positions.push_back((float)p.x);
                positions.push_back((float)p.y);
                positions.push_back((float)p.z);
            }

        color_buffer.assign(vertex_count * 4, 0.0f);
        for (int i = 0; i < vertex_count; i++) {
            color_buffer[i*4+0] = 0.1f;
            color_buffer[i*4+1] = 0.3f;
            color_buffer[i*4+2] = 1.0f;
            color_buffer[i*4+3] = 1.0f;
        }

        glBindBuffer(GL_ARRAY_BUFFER, vbo_position);
        glBufferData(GL_ARRAY_BUFFER, positions.size()*sizeof(float), positions.data(), GL_DYNAMIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
        glBufferData(GL_ARRAY_BUFFER, color_buffer.size()*sizeof(float), color_buffer.data(), GL_DYNAMIC_DRAW);
    }

    void highlight_non_delaunay()
    {
        int idx = 0;
        for (int fi = 0; fi < (int)mesh.faces.size(); ++fi)
        {
            bool bad = false;
            for (int ei = 0; ei < 3; ++ei)
                if (!mesh.is_locally_delaunay(fi, ei)) bad = true;

            float r = bad ? 0.7f : 0.1f;
            float g = bad ? 0.2f : 0.3f;
            float b = bad ? 0.2f : 1.0f;

            for (int i = 0; i < 3; i++) {
                color_buffer[idx*4+0] = r;
                color_buffer[idx*4+1] = g;
                color_buffer[idx*4+2] = b;
                color_buffer[idx*4+3] = 1.0f;
                idx++;
            }
        }
        glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
        glBufferData(GL_ARRAY_BUFFER, color_buffer.size()*sizeof(float), color_buffer.data(), GL_DYNAMIC_DRAW);
    }

    void rebuild_voronoi_buffer()
    {
        auto edges = mesh.build_voronoi_edges();
        std::cout << "Total faces: " << mesh.faces.size() 
                << " Voronoi edges: " << edges.size() << "\n";

        for (int fi = 0; fi < (int)mesh.faces.size(); ++fi)
        {
            Vec3 c = mesh.circumcenter(fi);
            std::cout << "Face " << fi << " circumcenter: (" 
                    << c.x << ", " << c.y << ")\n";
        }
        voronoi_vertex_count = (int)edges.size() * 2;

        std::vector<float> positions;
        positions.reserve(voronoi_vertex_count * 3);
        for (const auto& e : edges)
        {
            positions.push_back((float)e.first.x);
            positions.push_back((float)e.first.y);
            positions.push_back((float)e.first.z);
            positions.push_back((float)e.second.x);
            positions.push_back((float)e.second.y);
            positions.push_back((float)e.second.z);
        }

        glBindBuffer(GL_ARRAY_BUFFER, vbo_voronoi);
        glBufferData(GL_ARRAY_BUFFER,
                    positions.size() * sizeof(float),
                    positions.data(),
                    GL_DYNAMIC_DRAW);
    }

    MeshSewn mesh;
    std::vector<Vec3> pts;

    GLuint program;
    GLuint vao, vbo_position, vbo_color;
    GLuint vbo_voronoi;
    GLuint vao_voronoi;
    int vertex_count;
    int naive_index;
    int delaunay_index;
    bool delaunay_done;
    int point_index;
    int voronoi_vertex_count;
    bool show_voronoi;
    std::vector<float> color_buffer;
};

int main()
{
    int choice;
    std::cout << "1 - create meshes\n";
    std::cout << "2 - load mesh and open viewer (TP2-3)\n";
    std::cout << "3 - split operations (TP4)\n";
    std::cout << "4 - Delaunay and Voronoi\n";
    std::cin >> choice;

    if(choice == 1)
    {
        auto tet = make_tetra_sewn();
        auto pyr = make_pyramid_sewn();
        auto box = make_bbox_infinite_sewn();

        tet.print_neighbors("tetra");
        tet.check_sewing_verbose("tetra");

        pyr.print_neighbors("pyramid");
        pyr.check_sewing_verbose("pyramid");

        box.print_neighbors("bbox_infinite");
        box.check_sewing_verbose("bbox_infinite");

        tet.save_off("data/tetra.off");
        pyr.save_off("data/pyramid.off");
        box.save_off("data/bbox_infinite.off");

        {
            std::cout << "\n========== TEST split_triangle ==========\n";
            auto m1 = make_tetra_sewn();
            std::cout << "Initial tetrahedron:\n";
            m1.print_neighbors("before_split");
            int new_v = m1.split_triangle(0, {0.33, 0.33, 0.0});
            std::cout << "Created new vertex index: " << new_v << "\n";
            std::cout << "After split_triangle:\n";
            m1.check_sewing_verbose("split_triangle");
            std::cout << "Vertices count: " << m1.vertices.size() << " (was 4)\n";
            std::cout << "Faces count: " << m1.faces.size() << " (was 4)\n";
        }
        
        {
            std::cout << "\n========== TEST split_edge ==========\n";
            auto m2 = make_pyramid_sewn();
            int fi = 0, ei = 0;
            int a = m2.faces[fi].v[ei];
            int b = m2.faces[fi].v[(ei+1)%3];
            Vec3 pa = m2.vertices[a].p;
            Vec3 pb = m2.vertices[b].p;
            Vec3 mid = {(pa.x+pb.x)/2, (pa.y+pb.y)/2, (pa.z+pb.z)/2};
            m2.split_edge(fi, ei, mid);
            m2.check_sewing_verbose("split_edge");
            m2.save_off("data/split_edge.off");
        }                    
        
        {
            std::cout << "\n========== TEST flip_edge ==========\n";
            auto m3 = make_pyramid_sewn();
            std::cout << "Initial pyramid:\n";
            m3.print_neighbors("before_flip");
            
            std::cout << "Attempting to flip edge (face 0, edge 2)\n";
            std::cout << "Face 0 vertices: (" 
                      << m3.faces[0].v[0] << ","
                      << m3.faces[0].v[1] << ","
                      << m3.faces[0].v[2] << ")\n";
            std::cout << "Face 0 neighbor at edge 2: " << m3.faces[0].neigh[2] << "\n";
            
            m3.flip_edge(0, 2);
            
            std::cout << "\nAfter flip_edge:\n";
            m3.check_sewing_verbose("flip_edge");
            
            std::cout << "\nDetailed neighbor check after flip:\n";
            for (size_t i = 0; i < m3.faces.size(); i++) {
                std::cout << "Face " << i << " (";
                for (int j = 0; j < 3; j++) {
                    std::cout << m3.faces[i].v[j];
                    if (j < 2) std::cout << ",";
                }
                std::cout << ") neighbors: [";
                for (int j = 0; j < 3; j++) {
                    std::cout << m3.faces[i].neigh[j];
                    if (j < 2) std::cout << ",";
                }
                std::cout << "]\n";
            }
            
            std::cout << "\nSymmetry check:\n";
            bool ok = true;
            for (size_t i = 0; i < m3.faces.size(); i++) {
                for (int j = 0; j < 3; j++) {
                    int n = m3.faces[i].neigh[j];
                    if (n != -1) {
                        int ej = m3.find_edge_index(m3.faces[n], 
                                                    m3.faces[i].v[(j+1)%3], 
                                                    m3.faces[i].v[j]);
                        if (ej == -1 || m3.faces[n].neigh[ej] != (int)i) {
                            std::cout << "ASYMMETRY: Face " << i 
                                      << "edge " << j << " -> face " << n 
                                      << ", but no back reference\n";
                            ok = false;
                        }
                    }
                }
            }
            if (ok) {
                std::cout << "All neighbors are symmetric\n";
            }            
        }
        {
            MeshSewn m;

            Vec3 A{0, 0, 0};
            Vec3 B{1, 0, 0};
            Vec3 C{0, 1, 0};
            Vec3 D{0, -1, 0};
            Vec3 E{2, 0, 0};

            double r1 = m.orientation(A, B, C);
            double r2 = m.orientation(A, B, D);
            double r3 = m.orientation(A, B, E);

            std::cout << "\n== orientation test ==\n";
            std::cout << "A,B,C (CCW expected > 0): " << r1
                    << (r1 > 0 ? "  [OK]" : "  [FAIL]") << "\n";
            std::cout << "A,B,D (CW  expected < 0): " << r2
                    << (r2 < 0 ? "  [OK]" : "  [FAIL]") << "\n";
            std::cout << "A,B,E (collinear = 0)   : " << r3
                    << (r3 == 0 ? "  [OK]" : "  [FAIL]") << "\n";
        }

        {
            MeshSewn m;

            Vec3 A{0, 0, 0};
            Vec3 B{4, 0, 0};
            Vec3 C{0, 4, 0};

            Vec3 P_inside {1, 1, 0};
            Vec3 P_outside {3, 3, 0};
            Vec3 P_boundary {2, 0, 0};
            Vec3 P_vertex {0, 0, 0};

            double r1 = m.in_triangle(A, B, C, P_inside);
            double r2 = m.in_triangle(A, B, C, P_outside);
            double r3 = m.in_triangle(A, B, C, P_boundary);
            double r4 = m.in_triangle(A, B, C, P_vertex);

            std::cout << "\n== in_triangle test ==\n";
            std::cout << "P inside (expected > 0): " << r1
                    << (r1 > 0 ? "[OK]" : "[FAIL]") << "\n";
            std::cout << "P outside (expected < 0): " << r2
                    << (r2 < 0 ? "[OK]" : "[FAIL]") << "\n";
            std::cout << "P boundary (expected = 0): " << r3
                    << (r3 == 0 ? "[OK]" : "[FAIL]") << "\n";
            std::cout << "P vertex (expected = 0): " << r4
                    << (r4 == 0 ? "[OK]" : "[FAIL]") << "\n";
        }

        {
            MeshSewn m;

            // Points to insert
            std::vector<Vec3> pts = {
                {0.2, 0.2, 0}, {0.7, 0.2, 0}, {0.5, 0.7, 0},
                {0.3, 0.5, 0}, {0.6, 0.5, 0}
            };

            m.make_bounding_box(pts, 0.5);
            std::cout << "\n== naive insertion test ==\n";
            std::cout << "Initial: " << m.vertices.size()
                    << " vertices, " << m.faces.size() << " faces\n";

            for (int i = 0; i < (int)pts.size(); ++i) {
                int vi = m.insert_point(pts[i]);
                std::cout << "Insert point " << i
                        << " -> vertex " << vi
                        << " | faces: " << m.faces.size() << "\n";
            }

            bool ok = m.check_sewing_verbose("naive_insertion");
            // m.save_off("data/naive_insertion.off");
        }

        {
            std::cout << "\n== Delaunay test ==\n";

            std::vector<Vec3> pts = {
                {0.2, 0.2, 0}, {0.7, 0.1, 0}, {0.5, 0.7, 0},
                {0.3, 0.5, 0}, {0.6, 0.4, 0}, {0.1, 0.8, 0},
                {0.8, 0.6, 0}, {0.4, 0.3, 0}
            };

            // Test 1: make_delaunay on naive triangulation
            {
                MeshSewn m;
                m.make_bounding_box(pts, 0.5);
                for (const auto& p : pts)
                    m.insert_point(p);

                std::cout << "Before make_delaunay:\n";
                int bad = 0;
                for (int fi = 0; fi < (int)m.faces.size(); ++fi)
                    for (int ei = 0; ei < 3; ++ei)
                        if (!m.is_locally_delaunay(fi, ei)) ++bad;
                std::cout << "Non-Delaunay edges: " << bad << "\n";

                m.make_delaunay();

                bad = 0;
                for (int fi = 0; fi < (int)m.faces.size(); ++fi)
                    for (int ei = 0; ei < 3; ++ei)
                        if (!m.is_locally_delaunay(fi, ei)) ++bad;
                std::cout << "After make_delaunay:\n";
                std::cout << "Non-Delaunay edges: " << bad
                        << (bad == 0 ? "[OK]" : "[FAIL]") << "\n";
                m.check_sewing_verbose("make_delaunay");
                m.save_off("data/delaunay.off");
            }

            // Test 2: incremental insert_point_delaunay
            {
                MeshSewn m;
                m.make_bounding_box(pts, 0.5);

                for (int i = 0; i < (int)pts.size(); ++i)
                {
                    m.insert_point_delaunay(pts[i]);

                    // Check Delaunay property after each insertion
                    int bad = 0;
                    for (int fi = 0; fi < (int)m.faces.size(); ++fi)
                        for (int ei = 0; ei < 3; ++ei)
                            if (!m.is_locally_delaunay(fi, ei)) ++bad;

                    std::cout << "After insert " << i
                            << ": faces=" << m.faces.size()
                            << " bad_edges=" << bad
                            << (bad == 0 ? "[OK]" : "[FAIL]") << "\n";
                }

                m.check_sewing_verbose("insert_point_delaunay");
                m.save_off("data/delaunay_incremental.off");
            }
        }

        std::cout << "\nMeshes created and saved.\n";
        return 0;
    }

    if(choice == 2)
    {
        MeshSewn m;
        m.load_off("data/queen.off");

        Viewer viewer(m);
        viewer.run();
    }

    if(choice == 3)
    {
        MeshSewn m = make_pyramid_sewn();
        std::cout << "Controls:\n";
        std::cout << "T - split_triangle\n";
        std::cout << "E - split_edge\n";
        std::cout << "F - edge flip\n";
        std::cout << "R - reset to original mesh\n";

        Viewer viewer(m);
        viewer.run();
    }

    if(choice == 4)
    {
        std::cout << "Controls:\n";
        std::cout << "I - naive insert next point\n";
        std::cout << "D - insert next point (Delaunay)\n";
        std::cout << "Space - one Lawson flip step\n";
        std::cout << "R - reset\n";
        std::cout << "N - toggle Voronoi overlay\n";

        DelaunayViewer viewer;
        viewer.run();
    }

    return 0;
}