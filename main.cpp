#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <algorithm>

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
};

int main()
{
    int choice;
    std::cout << "1 - create test meshes only\n";
    std::cout << "2 - load queen and open viewer\n";
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

    return 0;
}