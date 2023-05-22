#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
using namespace std;
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/* define appropriate data types
(Color, RayType, Sphere) */

class Color {
public:
    float red;
    float green;
    float blue;
    float eta;
    void set_color(float, float, float);
    void set_eta(float);
};

void Color::set_color(float r, float g, float b) {
    red = r;
    green = g;
    blue = b;
}

void Color::set_eta(float eta) {
    this->eta = eta;
}

class MtlColor {
public:
    Color dc;
    Color sc;
    float ka;
    float kd;
    float ks;
    float n;
    float m_alpha;
    float m_eta;
    void set_mtl(Color, Color, float, float, float, float, float, float);
};

void MtlColor::set_mtl(Color dc, Color sc, float ka, float kd, float ks, float n, float m_alpha, float m_eta) {
    this->ka = ka;
    this->kd = kd;
    this->ks = ks;
    this->n = n;
    this->dc.set_color(dc.red, dc.green, dc.blue);
    this->sc.set_color(sc.red, sc.green, sc.blue);
    this->m_alpha = m_alpha;
    this->m_eta = m_eta;
}

class RayType {
public:
    float x;
    float y;
    float z;
    float dx;
    float dy;
    float dz;
    RayType(float, float, float);
    void set_dir(float, float, float);
};

RayType::RayType(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

void RayType::set_dir(float dx, float dy, float dz) {
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
}

class Light {
public:
    float x;
    float y;
    float z;
    float w;
    Color lc;
    float c1;
    float c2;
    float c3;
    bool lsa;
    Light(float, float, float, float, Color);
    Light(float, float, float, float, Color, float, float, float);
    Light();
};

Light::Light() {
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->w = 0.0;
    this->c1 = 0.0;
    this->c2 = 0.0;
    this->c3 = 0.0;
}

Light::Light(float x, float y, float z, float w, Color lc) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
    this->lsa = false;
    this->lc.set_color(lc.red, lc.green, lc.blue);
}

Light::Light(float x, float y, float z, float w, Color lc, float c1, float c2, float c3) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
    this->c1 = c1;
    this->c2 = c2;
    this->c3 = c3;
    this->lsa = true;
    this->lc.set_color(lc.red, lc.green, lc.blue);
}

class Sphere {
public:
    float center_x;
    float center_y;
    float center_z;
    float radius;
    int id;
    MtlColor m;
    vector<vector<vector<float> > > texture;
    void set_center(float, float, float, float, MtlColor, vector<vector<vector<float> > >, int);
};

void Sphere::set_center(float x, float y, float z, float r, MtlColor c, vector<vector<vector<float> > > texture, int id) {
    center_x = x;
    center_y = y;
    center_z = z;
    radius = r;
    this->id = id;
    this->m.set_mtl(c.dc, c.sc, c.ka, c.kd, c.ks, c.n, c.m_alpha, c.m_eta);
    this->texture = texture;
}

class Triangle {
public:
    vector<vector<float> > tri;
    vector<vector<float> > vn;
    vector<vector<float> > vt;
    vector<vector<vector<float> > > texture;
    MtlColor m;
    float alpha;
    float beta;
    float gamma;
    int id;
    vector<float> n;
    void set_face(vector<float>, vector<vector<float> >, MtlColor, vector<vector<vector<float> > >);
    void set_id(int);
    void set_vn(vector<float>, vector<vector<float> >);
    void set_vt(vector<float>, vector<vector<float> >);
    void set_n(vector<float>);
    void set_bcoordinates(float, float, float);
};


void Triangle::set_face(vector<float> face, vector<vector<float> > vertex, MtlColor c, vector<vector<vector<float> > > texture) {
    for (int i = 0; i < (int)face.size(); i++) {
        tri.push_back(vertex[face[i] - 1]);
    }
    this->m.set_mtl(c.dc, c.sc, c.ka, c.kd, c.ks, c.n, c.m_alpha, c.m_eta);
    this->texture = texture;
}

void Triangle::set_vn(vector<float> face, vector<vector<float> > vn) {
    for (int i = 0; i < (int)face.size(); i++) {
        this->vn.push_back(vn[face[i] - 1]);
    }
}

void Triangle::set_vt(vector<float> face, vector<vector<float> > vt) {
    for (int i = 0; i < (int)face.size(); i++) {
        this->vt.push_back(vt[face[i] - 1]);
    }
}

void Triangle::set_n(vector<float> n) {
    this->n = n;
}

void Triangle::set_id(int id) {
    this->id = id;
}

void Triangle::set_bcoordinates(float alpha, float beta, float gamma) {
    this->alpha = alpha;
    this->beta = beta;
    this->gamma = gamma;
}

// global variables
MtlColor mtlcolor;
vector<vector<vector<float> > > texture;
vector<Sphere> sphere_vector;
vector<Light> light_vector;
vector<Triangle> tri_vector;
vector<vector<float> > face_v_vector;
vector<vector<float> > face_vn_vector;
vector<vector<float> > face_vt_vector;
vector<vector<float> > vn_vector;
vector<vector<float> > vt_vector;
vector<vector<float> > vertex_vector;
vector<float> viewdir;
Color bkgcolor;
int width;
int height;
vector<float> eye_vector;
vector<float> updir;
vector<float> depth_cue;
bool depth = false;
int sphere_id = 0;
int triangle_id = 0;
int texture_width = 0;
int texture_height = 0;

#define INT_EYE 0
#define INT_SHADOW 1
#define INT_MIRROR 2
#define INT_TRANS_IN 3
#define INT_TRANS_OUT 4

// function for calculating the cross product
vector<float> cProduct(vector<float> a, vector<float> b) {
    vector<float> r(3, 0);
    r.at(0) = a.at(1) * b.at(2) - a.at(2) * b.at(1);
    r.at(1) = a.at(2) * b.at(0) - a.at(0) * b.at(2);
    r.at(2) = a.at(0) * b.at(1) - a.at(1) * b.at(0);
    return r;
}
// function for calculating the dot product
float dProduct(vector<float> a, vector<float> b) {
    return (a.at(0) * b.at(0) + a.at(1) * b.at(1) + a.at(2) * b.at(2));
}
// function for multiplying vectors with weight a
vector<float> vMultiply(float a, vector<float> b) {
    vector<float> r(3, 0);
    r.at(0) = a * b.at(0);
    r.at(1) = a * b.at(1);
    r.at(2) = a * b.at(2);
    return r;
}

// function for matrix multiplication
vector<float> vMul(vector<float> a, vector<float> b) {
    vector<float> r(3, 0);
    r.at(0) = a.at(0) * b.at(0);
    r.at(1) = a.at(1) * b.at(1);
    r.at(2) = a.at(2) * b.at(2);
    return r;
}

// function for adding vectors
vector<float> vAdd(vector<float> a, vector<float> b) {
    vector<float> r(3, 0);
    r.at(0) = a.at(0) + b.at(0);
    r.at(1) = a.at(1) + b.at(1);
    r.at(2) = a.at(2) + b.at(2);
    return r;
}
// function for subtracting vectors
vector<float> vSubtract(vector<float> a, vector<float> b) {
    vector<float> r(3, 0);
    r.at(0) = a.at(0) - b.at(0);
    r.at(1) = a.at(1) - b.at(1);
    r.at(2) = a.at(2) - b.at(2);
    return r;
}
// function for normalizing vectors
vector<float> normalize(vector<float> a) {
    vector<float> r(3, 0);
    float len = sqrt(pow(a.at(0), 2.0) + pow(a.at(1), 2.0) + pow(a.at(2), 2.0));
    r.at(0) = a.at(0) / len;
    r.at(1) = a.at(1) / len;
    r.at(2) = a.at(2) / len;
    return r;
}

Color TraceRay(RayType ray, int id_check, float* pshadow_level, float dist, int depth_limit, int cmd);

Color Shade_Ray(Sphere* psphere, Triangle* ptriangle, vector<float> i_point, vector<float> ray_origin, vector<float> ray_dir, int depth_level, int cmd) {
    // returns the intersected sphere's/triangle's color
    // The Phong Illumination Model
    // initializing important variables
    vector<float> sum(3, 0);
    vector<float> Os(3, 0);
    float kd;
    float ks;
    float n;
    vector<float> dir_surface;
    vector<float> ic_half1;
    vector<float> direction_I;
    vector<float> direction_I_p;
    float eta;
    float opacity;
    float tu;
    float tv;
    vector<float> Od(3, 0);
    vector<vector<vector<float> > > * ptexture;
    float eta_i = bkgcolor.eta;

    // SPHERE
    if (psphere->radius != 0) {
        // texture
        if (psphere->texture.at(0).at(0).at(0) != -1) {
            float nx = (i_point.at(0) - psphere->center_x) / psphere->radius;
            float ny = (i_point.at(1) - psphere->center_y) / psphere->radius;
            float nz = (i_point.at(2) - psphere->center_z) / psphere->radius;
            // cout << "i: " << i_point.at(2) << " nz: " << nz << endl;
            float phi = acos(nz);
            float theta = atan2(ny, nx);
            tu = 0.5 + (theta / (2 * M_PI));
            tv = phi / M_PI;
            // use the coordinate to look up corressponding color and change Od
            int i = round(tu * (texture_width - 1));
            int j = round(tv * (texture_height - 1));
            j = texture_height - 1 - j;
            if (j >= texture_height - 2) {
                j = texture_height - 2;
            }
            if (j < 0) {
                j = 0;
            }
            if (i >= texture_width - 2) {
                i = texture_width - 2;
            }
            if (i < 0) {
                i = 0;
            }
            Od = psphere->texture.at(j).at(i);
        }
        else {
            Od.at(0) = psphere->m.dc.red;
            Od.at(1) = psphere->m.dc.green;
            Od.at(2) = psphere->m.dc.blue;
        }
        vector<float> sphere_c(3, 0);
        sphere_c.at(0) = psphere->center_x;
        sphere_c.at(1) = psphere->center_y;
        sphere_c.at(2) = psphere->center_z;
        dir_surface = vMultiply((1.0 / psphere->radius), vSubtract(i_point, sphere_c));
        
        Os.at(0) = psphere->m.sc.red;
        Os.at(1) = psphere->m.sc.green;
        Os.at(2) = psphere->m.sc.blue;

        ic_half1 = vMultiply(psphere->m.ka, Od);
        kd = psphere->m.kd;
        ks = psphere->m.ks;
        n = psphere->m.n;
        eta = psphere->m.m_eta;
        opacity = psphere->m.m_alpha;
    }

    // TRIANGLE
    else {
        // texture
        if (ptriangle->vt.size() != 0) {
            tu = ptriangle->alpha * ptriangle->vt.at(0).at(0) + ptriangle->beta * ptriangle->vt.at(1).at(0) + ptriangle->gamma * ptriangle->vt.at(2).at(0);
            tv = ptriangle->alpha * ptriangle->vt.at(0).at(1) + ptriangle->beta * ptriangle->vt.at(1).at(1) + ptriangle->gamma * ptriangle->vt.at(2).at(1);
            // use the coordinate to look up corressponding color and change Od
            int i = round(tu * (texture_width - 1));
            int j = round(tv * (texture_height - 1));
            j = texture_height - 1 - j;
            if (j >= texture_height - 2) {
                j = texture_height - 2;
            }
            if (j < 0) {
                j = 0;
            }
            if (i >= texture_width - 2) {
                i = texture_width - 2;
            }
            if (i < 0) {
                i = 0;
            }
            Od = psphere->texture.at(j).at(i);
        }
        else {
            Od.at(0) = ptriangle->m.dc.red;
            Od.at(1) = ptriangle->m.dc.green;
            Od.at(2) = ptriangle->m.dc.blue;
        }
        dir_surface = normalize(ptriangle->n);
        if (ptriangle->vn.size() != 0) {
            vector<float> i_n = vAdd(vAdd(vMultiply(ptriangle->alpha, ptriangle->vn.at(0)), vMultiply(ptriangle->beta, ptriangle->vn.at(1))), vMultiply(ptriangle->gamma, ptriangle->vn.at(2)));
            dir_surface = normalize(i_n);
        }
        Os.at(0) = ptriangle->m.sc.red;
        Os.at(1) = ptriangle->m.sc.green;
        Os.at(2) = ptriangle->m.sc.blue;
        ic_half1 = vMultiply(ptriangle->m.ka, Od);
        kd = ptriangle->m.kd;
        ks = ptriangle->m.ks;
        n = ptriangle->m.n;
        eta = ptriangle->m.m_eta;
        opacity = ptriangle->m.m_alpha;
    }

    if (depth_level < 1) {
        Color r_color;
        r_color.set_color(Od.at(0), Od.at(1), Od.at(2));
        return r_color;
    }

    //direction_I_p = normalize(vSubtract(ray_origin, i_point));
    direction_I_p = normalize(vSubtract(i_point, ray_origin));
    direction_I = vMultiply(-1, direction_I_p);

    if (cmd == INT_TRANS_OUT) {
        // exit from transparent object
        eta_i = eta;
        eta = bkgcolor.eta;
        dir_surface = vMultiply(-1, dir_surface);
    }

    Light light;
    float shadow_level;
    // attentuation function decleration
    float f_att = 0.0;
    float dist = 0.0;
    RayType shadow_ray(i_point.at(0), i_point.at(1), i_point.at(2));
    // making calculations for each light source
    for (int i = 0; i < (int)light_vector.size(); i++) {
        // Li
        shadow_level = 1;
        light = light_vector.at(i);
        vector<float> l(3, 0);
        dist = sqrt(pow((light.x - i_point.at(0)), 2.0) + pow((light.y - i_point.at(1)), 2.0) + pow((light.z - i_point.at(2)), 2.0));
        // cout << dist << endl;
        if (light.lsa == true) {
            f_att = 1.0 / (light.c1 + light.c2 * dist + light.c3 * dist * dist);
        }
        else {
            f_att = 1;
        }
        // cout << "lred: " << light.lc.red << " lgreen: " << light.lc.green << " lblue: " << light.lc.blue << endl;
        // directional light
        if (light.w == 0) {
            l.at(0) = -1.0 * light.x;
            l.at(1) = -1.0 * light.y;
            l.at(2) = -1.0 * light.z;
            l = normalize(l);
        }
        // point light
        else {
            l.at(0) = light.x;
            l.at(1) = light.y;
            l.at(2) = light.z;
            l = vSubtract(l, i_point);
            l = normalize(l);
        }
        // recursively call TraceRay if there is no shadow
        shadow_ray.dx = l.at(0);
        shadow_ray.dy = l.at(1);
        shadow_ray.dz = l.at(2);

        if (psphere->radius != 0)
            TraceRay(shadow_ray, psphere->id, &shadow_level, dist, depth_level, INT_SHADOW);
        else
            TraceRay(shadow_ray, ptriangle->id, &shadow_level, dist, depth_level, INT_SHADOW);

        // ILi
        vector<float> ilc(3, 0);
        ilc.at(0) = light.lc.red;
        ilc.at(1) = light.lc.green;
        ilc.at(2) = light.lc.blue;
        // Hi
        vector<float> h;
        if (cmd != INT_EYE) {
            vector<float> v = vMultiply(-1, ray_dir);
            h = normalize(vAdd(l, v));
        }
        else {
            h = normalize(vAdd(l, direction_I));
        }
        // I
        vector<float> ic_mid = vMultiply(max((float)0, dProduct(dir_surface, l)), vMultiply(kd, Od));
        vector<float> ic_mid2 = vMultiply(pow(max((float)0, dProduct(dir_surface, h)), n), vMultiply(ks, Os));
        vector<float> ic_half2 = vAdd(ic_mid, ic_mid2);
        vector<float> ic2 = vMultiply(f_att, vMul(ilc, ic_half2));
        sum = vAdd(sum, ic2);
        sum = vMultiply(shadow_level, sum);
    }
    // color values for the final summation can't be > 1
    // adding the specular term
    sum = vAdd(ic_half1, sum);

    // calculating direction of the reflected ray
    // calculating fresnel reflectance
    float f0;
    float cos_theta_i = dProduct(direction_I, dir_surface);
    float fresnel_ref;
    float sin_theta_i = sqrt(1 - pow(cos_theta_i, 2.0));
    
    f0 = pow(((eta - eta_i) / (eta + eta_i)), 2.0);
    fresnel_ref = f0 + (1.0 - f0) * pow((1.0 - cos_theta_i), 5.0);

    if (cmd != INT_TRANS_OUT)
    if (ks != 0) {
        // R = 2(N.I)N - I
        Color mirror_color;
        //               R                                   N.I            N            I
        vector<float> s_ref = vSubtract( vMultiply( (2 * cos_theta_i), dir_surface), direction_I);
        // recursively trace a ray in the s_ref direction
        // make new mirror ray object

        RayType mirror_ray(i_point.at(0), i_point.at(1), i_point.at(2));

        mirror_ray.dx = s_ref.at(0);
        mirror_ray.dy = s_ref.at(1);
        mirror_ray.dz = s_ref.at(2);

        if (psphere->radius != 0)
            mirror_color = TraceRay(mirror_ray, psphere->id, &shadow_level, dist, depth_level, INT_MIRROR);
        else
            mirror_color = TraceRay(mirror_ray, ptriangle->id, &shadow_level, dist, depth_level, INT_MIRROR);

        vector<float> mirror_color_vector(3, 0);
        mirror_color_vector.at(0) = mirror_color.red;
        mirror_color_vector.at(1) = mirror_color.green;
        mirror_color_vector.at(2) = mirror_color.blue;

        // add mirror part
        sum = vAdd(sum, vMultiply(fresnel_ref, mirror_color_vector));
    }

    // refraction
    // TIR
    //

    if ((opacity < 1) && (eta >= 1) && (sin_theta_i < (eta / eta_i)) && ks != 0) {
        Color t_color;
        float t_dir_p1 = sqrt(1 - (pow((eta_i / eta), 2.0) * (1 - pow(cos_theta_i, 2.0))));
        vector<float> t_dir_p2 = vMultiply((eta_i / eta), vSubtract(vMultiply(cos_theta_i, dir_surface), direction_I));
        vector<float> t_dir = vAdd(vMultiply(t_dir_p1, vMultiply(-1, dir_surface)), t_dir_p2);
        // calculate color returned by the recursively traced transmitted ray
        // this can be done by recursively trace a ray from the ray/surface intersection point through the surface in the refracted direction of transmission
        RayType t_ray(i_point.at(0), i_point.at(1), i_point.at(2));
        t_ray.dx = t_dir.at(0);
        t_ray.dy = t_dir.at(1);
        t_ray.dz = t_dir.at(2);
        //EYE->TRANS_IN : first transparent object
        //SHADOW->TRANS_IN : ?
        //MIRROR->TRANS_IN : prev object mirror to transparent
        //TRANS_IN->TRANS_IN : impossible
        //TRANS_OUT->TRANS_IN : go on
        if (cmd == INT_TRANS_OUT) {
            if (psphere->radius != 0)
                t_color = TraceRay(t_ray, psphere->id, &shadow_level, dist, depth_level, INT_EYE);
            else
                t_color = TraceRay(t_ray, ptriangle->id, &shadow_level, dist, depth_level, INT_EYE);
        }
        else {
            if (psphere->radius != 0)
                t_color = TraceRay(t_ray, psphere->id, &shadow_level, dist, depth_level, INT_TRANS_IN);
            else
                t_color = TraceRay(t_ray, ptriangle->id, &shadow_level, dist, depth_level, INT_TRANS_IN);
        }

        vector<float> t_color_vector(3, 0);
        t_color_vector.at(0) = t_color.red;
        t_color_vector.at(1) = t_color.green;
        t_color_vector.at(2) = t_color.blue;
        // transparency
        vector<float> t_part = vMultiply((1 - fresnel_ref) * (1 - opacity), t_color_vector);
        // add transparency part
        sum = vAdd(sum, t_part);
    }

    if (sum.at(0) > 1) {
        sum.at(0) = 1;
    }
    if (sum.at(1) > 1) {
        sum.at(1) = 1;
    }
    if (sum.at(2) > 1) {
        sum.at(2) = 1;
    }

    // optional depth cueing (1=true)
    if (depth == true) {
        // do calculations
        float a_dc;
        vector<float> i_dc(3, 0);
        i_dc.at(0) = depth_cue.at(0);
        i_dc.at(1) = depth_cue.at(1);
        i_dc.at(2) = depth_cue.at(2);
        float d_obj = sqrt(pow((eye_vector.at(0) - i_point.at(0)), 2.0) + pow((eye_vector.at(1) - i_point.at(1)), 2.0) + pow((eye_vector.at(2) - i_point.at(2)), 2.0));
        if (d_obj <= depth_cue.at(6)) {
            a_dc = depth_cue.at(3);
        }
        else if (d_obj >= depth_cue.at(6) && d_obj <= depth_cue.at(5)) {
            a_dc = (depth_cue.at(5) - d_obj) / (depth_cue.at(5) - depth_cue.at(6));
        }
        else if (d_obj >= depth_cue.at(5)) {
            a_dc = depth_cue.at(4);
        }
        sum = vAdd(vMultiply(a_dc, sum), vMultiply((1.0 - a_dc), i_dc));
    }

    Color r_color;
    r_color.set_color(sum.at(0), sum.at(1), sum.at(2));
    return r_color;
}

Color TraceRay(RayType ray, int id_check, float* pshadow_level, float dist, int depth_level, int cmd) {
    /* if a valid intersection is found, call Shade_Ray() to  determine the pixel
    color; otherwise, return the background color */
    // B^2 - 4C
    if (depth_level-- < 1) {
        return bkgcolor;
    }

    Sphere *p_sph;
    Triangle no_triangle;
    for (int i = 0; i < (int)sphere_vector.size(); i++) {
        p_sph = &sphere_vector.at(i);
        if (cmd != INT_TRANS_IN) {
            //inside sphere check own
            if (p_sph->id == id_check) {
                continue;
            }
        }
        else {
            //Chech same sphere
            if (p_sph->id != id_check) {
                continue;
            }
        }

        float b = 2 * (ray.dx * (ray.x - p_sph->center_x) +
            ray.dy * (ray.y - p_sph->center_y) +
            ray.dz * (ray.z - p_sph->center_z));

        float c = pow(ray.x - p_sph->center_x, 2.0) +
            pow(ray.y - p_sph->center_y, 2.0) +
            pow(ray.z - p_sph->center_z, 2.0) - pow(p_sph->radius, 2.0);

        // calculate the discriminant
        float discriminant = b * b - 4 * c;
        float t;
        if (cmd == INT_TRANS_IN) {
            t = (b * -1 + sqrt(discriminant)) / 2;
        }
        else {
            t = (b * -1 - sqrt(discriminant)) / 2;
        }
        vector<float> intersection_point(3, 0);
        vector<float> dir(3, 0);
        dir.at(0) = ray.dx;
        dir.at(1) = ray.dy;
        dir.at(2) = ray.dz;
        vector<float> pt(3, 0);
        pt.at(0) = ray.x;
        pt.at(1) = ray.y;
        pt.at(2) = ray.z;
        intersection_point = vAdd(pt, vMultiply(t, dir));

        // intersection and shadow check
        if (discriminant > 0) {
            if (cmd == INT_MIRROR && t>0) {
                return Shade_Ray(p_sph, &no_triangle, intersection_point, pt, dir, depth_level, cmd);
            }
            else if (cmd == INT_MIRROR && p_sph->texture.at(0).at(0).at(0) != -1) {
                return Shade_Ray(p_sph, &no_triangle, intersection_point, pt, dir, depth_level, cmd);
            }
            else if (cmd == INT_TRANS_IN) {
                // entering the sphere
                return Shade_Ray(p_sph, &no_triangle, intersection_point, pt, dir, depth_level, INT_TRANS_OUT);
            }
            else {
                // if shadow ray
                if (id_check > 0) {
                    if (t > 0 && t < dist) {
                        // s = s * (1-opacity)
                        *pshadow_level = *pshadow_level * (1 - p_sph->m.m_alpha);
                    }
                }
                // intersection of object that has no reflectance
                else {
                    return Shade_Ray(p_sph, &no_triangle, intersection_point, pt, dir, depth_level, cmd);
                }
            }
        }
    }
    // triangle intersection
    Triangle *p_tri;
    Sphere no_sphere;
    no_sphere.radius = 0;
    for (int i = 0; i < (int)tri_vector.size(); i++) {
        p_tri = &tri_vector.at(i);
        if (cmd != INT_TRANS_IN) {
            //inside sphere check own
            if (p_tri->id == id_check) {
                continue;
            }
        }
        else {
            //Chech same sphere
            if (p_tri->id != id_check) {
                continue;
            }
        }

        vector<float> p0(3, 0);
        p0.at(0) = p_tri->tri.at(0).at(0);
        p0.at(1) = p_tri->tri.at(0).at(1);
        p0.at(2) = p_tri->tri.at(0).at(2);
        vector<float> p1(3, 0);
        p1.at(0) = p_tri->tri.at(1).at(0);
        p1.at(1) = p_tri->tri.at(1).at(1);
        p1.at(2) = p_tri->tri.at(1).at(2);
        vector<float> p2(3, 0);
        p2.at(0) = p_tri->tri.at(2).at(0);
        p2.at(1) = p_tri->tri.at(2).at(1);
        p2.at(2) = p_tri->tri.at(2).at(2);
        vector<float> e1 = vSubtract(p1, p0);
        vector<float> e2 = vSubtract(p2, p0);
        vector<float> n = cProduct(e1, e2);
        p_tri->set_n(n);
        float d = -1 * (n.at(0) * p0.at(0) + n.at(1) * p0.at(1) + n.at(2) * p0.at(2));

        // calculate the denominator
        float denominator = n.at(0) * ray.dx + n.at(1) * ray.dy + n.at(2) * ray.dz;
        if (denominator == 0) {
            continue;
        }
        float t = -1 * ((n.at(0) * ray.x + n.at(1) * ray.y + n.at(2) * ray.z + d) / denominator);
        vector<float> intersection_point(3, 0);
        vector<float> dir(3, 0);
        dir.at(0) = ray.dx;
        dir.at(1) = ray.dy;
        dir.at(2) = ray.dz;
        vector<float> pt(3, 0);
        pt.at(0) = ray.x;
        pt.at(1) = ray.y;
        pt.at(2) = ray.z;
        intersection_point = vAdd(pt, vMultiply(t, dir));
        // calculating the barycentric coordinates
        vector<float> p = vSubtract(intersection_point, p0);
        float d11 = dProduct(e1, e1);
        float d22 = dProduct(e2, e2);
        float d12 = dProduct(e1, e2);
        float d1p = dProduct(e1, p);
        float d2p = dProduct(e2, p);
        float det = d11 * d22 - d12 * d12;
        if (det == 0) {
            continue;
        }
        float beta = (d22 * d1p - d12 * d2p) / det;
        float gamma = (d11 * d2p - d12 * d1p) / det;
        float alpha = 1 - (beta + gamma);
        // set triangle's alpha beta gamma
        p_tri->set_bcoordinates(alpha, beta, gamma);
        if (!(alpha <= 1 && alpha >= 0 && beta <= 1 && beta >= 0 && gamma <= 1 && gamma >= 0)) {
            continue;
        }

        if (cmd == INT_MIRROR && t>0) {
            return Shade_Ray(&no_sphere, p_tri, intersection_point, pt, dir, depth_level, cmd);
        }
        else if (cmd == INT_MIRROR && p_tri->vt.size() != 0) {
            return Shade_Ray(&no_sphere, p_tri, intersection_point, pt, dir, depth_level, cmd);
        }
        else if (cmd == INT_TRANS_IN) {
            //
            Shade_Ray(&no_sphere, p_tri, intersection_point, pt, dir, depth_level, INT_TRANS_OUT);
        }
        else {
            if (id_check > 0) {
                if (t > 0 && t < dist) {
                    *pshadow_level = *pshadow_level * (1 - p_tri->m.m_alpha);
                    *pshadow_level = 0;
                }
            }
            else {
                return Shade_Ray(&no_sphere, p_tri, intersection_point, pt, dir, depth_level, cmd);
            }
        }
    }
    return bkgcolor;
}

// function for parsing lines
vector<float> parseFloat(vector<float> arr, string rest, string delimeter) {
    string token;
    size_t pos = 0;
    while ((pos = rest.find(delimeter)) != string::npos) {
        token = rest.substr(0, pos);
        arr.push_back(atof(token.c_str()));
        rest.erase(0, pos + delimeter.length());
    }
    pos = rest.find(delimeter);
    token = rest.substr(0, pos);
    arr.push_back(atof(token.c_str()));
    // for(int i=0; i<(int)arr.size(); i++){cout << arr.at(i) << endl;}
    return arr;
}

vector<string> parseFace(vector<string> arr, string rest, string delimeter) {
    string token;
    size_t pos = 0;
    while ((pos = rest.find(delimeter)) != string::npos) {
        token = rest.substr(0, pos);
        arr.push_back(token.c_str());
        rest.erase(0, pos + delimeter.length());
    }
    pos = rest.find(delimeter);
    token = rest.substr(0, pos);
    arr.push_back(token.c_str());
    // for(int i=0; i<(int)arr.size(); i++){cout << arr.at(i) << endl;}
    return arr;
}

// MAIN
int main(int argc, const char* argv[]) {
    if (argc < 2) {
        cout << "No command line arguments" << endl;
        return -1;
    }

    // defining important variables
    ifstream file;
    file.open(argv[1]);
    string line;
    string keys[16] = { "imsize", "eye", "viewdir", "updir", "hfov", "bkgcolor", "mtlcolor", "sphere", "light", "attlight", "depthcueing", "v", "f", "vn", "vt", "texture" };
    vector<float> bkg;
    float hfov;
    float dist = 1;
    vector<bool> checklist(8, false);

    // file operations and parsing lines
    if (file.is_open()) {
        while (getline(file, line)) {
            cout << line << endl;
            if (line.at(0) == '#') continue;

            size_t start = line.find_first_of(" ");
            string key = line.substr(0, start);
            string rest = line.substr(start + 1, line.length());

            if (strcmp(key.c_str(), keys[0].c_str()) == 0) {
                // imsize
                string str2 = line.substr(start + 1, line.length());
                size_t end = str2.find_first_of(" ");
                width = atoi(str2.substr(0, end).c_str());
                height = atoi(str2.substr(end + 1, str2.length()).c_str());
                cout << "height: " << height << " width: " << width << endl;
                if (height > 0 && width > 0) {
                    checklist.at(0) = true;
                }
            }
            else if (strcmp(key.c_str(), keys[1].c_str()) == 0) {
                // eye
                eye_vector = parseFloat(eye_vector, rest, " ");
                // error check
                if (eye_vector.size() != 3) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                checklist.at(1) = true;
            }
            else if (strcmp(key.c_str(), keys[2].c_str()) == 0) {
                // viewdir
                viewdir = parseFloat(viewdir, rest, " ");
                // error check
                if (viewdir.size() != 3) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                checklist.at(2) = true;
            }
            else if (strcmp(key.c_str(), keys[3].c_str()) == 0) {
                // updir
                updir = parseFloat(updir, rest, " ");
                // error check
                if (updir.size() != 3) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                checklist.at(3) = true;
            }
            else if (strcmp(key.c_str(), keys[4].c_str()) == 0) {
                // hfov
                string str2 = line.substr(start + 1, line.length());
                hfov = atoi(str2.substr(0, str2.length()).c_str());
                checklist.at(4) = true;
            }
            else if (strcmp(key.c_str(), keys[5].c_str()) == 0) {
                // bkgcolor
                bkg = parseFloat(bkg, rest, " ");
                // error check
                if (bkg.size() != 4) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                bkgcolor.set_color(bkg.at(0), bkg.at(1), bkg.at(2));
                bkgcolor.set_eta(bkg.at(3));
                checklist.at(5) = true;
            }
            else if (strcmp(key.c_str(), keys[10].c_str()) == 0) {
                // depthcueing
                depth_cue = parseFloat(depth_cue, rest, " ");
                // error check
                if (depth_cue.size() != 7) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                depth = true;
            }
            else if (strcmp(key.c_str(), keys[6].c_str()) == 0) {
                // mtlcolor
                vector<float> mtl;
                mtl = parseFloat(mtl, rest, " ");
                // error check
                if (mtl.size() != 12) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                Color dc;
                Color sc;
                dc.set_color(mtl.at(0), mtl.at(1), mtl.at(2));
                sc.set_color(mtl.at(3), mtl.at(4), mtl.at(5));
                mtlcolor.set_mtl(dc, sc, mtl.at(6), mtl.at(7), mtl.at(8), mtl.at(9), mtl.at(10), mtl.at(11));
                checklist.at(6) = true;
            }
            else if (strcmp(key.c_str(), keys[15].c_str()) == 0) {
                // texture file
                ifstream t_file;
                rest.erase(rest.size() - 1);
                // cout << rest.c_str() << endl;
                t_file.open(rest.c_str());
                string tline;
                vector<vector<float> > arr3;
                int count = 0;
                if (t_file.is_open()) {
                    while (getline(t_file, tline)) {
                        count++;
                        vector<string> arr;
                        arr = parseFace(arr, tline, "  ");
                        if (count == 1) {
                            vector<float> arr2;
                            arr2 = parseFloat(arr2, arr.at(0), " ");
                            texture_width = arr2.at(1);
                            texture_height = arr2.at(2);
                            cout << "height: " << texture_height << endl;
                            cout << "width: " << texture_width << endl;
                            continue;
                        }
                        for (int i = 0; i < (int)arr.size(); i++) {
                            vector<float> arr2;
                            if (arr.at(i) == "")
                                continue;

                            arr2 = parseFloat(arr2, arr.at(i), " ");
                            if (arr2.size() == 3) {
                                arr2.at(0) = arr2.at(0) / 255.0;
                                arr2.at(1) = arr2.at(1) / 255.0;
                                arr2.at(2) = arr2.at(2) / 255.0;
                                arr3.push_back(arr2);
                            }
                            else {
                                cout << "Texture Err (don't mind this, empty spaces end of every line): " << count << ", i: " << i << endl;
                            }
                        }
                        texture.push_back(arr3);
                        arr3.clear();
                    }
                    cout << "Texture Ok: " << count << endl;
                }
                else {
                    cout << "File isn't open" << endl;
                    return -1;
                }
                t_file.close();
            }
            else if (strcmp(key.c_str(), keys[7].c_str()) == 0) {
                // sphere
                vector<float> sph;
                sph = parseFloat(sph, rest, " ");
                // error check
                if (sph.size() != 4) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                Sphere sphere;
                if (texture.size() != 0) {
                    sphere.set_center(sph.at(0), sph.at(1), sph.at(2), sph.at(3), mtlcolor, texture, ++sphere_id);
                }
                else {
                    vector<vector<vector<float> > > no_texture(3, vector<vector<float> >(3, vector<float>(3, -1)));
                    sphere.set_center(sph.at(0), sph.at(1), sph.at(2), sph.at(3), mtlcolor, no_texture, ++sphere_id);
                }
                sphere_vector.push_back(sphere); // in case we have multiple entries
                checklist.at(7) = true;
            }
            else if (strcmp(key.c_str(), keys[8].c_str()) == 0 || strcmp(key.c_str(), keys[9].c_str()) == 0) {
                // light or attlight
                vector<float> lt;
                lt = parseFloat(lt, rest, " ");
                Color lc;
                // error checks
                if (strcmp(key.c_str(), keys[8].c_str()) == 0) {
                    if (lt.size() != 7) {
                        cout << "Incorrect input" << endl;
                        return -1;
                    }
                    lc.set_color(lt.at(4), lt.at(5), lt.at(6));
                    Light light(lt.at(0), lt.at(1), lt.at(2), lt.at(3), lc);
                    light_vector.push_back(light); // in case we have multiple entries
                }
                else {
                    if (lt.size() != 10) {
                        cout << "Incorrect input" << endl;
                        return -1;
                    }
                    lc.set_color(lt.at(4), lt.at(5), lt.at(6));
                    Light light(lt.at(0), lt.at(1), lt.at(2), lt.at(3), lc, lt.at(7), lt.at(8), lt.at(9));
                    light_vector.push_back(light); // in case we have multiple entries
                }
                if (lt.at(3) != 0.0 && lt.at(3) != 1.0) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                checklist.at(7) = true;
            }
            else if (strcmp(key.c_str(), keys[11].c_str()) == 0) {
                // triangle vertex positions
                vector<float> vertex;
                vertex = parseFloat(vertex, rest, " ");
                // error check
                if (vertex.size() != 3) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                vertex_vector.push_back(vertex); // in case we have multiple entries
            }
            else if (strcmp(key.c_str(), keys[13].c_str()) == 0) {
                // vertex normal vectors
                vector<float> vn;
                vn = parseFloat(vn, rest, " ");
                // error check
                if (vn.size() != 3) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                vn_vector.push_back(vn); // in case we have multiple entries
            }
            else if (strcmp(key.c_str(), keys[14].c_str()) == 0) {
                // vertex normal vectors
                vector<float> vt;
                vt = parseFloat(vt, rest, " ");
                // error check
                if (vt.size() != 2) {
                    cout << "Incorrect input" << endl;
                    return -1;
                }
                vt_vector.push_back(vt); // in case we have multiple entries
            }
            else if (strcmp(key.c_str(), keys[12].c_str()) == 0) {
                // triangle faces
                vector<float> face_v;
                vector<float> face_vn;
                vector<float> face_vt;
                vector<string> test;
                vector<string> test2;
                vector<string> test3;
                vector<string> test4;
                test = parseFace(test, rest, " ");
                size_t found = test.at(0).find("//");
                size_t found2 = test.at(0).find("/");
                // smooth shaded untextured triangle
                if (found != string::npos) {
                    test2 = parseFace(test2, test.at(0), "//");
                    face_v.push_back(atof(test2.at(0).c_str()));
                    face_vn.push_back(atof(test2.at(1).c_str()));
                    test3 = parseFace(test3, test.at(1), "//");
                    face_v.push_back(atof(test3.at(0).c_str()));
                    face_vn.push_back(atof(test3.at(1).c_str()));
                    test4 = parseFace(test4, test.at(2), "//");
                    face_v.push_back(atof(test4.at(0).c_str()));
                    face_vn.push_back(atof(test4.at(1).c_str()));
                }
                else {
                    if (found2 != string::npos) {
                        test2 = parseFace(test2, test.at(0), "/");
                        // smooth shaded textured triangle
                        if (test2.size() == 3) {
                            face_v.push_back(atof(test2.at(0).c_str()));
                            face_vn.push_back(atof(test2.at(1).c_str()));
                            face_vt.push_back(atof(test2.at(2).c_str()));
                            test3 = parseFace(test3, test.at(1), "/");
                            face_v.push_back(atof(test3.at(0).c_str()));
                            face_vn.push_back(atof(test3.at(1).c_str()));
                            face_vt.push_back(atof(test3.at(2).c_str()));
                            test4 = parseFace(test4, test.at(2), "/");
                            face_v.push_back(atof(test4.at(0).c_str()));
                            face_vn.push_back(atof(test4.at(1).c_str()));
                            face_vt.push_back(atof(test4.at(2).c_str()));
                        }
                        // without smooth shaded textured triangle
                        else if (test2.size() == 2) {
                            face_v.push_back(atof(test2.at(0).c_str()));
                            face_vt.push_back(atof(test2.at(1).c_str()));
                            test3 = parseFace(test3, test.at(1), "/");
                            face_v.push_back(atof(test3.at(0).c_str()));
                            face_vt.push_back(atof(test3.at(1).c_str()));
                            test4 = parseFace(test4, test.at(2), "/");
                            face_v.push_back(atof(test4.at(0).c_str()));
                            face_vt.push_back(atof(test4.at(1).c_str()));
                        }
                    }
                    // without smooth shaded and no texture triangle
                    else {
                        face_v.push_back(atof(test.at(0).c_str()));
                        face_v.push_back(atof(test.at(1).c_str()));
                        face_v.push_back(atof(test.at(2).c_str()));
                    }
                }
                // in case we have multiple entries
                face_v_vector.push_back(face_v);
                face_vt_vector.push_back(face_vt);
                face_vn_vector.push_back(face_vn);
                Triangle triangle;
                if (texture.size() != 0) {
                    triangle.set_face(face_v, vertex_vector, mtlcolor, texture);
                }
                else {
                    vector<vector<vector<float> > > no_texture(3, vector<vector<float> >(3, vector<float>(3, -1)));
                    triangle.set_face(face_v, vertex_vector, mtlcolor, no_texture);
                }
                tri_vector.push_back(triangle);
            }
            else if (strcmp(key.c_str(), "") == 0) {
                continue;
            }
            else {
                cout << "Incorrect input" << endl;
                return -1;
            }
        }
        // additional error checks regarding the input file
        for (int i = 0; i < (int)checklist.size(); i++) {
            if (checklist.at(i) == false) {
                cout << "Incorrect input" << endl;
                return -1;
            }
        }
        file.close();
    }
    else {
        cout << "File is not open" << endl;
        return -1;
    }

    /* initialize pixel array for output image */
    vector<Color> pixel_arr(width * height, Color());

    // initialize Triangle class using the face and vertex vectors
    ++triangle_id;
    for (int i = 0; i < (int)tri_vector.size(); i++) {
        tri_vector.at(i).set_id(triangle_id);
        if (face_vn_vector.size() != 0) {
            tri_vector.at(i).set_vn(face_vn_vector[i], vn_vector);
        }
        if (face_vt_vector.size() != 0) {
            tri_vector.at(i).set_vt(face_vt_vector[i], vt_vector);
        }
    }

    /* perform any required preliminary calculations */
    vector<float> u_prime = cProduct(viewdir, updir);
    vector<float> u = normalize(u_prime);
    vector<float> v_prime = cProduct(u, viewdir);
    vector<float> v = normalize(v_prime);

    float w = 2 * dist * tan(hfov * M_PI / 180 * 1 / 2);
    float aspect_ratio = (float)width / height;
    float h = w / aspect_ratio;

    vector<float> dn = vMultiply(dist, viewdir);
    vector<float> wu = vMultiply(w / 2, u);
    vector<float> hv = vMultiply(h / 2, v);

    vector<float> ul = vAdd(vSubtract(vAdd(eye_vector, dn), wu), hv);
    vector<float> ur = vAdd(vAdd(vAdd(eye_vector, dn), wu), hv);
    vector<float> ll = vSubtract(vSubtract(vAdd(eye_vector, dn), wu), hv);
    vector<float> lr = vSubtract(vAdd(vAdd(eye_vector, dn), wu), hv);

    vector<float> delta_h = vMultiply(1.0 / (width - 1), vSubtract(ur, ul));
    vector<float> delta_v = vMultiply(1.0 / (height - 1), vSubtract(ll, ul));

    RayType ray(eye_vector.at(0), eye_vector.at(1), eye_vector.at(2));
    vector<float> ray_d(3, 0);
    vector<float> ray_p(3, 0);
    float shadow_level;

    /* for each pixel in the image array: */
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // ray position
            ray_p = vAdd(vAdd(ul, vMultiply(x, delta_h)), vMultiply(y, delta_v));
            // ray direction
            ray_d.at(0) = ray_p.at(0) - eye_vector.at(0);
            ray_d.at(1) = ray_p.at(1) - eye_vector.at(1);
            ray_d.at(2) = ray_p.at(2) - eye_vector.at(2);
            ray_d = normalize(ray_d);
            // ray variable to pass in for the TraceRay function
            ray.set_dir(ray_d.at(0), ray_d.at(1), ray_d.at(2));
            /* call Trace_Ray() with appropriate parameters */
            /* use the value returned by Trace_Ray() to define the pixel color */
            shadow_level = 1;
            // cout << "before traceray" << endl;
            pixel_arr.at(x + y * width) = TraceRay(ray, 0, &shadow_level, 1.0, 10, INT_EYE);
            // cout << "x: " << x << endl;
        }
        cout << "y: " << y << endl;
    }
    // put this in a function instead of doing everything in main
    /* write the final image to an output file */
    string output = argv[1];
    size_t dot = output.find_last_of(".");
    output = output.substr(0, dot);
    output += ".ppm";
    cout << "name of file: " << output << endl;
    ofstream fout(output.c_str());
    if (fout.fail()) {
        cout << "Couldn't create output file" << endl;
        return -1;
    }
    // create headers
    fout << "P3\n";
    fout << width << " " << height << "\n";
    fout << "255\n";
    // write from pixel array for the body of the ppm file
    for (int i = 0; i < width * height; ++i) {
        fout << pixel_arr.at(i).red * 255 << " ";
        fout << pixel_arr.at(i).green * 255 << " ";
        fout << pixel_arr.at(i).blue * 255 << "\n";
    }
    fout.close();
    return 0; //yayy
}
