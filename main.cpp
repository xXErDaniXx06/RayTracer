#include <iostream>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm>
#include <vector>
#include <memory>
#include <omp.h> // <--- 1. IMPORTANTE: Librería OpenMP

// --- CONSTANTES ---
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// --- UTILIDADES ---
inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

inline double random_double(double min, double max) {
    return min + (max-min)*random_double();
}

// --- VECTORES ---
struct vec3 {
    double x, y, z;
    vec3 operator+(const vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    vec3 operator-(const vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    vec3 operator-() const { return {-x, -y, -z}; }
    vec3 operator*(double t) const { return {x * t, y * t, z * t}; }
    vec3 operator*(const vec3& v) const { return {x * v.x, y * v.y, z * v.z}; }
    vec3 operator/(double t) const { return {x / t, y / t, z / t}; }
    double dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    
    vec3 cross(const vec3& v) const {
        return {y * v.z - z * v.y,
                z * v.x - x * v.z,
                x * v.y - y * v.x};
    }

    double length() const { return std::sqrt(x*x + y*y + z*z); }
    double length_squared() const { return x*x + y*y + z*z; }
    vec3 normalize() const { return *this / length(); }
};

vec3 operator*(double t, const vec3& v) { return {t*v.x, t*v.y, t*v.z}; }

vec3 random_in_unit_sphere() {
    while (true) {
        vec3 p = {random_double(-1,1), random_double(-1,1), random_double(-1,1)};
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

vec3 random_in_unit_disk() {
    while (true) {
        vec3 p = {random_double(-1,1), random_double(-1,1), 0};
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2 * v.dot(n) * n;
}

vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = std::fmin((-uv).dot(n), 1.0);
    vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec3 r_out_parallel = -std::sqrt(std::fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

double reflectance(double cosine, double ref_idx) {
    auto r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*std::pow((1 - cosine),5);
}

// --- MATERIALES ---
enum MaterialType { LAMBERTIAN, METAL, DIELECTRIC };

struct Material {
    MaterialType type;
    vec3 albedo;      
    double fuzz;      
    double ref_idx;   
};

struct HitRecord {
    vec3 p;
    vec3 normal;
    double t;
    bool front_face;
    Material mat;

    void set_face_normal(const vec3& r_direction, const vec3& outward_normal) {
        front_face = r_direction.dot(outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

struct Sphere {
    vec3 center;
    double radius;
    Material mat;
};

// --- CÁMARA ---
struct Camera {
    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 u, v, w;
    double lens_radius;

    Camera(vec3 lookfrom, vec3 lookat, vec3 vup, double vfov, double aspect_ratio, double aperture, double focus_dist) {
        auto theta = vfov * pi / 180.0;
        auto h = tan(theta/2);
        auto viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;

        w = (lookfrom - lookat).normalize();
        u = vup.cross(w).normalize();
        v = w.cross(u);

        origin = lookfrom;
        horizontal = focus_dist * viewport_width * u;
        vertical = focus_dist * viewport_height * v;
        lower_left_corner = origin - horizontal/2 - vertical/2 - focus_dist*w;

        lens_radius = aperture / 2;
    }

    vec3 get_ray_direction(double s, double t) const {
        vec3 rd = lens_radius * random_in_unit_disk();
        vec3 offset = u * rd.x + v * rd.y;
        return lower_left_corner + s*horizontal + t*vertical - origin - offset;
    }
    
    vec3 get_ray_origin() const {
        vec3 rd = lens_radius * random_in_unit_disk();
        vec3 offset = u * rd.x + v * rd.y;
        return origin + offset;
    }
};

// --- FÍSICA ---
bool hit_sphere_obj(const Sphere& s, const vec3& r_origin, const vec3& r_dir, double t_min, double t_max, HitRecord& rec) {
    vec3 oc = r_origin - s.center;
    double a = r_dir.length_squared();
    double half_b = oc.dot(r_dir);
    double c = oc.length_squared() - s.radius*s.radius;
    double discriminant = half_b*half_b - a*c;

    if (discriminant < 0) return false;
    double sqrtd = std::sqrt(discriminant);

    double root = (-half_b - sqrtd) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || root > t_max)
            return false;
    }

    rec.t = root;
    rec.p = r_origin + r_dir * rec.t;
    vec3 outward_normal = (rec.p - s.center) / s.radius;
    rec.set_face_normal(r_dir, outward_normal);
    rec.mat = s.mat;
    return true;
}

bool hit_world(const std::vector<Sphere>& world, const vec3& r_origin, const vec3& r_dir, double t_min, double t_max, HitRecord& rec) {
    HitRecord temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;

    for (const auto& object : world) {
        if (hit_sphere_obj(object, r_origin, r_dir, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    return hit_anything;
}

vec3 ray_color(const vec3& r_origin, const vec3& r_dir, const std::vector<Sphere>& world, int depth) {
    HitRecord rec;
    if (depth <= 0) return {0,0,0};

    if (hit_world(world, r_origin, r_dir, 0.001, infinity, rec)) {
        vec3 scattered_origin, scattered_dir;
        vec3 attenuation;
        bool scatter = false;

        if (rec.mat.type == LAMBERTIAN) {
            vec3 scatter_direction = rec.normal + random_in_unit_sphere().normalize();
            if (scatter_direction.length_squared() < 1e-8) scatter_direction = rec.normal;
            scattered_origin = rec.p;
            scattered_dir = scatter_direction;
            attenuation = rec.mat.albedo;
            scatter = true;
        } else if (rec.mat.type == METAL) {
            vec3 reflected = reflect(r_dir.normalize(), rec.normal);
            scattered_origin = rec.p;
            scattered_dir = reflected + rec.mat.fuzz * random_in_unit_sphere();
            attenuation = rec.mat.albedo;
            scatter = (scattered_dir.dot(rec.normal) > 0);
        } else if (rec.mat.type == DIELECTRIC) {
            attenuation = {1.0, 1.0, 1.0};
            double refraction_ratio = rec.front_face ? (1.0/rec.mat.ref_idx) : rec.mat.ref_idx;
            vec3 unit_direction = r_dir.normalize();
            double cos_theta = std::fmin((-unit_direction).dot(rec.normal), 1.0);
            double sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);
            bool cannot_refract = refraction_ratio * sin_theta > 1.0;
            vec3 direction;
            if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
                direction = reflect(unit_direction, rec.normal);
            else
                direction = refract(unit_direction, rec.normal, refraction_ratio);
            scattered_origin = rec.p;
            scattered_dir = direction;
            scatter = true;
        }

        if (scatter) {
            return attenuation * ray_color(scattered_origin, scattered_dir, world, depth-1);
        }
        return {0,0,0};
    }

    vec3 unit_direction = r_dir.normalize();
    double t = 0.5*(unit_direction.y + 1.0);
    return (1.0-t)*vec3{1.0, 1.0, 1.0} + t*vec3{0.5, 0.7, 1.0};
}

// Función auxiliar para imprimir el color final
void write_color_buffer(std::ostream &out, const std::vector<vec3>& buffer, int width, int height, int samples_per_pixel) {
    for (const auto& pixel_color : buffer) {
        auto r = pixel_color.x;
        auto g = pixel_color.y;
        auto b = pixel_color.z;
        double scale = 1.0 / samples_per_pixel;
        r = std::sqrt(scale * r);
        g = std::sqrt(scale * g);
        b = std::sqrt(scale * b);
        out << static_cast<int>(256 * std::clamp(r, 0.0, 0.999)) << ' '
            << static_cast<int>(256 * std::clamp(g, 0.0, 0.999)) << ' '
            << static_cast<int>(256 * std::clamp(b, 0.0, 0.999)) << '\n';
    }
}

// --- ESCENA ---
std::vector<Sphere> random_scene() {
    std::vector<Sphere> world;
    world.push_back({ {0,-1000,0}, 1000, {LAMBERTIAN, {0.5, 0.5, 0.5}, 0, 0} });

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            double choose_mat = random_double();
            vec3 center = {a + 0.9*random_double(), 0.2, b + 0.9*random_double()};
            if ((center - vec3{4, 0.2, 0}).length() > 0.9) {
                if (choose_mat < 0.8) {
                    vec3 albedo = vec3{random_double(), random_double(), random_double()} * vec3{random_double(), random_double(), random_double()};
                    world.push_back({center, 0.2, {LAMBERTIAN, albedo, 0, 0}});
                } else if (choose_mat < 0.95) {
                    vec3 albedo = {random_double(0.5, 1), random_double(0.5, 1), random_double(0.5, 1)};
                    double fuzz = random_double(0, 0.5);
                    world.push_back({center, 0.2, {METAL, albedo, fuzz, 0}});
                } else {
                    world.push_back({center, 0.2, {DIELECTRIC, {1,1,1}, 0, 1.5}});
                }
            }
        }
    }
    world.push_back({ {0, 1, 0}, 1.0, {DIELECTRIC, {1,1,1}, 0, 1.5} });
    world.push_back({ {-4, 1, 0}, 1.0, {LAMBERTIAN, {0.4, 0.2, 0.1}, 0, 0} });
    world.push_back({ {4, 1, 0}, 1.0, {METAL, {0.7, 0.6, 0.5}, 0.0, 0} });
    return world;
}

// --- MAIN ---
int main() {
    // 2. CONFIGURACIÓN
    omp_set_num_threads(12);
    const int image_width = 1200; 
    const int image_height = 800;
    const int samples_per_pixel = 500; // Alta calidad
    const int max_depth = 50;

    auto world = random_scene();

    vec3 lookfrom = {13, 2, 3};
    vec3 lookat = {0, 0, 0};
    vec3 vup = {0, 1, 0};
    double dist_to_focus = 10.0;
    double aperture = 0.1;

    Camera cam(lookfrom, lookat, vup, 20, double(image_width)/image_height, aperture, dist_to_focus);

    // 3. BUFFER DE MEMORIA (Aquí guardamos la imagen antes de imprimir)
    // El tamaño es ancho * alto píxeles
    std::vector<vec3> image_buffer(image_width * image_height);

    std::clog << "Renderizando en paralelo... (Esto usará 100% CPU)\n";

    // 4. PARALELIZACIÓN REAL
    // Quitamos la impresión de dentro del bucle
    #pragma omp parallel for schedule(dynamic, 50)
    for (int j = image_height-1; j >= 0; --j) {
        // No imprimimos 'Lineas restantes' porque explotaría la consola
        
        for (int i = 0; i < image_width; ++i) {
            vec3 pixel_color = {0, 0, 0};
            for (int s = 0; s < samples_per_pixel; ++s) {
                double u = (i + random_double()) / (image_width-1);
                double v = (j + random_double()) / (image_height-1);
                
                vec3 ray_o = cam.get_ray_origin();
                vec3 ray_d = cam.get_ray_direction(u, v);
                
                pixel_color = pixel_color + ray_color(ray_o, ray_d, world, max_depth);
            }
            
            // GUARDAMOS EN MEMORIA EN LUGAR DE IMPRIMIR
            // Calculamos el índice lineal del vector
            // La imagen se guarda de arriba a abajo en el archivo, así que invertimos J para el índice
            int index = (image_height - 1 - j) * image_width + i;
            image_buffer[index] = pixel_color;
        }
    }

    std::clog << "Renderizado terminado. Guardando archivo...\n";

    // 5. IMPRIMIR TODO DE GOLPE (Serializado)
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";
    write_color_buffer(std::cout, image_buffer, image_width, image_height, samples_per_pixel);

    std::clog << "Hecho.\n";
}