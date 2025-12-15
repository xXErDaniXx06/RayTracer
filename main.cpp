#include <iostream>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm> // ¡Importante para que no te de error!

// --- CONSTANTES Y UTILIDADES ---
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

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
    vec3 operator*(double t) const { return {x * t, y * t, z * t}; }
    vec3 operator/(double t) const { return {x / t, y / t, z / t}; }
    double dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    double length() const { return std::sqrt(x*x + y*y + z*z); }
    
    vec3 normalize() const { 
        double len = length();
        return {x/len, y/len, z/len}; 
    }
};

vec3 operator*(double t, const vec3& v) { return {t*v.x, t*v.y, t*v.z}; }

// Función para rebotar la luz como un espejo (METAL)
vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2 * v.dot(n) * n;
}

// Función para rebotar la luz aleatoriamente (MATE)
vec3 random_in_unit_sphere() {
    while (true) {
        vec3 p = {random_double(-1,1), random_double(-1,1), random_double(-1,1)};
        if (p.dot(p) >= 1) continue;
        return p;
    }
}

// --- COLOR ---
void write_color(std::ostream &out, vec3 pixel_color, int samples_per_pixel) {
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

// --- ESCENA ---
// Calcula si golpeamos una esfera
double hit_sphere(const vec3& center, double radius, const vec3& origin, const vec3& direction) {
    vec3 oc = origin - center;
    double a = direction.dot(direction);
    double b = 2.0 * oc.dot(direction);
    double c = oc.dot(oc) - radius*radius;
    double discriminant = b*b - 4*a*c;
    if (discriminant < 0) return -1.0;
    return (-b - std::sqrt(discriminant)) / (2.0*a);
}

// --- MOTOR DE LUZ ---
vec3 ray_color(const vec3& origin, const vec3& direction, int depth) {
    if (depth <= 0) return {0,0,0};

    // Definimos 2 esferas en el mundo
    vec3 esfera_mate_pos = {-0.6, 0, -1};  // Izquierda
    vec3 esfera_metal_pos = {0.6, 0, -1};  // Derecha
    double radio = 0.5;

    // Comprobamos golpes
    double t_mate = hit_sphere(esfera_mate_pos, radio, origin, direction);
    double t_metal = hit_sphere(esfera_metal_pos, radio, origin, direction);

    // Lógica simple para saber cuál está más cerca (Z-Buffer manual)
    double t_final = infinity;
    int tipo_material = 0; // 0=nada, 1=mate, 2=metal
    vec3 pos_esfera_final;

    if (t_mate > 0.001 && t_mate < t_final) {
        t_final = t_mate;
        tipo_material = 1;
        pos_esfera_final = esfera_mate_pos;
    }
    if (t_metal > 0.001 && t_metal < t_final) {
        t_final = t_metal;
        tipo_material = 2;
        pos_esfera_final = esfera_metal_pos;
    }

    // SI GOLPEAMOS ALGO
    if (t_final < infinity) {
        vec3 hit_point = origin + (direction * t_final);
        vec3 N = (hit_point - pos_esfera_final).normalize();
        vec3 target, attenuation;

        if (tipo_material == 1) { 
            // --- MATERIAL MATE (Azul/Gris) ---
            target = hit_point + N + random_in_unit_sphere();
            attenuation = {0.1, 0.2, 0.5}; // Color de la bola
            return {attenuation.x * ray_color(hit_point, target - hit_point, depth-1).x,
                    attenuation.y * ray_color(hit_point, target - hit_point, depth-1).y,
                    attenuation.z * ray_color(hit_point, target - hit_point, depth-1).z};
        } 
        else if (tipo_material == 2) {
            // --- MATERIAL METAL (Dorado) ---
            vec3 reflected = reflect(direction.normalize(), N);
            attenuation = {0.8, 0.6, 0.2}; // Color Oro
            // Añadimos un poco de "fuzz" (borrosidad) multiplicando random por 0.1
            vec3 scattered = reflected + (random_in_unit_sphere() * 0.1); 
            
            if (scattered.dot(N) > 0) {
                 vec3 col = ray_color(hit_point, scattered, depth-1);
                 return {attenuation.x * col.x, attenuation.y * col.y, attenuation.z * col.z};
            }
            return {0,0,0};
        }
    }

    // Fondo (Cielo)
    vec3 unit_direction = direction.normalize();
    double t = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - t)*vec3{1.0, 1.0, 1.0} + t*vec3{0.5, 0.7, 1.0};
}

// --- MAIN ---
int main() {
    const int image_width = 400;
    const int image_height = 225;
    const int samples_per_pixel = 50; 
    const int max_depth = 50;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    vec3 origin = {0, 0, 0};
    vec3 horizontal = {4.0, 0, 0};
    vec3 vertical = {0, 2.25, 0};
    vec3 lower_left_corner = origin - horizontal/2 - vertical/2 - vec3{0, 0, 1};

    for (int j = image_height-1; j >= 0; --j) {
        std::clog << "\rLineas restantes: " << j << " " << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 pixel_color = {0, 0, 0};
            for (int s = 0; s < samples_per_pixel; ++s) {
                double u = (i + random_double()) / (image_width-1);
                double v = (j + random_double()) / (image_height-1);
                vec3 r = lower_left_corner + u*horizontal + v*vertical - origin;
                pixel_color = pixel_color + ray_color(origin, r, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }
    std::clog << "\rHecho.\n";
}