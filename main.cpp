#include <iostream>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm>

// --- UTILIDADES ---
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Generador de números aleatorios
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
    double length_squared() const { return x*x + y*y + z*z; }
    
    vec3 normalize() const { 
        double len = std::sqrt(length_squared());
        return {x/len, y/len, z/len}; 
    }
};

// Funciones extra para vectores
vec3 operator*(double t, const vec3& v) { return {t*v.x, t*v.y, t*v.z}; }

// Generar un vector aleatorio dentro de una esfera unitaria (para el rebote mate)
vec3 random_in_unit_sphere() {
    while (true) {
        vec3 p = {random_double(-1,1), random_double(-1,1), random_double(-1,1)};
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

// --- COLOR ---
void write_color(std::ostream &out, vec3 pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x;
    auto g = pixel_color.y;
    auto b = pixel_color.z;

    // Dividimos el color por el número de muestras (promedio) y aplicamos Gamma 2.0 (raíz cuadrada)
    // Esto hace que la imagen se vea más brillante y realista
    double scale = 1.0 / samples_per_pixel;
    r = std::sqrt(scale * r);
    g = std::sqrt(scale * g);
    b = std::sqrt(scale * b);

    out << static_cast<int>(256 * std::clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * std::clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * std::clamp(b, 0.0, 0.999)) << '\n';
}

// --- FÍSICA ---
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
// Ahora la función es recursiva: depth controla cuántas veces rebota el rayo
vec3 ray_color(const vec3& origin, const vec3& direction, int depth) {
    // Si el rayo ha rebotado muchas veces, deja de calcular luz (se vuelve negro)
    if (depth <= 0) return {0,0,0};

    vec3 sphere_center = {0, 0, -1};
    double t = hit_sphere(sphere_center, 0.5, origin, direction);
    
    if (t > 0.001) { // 0.001 para arreglar el "acné de sombra"
        vec3 hit_point = origin + (direction * t);
        vec3 N = (hit_point - sphere_center).normalize();
        
        // MATERIAL MATE: El rayo rebota en una dirección aleatoria
        vec3 target = hit_point + N + random_in_unit_sphere();
        
        // Lanzamos un nuevo rayo desde el punto de choque hacia el objetivo aleatorio
        // Y absorbemos el 50% de la luz (0.5) en cada rebote
        return 0.5 * ray_color(hit_point, target - hit_point, depth - 1);
    }
    
    // Fondo (Cielo)
    vec3 unit_direction = direction.normalize();
    t = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - t)*vec3{1.0, 1.0, 1.0} + t*vec3{0.5, 0.7, 1.0};
}

int main() {
    const int image_width = 400;
    const int image_height = 225;
    const int samples_per_pixel = 50; // ¡Antialiasing! Lanzamos 50 rayos por cada píxel
    const int max_depth = 50; // Máximo 50 rebotes de luz

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    vec3 origin = {0, 0, 0};
    vec3 horizontal = {4.0, 0, 0};
    vec3 vertical = {0, 2.25, 0};
    vec3 lower_left_corner = origin - horizontal/2 - vertical/2 - vec3{0, 0, 1};

    for (int j = image_height-1; j >= 0; --j) {
        std::clog << "\rLineas restantes: " << j << " " << std::flush;
        for (int i = 0; i < image_width; ++i) {
            vec3 pixel_color = {0, 0, 0};
            
            // Antialiasing: Promediamos varios rayos vecinos
            for (int s = 0; s < samples_per_pixel; ++s) {
                double u = (i + random_double()) / (image_width-1);
                double v = (j + random_double()) / (image_height-1);
                vec3 r = lower_left_corner + u*horizontal + v*vertical - origin;
                pixel_color = pixel_color + ray_color(origin, r, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }
    std::clog << "\rHecho.                 \n";
}