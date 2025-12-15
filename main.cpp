#include <iostream>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm>

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

struct vec3 {
    double x, y, z;
    
    // --- OPERADORES BÁSICOS ---
    // Suma
    vec3 operator+(const vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    // Resta
    vec3 operator-(const vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    // Multiplicación por número (escalar)
    vec3 operator*(double t) const { return {x * t, y * t, z * t}; }
    // División por número
    vec3 operator/(double t) const { return {x / t, y / t, z / t}; }

    // --- ¡ESTOS ERAN LOS QUE FALTABAN! ---
    // 1. Negación (Menos unario): Para hacer "-v"
    vec3 operator-() const { return {-x, -y, -z}; }
    
    // 2. Multiplicación de Vector * Vector (para mezclar colores)
    vec3 operator*(const vec3& v) const { return {x * v.x, y * v.y, z * v.z}; }
    // -------------------------------------

    double dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    double length() const { return std::sqrt(x*x + y*y + z*z); }
    vec3 normalize() const { return *this / length(); }
};

// Operador auxiliar para permitir "double * vector" (ej: 0.5 * v)
vec3 operator*(double t, const vec3& v) { return {t*v.x, t*v.y, t*v.z}; }

// --- FÍSICA DE LA LUZ ---
vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2 * v.dot(n) * n;
}

// Refracción (Ley de Snell)
vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = std::fmin((-uv).dot(n), 1.0);
    vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec3 r_out_parallel = -std::sqrt(std::fabs(1.0 - r_out_perp.dot(r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}

// Aproximación de Schlick (Reflectancia en bordes de cristal)
double reflectance(double cosine, double ref_idx) {
    auto r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*std::pow((1 - cosine),5);
}

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

// --- GEOMETRÍA ---
double hit_sphere(const vec3& center, double radius, const vec3& origin, const vec3& direction) {
    vec3 oc = origin - center;
    double a = direction.dot(direction);
    double b = 2.0 * oc.dot(direction);
    double c = oc.dot(oc) - radius*radius;
    double discriminant = b*b - 4*a*c;
    if (discriminant < 0) return -1.0;
    return (-b - std::sqrt(discriminant)) / (2.0*a);
}

// --- MOTOR DE RENDERIZADO ---
vec3 ray_color(const vec3& origin, const vec3& direction, int depth) {
    if (depth <= 0) return {0,0,0};

    // ESCENA: 3 Esferas
    vec3 pos_izq = {-1.0, 0, -1};
    vec3 pos_cen = {0.0, 0, -1};
    vec3 pos_der = {1.0, 0, -1};
    double radio = 0.5;

    double t_final = infinity;
    int material = 0; // 1=Mate, 2=Metal, 3=Cristal
    vec3 pos_final;

    // Buscamos el choque más cercano
    double t;
    if ((t = hit_sphere(pos_izq, radio, origin, direction)) > 0.001 && t < t_final) { t_final = t; material = 1; pos_final = pos_izq; }
    if ((t = hit_sphere(pos_cen, radio, origin, direction)) > 0.001 && t < t_final) { t_final = t; material = 3; pos_final = pos_cen; }
    if ((t = hit_sphere(pos_der, radio, origin, direction)) > 0.001 && t < t_final) { t_final = t; material = 2; pos_final = pos_der; }

    if (t_final < infinity) {
        vec3 hit_point = origin + (direction * t_final);
        vec3 N = (hit_point - pos_final).normalize();
        
        // --- MATE (Rojo Arcilla) ---
        if (material == 1) {
            vec3 target = hit_point + N + random_in_unit_sphere();
            // Ahora sí funcionará la multiplicación vec3 * vec3
            return 0.5 * vec3{0.8, 0.3, 0.3} * ray_color(hit_point, target - hit_point, depth-1); 
        }
        
        // --- METAL (Plateado Brillante) ---
        if (material == 2) {
            vec3 reflected = reflect(direction.normalize(), N);
            vec3 scattered = reflected + (random_in_unit_sphere() * 0.0);
            if (scattered.dot(N) > 0)
                 return vec3{0.8, 0.8, 0.8} * ray_color(hit_point, scattered, depth-1);
            return {0,0,0};
        }

        // --- CRISTAL (Vidrio) ---
        if (material == 3) {
            double ir = 1.5; // Índice de refracción
            double refraction_ratio = direction.dot(N) < 0 ? (1.0/ir) : ir;
            vec3 unit_dir = direction.normalize();
            // Aquí usamos el operador unario (-) que arreglamos
            vec3 normal_correcta = direction.dot(N) < 0 ? N : -N; 

            double cos_theta = std::fmin((-unit_dir).dot(normal_correcta), 1.0);
            double sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);

            bool cannot_refract = refraction_ratio * sin_theta > 1.0;
            vec3 direction_out;

            if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
                direction_out = reflect(unit_dir, normal_correcta);
            else
                direction_out = refract(unit_dir, normal_correcta, refraction_ratio);

            return vec3{1.0, 1.0, 1.0} * ray_color(hit_point, direction_out, depth-1);
        }
    }

    // Fondo Cielo
    vec3 unit_direction = direction.normalize();
    double a = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - a)*vec3{1.0, 1.0, 1.0} + a*vec3{0.5, 0.7, 1.0};
}

int main() {
    // CAMBIA ESTO:
const int image_width = 1200; // Antes 400 -> Más resolución
const int image_height = 675; // Antes 225
const int samples_per_pixel = 200; // Antes 50 -> Menos ruido (pero tarda x4 más)
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
    std::clog << "\rHecho.                 \n";
}