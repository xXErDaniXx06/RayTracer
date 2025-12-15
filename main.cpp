#include <iostream>
#include <cmath>

// --- PARTE COMPLEJA: Estructura de Vector 3D ---
// Esto nos permite sumar puntos y vectores como si fueran números simples
struct vec3 {
    double x, y, z;

    vec3 operator+(const vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    vec3 operator-(const vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    vec3 operator*(double t) const { return {x * t, y * t, z * t}; }
    vec3 operator/(double t) const { return {x / t, y / t, z / t}; }
    
    // Producto punto (básico para calcular ángulos y luz)
    double dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    
    // Longitud del vector
    double length() const { return std::sqrt(dot(*this)); }
    
    // Convertir a vector unitario (longitud 1)
    vec3 normalize() const { return *this / length(); }
};

// Función para escribir colores en el archivo
void write_color(std::ostream &out, vec3 pixel_color) {
    out << static_cast<int>(255.999 * pixel_color.x) << ' '
        << static_cast<int>(255.999 * pixel_color.y) << ' '
        << static_cast<int>(255.999 * pixel_color.z) << '\n';
}

// --- FÍSICA DEL RAYO ---

// Calcula si golpeamos la esfera y A QUÉ DISTANCIA (t)
double hit_sphere(const vec3& center, double radius, const vec3& origin, const vec3& direction) {
    vec3 oc = origin - center;
    
    // Ecuación cuadrática simplificada con vectores
    double a = direction.dot(direction);
    double b = 2.0 * oc.dot(direction);
    double c = oc.dot(oc) - radius * radius;
    double discriminant = b*b - 4*a*c;
    
    if (discriminant < 0) {
        return -1.0; // No golpea
    } else {
        return (-b - std::sqrt(discriminant)) / (2.0*a); // Devuelve la distancia 't' más cercana
    }
}

// Función principal que decide el color de cada píxel
vec3 ray_color(const vec3& origin, const vec3& direction) {
    // Definimos una esfera en (0,0,-1) con radio 0.5
    vec3 sphere_center = {0, 0, -1};
    double t = hit_sphere(sphere_center, 0.5, origin, direction);
    
    if (t > 0.0) {
        // ¡GOLPE!
        // Calculamos el punto exacto donde chocó el rayo: P = Origen + t*Dirección
        vec3 hit_point = origin + (direction * t);
        
        // Calculamos la "Normal": vector que apunta hacia afuera de la superficie
        // N = (PuntoGolpe - Centro) normalizado
        vec3 N = (hit_point - sphere_center).normalize();
        
        // Mapeamos las coordenadas de la normal (-1 a 1) a colores (0 a 1)
        return vec3{N.x+1, N.y+1, N.z+1} * 0.5;
    }
    
    // SI NO GOLPEA: Fondo cielo (degradado)
    vec3 unit_direction = direction.normalize();
    t = 0.5 * (unit_direction.y + 1.0);
    // Interpolación lineal (Lerp) entre blanco y azul
    return vec3{1.0, 1.0, 1.0} * (1.0 - t) + vec3{0.5, 0.7, 1.0} * t;
}

// --- MAIN ---
int main() {
    // Configuración de imagen
    const double aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // Configuración de cámara y viewport
    double viewport_height = 2.0;
    double viewport_width = aspect_ratio * viewport_height;
    double focal_length = 1.0;

    vec3 origin = {0, 0, 0};
    vec3 horizontal = {viewport_width, 0, 0};
    vec3 vertical = {0, viewport_height, 0};
    // Esquina inferior izquierda del "marco" virtual
    vec3 lower_left_corner = origin - (horizontal/2) - (vertical/2) - vec3{0, 0, focal_length};

    // Renderizado
    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j) {
        std::clog << "\rLineas restantes: " << j << " " << std::flush;
        for (int i = 0; i < image_width; ++i) {
            double u = double(i) / (image_width - 1);
            double v = double(j) / (image_height - 1);
            
            // Creamos el rayo
            vec3 direction = lower_left_corner + (horizontal * u) + (vertical * v) - origin;
            
            vec3 pixel_color = ray_color(origin, direction);
            write_color(std::cout, pixel_color);
        }
    }
    std::clog << "\rHecho.                 \n";
}