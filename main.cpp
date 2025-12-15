#include <iostream>

int main() {
    // 1. Definimos el tamaño de la imagen
    const int image_width = 256;
    const int image_height = 256;

    // 2. Imprimimos la cabecera del formato PPM (P3 = colores ASCII)
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    // 3. Los píxeles se escriben de arriba a abajo, de izquierda a derecha
    for (int j = 0; j < image_height; ++j) {
        // Imprimimos un log de progreso en la consola de error (cerr) para no ensuciar la imagen
        std::clog << "\rLíneas restantes: " << (image_height - j) << ' ' << std::flush;
        
        for (int i = 0; i < image_width; ++i) {
            // Normalizamos las coordenadas (de 0 a 1)
            double r = double(i) / (image_width - 1);
            double g = double(j) / (image_height - 1);
            double b = 0.25;

            // Convertimos de 0..1 a 0..255 (entero)
            int ir = static_cast<int>(255.999 * r);
            int ig = static_cast<int>(255.999 * g);
            int ib = static_cast<int>(255.999 * b);

            // Escribimos el pixel
            std::cout << ir << ' ' << ig << ' ' << ib << '\n';
        }
    }

    std::clog << "\rHecho.                 \n";
}