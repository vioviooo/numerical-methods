#include <iostream>
#include <cmath>

double f(double x) {
    return x / pow(3 * x + 4, 3);
}

double rectangleApproximationMethod(double (*foo)(double), double a, double b, double h) {
    double sum = 0.0;
    int n = round((b - a) / h);

    for (int i = 1; i <= n; ++i) {
        double mid = a + h * (i - 0.5);
        sum += foo(mid);
    }

    return sum * h; 
}

double trapezoidApproximationMethod(double (*foo)(double), double a, double b, double h) {
    double sum = 0.0;
    int n = (b - a) / h; 

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += foo(x);
    }

    sum += (foo(a) + foo(b)) / 2.0;
    return sum * h;
}

double simpsonsApproximationMethod(double (*foo)(double), double a, double b, double h) {
    int n = (b - a) / h;
    if (n % 2 != 0) {
        std::cerr << "Error: Number of intervals must be even for Метод Симпсона." << '\n';
        return -1.0;
    }

    double sum = foo(a) + foo(b);

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        if (i % 2 == 0) {
            sum += 2 * foo(x);
        } else {
            sum += 4 * foo(x);
        }
    }

    return (h / 3.0) * sum;
}

double rungeRombergError(double I_hk, double I_h, int p) {
    return (I_h + (I_h - I_hk) / (pow(2, p) - 1));
}

int main() {
    double a = -1.0;
    double b = 1.0;

    double h = 1e-6;
    double true_approximation = rectangleApproximationMethod(f, a, b, h);
    std::cout << "> True approximation: " << true_approximation << '\n'; 

    // * метод прямоугольиков
    double h1 = 0.5, h2 = 0.25;
    double I1_rect = rectangleApproximationMethod(f, a, b, h1);
    double I2_rect = rectangleApproximationMethod(f, a, b, h2);
    double error_rect = rungeRombergError(I1_rect, I2_rect, 2); 
    
    std::cout << "Метод прямоугольников (h1 = 0.5): " << I1_rect << '\n';
    std::cout << "Метод прямоугольников (h2 = 0.25): " << I2_rect << '\n';
    std::cout << "Рунге-Ромберга для Метода прямоугольников: " << error_rect << '\n';
    std::cout << "Метод прямоугольников error diff: " << abs(true_approximation - error_rect) << '\n';

    std::cout << '\n';

    // * метод трапеций
    double I1_trap = trapezoidApproximationMethod(f, a, b, h1);
    double I2_trap = trapezoidApproximationMethod(f, a, b, h2);
    double error_trap = rungeRombergError(I1_trap, I2_trap, 2); 
    
    std::cout << "Метод трапеций (h1 = 0.5): " << I1_trap << '\n';
    std::cout << "Метод трапеций (h2 = 0.25): " << I2_trap << '\n';
    std::cout << "Рунге-Ромберга для Метода трапеций: " << error_trap << '\n';
    std::cout << "Метод трапеций error diff: " << abs(true_approximation - error_trap) << '\n';

    std::cout << '\n';


    // * метод Симпсона
    if (static_cast<int>((b - a) / h1) % 2 != 0) {
        h1 = (b - a) / (static_cast<int>((b - a) / h1) + 1);
    }
    double I1_simp = simpsonsApproximationMethod(f, a, b, h1);
    double I2_simp = simpsonsApproximationMethod(f, a, b, h2);
    double error_simp = rungeRombergError(I1_simp, I2_simp, 4); 
    
    std::cout << "Метод Симпсона (h1 = 0.5): " << I1_simp << '\n';
    std::cout << "Метод Симпсона (h2 = 0.25): " << I2_simp << '\n';
    std::cout << "Рунге-Ромберга для Метод Симпсона: " << error_simp << '\n';
    std::cout << "Метод Симпсона error diff: " << abs(true_approximation - error_simp) << '\n';


    return 0;
}
