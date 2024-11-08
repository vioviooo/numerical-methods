#include <iostream>
#include <cmath>
#include <functional>
#include <stdexcept>

using namespace std;

// f(x)
double f(double x) {
    return sqrt(1 - x * x) - exp(x) + 0.1;
}

// производная f(x)
double f_prime(double x) {
    return -x / sqrt(1 - x * x) - exp(x);
}

// вторая производная
double f_double_prime(double x) {
    return -exp(x) - (1 / sqrt(1 - x * x)) - (x * x) / pow(1 - x * x, 1.5);
}

// функция для метода простых итераций, чтобы метод сходился
double g(double x) {
    return log(0.1 + sqrt(1 - x * x));
}

double bisection(double a, double b, double epsilon) {
    double mid;
    while ((b - a) / 2 > epsilon) {
        mid = (a + b) / 2;
        if (f(mid) == 0.0)
            return mid;
        else if (f(a) * f(mid) < 0)
            b = mid;
        else
            a = mid;
    }
    return (a + b) / 2;
}

double g_prime(double x) {
    return -x / (sqrt(1 - x * x) * (0.1 + sqrt(1 - x * x)));
}

double compute_lipschitz_constant(double a, double b) {
    const int N = 1000;
    double max_q = 0.0;
    for (int i = 0; i <= N; ++i) {
        double x = a + i * (b - a) / N;
        double q = fabs(g_prime(x));
        if (q > max_q) {
            max_q = q;
        }
    }
    return max_q;
}

double simple_iteration(double x0, double epsilon) {
    double q = compute_lipschitz_constant(0.0, 0.1);
    std::cout << "> q = " << q << '\n';
    if (q >= 1) {
        throw std::runtime_error("Ошибка! Проверка Липшица не выполнена");
    }
    double x = x0, x_prev;
    do {
        x_prev = x;
        x = g(x_prev); // * (q / (1 - q))
    } while ((q / (1 - q)) * fabs(x - x_prev) > epsilon);
    return x;
}

double newton(double x0, double epsilon) {
    double x = x0;
    while (fabs(f(x)) > epsilon) {
        if (f(x) / f_double_prime(x) <= 0) {
            cout << "Ошибка! Условие f(x) / f''(x) > 0 не выполнено\n";
        }
        
        x = x - f(x) / f_prime(x);
    }
    return x;
}

double secant(double x0, double x1, double epsilon) {
    double x2;
    while (fabs(x1 - x0) > epsilon) {
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        x0 = x1;
        x1 = x2;
    }
    return x1;
}


int main() {
    double epsilon = 1e-6;

    cout << "> Метод дихотомии: " << bisection(0.0, 0.1, epsilon) << '\n';

    double res = simple_iteration(0.1, epsilon);
    cout << "> Метод простой итерации: " << res << '\n';

    cout << "> Метод Ньютона: " << newton(0.1, epsilon) << '\n';

    cout << "> Метод секущих: " << secant(0.1, 0.2, epsilon) << '\n';

    return 0;
}