#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <tuple>

// функции для вычисления f1(x1, x2) и f2(x1, x2)
double f1(double x1, double x2, double a) {
    return (x1 * x1 + a * a) * x2 - a * a * a;
}

double f2(double x1, double x2, double a) {
    return (x1 - a / 2) * (x1 - a / 2) + (x2 - a / 2) * (x2 - a / 2) - a * a;
}

// функции для вычисления частных производных (для метода ньютона)
double df1_dx1(double x1, double x2, double a) {
    return 2 * x1 * x2;
}

double df1_dx2(double x1, double x2, double a) {
    return x1 * x1 + a * a;
}

double df2_dx1(double x1, double x2, double a) {
    return 2 * (x1 - a / 2);
}

double df2_dx2(double x1, double x2, double a) {
    return 2 * (x2 - a / 2);
}

double dg1_dx1(double x1, double a) {
    return -2 * a * a * a * x1 / std::pow(x1 * x1 + a * a, 2);
}

double dg2_dx2(double x2, double a) {
    double denom = std::sqrt(a * a - (x2 - a / 2) * (x2 - a / 2));
    return (x2 - a / 2) / denom;
}

double compute_lipschitz_constant(double a, double x1, double x2) {
    double q1 = std::abs(dg1_dx1(x1, a));
    double q2 = std::abs(dg2_dx2(x2, a));
    return std::max(q1, q2);
}

// функция для проверки условия Липшица с константой q
bool check_lipschitz_condition(double q) {
    if (q >= 1.0) {
        std::cerr << "Ошибка! Условие Липшица не выполнено: q >= 1\n";
        return false;
    }
    return true;
}

// функция для вычисления определителя матрицы 2x2
double determinant(double a11, double a12, double a21, double a22) {
    return a11 * a22 - a12 * a21;
}

std::tuple<double, double> simple_iteration(double epsilon, double a, double &x1, double &x2, int &iterCount) {
    double q = compute_lipschitz_constant(a, x1, x2);
    
    std::cout << "> q = " << q << '\n';

    if (!check_lipschitz_condition(q)) {
        throw std::runtime_error("Ошибка! Условие Липшица не выполнено\n");
    }
    
    // std::cout << "\n> Метод простых итераций ----\n";
    // std::cout << "\nНачальное приближение \n--->\n";
    // std::cout << "x1 = " << x1 << '\n';
    // std::cout << "x2 = " << x2 << '\n';

    std::vector<double> prev_x = {x1, x2};
    do {
        // std::cout << "\nИтерация " << iterCount + 1 << "\n";

        prev_x = {x1, x2};

        x2 = a * a * a / (x1 * x1 + a * a);                          // Используем f1(x1, x2) = 0
        x1 = a / 2 + std::sqrt(a * a - (x2 - a / 2) * (x2 - a / 2)); // Используем f2(x1, x2) = 0

        // std::cout << "x1 = " << x1 << '\n';
        // std::cout << "x2 = " << x2 << '\n';

        ++iterCount;
    } while (std::max(std::abs(x1 - prev_x[0]), std::abs(x2 - prev_x[1])) > epsilon);
    
    return {x1, x2};
}

std::tuple<double, double> newton_method(double epsilon, double a, double &x1, double &x2, int &iterCount) {
    // std::cout << "\n---- Метод Ньютона ----\n";
    // std::cout << "\nНачальное приближение \n--->\n";
    // std::cout << "x1 = " << x1 << '\n';
    // std::cout << "x2 = " << x2 << '\n';
    std::vector<double> prev_x = {x1, x2};
    do {
        // std::cout << "\nИтерация " << iterCount + 1 << "\n";

        prev_x = {x1, x2};

        // определитель якобиана
        double det = determinant(df1_dx1(x1, x2, a), df1_dx2(x1, x2, a), df2_dx1(x1, x2, a), df2_dx2(x1, x2, a));
        if (std::abs(det) < epsilon) {
            throw std::runtime_error("Определитель якобиана близок к нулю!");
        }

        x1 = x1 - determinant(f1(x1, x2, a), df1_dx2(x1, x2, a), f2(x1, x2, a), df2_dx2(x1, x2, a)) / det;
        x2 = x2 - determinant(df1_dx1(x1, x2, a), f1(x1, x2, a), df2_dx1(x1, x2, a), f2(x1, x2, a)) / det;
        // std::cout << "x1 = " << x1 << '\n';
        // std::cout << "x2 = " << x2 << '\n';
        ++iterCount;
    } while (std::max(std::abs(x1 - prev_x[0]), std::abs(x2 - prev_x[1])) > epsilon);
    
    return {x1, x2};
}

std::pair<double, double> find_x0(double a1, double b1, double a2, double b2, double a) {
    double x1_mid = (a1 + b1) / 2.0;
    double x2_mid = (a2 + b2) / 2.0;

    if (std::abs(df1_dx1(a1, x2_mid, a)) < 1.0) {
        return {a1, x2_mid};
    } else if (std::abs(df1_dx1(b1, x2_mid, a)) < 1.0) {
        return {b1, x2_mid};
    }

    if (std::abs(df2_dx2(x1_mid, a2, a)) < 1.0) {
        return {x1_mid, a2};
    } else if (std::abs(df2_dx2(x1_mid, b2, a)) < 1.0) {
        return {x1_mid, b2};
    }

    return {x1_mid, x2_mid};
}


int main() {
    double epsilon = 0.0001; // Задаем точность
    double a = 4;            // Значение a

    double a1 = 4, a2 = 0;
    double b1 = 6, b2 = 2;

    // double x1 = (a1 + b1) * 0.5;
    // double x2 = (a2 + b2) * 0.5;

    auto [x1, x2] = find_x0(a1, b1, a2, b2, a);

    double maxF1 = std::max(std::abs(f1(a1, a2, a)), std::abs(f1(b1, b2, a)));
    double maxF2 = std::max(std::abs(f2(a1, a2, a)), std::abs(f2(b1, b2, a)));

    std::cout << "\n> Максимумы по модулю:\n";
    std::cout << "> |f1(x1, x2)| max = " << maxF1 << '\n';
    std::cout << "> |f2(x1, x2)| max = " << maxF2 << '\n';

    int iterCount = 0;

    // * подставить значения области и найти макс по модулю

    try {
        auto [simpleX1, simpleX2] = simple_iteration(epsilon, a, x1, x2, iterCount);
        std::cout << "\n> Результат решения (метод простых итераций): \n\n";
        std::cout << "> x1 = " << simpleX1 << '\n';
        std::cout << "> x2 = " << simpleX2 << '\n';
        std::cout << "> Количество итераций: " << iterCount << "\n\n";

        iterCount = 0; 
        x1 = (a1 + b1) * 0.5;
        x2 = (a2 + b2) * 0.5;
        auto [newtonX1, newtonX2] = newton_method(epsilon, a, x1, x2, iterCount);
        std::cout << "\n> Результат решения (метод Ньютона): \n\n";
        std::cout << "> x1 = " << newtonX1 << '\n';
        std::cout << "> x2 = " << newtonX2 << '\n';
        std::cout << "> Количество итераций: " << iterCount << "\n\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    return 0;
}
