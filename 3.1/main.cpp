#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <sstream>

double f(double x) {
    if (std::abs(x) < 1e-6) {
        return 0; // * tan(0) = 0
    } else {
        return -tan(x);
    }
}

double d4f(double x) {
    if (std::abs(x) < 1e-6) {
        return 0.0;
    } else {
        return cos(x) * (6 * sin(x) * pow(cos(x), -4) + 12 * cos(x) * sin(x) * pow(cos(x), -5) + 30 * pow(sin(x), 3) * pow(cos(x), -6)) - sin(x) * (2 * pow(cos(x), -3) + 6 * pow(sin(x), 2) * pow(cos(x), -5));
    }
}

double divided_difference(const std::vector<double> &x, const std::vector<double> &y, int i, int j) {
    if (i == j) {
        return y[i];
    } else {
        return (divided_difference(x, y, i + 1, j) - divided_difference(x, y, i, j - 1)) / (x[j] - x[i]);
    }
}
// * результ -разделенная разность умножается на соответствующее произведение.
double newton_polynomial(const std::vector<double> &x, const std::vector<double> &y, double x_star) {
    double result = y[0];
    double product = 1.0;
    for (size_t i = 1; i < x.size(); ++i) {
        product *= (x_star - x[i - 1]);
        result += divided_difference(x, y, 0, i) * product;
    }
    return result;
}

double lagrange_polynomial(const std::vector<double> &x, const std::vector<double> &y, double x_star) {
    double result = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double term = y[i];
        for (size_t j = 0; j < x.size(); ++j) {
            if (j != i)
                term *= (x_star - x[j]) / (x[i] - x[j]);
        }
        result += term;
    }
    return result;
}

double find_max_fourth_derivative(const std::vector<double> &x, double min_x, double max_x) {
    double max_value = 0.0;
    int num_points = 10000;
    double step = (max_x - min_x) / num_points;

    for (double cur_x = min_x; cur_x <= max_x; cur_x += step) {
        if (std::abs(cur_x) < 1e-6)
            continue;
        double value = std::abs(d4f(cur_x));
        max_value = std::max(max_value, value);
    }
    return max_value;
}

std::string format_double(double value, int precision = 6) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

// * начальное значение на знаменатель базисного многочлена
std::string build_lagrange_polynomial(const std::vector<double> &x, const std::vector<double> &y) {
    std::string result;
    for (size_t i = 0; i < x.size(); ++i) {
        double term = y[i];
        double den = 1.0;
        for (size_t j = 0; j < x.size(); ++j) {
            if (i != j) {
                den *= (x[i] - x[j]);
            }
        }
        double coeff = term / den;
        if (i > 0) {
            result += (coeff >= 0 ? " + " : " - ");
        } else if (coeff < 0) {
            result += "-";
        }
        result += format_double(std::abs(coeff));

        for (size_t j = 0; j < x.size(); ++j) {
            if (i != j) {
                result += " * (x - ";
                result += format_double(x[j]);
                result += ")";
            }
        }
    }
    return result;
}

std::string build_newton_polynomial(const std::vector<double> &x, const std::vector<double> &y) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(6);
    std::vector<double> diff;
    for (size_t i = 0; i < y.size(); ++i) {
        diff.push_back(divided_difference(x, y, 0, i));
    }
    ss << diff[0];
    for (size_t i = 1; i < diff.size(); ++i) {
        ss << (diff[i] >= 0 ? " + " : " - ");
        ss << std::abs(diff[i]);
        for (size_t j = 0; j < i; ++j)
            ss << " * (x - " << x[j] << ")";
    }
    return ss.str();
}

void solve_interpolation(const std::vector<double> &x, double x_star) {
    std::vector<double> y;

    for (double xi : x) {
        y.push_back(f(xi));
    }

    double lagrange_val = lagrange_polynomial(x, y, x_star);
    double newton_val = newton_polynomial(x, y, x_star);
    double actual_val = f(x_star);

    double min_x = *std::min_element(x.begin(), x.end());
    double max_x = *std::max_element(x.begin(), x.end());
    double max_deriv = find_max_fourth_derivative(x, min_x, max_x);

    std::string lagrange_poly = build_lagrange_polynomial(x, y);
    std::string newton_poly = build_newton_polynomial(x, y);

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Интерполяция для точек: ";
    for (double val : x) {
        std::cout << std::to_string(val) << " ";
    }
    std::cout << "\n";
    std::cout << "Значение X_star = " << std::to_string(x_star) << "\n\n";

    std::cout << "Многочлен Лагранжа:\n";
    std::cout << "> Полином:\n" << lagrange_poly << "\n\n";
    std::cout << "> Значение:\n" << std::to_string(lagrange_val) << "\n";
    std::cout << "> Погрешность:\n" << std::to_string(std::abs(actual_val - lagrange_val)) << "\n\n";

    std::cout << "Многочлен Ньютона:\n";
    std::cout << "> Полином:\n" << newton_poly << "\n\n";
    std::cout << "> Значение:\n" << std::to_string(newton_val) << "\n";
    std::cout << "> Погрешность:\n" << std::to_string(std::abs(actual_val - newton_val)) << "\n\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    double x_star = 3.0 * M_PI / 16.0;

    try {
        solve_interpolation({0, M_PI / 8, 2 * M_PI / 8, 3 * M_PI / 8}, x_star);
        std::cout << '\n' << '\n' << '\n';
        solve_interpolation({0, M_PI / 8, M_PI / 3, 3 * M_PI / 8}, x_star);
    } catch (const std::exception &e) {
        std::cerr << "* Ошибка: " << e.what() << '\n';
        return 1;
    }

    return 0;
}

// начальное значение члена многочлена.
//   * den = 1.0; — инициализация знаменателя базисного многочлена Лагранжа.
// коэф = начальное значение на знаменатель базисного многочлена
// Цикл  for (size_t i = 0; i < y.size(); ++i) вычисляет разделенные разности  с помощью функции divided_difference(x, y, 0, i) и сохраняет их в векторе diff.
