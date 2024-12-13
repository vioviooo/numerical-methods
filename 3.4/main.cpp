#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept> // * производные
#include <limits> 

using namespace std;

int nearestIndex(double eps, double x_star, const vector<double> &x) {
    int n = x.size();
    int best_index = 0;
    double min_diff = abs(x_star - x[0]);
    for (int i = 1; i < n; ++i) {
        double diff = abs(x_star - x[i]);
        if (diff < min_diff) {
            min_diff = diff;
            best_index = i;
        }
    }

    if (min_diff > eps) {
        throw runtime_error("* Ошибка: Значение x_star не найдено в заданном допуске");
    }
    return best_index;
}

double findDerivative(const vector<double> &x, const vector<double> &y, double target_x, int deriv_type) {
    int n = x.size();
    int index = nearestIndex(1e-6, target_x, x);

    if (deriv_type == 1 or deriv_type == 2 or deriv_type == 3 or deriv_type == 4) {
        // * проверка на невозможность вычисления производной на границах
        if (index == 0 and (deriv_type == 1 or deriv_type == 3 or deriv_type == 4)) {
            throw runtime_error("* Ошибка: невозможно вычислить производную в начале интервала");
        }
        if (index == n - 1 and (deriv_type == 2 or deriv_type == 3 or deriv_type == 4)) {
            throw runtime_error("* Ошибка: невозможно вычислить производную в конце интервала");
        }
        if (index == n - 2 and deriv_type == 4) {
            throw runtime_error("* Ошибка: невозможно вычислить вторую производную из-за недостатка данных");
        }

        if (deriv_type == 1) {
            return (y[index] - y[index - 1]) / (x[index] - x[index - 1]);
        } else if (deriv_type == 2) {
            return (y[index + 1] - y[index]) / (x[index + 1] - x[index]);
        } else if (deriv_type == 3) {
            return (y[index + 1] - y[index - 1]) / (x[index + 1] - x[index - 1]);
        } else if (deriv_type == 4) { // * 2-я производная
            double h = x[index + 1] - x[index];
            if (abs(h) < numeric_limits<double>::epsilon()) {
                throw runtime_error("* Ошибка: Шаг h слишком мал");
            }
            return (y[index + 1] - 2 * y[index] + y[index - 1]) / (h * h);
        }
    } else {
        throw runtime_error("* Ошибка: Неверный тип производной");
    }
    return 0;
}


int main() {
    vector<double> x = {1.0, 1.5, 2.0, 2.5, 3.0};
    vector<double> y = {0.0, 0.40547, 0.69315, 0.91629, 1.0986};
    double target_x = 2.0;

    std::cout << "> Точка x = " << target_x << '\n';

    try {
        std::cout << "> Левая производная: " << findDerivative(x, y, target_x, 1) << '\n';
        std::cout << "> Правая производная: " << findDerivative(x, y, target_x, 2) << '\n';
        std::cout << "> Центральная производная: " << findDerivative(x, y, target_x, 3) << '\n';
        std::cout << "> Вторая производная: " << findDerivative(x, y, target_x, 4) << '\n';

        std::cout << "> Центральная производная через левую и правую производную: " << (findDerivative(x, y, target_x, 1) + findDerivative(x, y, target_x, 2)) / 2.0 << '\n';
    } catch (const runtime_error &error) {
        std::cerr << "* Ошибка: " << error.what() << '\n';
    }

    return 0;
}
