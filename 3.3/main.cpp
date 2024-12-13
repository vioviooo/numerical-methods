#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// * Функция для решения нормальной системы МНК
std::vector<double> solveNormalSystemMNK(const std::vector<double> &x, const std::vector<double> &y, int degree) {
    int n = x.size();
    int m = degree + 1;
    std::vector<std::vector<double>> vec1(m, std::vector<double>(m, 0.0));
    std::vector<double> vec2(m, 0.0);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < n; ++k) {
                vec1[i][j] += std::pow(x[k], i + j);
            }
        }
        for (int k = 0; k < n; ++k) {
            vec2[i] += std::pow(x[k], i) * y[k];
        }
    }

    // * Метод Гаусса
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            double factor = vec1[j][i] / vec1[i][i];
            for (int k = i; k < m; ++k) {
                vec1[j][k] -= factor * vec1[i][k];
            }
            vec2[j] -= factor * vec2[i];
        }
    }

    std::vector<double> coeffs(m);
    for (int i = m - 1; i >= 0; --i) {
        coeffs[i] = vec2[i];
        for (int j = i + 1; j < m; ++j) {
            coeffs[i] -= vec1[i][j] * coeffs[j];
        }
        coeffs[i] /= vec1[i][i];
    }

    return coeffs;
}

// * Функция для вычисления значений многочлена
double evaluatePolynomial(const std::vector<double> &coeffs, double x) {
    double result = 0.0;
    double power = 1.0;
    for (const auto &coeff : coeffs) {
        result += coeff * power;
        power *= x;
    }
    return result;
}

// * Функция для вычисления суммы квадратов ошибок
double calculateErrorSum(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &coeffs) {
    double errorSum = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double predicted = evaluatePolynomial(coeffs, x[i]);
        errorSum += std::pow(predicted - y[i], 2);
    }
    return errorSum;
}

int main() {
    std::vector<double> x = {-0.9, 0.0, 0.9, 1.8, 2.7, 3.6};
    std::vector<double> y = {-0.36892, 0.0, 0.36892, 0.85408, 1.7856, 6.3138};

    // * Аппроксимация 1-ой степени
    auto coeffs1 = solveNormalSystemMNK(x, y, 1);
    double errSum1 = calculateErrorSum(x, y, coeffs1);

    std::cout << "> Коэффициенты многочлена 1-й степени:\n";
    for (size_t i = 0; i < coeffs1.size(); ++i) {
        std::cout << "coef_" << i << " = " << coeffs1[i] << "\n";
    }
    std::cout << "> Сумма квадратов ошибок 1-й степени: " << errSum1 << "\n\n";

    // * Аппроксимация 2-ой степени
    auto coeffs2 = solveNormalSystemMNK(x, y, 2);
    double errSum2 = calculateErrorSum(x, y, coeffs2);

    std::cout << "> Коэффициенты многочлена 2-й степени:\n";
    for (size_t i = 0; i < coeffs2.size(); ++i) {
        std::cout << "coef_" << i << " = " << coeffs2[i] << "\n";
    }

    std::cout << "> Сумма квадратов ошибок 2-й степени: " << errSum2 << "\n";

    return 0;
}
