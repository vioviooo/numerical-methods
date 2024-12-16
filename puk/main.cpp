#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>

// * Функция для вычисления базисных функций
std::vector<double> basisFunctions(double x) {
    return {1.0, std::log(x), 1.0 / (x * x), 1.0 / x, x, x * x, x * x * x};
}

// * Решение нормальной системы для метода наименьших квадратов
std::vector<double> solveNormalSystemMNK(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    int m = basisFunctions(x[0]).size();
    std::vector<std::vector<double>> A(m, std::vector<double>(m, 0.0));
    std::vector<double> B(m, 0.0);

    // * Заполнение матрицы нормального уравнения и вектора
    for (int i = 0; i < n; ++i) {
        auto basis = basisFunctions(x[i]);
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < m; ++k) {
                A[j][k] += basis[j] * basis[k];
            }
            B[j] += basis[j] * y[i];
        }
    }

    // * Решение нормальной системы методом Гаусса
    for (int i = 0; i < m; ++i) {
        // * Масштабирование опорной строки
        double pivot = A[i][i];
        for (int j = i; j < m; ++j) {
            A[i][j] /= pivot;
        }
        B[i] /= pivot;

        // * Элиминция под опорным элементом
        for (int j = i + 1; j < m; ++j) {
            double factor = A[j][i];
            for (int k = i; k < m; ++k) {
                A[j][k] -= factor * A[i][k];
            }
            B[j] -= factor * B[i];
        }
    }

    // * Обратная подстановка
    std::vector<double> coeffs(m, 0.0);
    for (int i = m - 1; i >= 0; --i) {
        coeffs[i] = B[i];
        for (int j = i + 1; j < m; ++j) {
            coeffs[i] -= A[i][j] * coeffs[j];
        }
    }

    return coeffs;
}

// * Оценка аппроксимации
double evaluateApproximation(const std::vector<double>& coeffs, double x) {
    auto basis = basisFunctions(x);
    double result = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * basis[i];
    }
    return result;
}

// * Вычисление суммы квадратов ошибок
double calculateErrorSum(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& coeffs) {
    double errorSum = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double predicted = evaluateApproximation(coeffs, x[i]);
        errorSum += std::pow(predicted - y[i], 2);
    }
    return errorSum;
}

// * Чтение данных из файла
bool readDataFromFile(const std::string& filename, std::vector<double>& x, std::vector<double>& y) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Ошибка: Не удалось открыть файл " << filename << "\n";
        return false;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        double xi, yi;
        if (!(iss >> xi >> yi)) {
            std::cerr << "Ошибка: Неверный формат данных в файле.\n";
            return false;
        }
        x.push_back(xi); // * Преобразование значений x
        y.push_back(yi);
    }

    inputFile.close();
    return true;
}

int main() {
    std::vector<double> x, y;
    
    std::string filename = "water_data.txt";
    if (!readDataFromFile(filename, x, y)) {
        return 1;
    }
    
    auto coeffs = solveNormalSystemMNK(x, y);
    double errorSum = calculateErrorSum(x, y, coeffs);

    std::cout << "> Коэффициенты аппроксимации:\n";
    for (size_t i = 0; i < coeffs.size(); ++i) {
        std::cout << std::fixed << std::setprecision(20) << "coef_" << i << " = " << coeffs[i] << "\n";
    }
    std::cout << "> Сумма квадратов ошибок: " << errorSum << "\n";

    return 0;
}
