#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>

void print(const std::vector<double>& vec) {
    for (double val : vec) {
        std::cout << std::fixed << std::setprecision(7) << val << ' ';
    }
    std::cout << '\n';
}

void print(const std::vector<std::vector<double>>& vec) {
    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[i].size(); ++j) {
            std::cout << std::fixed << std::setprecision(7) << vec[i][j] << ' ';
        }
        std::cout << '\n';
    }
}

std::tuple<double, std::vector<double>> gaussian_elimination(std::vector<std::vector<double>> matrix, std::vector<double> coef) {
    int n = matrix.size();
    double det = 1;

    for (int i = 0; i < n; ++i) {
        int swap_row = i;
        // * выборка главного элемента для перестановки
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(matrix[k][i]) > std::abs(matrix[swap_row][i])) {
                swap_row = k;
            }
        }

        // * перестановка строк
        if (swap_row != i) {
            det *= -1; // * при перестановке строк знак определителя меняется
            std::swap(matrix[i], matrix[swap_row]);
            std::swap(coef[i], coef[swap_row]);
        }

        det *= matrix[i][i];

        // * прямой ход
        for (int k = i + 1; k < n; ++k) {
            if (matrix[i][i] == 0) {
                throw std::logic_error("* Деление на ноль\n");
            }
            double factor = matrix[k][i] / matrix[i][i];
            coef[k] -= factor * coef[i];
            for (int j = i; j < n; ++j) {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
        std::cout << "Приводим матрицу к верхнетреугольной, итерация " << i + 1 << ": " << '\n';
        print(matrix);
        std::cout << "\n\n";
    }

    // * обратный ход
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = coef[i] / matrix[i][i];
        for (int k = i - 1; k >= 0; --k) {
            coef[k] -= matrix[k][i] * x[i];
        }
    }

    return {det, x};
}

std::vector<std::vector<double>> gaussian_elimination_inverse(std::vector<std::vector<double>> matrix) {
    int n = matrix.size();

    // * создаем расширенную матрицу% слева - исходная справа - единичная
    std::vector<std::vector<double>> augmented_matrix(n, std::vector<double>(2 * n, 0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented_matrix[i][j] = matrix[i][j];
        }
        augmented_matrix[i][n + i] = 1.0;
    }
    // std::cout << '\n';
    // print(augmented_matrix);
    // std::cout << '\n';

    // * прямой ход
    for (int i = 0; i < n; ++i) {
        int swap_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(augmented_matrix[k][i]) > std::abs(augmented_matrix[swap_row][i])) {
                swap_row = k;
            }
        }

        std::swap(augmented_matrix[i], augmented_matrix[swap_row]);

        if (augmented_matrix[i][i] == 0) {
            throw std::logic_error("* Деление на ноль\n");
        }

        // * нормализация строки
        double divisor = augmented_matrix[i][i];
        for (int j = 0; j < 2 * n; ++j) {
            augmented_matrix[i][j] /= divisor;
        }

        // * вычитаем текущую строку из следующих строк
        for (int k = i + 1; k < n; ++k) {
            double factor = augmented_matrix[k][i];
            for (int j = 0; j < 2 * n; ++j) {
                augmented_matrix[k][j] -= factor * augmented_matrix[i][j];
            }
        }
        std::cout << '\n';
        print(augmented_matrix);
        std::cout << '\n';
    }

    // * обратный ход
    for (int i = n - 1; i >= 0; --i) {
        for (int k = i - 1; k >= 0; --k) {
            double factor = augmented_matrix[k][i];
            for (int j = 0; j < 2 * n; ++j) {
                augmented_matrix[k][j] -= factor * augmented_matrix[i][j];
            }
        }
    }

    std::vector<std::vector<double>> inverse_matrix(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverse_matrix[i][j] = augmented_matrix[i][n + j];
        }
    }

    return inverse_matrix;
}
template<typename T>
std::vector<std::vector<T>> multiply_matrices(const std::vector<std::vector<T>>& a, const std::vector<std::vector<T>>& b) {
    int a_rows = a.size();
    int a_cols = a[0].size();
    int b_rows = b.size();
    int b_cols = b[0].size();

    if (a_cols != b_rows) {
        throw std::invalid_argument("* Количество столбцов первой матрицы должно совпадать с количеством строк второй матрицы.");
    }

    std::vector<std::vector<T>> res(a_rows, std::vector<T>(b_cols, 0));

    for (int i = 0; i < a_rows; ++i) {
        for (int j = 0; j < b_cols; ++j) {
            for (int k = 0; k < a_cols; ++k) {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return res;
}


int main() {

    std::ifstream fin("matrix.txt");

    if (!fin) {
        std::cout << "* Файл не найден\n";
        return 1;
    }

    int n;
    fin >> n;

    // * матрица
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fin >> matrix[i][j];
        }
    }

    // * вектор свободных членов
    std::vector<double> coef(n);

    for (int i = 0; i < n; ++i) {
        fin >> coef[i];
    }

    auto [det, x] = std::make_tuple(0.0, std::vector<double>{});

    try {
        std::tie(det, x) = gaussian_elimination(matrix, coef);
    } catch (std::exception& e) {
        std::cout << e.what();
        return 1;
    }

    std::vector<std::vector<double>> inverse;
    try {
        inverse = gaussian_elimination_inverse(matrix);
    } catch (std::exception& e) {
        std::cout << e.what();
        return 1;
    }

    std::cout << "\n> Итоговое решение:\n";
    print(x);

    std::cout << "\n> Проверка a * x:\n";

    std::vector<std::vector<double>> x_matrix(x.size(), std::vector<double>(1));

    for (int i = 0; i < x.size(); ++i) {
        x_matrix[i][0] = x[i];
    }

    print(multiply_matrices(matrix, x_matrix));

    std::cout << "\n> Определитель:\ndet = " << det << '\n';

    std::cout << "\n> Обратная матрица:\n";
    print(inverse);

    std::cout << "\n> Проверка a * a^(-1):\n";
    print(multiply_matrices(matrix, inverse));

    return 0;
}

// * add нахождение обратнйо матрицы
// * calculate denominator