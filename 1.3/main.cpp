#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>

void print_matrix(const std::vector<double> &x) {
    for (auto elem : x) {
        std::cout << std::fixed << std::setprecision(10) << elem << ' ';
    }
    std::cout << '\n';
}

void print_matrix(const std::vector<std::vector<double>> &matrix) {
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[0].size(); ++j) {
            std::cout << std::fixed << std::setprecision(10) << matrix[i][j] << ' ';
        }
    }
    std::cout << '\n';
}

bool is_diagonally_dominant(const std::vector<std::vector<double>> &matrix) {
    size_t n = matrix.size();
    for (size_t i = 0; i < n; ++i) {
        double sum_d = std::abs(matrix[i][i]);
        double sum = 0.0;

        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                sum += std::abs(matrix[i][j]);
            }
        }

        if (sum_d <= sum) {
            return false;
        }
    }
    return true;
}

// * функция подсчета нормы разности двух векторов (L-infinity norm)
double calculate_norm(const std::vector<double> &x, const std::vector<double> &prev_x) {
    double norm = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        norm = std::max(norm, std::abs(x[i] - prev_x[i]));
    }
    return norm;
}

// * метод простых итераций
std::vector<double> simple_iterations_algorithm(const std::vector<std::vector<double>> &matrix, const std::vector<double> &b, double EPS, int &iter) {
    size_t n = matrix.size();
    // * начальное приближение
    std::vector<double> x(n, 0.0);
    std::vector<double> prev_x(n, 0.0);

    do {
        prev_x = x;
        for (size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    sum += matrix[i][j] * prev_x[j];
                }
            }
            if (matrix[i][i] == 0) {
                throw std::logic_error("* Деление на ноль\n");
            }
            x[i] = (b[i] - sum) / matrix[i][i];
        }
        ++iter;
        double norm = calculate_norm(x, prev_x);
        std::cout << std::fixed << std::setprecision(10)<< "Simple iterations norm: " << norm << '\n';
    } while (calculate_norm(x, prev_x) > EPS);
    
    return x;
}

std::vector<double> seidel_algorithm(const std::vector<std::vector<double>> &matrix, const std::vector<double> &b, double EPS, int &iter) {
    size_t n = matrix.size();

    // * начальное приближение
    std::vector<double> x(n, 0.0);

    do {
        std::vector<double> prev_x = x;
        for (size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < i; ++j) {
                sum += matrix[i][j] * x[j];
            }
            for (size_t j = i + 1; j < n; ++j) {
                sum += matrix[i][j] * prev_x[j];
            }
            if (matrix[i][i] == 0) {
                throw std::logic_error("* Деление на ноль\n");
            }
            x[i] = (b[i] - sum) / matrix[i][i];
        }
        ++iter;

        double norm = calculate_norm(x, prev_x);
        std::cout << std::fixed << std::setprecision(10) << "Seidel algorithm norm: " << norm << '\n';
        
        if (norm <= EPS) {
            break;
        }
    } while (true);
    
    return x;
}

template<typename T>
std::vector<std::vector<T>> multiply_matrices(const std::vector<std::vector<T>> &a, const std::vector<std::vector<T>> &b) {
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

template<typename T>
std::vector<std::vector<T>> multiply_matrices(const std::vector<std::vector<T>> &a, const std::vector<T> &b) {
    int a_rows = a.size();
    int a_cols = a[0].size();
    int b_rows = b.size();
    int b_cols = 1;

    if (a_cols != b_rows) {
        throw std::invalid_argument("* Количество столбцов первой матрицы должно совпадать с количеством строк второй матрицы.");
    }

    std::vector<std::vector<T>> res(a_rows, std::vector<T>(b_cols, 0));

    for (int i = 0; i < a_rows; ++i) {
        for (int j = 0; j < b_cols; ++j) {
            for (int k = 0; k < a_cols; ++k) {
                res[i][j] += a[i][k] * b[k];
            }
        }
    }

    return res;
}

int main() {

    std::ifstream fin("matrix.txt");

    int n;
    fin >> n;

    std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fin >> matrix[i][j];
        }
    }

    if (!is_diagonally_dominant(matrix)) {
        std::cout << "* Матрица не является диагонально доминирующей.\n";
        return 1;
    }

    std::vector<double> b(n, 0.0);
    for (int i = 0; i < n; ++i) {
        fin >> b[i];
    }

    double EPS = 1e-6;

    fin >> EPS;

    int simple_iter = 0;
    int seidel_iter = 0;

    std::vector<double> simple_iter_result;
    try {
       simple_iter_result = simple_iterations_algorithm(matrix, b, EPS, simple_iter);
    } catch(std::exception &e) {
        std::cout << e.what() << '\n';
        return 1;
    }

    std::cout << "> Решение методом простых итераций:" << '\n';
    print_matrix(simple_iter_result);
    std::cout << "> Количество итераций: " << simple_iter << '\n';

    std::vector<double> gauss_seidel_result;
    try {
        gauss_seidel_result = seidel_algorithm(matrix, b, EPS, seidel_iter);
    } catch(std::exception &e) {
        std::cout << e.what() << '\n';
        return 1;
    }
    std::cout << "\n> Решение методом Зейделя:" << '\n';
    print_matrix(gauss_seidel_result);
    std::cout << "> Количество итераций: " << seidel_iter << '\n';

    std::vector<std::vector<double>> check_simple_iter;
    try {
        check_simple_iter = multiply_matrices(matrix, simple_iter_result);
    } catch(std::exception &e) {
        std::cout << e.what() << '\n';
    }

    std::cout << "\n> Проверка решения методом простых итераций:\n";
    print_matrix(check_simple_iter);

    std::vector<std::vector<double>> check_seidel;
    try {
        check_seidel = multiply_matrices(matrix, gauss_seidel_result);
    } catch(std::exception &e) {
        std::cout << e.what() << '\n';
        return 1;
    }

    std::cout << "\n> Проверка решения методом Зейделя:\n";
    print_matrix(check_seidel);

    return 0;
}