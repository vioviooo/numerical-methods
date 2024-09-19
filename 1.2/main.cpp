#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

const double EPS = 1e-6;

tuple<vector<double>, vector<double>, vector<double>> solveThomasAlgorithm(const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &d) {
    int n = b.size();
    vector<double> p(n), q(n), x(n);

    if (abs(b[0]) < EPS) {
        throw logic_error("* Ошибка: Деление на ноль");
    }

    // * прямой ход
    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * p[i - 1];
        if (abs(denom) < EPS) {
            throw logic_error("* Ошибка: Деление на ноль");
        }
        p[i] = -c[i] / denom;
        q[i] = (d[i] - a[i] * q[i - 1]) / denom;
    }

    // * обратный ход
    x[n - 1] = q[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    return {p, q, x};
}

double findDeterminant(const vector<double> &a, const vector<double> &b, const vector<double> &p) {
    int n = b.size();
    double res = b[0] + a[0];

    for (int i = 1; i < n; ++i) {
        res *= b[i] + a[i] * p[i - 1];
    }

    return res;
}

bool isStable(const vector<double> &p) {
    int n = p.size();

    for (int i = 0; i < n; ++i) {
        if (abs(p[i]) > 1.0) return false;
    }

    return true;
}

bool isStableABC(const vector<double> &a, const vector<double> &b, const vector<double> &c) {
    int n = b.size();

    for (int i = 0; i < n; ++i) {
        if (a[i] != 0 and c[i] != 0) {
            if (abs(b[i]) <= abs(a[i]) + abs(c[i])) return false;
        } else {
            ++i;
        }
    }

    return true;
}

bool isTridiagonal(const vector<double> &a, const vector<double> &b, const vector<double> &c) {

    for (int i = 1; i < a.size(); ++i) {
        if (abs(a[i]) < EPS) return false;
    }

    for (int i = 0; i < b.size(); ++i) {
        if (abs(b[i]) < EPS) return false;
    }

    for (int i = 0; i < c.size() - 1; ++i) {
        if (abs(c[i]) < EPS) return false;
    }

    return true;
}

bool checkSolution(const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &x, const vector<double> &d) {
    int n = b.size();
    vector<vector<double>> matrix(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        matrix[i][i] = b[i];
        if (i > 0) matrix[i][i - 1] = a[i];
        if (i < n - 1) matrix[i][i + 1] = c[i];
    }

    vector<double> m(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            m[i] += matrix[i][j] * x[j];
        }
        if (abs(m[i] - d[i]) > EPS) return false;
    }
    return true;
}

int main() {

    ifstream fin("test.txt");

    int n;
    
    fin >> n;

    if (n < 1) {
        std::cout << "> Размерность n должна быть больше 0\n";
        return 1;
    }

    vector <double> a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);

    for (int i = 0; i < n - 1; ++i) {
        fin >> c[i];
    }

    for (int i = 0; i < n; ++i) {
        fin >> b[i];
    }

    for (int i = 1; i < n; ++i) {
        fin >> a[i];
    }

    for (int i = 0; i < n; ++i) {
        fin >> d[i];
    }

    std::cout << ((isStableABC(a,b,c)) ? "> Метод прогонки устойчив.\n" : "Метод прогонки НЕ устойчив.\n");
    
    auto [p, q, x] = std::make_tuple(vector<double>{}, vector<double>{}, vector<double>{});
    try {
        std::tie(p, q, x) = solveThomasAlgorithm(a, b, c, d);
    } catch(std::exception &e) {
        std::cout << e.what() << '\n';
        return 1;
    }

    std::cout << ((isStable(p)) ? "> Метод прогонки устойчив.\n" : "Метод прогонки НЕ устойчив.\n");

    std::cout << "> Итоговое решение: ";
    for (double elem : x) {
        std::cout << setprecision(6) << fixed << elem << ' ';
    }
    std::cout << '\n';

    if (checkSolution(a, b, c, x, d)) {
        std::cout << "> Решение корректно.\n";
    } else {
        std::cout << "> Решение НЕ корректно.\n";
    }

    double det = findDeterminant(a, b, p);
    std::cout << "> Определитель det = " << setprecision(6) << fixed << det << '\n';

    return 0;
}