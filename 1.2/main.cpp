#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

vector <double> solveThomasAlgorithm(const vector <double> &a, const vector <double> &b, const vector <double> &c, const vector <double> &d) {
    int n = b.size();
    vector <double> p(n), q(n), x(n);

    // прямой ход
    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * p[i - 1];
        p[i] = -c[i] / denom;
        q[i] = (d[i] - a[i] * q[i - 1]) / denom;
    }

    // обратный ход
    x[n - 1] = q[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    return x;
}

bool isTridiagonal(const vector <vector <double>> &matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (abs(i - j) > 1 && matrix[i][j] != 0) {
                return false;
            }
        }
    }
    return true;
}

double findDeterminant(const vector <double> &a, const vector <double> &b, const vector <double> &c) {
    int n = b.size();
    vector <double> p(n), q(n), x(n);

    p[0] = -c[0] / b[0];

    double res = b[0] + a[0];
    for (int i = 1; i < n; ++i) {
        res *= b[i] + a[i] * p[i - 1];
        double denom = b[i] + a[i] * p[i - 1];
        p[i] = -c[i] / denom;
    }

    return res;
}

bool checkSolution(const vector<vector<double>> &matrix, const vector<double> &x, const vector<double> &d) {
    int n = matrix.size();
    vector<double> m(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            m[i] += matrix[i][j] * x[j];
        }
        if (fabs(m[i] - d[i]) > 1e-6) return false;
    }
    return true;
}


bool isStable(const vector <double> &a, const vector <double> &b, const vector <double> &c, const vector <double> &d) {
    int n = b.size();
    vector <double> p(n);

    p[0] = -c[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * p[i - 1];
        p[i] = -c[i] / denom;
    }

    for (int i = 0; i < n; ++i) {
        if (fabs(p[i]) > 1) return false;
    }

    return true;
}

int main() {
    int n;
    
    cout << "> Введите n - размерность матрицы: ";
    cin >> n;

    vector <vector <double>> matrix(n, vector <double>(n + 1));

    cout << "> Введите матрицу вместе со столбцом свободных членов:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            cin >> matrix[i][j];
        }
    }

    if (!isTridiagonal(matrix)) {
        cout << "> Матрица не является тридиагональной.\n";
        return 0;
    }

    vector <double> a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);

    for (int i = 0; i < n; ++i) {
        if (i > 0) a[i] = matrix[i][i - 1]; // поддиагональ
        b[i] = matrix[i][i]; // главная диагональ
        if (i < n - 1) c[i] = matrix[i][i + 1]; // наддиагональ
        d[i] = matrix[i][n]; // вектор свободных членов 
    }

    bool flag = isStable(a, b, c, d);

    cout << ((flag == true) ? "> Метод прогонки устойчив.\n" : "Метод прогонки НЕ устойчив.\n");

    vector <double> x = solveThomasAlgorithm(a, b, c, d);

    cout << "> Итоговое решение: ";
    for (double elem : x) {
        cout << setprecision(6) << fixed << elem << ' ';
    }
    cout << '\n';

    if (checkSolution(matrix, x, d)) {
        cout << "> Решение корректно.\n";
    } else {
        cout << "> Решение НЕ корректно.\n";
    }

    double det = findDeterminant(a, b, c);
    cout << "> Определитель det = " << setprecision(6) << fixed << det << '\n';



    return 0;
}
// * add checking the correctness of the form of the input matrix
// * add det finding function
// * add checking the устойчивость
// * add function that checks the correctness of the solution