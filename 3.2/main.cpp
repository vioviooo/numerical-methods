#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// * cтруктура для хранения коэффициентов сплайна
struct Spline {
    double a, b, c, d; // * коэффициенты для каждого интервала
};

// * функция для решения системы уравнений для коэффициентов сплайна
vector<Spline> cubicSplineInterpolation(const vector<double>& x, const vector<double>& f) {
    int n = x.size();
    vector<Spline> splines(n - 1);

    // * математические коэффициенты для решения системы
    vector<double> h(n - 1), alpha(n - 1), l(n), mu(n), z(n);
    for (int i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    // * вычисление alpha, l, mu, z

    for (int i = 1; i < n - 1; ++i) {
        alpha[i] = (3.0 / h[i]) * (f[i + 1] - f[i]) - (3.0 / h[i - 1]) * (f[i] - f[i - 1]);
    }

    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    // * прямой ход
    for (int i = 1; i < n - 1; ++i) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1.0;
    z[n - 1] = 0.0;

    vector<double> c(n), b(n - 1), d(n - 1), a(n - 1);
    
    c[n - 1] = 0.0;

    // * вычисление a, b, c, d = обратный ход
    for (int j = n - 2; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (f[j + 1] - f[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        a[j] = f[j];
    }

    // * составляем сплайн
    for (int i = 0; i < n - 1; ++i) {
        splines[i].a = a[i];
        splines[i].b = b[i];
        splines[i].c = c[i];
        splines[i].d = d[i];
    }

    return splines;
}

// * функция для вычисления значения сплайна в точке x
double splineValue(const vector<Spline>& splines, const vector<double>& x, double x_value) {
    int n = x.size();

    // * находим интервал, к которому принадлежит точка x_value
    int i = 0;
    while (i < n - 1 and x_value > x[i + 1]) {
        ++i;
    }

    double dx = x_value - x[i];

    return splines[i].a + splines[i].b * dx + splines[i].c * dx * dx + splines[i].d * dx * dx * dx;
}

int main() {
    vector<double> x = {0.0, 0.9, 1.8, 2.7, 3.6};
    vector<double> f = {0.0, 0.36892, 0.85408, 1.7856, 6.3138};

    vector<Spline> splines = cubicSplineInterpolation(x, f);

    std::cout << "> Вывод всех коэф-ов: " << '\n';
    for (int i = 0; i < splines.size(); ++i) {
        cout << "Interval: " << fixed << setprecision(5) << x[i] << ' ' << x[i + 1] << "\nSpline values: " << splines[i].a << ' ' << splines[i].b << ' ' << splines[i].c << ' ' << splines[i].d << "\n\n";
    }

    double K = 1.5;
    double result = splineValue(splines, x, K);

    cout << "> Значение функции в точке x = " << K << " : " << result << endl;

    return 0;
}
