#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

double exactSolution(double x) {
    return exp(x * x) + exp(x * sqrt(2)) + exp(-x * sqrt(2));
}

// * define the first-order system: 
// * y1 = y, y2 = y' => y2' = 2 * y1 + 4 * x^2 * exp(x^2)
std::pair<double, double> systemEquationsPair(double x, double y1, double y2) {
    double dy1_dx = y2;
    double dy2_dx = 2 * y1 + 4 * x * x * exp(x * x);
    return {dy1_dx, dy2_dx};
}

void eulerMethod(double h, double x0, double y1_0, double y2_0, int steps, std::vector<double> &res) {
    double x = x0, y1 = y1_0, y2 = y2_0;
    for (int i = 0; i <= steps; ++i) {
        res.push_back(y1);

        auto [dy1, dy2] = systemEquationsPair(x, y1, y2);
        double y1_pred = y1 + h * dy1;
        double y2_pred = y2 + h * dy2;

        auto [dy1_corr, dy2_corr] = systemEquationsPair(x + h, y1_pred, y2_pred);
        y1 += h * (dy1 + dy1_corr) / 2;
        y2 += h * (dy2 + dy2_corr) / 2;
        x += h;
    }
}

void rungeKutta4(double h, double x0, double y1_0, double y2_0, int steps, std::vector<double> &res) {
    double x = x0, y1 = y1_0, y2 = y2_0;
    
    for (int i = 0; i <= steps; ++i) {
        res.push_back(y1);

        auto [k1_y1, k1_y2] = systemEquationsPair(x, y1, y2);
        auto [k2_y1, k2_y2] = systemEquationsPair(x + h / 2, y1 + h * k1_y1 / 2, y2 + h * k1_y2 / 2);
        auto [k3_y1, k3_y2] = systemEquationsPair(x + h / 2, y1 + h * k2_y1 / 2, y2 + h * k2_y2 / 2);
        auto [k4_y1, k4_y2] = systemEquationsPair(x + h, y1 + h * k3_y1, y2 + h * k3_y2);

        y1 += h * (k1_y1 + 2 * k2_y1 + 2 * k3_y1 + k4_y1) / 6;
        y2 += h * (k1_y2 + 2 * k2_y2 + 2 * k3_y2 + k4_y2) / 6;
        x += h;
    }
}

void adamsMethod(double h, double x0, double y1_0, double y2_0, int steps, std::vector<std::pair<double, double>> &res) {
    res.clear();
    res.push_back({y1_0, y2_0});
    
    std::vector<std::pair<double, double>> derivatives;

    double x = x0;

    for (int i = 0; i < 3; ++i) {
        auto [k1_y1, k1_y2] = systemEquationsPair(x, res[i].first, res[i].second);
        auto [k2_y1, k2_y2] = systemEquationsPair(x + h / 2, res[i].first + h / 2 * k1_y1, res[i].second + h / 2 * k1_y2);
        auto [k3_y1, k3_y2] = systemEquationsPair(x + h / 2, res[i].first + h / 2 * k2_y1, res[i].second + h / 2 * k2_y2);
        auto [k4_y1, k4_y2] = systemEquationsPair(x + h, res[i].first + h * k3_y1, res[i].second + h * k3_y2);

        double y1_next = res[i].first + h / 6 * (k1_y1 + 2 * k2_y1 + 2 * k3_y1 + k4_y1);
        double y2_next = res[i].second + h / 6 * (k1_y2 + 2 * k2_y2 + 2 * k3_y2 + k4_y2);

        res.push_back({y1_next, y2_next});

        derivatives.push_back({k1_y1, k1_y2});
        x += h;
    }

    auto [last_dy1, last_dy2] = systemEquationsPair(x, res.back().first, res.back().second);

    derivatives.push_back({last_dy1, last_dy2});

    for (int i = 3; i < steps; ++i) {
        auto [f_n, g_n] = derivatives[i];
        auto [f_n_1, g_n_1] = derivatives[i - 1];
        auto [f_n_2, g_n_2] = derivatives[i - 2];
        auto [f_n_3, g_n_3] = derivatives[i - 3];

        double y1_next = res[i].first + h / 24 * (55 * f_n - 59 * f_n_1 + 37 * f_n_2 - 9 * f_n_3);
        double y2_next = res[i].second + h / 24 * (55 * g_n - 59 * g_n_1 + 37 * g_n_2 - 9 * g_n_3);

        res.push_back({y1_next, y2_next});

        x += h;
        derivatives.push_back(systemEquationsPair(x, y1_next, y2_next));
    }
}

double rungeRombergError(double I_hk, double I_h, int p) {
    return (I_h + (I_h - I_hk) / (pow(2, p) - 1));
}

int main() {
    double h = 0.1, x0 = 0.0, y1_0 = 3.0, y2_0 = 0.0;
    double left = 0.0, right = 1.0;
    
    int steps = (right  - left) / h;

    std::vector<double> eulerResults, rkResults;
    std::vector<std::pair<double, double>> adamsResults;

    eulerMethod(h, x0, y1_0, y2_0, steps, eulerResults);
    rungeKutta4(h, x0, y1_0, y2_0, steps, rkResults);
    adamsMethod(h, x0, y1_0, y2_0, steps, adamsResults);

    std::cout << "\n> Полученные решения и отличие относительно точного: \n\n";

    std::cout << "x\t\tExact\t\tEuler\t\tEuler diff\tRunge-Kutt\tRK diff\t\tAdams\t\tAdams diff\n";
    for (int i = 0; i <= steps; ++i) {
        double x = x0 + i * h;
        double exact = exactSolution(x);
        std::cout << std::fixed << std::setprecision(8) << x << "\t" << exact << "\t" 
                  << eulerResults[i] << "\t"  << abs(exact - eulerResults[i]) << "\t" 
                  << rkResults[i] << "\t" << abs(exact - rkResults[i]) << "\t" 
                  << adamsResults[i].first << "\t" << abs(exact - adamsResults[i].first) << "\n";
    }

    std::vector<double> eulerResultsHalfStep, rkResultsHalfStep;
    std::vector<std::pair<double, double>> adamsResultsHalfStep;

    double hHalf = h / 2;
    int stepsHalf = static_cast<int>((right - left) / hHalf);
    eulerMethod(hHalf, x0, y1_0, y2_0, stepsHalf, eulerResultsHalfStep);
    rungeKutta4(hHalf, x0, y1_0, y2_0, stepsHalf, rkResultsHalfStep);
    adamsMethod(hHalf, x0, y1_0, y2_0, stepsHalf, adamsResultsHalfStep);

    std::cout << "\n> Улучшенные с помощью метода Рунге-Ромберга результаты и отличие от точного решения: \n\n";

    std::cout << "x\t\tExact\t\tEuler\t\tEuler RR\tRK4\t\tRK4 RR\t\tAdams\t\tAdams RR\n";

    for (int i = 0; i <= steps; ++i) {
        double x = x0 + i * h;
        double exact = exactSolution(x);

        double eulerHalfPoint = eulerResultsHalfStep[2 * i];
        double rkHalfPoint = rkResultsHalfStep[2 * i];
        double adamsHalfPoint = adamsResultsHalfStep[2 * i].first;

        double eulerRR = rungeRombergError(eulerResults[i], eulerHalfPoint, 2);
        double rkRR = rungeRombergError(rkResults[i], rkHalfPoint, 4);
        double adamsRR = rungeRombergError(adamsResults[i].first, adamsHalfPoint, 4);

        std::cout << std::fixed << std::setprecision(8) << x << "\t" << exact << "\t" 
                  << eulerRR << "\t" << abs(exact - eulerRR) << "\t"
                  << rkRR << "\t" << abs(exact - rkRR) << "\t"
                  << adamsRR << "\t" << abs(exact - adamsRR) << "\n";
    }

    return 0;
}