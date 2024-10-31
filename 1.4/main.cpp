#include <iostream>
#include <cmath>
#include <vector>

const double EPSILON = 1e-10;
const int MAX_ITER = 1e3;

void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double elem : row) {
            std::cout << elem << ' ';
        }
        std::cout << '\n';
    }
}

std::vector<double> matVecMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& x) {
    int n = matrix.size();
    std::vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += matrix[i][j] * x[j];
        }
    }
    return result;
}

void normalize(std::vector<double>& v) {
    double norm = 0.0;
    for (double elem : v) {
        norm += elem * elem;
    }
    norm = std::sqrt(norm);
    for (double& elem : v) {
        elem /= norm;
    }
}

void rotate(std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& V, int p, int q) {
    if (matrix[p][q] == 0) {
        return;
    }

    int n = matrix.size();

    double theta = 0.5 * atan2(2 * matrix[p][q], matrix[q][q] - matrix[p][p]);
    
    if (std::abs(matrix[p][p] - matrix[q][q]) < EPSILON) {
        theta = M_PI / 4;
    }
    
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);

    double app = matrix[p][p], aqq = matrix[q][q], apq = matrix[p][q];
    matrix[p][p] = cosTheta * cosTheta * app + sinTheta * sinTheta * aqq - 2 * sinTheta * cosTheta * apq;
    matrix[q][q] = sinTheta * sinTheta * app + cosTheta * cosTheta * aqq + 2 * sinTheta * cosTheta * apq;
    matrix[p][q] = matrix[q][p] = 0;

    for (int i = 0; i < n; i++) {
        if (i != p and i != q) {
            double aip = matrix[i][p], aiq = matrix[i][q];
            matrix[i][p] = matrix[p][i] = cosTheta * aip - sinTheta * aiq;
            matrix[i][q] = matrix[q][i] = sinTheta * aip + cosTheta * aiq;
        }
    }

    
    for (int i = 0; i < n; i++) {
        double vip = V[i][p], viq = V[i][q];
        V[i][p] = cosTheta * vip - sinTheta * viq;
        V[i][q] = sinTheta * vip + cosTheta * viq;
    }
}


void jacobiEigenvalMethod(std::vector<std::vector<double>>& matrix, std::vector<double>& eigenvals, std::vector<std::vector<double>>& eigenvecs) {
    int n = matrix.size();
    eigenvecs.assign(n, std::vector<double>(n, 0.0));

    
    for (int i = 0; i < n; i++) {
        eigenvecs[i][i] = 1.0;
    }

    for (int k = 0; k < MAX_ITER; k++) {
        int p = 0, q = 1;
        double max = 0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (std::abs(matrix[i][j]) > max) {
                    max = std::abs(matrix[i][j]);
                    p = i;
                    q = j;
                }
            }
        }
        if (max < EPSILON) {
            break; 
        }
        rotate(matrix, eigenvecs, p, q);
    }

    for (int i = 0; i < n; i++) {
        eigenvals[i] = matrix[i][i];
    }
}

std::pair<double, std::vector<double>> powerMethod(const std::vector<std::vector<double>>& matrix, int maxIterations = 1000, double eps = EPSILON) {
    int n = matrix.size();
    std::vector<double> x(n, 1.0);  
    double eigenval = 0.0;

    for (int k = 0; k < maxIterations; k++) {
        std::vector<double> y = matVecMultiply(matrix, x);
        
        double tmpEigenval = y[0] / x[0];
        
        normalize(y);
        
        if (fabs(tmpEigenval - eigenval) < eps) {
            return {tmpEigenval, y};
        }

        x = y;

        eigenval = tmpEigenval;
    }

    return {eigenval, x};
}

int main() {
    const std::vector<std::vector<double>> constA = {
        {5, 5, 3},
        {5, -4, 1},
        {3, 1, 2}
    };

    std::vector<std::vector<double>> matrix = constA;

    int n = matrix.size();
    std::vector<double> eigenvals(n);
    std::vector<std::vector<double>> eigenvecs;

    jacobiEigenvalMethod(matrix, eigenvals, eigenvecs);

    std::cout << "> Eigenvalues (Jacobi):\n";
    for (double elem : eigenvals) {
        std::cout << elem << ' ';
    }
    std::cout << "\n\n> Eigenvectors (Jacobi):\n";
    printMatrix(eigenvecs);

    auto [eigenval, eigenvector] = powerMethod(constA);

    std::cout << "\n> Dominant Eigenvalue (Power): " << eigenval << '\n';
    std::cout << "\n> Corresponding Eigenvector (Power): ";
    for (auto elem : eigenvector) {
        std::cout << elem << ' ';
    }
    std::cout << '\n';

    return 0;
}