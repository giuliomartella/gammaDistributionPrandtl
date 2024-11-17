#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

// functions list
// 1. theta vector
// 2. chord vector
// 3. A matrix
// 4. c vector
// 5. linear solver
// 6. coefficients
// 7. visualizations
std::vector<double> buildTheta(const int N);
std::vector<double> buildClAlpha(const int N);
std::vector<double> buildChord(const int N, double &AR);
Eigen::MatrixXd buildA(const int N, std::vector<double> &theta, std::vector<double> &chord, std::vector<double> &clAlpha);
std::vector<double> buildAlphaG(const int N);
std::vector<double> buildAlpha0(const int N);
Eigen::VectorXd buildCVector(const int N, std::vector<double> &clAlpha, std::vector<double> &alphaG, std::vector<double> &alpha0);
Eigen::VectorXd linearSolverEigen(const int N, Eigen::VectorXd& cVector, Eigen::MatrixXd &A);
double computeCl(double& AR, Eigen::VectorXd& bVector);
double computeCd(double& AR, Eigen::VectorXd& bVector, double cl);

int main()
{
    const int N = 10;
    double AR = 10.0;
    // const b = 1.0; implicit

    std::vector<double> theta = buildTheta(N);
    std::vector<double> clAlpha = buildClAlpha(N);
    std::vector<double> chord = buildChord(N, AR);
    Eigen::MatrixXd A = buildA(N, theta, chord, clAlpha);
    std::vector<double> alphaG = buildAlphaG(N);
    std::vector<double> alpha0 = buildAlpha0(N);
    Eigen::VectorXd cVector = buildCVector(N, clAlpha, alphaG, alpha0);
    Eigen::VectorXd bVector = linearSolverEigen(N, cVector, A);
    double cl = computeCl(AR, bVector);
    double cd = computeCd(AR, bVector, cl);
    double efficiency = cl / cd;


    return 0;
}

std::vector<double> buildTheta(const int N)
{
    std::vector<double> theta;
    for (int i = 0; i < N; i++)
    {
        theta.push_back(i * M_PI / N + M_PI / (2 * N));
    }
    return theta;
}

std::vector<double> buildClAlpha(const int N)
{
    std::vector<double> clAlpha;
    for (int i = 0; i < N; i++)
    {
        clAlpha.push_back(2.0 * M_PI);
    }
    return clAlpha;
}

std::vector<double> buildChord(const int N, double &AR)
{
    std::vector<double> chord;
    for (int i = 0; i < N; i++)
    {
        chord.push_back(1.0 / AR);
    }
    return chord;
}

Eigen::MatrixXd buildA(const int N, std::vector<double> &theta, std::vector<double> &chord, std::vector<double> &clAlpha)
{
    Eigen::MatrixXd A(N, N);
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            A(i,k) = -4.0 * sin((k + 1) * theta[i]) / chord[i] - clAlpha[i] * (k + 1) * sin((k + 1) * theta[i]) / sin(theta[i]);
        }
    }
    return A;
}

std::vector<double> buildAlphaG(const int N)
{
    std::vector<double> alphaG;
    for (int i = 0; i < N; i++)
    {
        alphaG.push_back(M_PI / 180.0);
    }
    return alphaG;
}

std::vector<double> buildAlpha0(const int N)
{
    std::vector<double> alpha0;
    for (int i = 0; i < N; i++)
    {
        alpha0.push_back(0.0);
    }
    return alpha0;
}

Eigen::VectorXd buildCVector(const int N, std::vector<double> &clAlpha, std::vector<double> &alphaG, std::vector<double> &alpha0)
{
    Eigen::VectorXd cVector(N);
    for (int i = 0; i < N; i++)
    {
        cVector(i) = (clAlpha[i] * alphaG[i] - alpha0[i]);
    }
    return cVector;
}

Eigen::VectorXd linearSolverEigen(const int N, Eigen::VectorXd& cVector, Eigen::MatrixXd& A)
{

    //  std::cout << "Matrix A:" << std::endl;
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         std::cout << std::setw(10) << A[i][j] << " "; // Formattazione per una stampa più leggibile
    //     }
    //     std::cout << std::endl;
    // }

    // Perform LU decomposition using PartialPivLU
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_decomp(A);

    // Solve the system using the LU decomposition
    Eigen::VectorXd bVector = lu_decomp.solve(cVector);


    std::cout << "Vector:" << std::endl;
    for (double val : bVector)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return bVector;
}


double computeCl(double& AR, Eigen::VectorXd& bVector)
{

    double cl = - M_PI * AR * bVector(0);
    return cl;
}

double computeCd(double& AR, Eigen::VectorXd& bVector, double cl)
{

    double delta = 0.0;
    for (int i = 1; i < bVector.size(); i++)
    {
        delta += (bVector(i) * bVector(0)) * (bVector(i) * bVector(0));
    }
    double cd = cl * cl / (M_PI * AR) * (1 + delta);
    return cd;
}
