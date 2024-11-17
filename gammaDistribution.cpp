#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
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
std::vector<double> buildChord(const int N, double &AR, double &TR, std::vector<double> &theta);
Eigen::MatrixXd buildA(const int N, std::vector<double> &theta, std::vector<double> &chord, std::vector<double> &clAlpha);
std::vector<double> buildAlphaG(const int N, double &twistPrime, std::vector<double> &theta);
std::vector<double> buildAlpha0(const int N);
Eigen::VectorXd buildCVector(const int N, std::vector<double> &clAlpha, std::vector<double> &alphaG, std::vector<double> &alpha0);
Eigen::VectorXd linearSolverEigen(const int N, Eigen::VectorXd &cVector, Eigen::MatrixXd &A);
double computeCl(double &AR, Eigen::VectorXd &bVector);
double computeCd(double &AR, Eigen::VectorXd &bVector, double cl);
void computeResults(double& TR, double& AR, double& twistPrime, double &cl, double &cd, double &efficiency);

int main()
{
   // Define the ranges for the degrees of freedom (TR, AR, twistPrime)
    std::vector<double> TR_values = {0.3, 0.5, 0.7};
    std::vector<double> AR_values = {8.0, 10.0, 12.0};
    std::vector<double> twistPrime_values = {0.0, M_PI/90.0, M_PI/45.0};

    // Open the output file to save the results
    std::ofstream outputFile("output.txt");

    // Write table header to file
    outputFile << std::setw(10) << "TR" << std::setw(10) << "AR" << std::setw(15) << "Twist Prime" 
               << std::setw(20) << "Lift Coefficient (Cl)" << std::setw(20) << "Induced Drag (Cd)" 
               << std::setw(20) << "Efficiency (Cl/Cd)" << std::endl;

    // Loop over all combinations of TR, AR, and twistPrime values
    for (double TR : TR_values)
    {
        for (double AR : AR_values)
        {
            for (double twistPrime : twistPrime_values)
            {
                // Declare variables for Cl, Cd, and Efficiency
                double cl, cd, efficiency;

                // Call computeResults to calculate Cl, Cd, and Efficiency for this set of parameters
                computeResults(TR, AR, twistPrime, cl, cd, efficiency);

                // Write the results into the output file in a well-formatted table
                outputFile << std::setw(10) << TR << std::setw(10) << AR << std::setw(15) << twistPrime
                           << std::setw(20) << cl << std::setw(20) << cd << std::setw(20) << efficiency << std::endl;
            }
        }
    }

    // Close the output file after writing the results
    outputFile.close();

    std::cout << "Results have been saved to output.txt" << std::endl;

    

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
   

std::vector<double> buildChord(const int N, double &AR, double &TR, std::vector<double> &theta)
{
    std::vector<double> chord;
    double z;
    for (int i = 0; i < N; i++)
    {
        z = 0.5 * cos(theta[i]);
        chord.push_back(1.0 / AR * TR * abs(z));
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
            A(i, k) = -4.0 * sin((k + 1) * theta[i]) / chord[i] - clAlpha[i] * (k + 1) * sin((k + 1) * theta[i]) / sin(theta[i]);
        }
    }
    return A;
}

std::vector<double> buildAlphaG(const int N, double &twistPrime, std::vector<double> &theta)
{
    std::vector<double> alphaG;
    double z;
    for (int i = 0; i < N; i++)
    {
        z = 0.5 * cos(theta[i]);
        alphaG.push_back(twistPrime * abs(z));
    }
    return alphaG;
}


std::vector<double> buildAlpha0(const int N)
{
    std::vector<double> alpha0;
    for (int i = 0; i < N; i++)
    {
        alpha0.push_back(-2.0);
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

Eigen::VectorXd linearSolverEigen(const int N, Eigen::VectorXd &cVector, Eigen::MatrixXd &A)
{

    // Perform LU decomposition using PartialPivLU
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_decomp(A);

    // Solve the system using the LU decomposition
    Eigen::VectorXd bVector = lu_decomp.solve(cVector);

    // std::cout << "Vector:" << std::endl;
    // for (double val : bVector)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    return bVector;
}

double computeCl(double &AR, Eigen::VectorXd &bVector)
{

    double cl = -M_PI * AR * bVector(0);
    return cl;
}

double computeCd(double &AR, Eigen::VectorXd &bVector, double cl)
{

    double delta = 0.0;
    for (int i = 1; i < bVector.size(); i++)
    {
        delta += (bVector(i) * bVector(0)) * (bVector(i) * bVector(0));
    }
    double cd = cl * cl / (M_PI * AR) * (1 + delta);
    return cd;
}

void computeResults(double& TR, double& AR, double& twistPrime, double &cl, double &cd, double &efficiency)
{
    const int N = 10;
    // const b = 1.0; implicit

    std::vector<double> theta = buildTheta(N);
    std::vector<double> clAlpha = buildClAlpha(N);
    std::vector<double> chord = buildChord(N, AR, TR, theta);
    Eigen::MatrixXd A = buildA(N, theta, chord, clAlpha);
    std::vector<double> alphaG = buildAlphaG(N, twistPrime, theta);
    std::vector<double> alpha0 = buildAlpha0(N);
    Eigen::VectorXd cVector = buildCVector(N, clAlpha, alphaG, alpha0);
    Eigen::VectorXd bVector = linearSolverEigen(N, cVector, A);
    cl = computeCl(AR, bVector);
    cd = computeCd(AR, bVector, cl);
    efficiency = cl / cd;

    

    // compute cl, cd ed efficiency
    std::cout << "Lift coefficient (Cl): " << cl << std::endl;
    std::cout << "Induced Drag Coefficient (Cd): " << cd << std::endl;
    std::cout << "Efficiency (Cl/Cd): " << efficiency << std::endl;

}
