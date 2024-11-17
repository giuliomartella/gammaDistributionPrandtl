

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

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
double **buildA(const int N, std::vector<double> &theta, std::vector<double> &chord, std::vector<double> &clAlpha);
void freeA(double **A, int rows);
std::vector<double> buildAlphaG(const int N);
std::vector<double> buildAlpha0(const int N);
std::vector<double> buildCVector(const int N, std::vector<double> &clAlpha, std::vector<double> &alphaG, std::vector<double> &alpha0);
std::vector<double> linearSolverJacobi(const int N, std::vector<double> &cVector, double **A);


int main()
{
    const int N = 10;
    double AR = 10.0;
    // const b = 1.0; implicit

    std::vector<double> theta = buildTheta(N);
    std::vector<double> clAlpha = buildClAlpha(N);
    std::vector<double> chord = buildChord(N, AR);
    double **A = buildA(N, theta, chord, clAlpha);
    std::vector<double> alphaG = buildAlphaG(N);
    std::vector<double> alpha0 = buildAlpha0(N);
    std::vector<double> cVector = buildCVector(N, clAlpha, alphaG, alpha0);
    std::vector<double> bVector = linearSolver(N, cVector, A);



    freeA(A, N);
  

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

double **buildA(const int N, std::vector<double> &theta, std::vector<double> &chord, std::vector<double> &clAlpha)
{
    double **A;
    A = new double *[N];

    for (int i = 0; i < N; i++)
    {
        A[i] = new double[N];
    }

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            A[i][k] = -4.0 * sin((k + 1) * theta[i]) / chord[i] - clAlpha[i] * (k + 1) * sin((k + 1) * theta[i]) / sin(theta[i]);
            
        }
    }
    return A;
}

void freeA(double **A, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        delete[] A[i];
    }
    delete[] A;
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

std::vector<double> buildCVector(const int N, std::vector<double> &clAlpha, std::vector<double> &alphaG, std::vector<double> &alpha0)
{
    std::vector<double> cVector;
    for (int i = 0; i < N; i++)
    {
        cVector.push_back(clAlpha[i] * alphaG[i] - alpha0[i]);
    }
    return cVector;
}

std::vector<double> linearSolverJacobi(const int N, std::vector<double> &cVector, double **A)
{
    // Jacobi method
    std::vector<double> bVector(N, 0.0);
    std::vector<double> bVectorNext(N, 0.0);

    // std::cout << "Matrix A:" << std::endl;
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         std::cout << std::setw(10) << A[i][j] << " "; // Formattazione per una stampa più leggibile
    //     }
    //     std::cout << std::endl;
    // }
    //    std::cout << "Vector:" << std::endl;
    // for (double val : cVector)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;
   
    bool end = false;
    int k = 0;
    double sum;
    double delta;
    while (end == false)
    {std::cout << "Matrix A:" << std::endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            std::cout << std::setw(10) << A[i][j] << " "; // Formattazione per una stampa più leggibile
        }
        std::cout << std::endl;
    }
       std::cout << "Vector:" << std::endl;
    for (double val : cVector)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;
        for (int i = 0; i < N; i++)
        {
            sum = 0.0;
            for (int j = 0; j < N; j++)
            {
                if (i != j)
                {
                    sum += A[i][j] * bVector[j];
                }
            }

            bVectorNext[i] = (cVector[i] - sum) / A[i][i];
        }

        k++;
        delta = 0.0;
        for (int i = 0; i < N; i++) {
            delta += (bVectorNext[i] - bVector[i]) * (bVectorNext[i] - bVector[i]);
        }
        delta = std::sqrt(delta);
        end = (k > 1000 || delta < 1e-8);
        for (int i = 0; i < N; i++)
        {
            bVector[i] = bVectorNext[i];
        }
        
    }

    return bVector;
}
