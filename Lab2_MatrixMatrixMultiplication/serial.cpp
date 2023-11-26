#include <stdio.h>
#include <stdlib.h>
// #include <conio.h>
#include <time.h>
#include <fstream>
#include <chrono>

// Function for simple initialization of matrix elements
void DummyDataInitialization(double *pAMatrix, double *pBMatrix, int Size)
{
    int i, j; // Loop variables
    for (i = 0; i < Size; i++)
        for (j = 0; j < Size; j++)
        {
            pAMatrix[i * Size + j] = 1;
            pBMatrix[i * Size + j] = 1;
        }
}
// Function for initialization of matrix elements from files
void ReadDataInitialization(double *pAMatrix, double *pBMatrix, int Size)
{
    std::ifstream fin("AMatrix" + std::to_string(Size) + ".txt");
    int i, j;
    for (i = 0; i < Size; i++)
    {
        for (j = 0; j < Size; j++)
        {
            fin >> pAMatrix[i * Size + j];
        }
    }

    std::ifstream fin2("BMatrix" + std::to_string(Size) + ".txt");

    for (i = 0; i < Size; i++)
    {
        for (j = 0; j < Size; j++)
        {
            fin2 >> pBMatrix[i * Size + j];
        }
    }
}
// Function for random initialization of matrix elements
void RandomDataInitialization(double *pAMatrix, double *pBMatrix,
                              int Size)
{
    int i, j; // Loop variables
    srand(unsigned(clock()));
    for (i = 0; i < Size; i++)
        for (j = 0; j < Size; j++)
        {
            pAMatrix[i * Size + j] = rand() / double(1000);
            pBMatrix[i * Size + j] = rand() / double(1000);
        }
}
// Function for memory allocation and initialization of matrix elements
void ProcessInitialization(double *&pAMatrix, double *&pBMatrix,
                           double *&pCMatrix, int &Size)
{
    printf("\nMatrices' size = %d\n", Size);
    // Memory allocation
    pAMatrix = new double[Size * Size];
    pBMatrix = new double[Size * Size];
    pCMatrix = new double[Size * Size];
    // Initialization of matrix elements
    ReadDataInitialization(pAMatrix, pBMatrix, Size);
    for (int i = 0; i < Size * Size; i++)
    {
        pCMatrix[i] = 0;
    }
}
// Function for formatted matrix output
void PrintMatrix(double *pMatrix, int RowCount, int ColCount)
{
    int i, j; // Loop variables
    for (i = 0; i < RowCount; i++)
    {
        for (j = 0; j < ColCount; j++)
            printf("%7.4f ", pMatrix[i * RowCount + j]);
        printf("\n");
    }
}
// Function for matrix multiplication
void SerialResultCalculation(double *pAMatrix, double *pBMatrix,
                             double *pCMatrix, int Size)
{
    int i, j, k; // Loop variables
    for (i = 0; i < Size; i++)
    {
        for (j = 0; j < Size; j++)
            for (k = 0; k < Size; k++)
                pCMatrix[i * Size + j] += pAMatrix[i * Size + k] * pBMatrix[k * Size + j];
    }
}
// Function for computational process termination
void ProcessTermination(double *pAMatrix, double *pBMatrix,
                        double *pCMatrix)
{
    delete[] pAMatrix;
    delete[] pBMatrix;
    delete[] pCMatrix;
}
int main(int argc, char *argv[])
{
    double *pAMatrix; // First argument of matrix multiplication
    double *pBMatrix; // Second argument of matrix multiplication
    double *pCMatrix; // Result matrix
    int Size;         // Size of matrices
    Size = atoi(argv[1]);
    double duration;
    printf("Serial matrix multiplication program\n");
    // Memory allocation and initialization of matrix elements
    ProcessInitialization(pAMatrix, pBMatrix, pCMatrix, Size);
    // Matrix output
    // printf("Initial A Matrix \n");
    // PrintMatrix(pAMatrix, Size, Size);
    // printf("Initial B Matrix \n");
    // PrintMatrix(pBMatrix, Size, Size);
    // Matrix multiplication
    auto start = std::chrono::system_clock::now();
    SerialResultCalculation(pAMatrix, pBMatrix, pCMatrix, Size);
    
    auto end = std::chrono::system_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Printing the result matrix
    // printf("\n Result Matrix: \n");
    // PrintMatrix(pCMatrix, Size, Size);
    // Printing the time spent by matrix multiplication
    printf("Time of serial execution  = %ld microseconds\n\n", elapsed.count());
    // Computational process termination
    ProcessTermination(pAMatrix, pBMatrix, pCMatrix);
}

// #include <stdio.h>
// #include <stdlib.h>
// #include <fstream>
// #include <iostream>
// #include <chrono>

// void ReadDataInitialization(double *pAMatrix, double *pBMatrix, int Size)
// {
//     std::ifstream fin("AMatrix" + std::to_string(Size) + ".txt");
//     int i, j;
//     for (i = 0; i < Size; i++)
//     {
//         for (j = 0; j < Size; j++)
//         {
//             fin >> pAMatrix[i * Size + j];
//         }
//     }

//     std::ifstream fin2("BMatrix" + std::to_string(Size) + ".txt");

//     for (i = 0; i < Size; i++)
//     {
//         for (j = 0; j < Size; j++)
//         {
//             fin >> pBMatrix[i * Size + j];
//         }
//     }
// }
// // Function for random setting the matrix and vector elements
// void RandomDataInitialization(double *pMatrix, double *pVector, int Size)
// {
//     int i, j; // Loop variables
//     srand(unsigned(clock()));
//     for (i = 0; i < Size; i++)
//     {
//         pVector[i] = rand() / double(1000);
//         for (j = 0; j < Size; j++)
//             pMatrix[i * Size + j] = rand() / double(1000);
//     }
// }
// // Function for memory allocation and data initialization
// void ProcessInitialization(double *&pMatrix, double *&pVector,
//                            double *&pResult, int &Size)
// {
//     // Memory allocation
//     pMatrix = new double[Size * Size];
//     pVector = new double[Size];
//     pResult = new double[Size];
//     // Setting the values of the matrix and vector elements
//     ReadDataInitialization(pMatrix, pVector, Size);
// }
// // Function for formatted matrix output
// void PrintMatrix(double *pMatrix, int RowCount, int ColCount)
// {
//     int i, j; // Loop variables
//     for (i = 0; i < RowCount; i++)
//     {
//         for (j = 0; j < ColCount; j++)
//             printf("%7.4f ", pMatrix[i * RowCount + j]);
//         printf("\n");
//     }
// }
// // Function for formatted vector output
// void PrintVector(double *pVector, int Size)
// {
//     int i;
//     for (i = 0; i < Size; i++)
//         printf("%7.4f ", pVector[i]);
// }
// // Function for matrix-vector multiplication
// void ResultCalculation(double *pMatrix, double *pVector, double *pResult,
//                        int Size)
// {
//     int i, j; // Loop variables
//     for (i = 0; i < Size; i++)
//     {
//         pResult[i] = 0;
//         for (j = 0; j < Size; j++)
//             pResult[i] += pMatrix[i * Size + j] * pVector[j];
//     }
// }
// // Function for computational process termination
// void ProcessTermination(double *pMatrix, double *pVector, double *pResult)
// {
//     delete[] pMatrix;
//     delete[] pVector;
//     delete[] pResult;
// }
// int main(int argc, char *argv[])
// {
//     double *pMatrix; // First argument - initial matrix
//     double *pVector; // Second argument - initial vector
//     double *pResult; // Result vector for matrix-vector multiplication
//     int Size;        // Sizes of initial matrix and vector

//     double duration;
//     // printf("Serial matrix-vector multiplication program\n");
//     // Memory allocation and data initialization

//     Size = atoi(argv[1]);
//     // printf("Matrix size: %d\n", Size);

//     ProcessInitialization(pMatrix, pVector, pResult, Size);
//     // Matrix and vector output
//     // printf("Initial Matrix \n");
//     // PrintMatrix(pMatrix, Size, Size);
//     // printf("Initial Vector \n");
//     // PrintVector(pVector, Size);
//     // Matrix-vector multiplication
//     auto start = std::chrono::system_clock::now();

//     ResultCalculation(pMatrix, pVector, pResult, Size);

//     auto end = std::chrono::system_clock::now();
//     // Printing the result vector
//     // printf("\n Result Vector: \n");
//     // PrintVector(pResult, Size);
//     // Printing the time spent by matrix-vector multiplication

//     auto elapsed =
//         std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//         // std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     std::cout << "Time of serial execution = " << elapsed.count() << " microseconds" << std::endl << std::endl;
//     // Computational process termination
//     ProcessTermination(pMatrix, pVector, pResult);
//     // getch();
// }