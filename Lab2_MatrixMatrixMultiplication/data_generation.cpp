#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <chrono>

// Function for random definition of matrix and vector elements
void RandomDataInitialization(double *pAMatrix, double *pBMatrix, int Size)
{
    int i, j; // Loop variables
    srand(unsigned(clock()));
    for (i = 0; i < Size; i++)
        for (j = 0; j < Size; j++){
            pAMatrix[i * Size + j] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/10));
            pBMatrix[i * Size + j] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/10));
            // float(rand() % 10);
        }
}

// Function for formatted matrix output
void PrintMatrix(double *pMatrix, int RowCount, int ColCount)
{
    int i, j; // Loop variables
    for (i = 0; i < RowCount; i++)
    {
        for (j = 0; j < ColCount; j++)
            printf("%7.4f ", pMatrix[i * ColCount + j]);
        printf("\n");
    }
}

void SaveMatrix(double *pMatrix, int RowCount, int ColCount, std::string filename)
{
    int i, j;
    std::ofstream fout(filename + std::to_string(ColCount) + ".txt");
    for (i = 0; i < RowCount; i++)
    {
        for (j = 0; j < ColCount; j++)
        {
            fout << std::fixed << std::setprecision(4) << pMatrix[i * ColCount + j] << " ";
        }
    }
}

int main(int argc, char *argv[])
{

    int Size = atoi(argv[1]);
    double *pAMatrix = new double[Size * Size];
    double *pBMatrix = new double[Size * Size];
    
    auto start = std::chrono::system_clock::now();
    
    RandomDataInitialization(pAMatrix, pBMatrix, Size);
    SaveMatrix(pAMatrix, Size, Size, "AMatrix");
    SaveMatrix(pBMatrix, Size, Size, "BMatrix");
    
    auto end = std::chrono::system_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << elapsed.count() << " microseconds" << std::endl;
    printf("\n");

    return 0;
}

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
            fin >> pBMatrix[i * Size + j];
        }
    }
}