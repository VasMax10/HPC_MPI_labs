#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <chrono>

// Function for random definition of matrix and vector elements
void RandomDataInitialization(double *pData, int Size)
{
    int i; // Loop variable
    srand(unsigned(clock()));
    for (i = 0; i < Size; i++)
    {
        // pVector[i] = float(rand() % 10);
        pData[i] = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / 10));
    }
}

// Function for formatted matrix output
void PrintMatrix(double *pMatrix, int RowCount, int ColCount)
{
    int i; // Loop variables
    for (i = 0; i < RowCount; i++)
    {
        printf("%7.4f ", pMatrix[i]);
    }
}

void SaveMatrix(double *pData, int RowCount, int ColCount, std::string filename)
{
    int i;
    std::ofstream fout(filename + std::to_string(ColCount) + ".txt");
    for (i = 0; i < RowCount; i++)
    {
        fout << std::fixed << std::setprecision(4) << pData[i] << " ";
    }
}

int main(int argc, char *argv[])
{

    int Size;
    Size = atoi(argv[1]);
    printf("%d", Size);
    double *pData = new double[Size];

    auto start = std::chrono::system_clock::now();

    RandomDataInitialization(pData, Size);
    SaveMatrix(pData, Size, Size, "data");

    auto end = std::chrono::system_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << elapsed.count() << " microseconds" << std::endl;
    printf("\n");

    return 0;
}