#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <chrono>

const double InfinitiesPercent = 50.0;

// Function for random definition of matrix and vector elements
// void RandomDataInitialization(double *pData, int Size)
// {
//     int i; // Loop variable
//     srand(unsigned(clock()));
//     for (i = 0; i < Size; i++)
//     {
//         // pVector[i] = float(rand() % 10);
//         pData[i] = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / 10));
//     }
// }

void RandomDataInitialization(int *pMatrix, int Size)
{
    srand((unsigned)time(0));
    for (int i = 0; i < Size; i++)
        for (int j = 0; j < Size; j++)
            if (i != j)
            {
                if ((rand() % 100) < InfinitiesPercent)
                    pMatrix[i * Size + j] = -1;
                else
                    pMatrix[i * Size + j] = (rand() + 1) % 1000;
            }
            else
                pMatrix[i * Size + j] = 0;
}


// Function for formatted matrix output
void PrintMatrix(int *pMatrix, int RowCount, int ColCount)
{
    int i, j; // Loop variables
    for (i = 0; i < RowCount; i++)
    {
        for (j = 0; j < ColCount; j++)
            printf("%d ", pMatrix[i * ColCount + j]);
        printf("\n");
    }
}

void SaveMatrix(int *pData, int RowCount, int ColCount, std::string filename)
{
    std::ofstream fout(filename + std::to_string(RowCount) + ".txt");
    for (int i = 0; i < RowCount; i++)
        for (int j = 0; j < ColCount; j++)
            fout << pData[i * ColCount + j] << " " ;
    fout.close();
}

int main(int argc, char *argv[])
{

    int Size;
    Size = atoi(argv[1]);

    printf("%d\n", Size);

    int *pData = new int[Size * Size];

    auto start = std::chrono::system_clock::now();

    RandomDataInitialization(pData, Size);
    // PrintMatrix(pData, Size, Size);
    SaveMatrix(pData, Size, Size, "graph");

    auto end = std::chrono::system_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << elapsed.count() << " microseconds" << std::endl;
    printf("\n");

    return 0;
}