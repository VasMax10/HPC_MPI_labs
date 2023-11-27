#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <chrono>
#include <fstream>

using namespace std;

const double InfinitiesPercent = 50.0;
const double RandomDataMultiplier = 10;

// Function for simple setting the initial data
void DummyDataInitialization(int *pMatrix, int Size)
{
    for (int i = 0; i < Size; i++)
        for (int j = i; j < Size; j++)
        {
            if (i == j)
                pMatrix[i * Size + j] = 0;
            else if (i == 0)
                pMatrix[i * Size + j] = j;
            else
                pMatrix[i * Size + j] = -1;
            pMatrix[j * Size + i] = pMatrix[i * Size + j];
        }
}
// Function for initializing the data by the random generator
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
                    pMatrix[i * Size + j] = rand() + 1;
            }
            else
                pMatrix[i * Size + j] = 0;
}

void ReadDataInitialization(int *pData, int Size)
{
    std::ifstream fin("graph" + std::to_string(Size) + ".txt");
    int i;
    for (i = 0; i < Size; i++)
    {
        for (int j = 0; j < Size; j++)
            fin >> pData[i * Size + j];
    }
    fin.close();
}

int Min(int A, int B)
{
    int Result = (A < B) ? A : B;
    if ((A < 0) && (B >= 0))
        Result = B;
    if ((B < 0) && (A >= 0))
        Result = A;
    if ((A < 0) && (B < 0))
        Result = -1;
    return Result;
}

// Function for formatted matrix output
void PrintMatrix(int *pMatrix, int RowCount, int ColCount)
{
    for (int i = 0; i < RowCount; i++)
    {
        for (int j = 0; j < ColCount; j++)
            printf("%7d", pMatrix[i * ColCount + j]);
        printf("\n");
    }
}

// Function for allocating the memory and setting the initial values
void ProcessInitialization(int *&pMatrix, int &Size)
{
    printf("Graph with %d vertices\n", Size);
    // Allocate memory for the adjacency matrix
    pMatrix = new int[Size * Size];
    // Data initalization
    // DummyDataInitialization(pMatrix, Size);
    // RandomDataInitialization(pMatrix, Size);
    ReadDataInitialization(pMatrix, Size);
}
// Function for computational process termination
void ProcessTermination(int *pMatrix)
{
    delete[] pMatrix;
}

// Function for the serial Floyd algorithm
void SerialFloyd(int *pMatrix, int Size)
{
    int t1, t2;
    for (int k = 0; k < Size; k++)
        for (int i = 0; i < Size; i++)
            for (int j = 0; j < Size; j++)
                if ((pMatrix[i * Size + k] != -1) &&
                    (pMatrix[k * Size + j] != -1))
                {
                    t1 = pMatrix[i * Size + j];
                    t2 = pMatrix[i * Size + k] + pMatrix[k * Size + j];
                    pMatrix[i * Size + j] = Min(t1, t2);
                }
}

int main(int argc, char *argv[])
{
    int *pMatrix; // Adjacency matrix
    int Size;     // Size of adjacency matrix

    Size = atoi(argv[1]);

    printf("Serial Floyd algorithm\n");
    // Process initialization
    ProcessInitialization(pMatrix, Size);
    // printf("The matrix before Floyd algorithm\n");
    // PrintMatrix(pMatrix, Size, Size);
    auto start = std::chrono::system_clock::now();
    // Parallel Floyd algorithm
    SerialFloyd(pMatrix, Size);
    auto end = std::chrono::system_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    printf("Time of serial execution = %ld microseconds\n\n", elapsed.count());
    // printf("The matrix after Floyd algorithm\n");
    // PrintMatrix(pMatrix, Size, Size);
    // duration = (finish - start) / double(CLOCKS_PER_SEC);
    // printf("Time of execution: %f\n", duration);
    // Process termination
    ProcessTermination(pMatrix);
    return 0;
}