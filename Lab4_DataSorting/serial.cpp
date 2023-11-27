#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <chrono>
#include <fstream>

using namespace std;

const double RandomDataMultiplier = 1000.0;

// Function for formatted data output
void PrintData(double *pData, int DataSize)
{
    for (int i = 0; i < DataSize; i++)
        printf("%7.4f ", pData[i]);
    printf("\n");
}
// Sorting by the standard library algorithm
void SerialStdSort(double *pData, int DataSize)
{
    sort(pData, pData + DataSize);
}

// Function for simple setting the initial data
void DummyDataInitialization(double *&pData, int &DataSize)
{
    for (int i = 0; i < DataSize; i++)
        pData[i] = DataSize - i;
}

void ReadDataInitialization(double *pData, int Size)
{
    std::ifstream fin("data" + std::to_string(Size) + ".txt");
    int i;
    for (i = 0; i < Size; i++)
    {
        fin >> pData[i];
    }
}
// Function for initializing the data by the random generator
void RandomDataInitialization(double *&pData, int &DataSize)
{
    srand((unsigned)time(0));
    for (int i = 0; i < DataSize; i++)
        pData[i] = double(rand()) / RAND_MAX * RandomDataMultiplier;
}

// Function for allocating the memory and setting the initial values
void ProcessInitialization(double *&pData, int &DataSize)
{
    // do
    // {
    //     printf("Enter the size of data to be sorted: ");
    //     scanf("%d", &DataSize);
    //     if (DataSize <= 0)
    //         printf("Data size should be greater than zero\n");
    // } while (DataSize <= 0);
    // printf("Sorting %d data items\n", DataSize);
    pData = new double[DataSize];

    // Simple setting the data
    // DummyDataInitialization(pData, DataSize);

    // Setting the data by the random generator
    // RandomDataInitialization(pData, DataSize);
    ReadDataInitialization(pData, DataSize);
}
// Function for computational process termination
void ProcessTermination(double *pData)
{
    delete[] pData;
}

// Function for the serial bubble sort algorithm
void SerialBubble(double *pData, int DataSize)
{
    double Tmp;
    for (int i = 1; i < DataSize; i++)
        for (int j = 0; j < DataSize - i; j++)
            if (pData[j] > pData[j + 1])
            {
                Tmp = pData[j];
                pData[j] = pData[j + 1];
                pData[j + 1] = Tmp;
            }
}

int main(int argc, char *argv[])
{
    double *pData = 0;
    int DataSize = atoi(argv[1]);
    // time_t start, finish;
    double duration = 0.0;
    printf("Serial bubble sort program\n");
    // Process initialization
    ProcessInitialization(pData, DataSize);
    // printf("Data before sorting\n");
    // PrintData(pData, DataSize);
    // Serial bubble sor
    auto start = std::chrono::system_clock::now();
    SerialBubble(pData, DataSize);
    auto end = std::chrono::system_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    printf("Time of serial execution = %ld microseconds\n\n", elapsed.count());
    // Process termination
    ProcessTermination(pData);
    return 0;
}