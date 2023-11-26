#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <mpi.h>
#include <chrono>

int ProcNum = 0;  // Number of available processes
int ProcRank = 0; // Rank of current process

void ReadDataInitialization(double *pMatrix, double *pVector, int Size)
{
    std::ifstream fin("matrix" + std::to_string(Size) + ".txt");
    int i, j;
    for (i = 0; i < Size; i++)
    {
        for (j = 0; j < Size; j++)
        {
            fin >> pMatrix[i * Size + j];
        }
    }

    std::ifstream fin2("vector" + std::to_string(Size) + ".txt");
    for (i = 0; i < Size; i++)
    {
        fin2 >> pVector[i];
    }
}
// Function for memory allocation and data initialization
void ProcessInitialization(double *&pMatrix, double *&pVector,
                           double *&pResult, double *&pProcRows, double *&pProcResult,
                           int &Size, int &RowNum)
{
    int RestRows; // Number of rows, that haven�t been distributed yet
    int i;        // Loop variable
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    RestRows = Size;
    for (i = 0; i < ProcRank; i++)
        RestRows = RestRows - RestRows / (ProcNum - i);
    RowNum = RestRows / (ProcNum - ProcRank);
    pVector = new double[Size];
    pResult = new double[Size];
    pProcRows = new double[RowNum * Size];
    pProcResult = new double[RowNum];
    if (ProcRank == 0)
    {
        pMatrix = new double[Size * Size];
        ReadDataInitialization(pMatrix, pVector, Size);
    }
}
// Data distribution among the processes
void DataDistribution(double *pMatrix, double *pProcRows, double *pVector,
                      int Size, int RowNum, int64_t &scatterTime)
{
    int *pSendNum;       // the number of elements sent to the process
    int *pSendInd;       // the index of the first data element sent to the process
    int RestRows = Size; // Number of rows, that haven�t been distributed yet
    MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Alloc memory for temporary objects
    pSendInd = new int[ProcNum];
    pSendNum = new int[ProcNum];
    // Define the disposition of the matrix rows for current process
    RowNum = (Size / ProcNum);
    pSendNum[0] = RowNum * Size;
    pSendInd[0] = 0;
    for (int i = 1; i < ProcNum; i++)
    {
        RestRows -= RowNum;
        RowNum = RestRows / (ProcNum - i);
        pSendNum[i] = RowNum * Size;
        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
    }
    std::chrono::time_point<std::chrono::system_clock> start;
    if (ProcRank == 0)
        start = std::chrono::system_clock::now();
    // Scatter the rows
    // MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows, pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int sendCount = pSendNum[ProcRank]; // Number of elements to send to each process

    MPI_Scatter(pMatrix, sendCount, MPI_DOUBLE, pProcRows, sendCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // MPI_Scatter(pMatrix, Size * Size / ProcNum, MPI_DOUBLE, pProcRows,
    //             Size * Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (ProcRank == 0)
    {
        auto end = std::chrono::system_clock::now();
        auto elapsed =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        // std::cout << "scatterv time: " << elapsed.count() << " microseconds" << std::endl;
        scatterTime = elapsed.count();
    }
    // Free the memory
    delete[] pSendNum;
    delete[] pSendInd;
}
// Function for gathering the result vector
void ResultReplication(double *pProcResult, double *pResult, int Size,
                       int RowNum)
{
    int i;               // Loop variable
    int *pReceiveNum;    // Number of elements, that current process sends
    int *pReceiveInd;    /* Index of the first element from current process
       in result vector */
    int RestRows = Size; // Number of rows, that haven�t been distributed yet
    // Alloc memory for temporary objects
    pReceiveNum = new int[ProcNum];
    pReceiveInd = new int[ProcNum];
    // Define the disposition of the result vector block of current processor
    pReceiveInd[0] = 0;
    pReceiveNum[0] = Size / ProcNum;
    for (i = 1; i < ProcNum; i++)
    {
        RestRows -= pReceiveNum[i - 1];
        pReceiveNum[i] = RestRows / (ProcNum - i);
        pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
    }
    // Gather the whole result vector on every processor
    MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_DOUBLE, pResult,
                   pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
    
    // Free the memory
    delete[] pReceiveNum;
    delete[] pReceiveInd;
}
// Function for sequential matrix-vector multiplication
void SerialResultCalculation(double *pMatrix, double *pVector, double *pResult, int Size)
{
    int i, j; // Loop variables
    for (i = 0; i < Size; i++)
    {
        pResult[i] = 0;
        for (j = 0; j < Size; j++)
            pResult[i] += pMatrix[i * Size + j] * pVector[j];
    }
}
// Function for calculating partial matrix-vector multiplication
void ParallelResultCalculation(double *pProcRows, double *pVector, double *pProcResult, int Size, int RowNum)
{
    int i, j; // Loop variables
    for (i = 0; i < RowNum; i++)
    {
        pProcResult[i] = 0;
        for (j = 0; j < Size; j++)
            pProcResult[i] += pProcRows[i * Size + j] * pVector[j];
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
// Function for formatted vector output
void PrintVector(double *pVector, int Size)
{
    int i;
    for (i = 0; i < Size; i++)
        printf("%7.4f ", pVector[i]);
}
void TestDistribution(double *pMatrix, double *pVector, double *pProcRows,
                      int Size, int RowNum)
{
    if (ProcRank == 0)
    {
        printf("Initial Matrix: \n");
        PrintMatrix(pMatrix, Size, Size);
        printf("Initial Vector: \n");
        PrintVector(pVector, Size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < ProcNum; i++)
    {
        if (ProcRank == i)
        {
            printf("\nProcRank = %d \n", ProcRank);
            printf(" Matrix Stripe:\n");
            PrintMatrix(pProcRows, RowNum, Size);
            printf(" Vector: \n");
            PrintVector(pVector, Size);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
void TestPartialResults(double *pProcResult, int RowNum)
{
    int i; // Loop variables
    for (i = 0; i < ProcNum; i++)
    {
        if (ProcRank == i)
        {
            printf("\nProcRank = %d \n Part of result vector: \n", ProcRank);
            PrintVector(pProcResult, RowNum);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
void TestResult(double *pMatrix, double *pVector, double *pResult,
                int Size)
{
    // Buffer for storing the result of serial matrix-vector multiplication
    double *pSerialResult;
    // Flag, that shows wheather the vectors are identical or not
    int equal = 0;
    int i; // Loop variable
    if (ProcRank == 0)
    {
        pSerialResult = new double[Size];
        SerialResultCalculation(pMatrix, pVector, pSerialResult, Size);
        for (i = 0; i < Size; i++)
        {
            if (pResult[i] != pSerialResult[i])
                equal = 1;
        }
        if (equal == 1)
            printf("The results of serial and parallel algorithms "
                   "are NOT identical. Check your code.");
        else
            printf("The results of serial and parallel algorithms "
                   "are identical.");
    }
}
// Function for computational process termination
void ProcessTermination(double *pMatrix, double *pVector, double *pResult,
                        double *pProcRows, double *pProcResult)
{
    if (ProcRank == 0)
        delete[] pMatrix;
    delete[] pVector;
    delete[] pResult;
    delete[] pProcRows;
    delete[] pProcResult;
}
int main(int argc, char *argv[])
{
    double *pMatrix; // The first argument - initial matrix
    double *pVector; // The second argument - initial vector
    double *pResult; // Result vector for matrix-vector multiplication
    int Size;        // Sizes of initial matrix and vector
    double *pProcRows;
    double *pProcResult;
    int RowNum;
    double Start, Finish, DurationParallel, DurationSerial;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0)
    {
        std::cout << "Number of processes = " << ProcNum << std::endl;
        // Size = ;
        Size = atoi(argv[1]);
        // printf("Matrix size: %d\n", Size);
    }

    ProcessInitialization(pMatrix, pVector, pResult, pProcRows, pProcResult, Size, RowNum);

    auto start = std::chrono::system_clock::now();
    int64_t scatterTime;

    DataDistribution(pMatrix, pProcRows, pVector, Size, RowNum, scatterTime);
    ParallelResultCalculation(pProcRows, pVector, pProcResult, Size, RowNum);
    ResultReplication(pProcResult, pResult, Size, RowNum);

    auto end = std::chrono::system_clock::now();
    DurationParallel = Finish - Start;

    // Start = MPI_Wtime();
    TestResult(pMatrix, pVector, pResult, Size);
    auto elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // Finish = MPI_Wtime();

    // DurationSerial = Finish - Start;
    if (ProcRank == 0)
    {
        printf("\n");
        printf("Time of parallel execution with scatterv issue = %ld microseconds\n\n", elapsed.count());
        printf("Time of parallel clean execution = %ld microseconds\n\n", elapsed.count() - scatterTime);

        // printf("Time of serial execution = %f\n", DurationSerial);
    }
    ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcResult);
    MPI_Finalize();
}