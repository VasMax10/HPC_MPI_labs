#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <chrono>
#include <fstream>

int ProcNum = 0;   // Number of available processes
int ProcRank = 0;  // Rank of current process
int GridSize;      // Size of virtual processor grid
int GridCoords[2]; // Coordinates of current processor in grid
MPI_Comm GridComm; // Grid communicator
MPI_Comm ColComm;  // Column communicator
MPI_Comm RowComm;  // Row communicator
/// Function for simple initialization of matrix elements
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
// Function for block multiplication
void BlockMultiplication(double *pAblock, double *pBblock,
                         double *pCblock, int Size)
{
    SerialResultCalculation(pAblock, pBblock, pCblock, Size);
}
// Function for creating the two-dimensional grid communicator
// and communicators for each row and each column of the grid
void CreateGridCommunicators()
{
    int DimSize[2];  // Number of processes in each dimension of the grid
    int Periodic[2]; // =1, if the grid dimension should be periodic
    int Subdims[2];  // =1, if the grid dimension should be fixed
    DimSize[0] = GridSize;
    DimSize[1] = GridSize;
    Periodic[0] = 0;
    Periodic[1] = 0;
    // Creation of the Cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);
    // Determination of the cartesian coordinates for every process
    MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);
    // Creating communicators for rows
    Subdims[0] = 0; // Dimensionality fixing
    Subdims[1] = 1; // The presence of the given dimension in the subgrid
    MPI_Cart_sub(GridComm, Subdims, &RowComm);
    // Creating communicators for columns
    Subdims[0] = 1;
    Subdims[1] = 0;
    MPI_Cart_sub(GridComm, Subdims, &ColComm);
}
// Function for memory allocation and data initialization
void ProcessInitialization(double *&pAMatrix, double *&pBMatrix,
                           double *&pCMatrix, double *&pAblock, double *&pBblock, double *&pCblock,
                           double *&pTemporaryAblock, int &Size, int &BlockSize)
{
    if (ProcRank == 0)
    {
        if (Size % GridSize != 0)
        {
            printf("Size of matrices must be divisible by the grid size!\n");
        }
    }
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    BlockSize = Size / GridSize;
    pAblock = new double[BlockSize * BlockSize];
    pBblock = new double[BlockSize * BlockSize];
    pCblock = new double[BlockSize * BlockSize];
    pTemporaryAblock = new double[BlockSize * BlockSize];
    for (int i = 0; i < BlockSize * BlockSize; i++)
    {
        pCblock[i] = 0;
    }
    if (ProcRank == 0)
    {
        pAMatrix = new double[Size * Size];
        pBMatrix = new double[Size * Size];
        pCMatrix = new double[Size * Size];
        ReadDataInitialization(pAMatrix, pBMatrix, Size);
        // RandomDataInitialization(pAMatrix, pBMatrix, Size);
    }
}
// Function for checkerboard matrix decomposition
void CheckerboardMatrixScatter(double *pMatrix, double *pMatrixBlock,
                               int Size, int BlockSize)
{
    double *MatrixRow = new double[BlockSize * Size];
    if (GridCoords[1] == 0)
    {
        MPI_Scatter(pMatrix, BlockSize * Size, MPI_DOUBLE, MatrixRow,
                    BlockSize * Size, MPI_DOUBLE, 0, ColComm);
    }
    for (int i = 0; i < BlockSize; i++)
    {
        MPI_Scatter(&MatrixRow[i * Size], BlockSize, MPI_DOUBLE,
                    &(pMatrixBlock[i * BlockSize]), BlockSize, MPI_DOUBLE, 0, RowComm);
    }
    delete[] MatrixRow;
}
// Data distribution among the processes
void DataDistribution(double *pAMatrix, double *pBMatrix, double *pMatrixAblock, double *pBblock, int Size, int BlockSize, int64_t &scatterTime)
{
    auto start = std::chrono::system_clock::now();
    // Scatter the matrix among the processes of the first grid column
    CheckerboardMatrixScatter(pAMatrix, pMatrixAblock, Size, BlockSize);
    CheckerboardMatrixScatter(pBMatrix, pBblock, Size, BlockSize);
    auto end = std::chrono::system_clock::now();
    scatterTime =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}
// Function for gathering the result matrix
void ResultCollection(double *pCMatrix, double *pCblock, int Size,
                      int BlockSize)
{
    double *pResultRow = new double[Size * BlockSize];
    for (int i = 0; i < BlockSize; i++)
    {
        MPI_Gather(&pCblock[i * BlockSize], BlockSize, MPI_DOUBLE,
                   &pResultRow[i * Size], BlockSize, MPI_DOUBLE, 0, RowComm);
    }
    if (GridCoords[1] == 0)
    {
        MPI_Gather(pResultRow, BlockSize * Size, MPI_DOUBLE, pCMatrix,
                   BlockSize * Size, MPI_DOUBLE, 0, ColComm);
    }
    delete[] pResultRow;
}
// Broadcasting blocks of the matrix A to process grid rows
void ABlockCommunication(int iter, double *pAblock, double *pMatrixAblock,
                         int BlockSize)
{
    // Defining the leading process of the process grid row
    int Pivot = (GridCoords[0] + iter) % GridSize;
    // Copying the transmitted block in a separate memory buffer
    if (GridCoords[1] == Pivot)
    {
        for (int i = 0; i < BlockSize * BlockSize; i++)
            pAblock[i] = pMatrixAblock[i];
    }
    // Block broadcasting
    MPI_Bcast(pAblock, BlockSize * BlockSize, MPI_DOUBLE, Pivot, RowComm);
}
// Function for cyclic shifting the blocks of the matrix B
void BblockCommunication(double *pBblock, int BlockSize)
{
    MPI_Status Status;
    int NextProc = GridCoords[0] + 1;
    if (GridCoords[0] == GridSize - 1)
        NextProc = 0;
    int PrevProc = GridCoords[0] - 1;
    if (GridCoords[0] == 0)
        PrevProc = GridSize - 1;
    MPI_Sendrecv_replace(pBblock, BlockSize * BlockSize, MPI_DOUBLE,
                         NextProc, 0, PrevProc, 0, ColComm, &Status);
}
// Function for parallel execution of the Fox method
void ParallelResultCalculation(double *pAblock, double *pMatrixAblock,
                               double *pBblock, double *pCblock, int BlockSize)
{
    for (int iter = 0; iter < GridSize; iter++)
    {
        // Sending blocks of matrix A to the process grid rows
        ABlockCommunication(iter, pAblock, pMatrixAblock, BlockSize);
        // Block multiplication
        BlockMultiplication(pAblock, pBblock, pCblock, BlockSize);
        // Cyclic shift of blocks of matrix B in process grid columns
        BblockCommunication(pBblock, BlockSize);
    }
}
// Test printing of the matrix block
void TestBlocks(double *pBlock, int BlockSize, char str[])
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcRank == 0)
    {
        printf("%s \n", str);
    }
    for (int i = 0; i < ProcNum; i++)
    {
        if (ProcRank == i)
        {
            printf("ProcRank = %d \n", ProcRank);
            PrintMatrix(pBlock, BlockSize, BlockSize);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
// Function for testing the matrix multiplication result
void TestResult(double *pAMatrix, double *pBMatrix, double *pCMatrix,
                int Size)
{
    double *pSerialResult;   // Result matrix of serial multiplication
    double Accuracy = 1.e-6; // Comparison accuracy
    int equal = 0;           // =1, if the matrices are not equal
    int i;                   // Loop variable
    if (ProcRank == 0)
    {
        pSerialResult = new double[Size * Size];
        for (i = 0; i < Size * Size; i++)
        {
            pSerialResult[i] = 0;
        }
        BlockMultiplication(pAMatrix, pBMatrix, pSerialResult, Size);
        for (i = 0; i < Size * Size; i++)
        {
            if (fabs(pSerialResult[i] - pCMatrix[i]) >= Accuracy)
                equal = 1;
        }
        if (equal == 1)
            printf("The results of serial and parallel algorithms are NOT"
                   "identical. Check your code.");
        else
            printf("The results of serial and parallel algorithms are "
                   "identical. ");
    }
}
// Function for computational process termination
void ProcessTermination(double *pAMatrix, double *pBMatrix,
                        double *pCMatrix, double *pAblock, double *pBblock, double *pCblock,
                        double *pMatrixAblock)
{
    if (ProcRank == 0)
    {
        delete[] pAMatrix;
        delete[] pBMatrix;
        delete[] pCMatrix;
    }
    delete[] pAblock;
    delete[] pBblock;
    delete[] pCblock;
    delete[] pMatrixAblock;
}
int main(int argc, char *argv[])
{
    double *pAMatrix; // First argument of matrix multiplication
    double *pBMatrix; // Second argument of matrix multiplication
    double *pCMatrix; // Result matrix
    int Size;         // Size of matrices
    int BlockSize;    // Sizes of matrix blocks
    double *pAblock;  // Initial block of matrix A
    double *pBblock;  // Initial block of matrix B
    double *pCblock;  // Block of result matrix C
    double *pMatrixAblock;
    double Start, Finish, Duration;
    setvbuf(stdout, 0, _IONBF, 0);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0)
    {
        Size = atoi(argv[1]);
        printf("Matrices' size = %d\n", Size);
    }

    GridSize = sqrt((double)ProcNum);
    if (ProcNum != GridSize * GridSize)
    {
        if (ProcRank == 0)
        {
            printf("Number of processes must be a perfect square \n");
        }
    }
    else
    {
        if (ProcRank == 0)
            printf("Parallel matrix multiplication program\n");

        // Creating the cartesian grid, row and column communcators
        CreateGridCommunicators();
        // Memory allocation and initialization of matrix elements
        ProcessInitialization(pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock,
                              pCblock, pMatrixAblock, Size, BlockSize);

        auto start = std::chrono::system_clock::now();
        int64_t scatterTime;
        DataDistribution(pAMatrix, pBMatrix, pMatrixAblock, pBblock, Size,
                         BlockSize, scatterTime);
        // Execution of the Fox method
        ParallelResultCalculation(pAblock, pMatrixAblock, pBblock,
                                  pCblock, BlockSize);
        // Gathering the result matrix
        ResultCollection(pCMatrix, pCblock, Size, BlockSize);
        TestResult(pAMatrix, pBMatrix, pCMatrix, Size);

        auto end = std::chrono::system_clock::now();
        auto elapsed =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        if (ProcRank == 0)
        {
            printf("\n");
            printf("Time of parallel execution with scatterv issue = %ld microseconds\n\n", elapsed.count());
            printf("Time of parallel clean execution = %ld microseconds\n\n", elapsed.count() - scatterTime);
        }

        // Process Termination
        ProcessTermination(pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock,
                           pCblock, pMatrixAblock);
    }
    MPI_Finalize();
}

// #include <stdio.h>
// #include <stdlib.h>
// #include <iostream>
// #include <fstream>
// #include <time.h>
// #include <mpi.h>
// #include <chrono>

// int ProcNum = 0;  // Number of available processes
// int ProcRank = 0; // Rank of current process

// void ReadDataInitialization(double *pMatrix, double *pVector, int Size)
// {
//     std::ifstream fin("matrix" + std::to_string(Size) + ".txt");
//     int i, j;
//     for (i = 0; i < Size; i++)
//     {
//         for (j = 0; j < Size; j++)
//         {
//             fin >> pMatrix[i * Size + j];
//         }
//     }

//     std::ifstream fin2("vector" + std::to_string(Size) + ".txt");
//     for (i = 0; i < Size; i++)
//     {
//         fin2 >> pVector[i];
//     }
// }
// // Function for memory allocation and data initialization
// void ProcessInitialization(double *&pMatrix, double *&pVector,
//                            double *&pResult, double *&pProcRows, double *&pProcResult,
//                            int &Size, int &RowNum)
// {
//     int RestRows; // Number of rows, that haven�t been distributed yet
//     int i;        // Loop variable
//     MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//     RestRows = Size;
//     for (i = 0; i < ProcRank; i++)
//         RestRows = RestRows - RestRows / (ProcNum - i);
//     RowNum = RestRows / (ProcNum - ProcRank);
//     pVector = new double[Size];
//     pResult = new double[Size];
//     pProcRows = new double[RowNum * Size];
//     pProcResult = new double[RowNum];
//     if (ProcRank == 0)
//     {
//         pMatrix = new double[Size * Size];
//         ReadDataInitialization(pMatrix, pVector, Size);
//     }
// }
// // Data distribution among the processes
// void DataDistribution(double *pMatrix, double *pProcRows, double *pVector,
//                       int Size, int RowNum, int64_t &scatterTime)
// {
//     int *pSendNum;       // the number of elements sent to the process
//     int *pSendInd;       // the index of the first data element sent to the process
//     int RestRows = Size; // Number of rows, that haven�t been distributed yet
//     MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     // Alloc memory for temporary objects
//     pSendInd = new int[ProcNum];
//     pSendNum = new int[ProcNum];
//     // Define the disposition of the matrix rows for current process
//     RowNum = (Size / ProcNum);
//     pSendNum[0] = RowNum * Size;
//     pSendInd[0] = 0;
//     for (int i = 1; i < ProcNum; i++)
//     {
//         RestRows -= RowNum;
//         RowNum = RestRows / (ProcNum - i);
//         pSendNum[i] = RowNum * Size;
//         pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
//     }
//     std::chrono::time_point<std::chrono::system_clock> start;
//     if (ProcRank == 0)
//         start = std::chrono::system_clock::now();
//     // Scatter the rows
//     MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows,
//                  pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     if (ProcRank == 0)
//     {
//         auto end = std::chrono::system_clock::now();
//         auto elapsed =
//             std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//         // std::cout << "scatterv time: " << elapsed.count() << " microseconds" << std::endl;
//         scatterTime = elapsed.count();
//     }
//     // Free the memory
//     delete[] pSendNum;
//     delete[] pSendInd;
// }
// // Function for gathering the result vector
// void ResultReplication(double *pProcResult, double *pResult, int Size,
//                        int RowNum)
// {
//     int i;               // Loop variable
//     int *pReceiveNum;    // Number of elements, that current process sends
//     int *pReceiveInd;    /* Index of the first element from current process
//        in result vector */
//     int RestRows = Size; // Number of rows, that haven�t been distributed yet
//     // Alloc memory for temporary objects
//     pReceiveNum = new int[ProcNum];
//     pReceiveInd = new int[ProcNum];
//     // Define the disposition of the result vector block of current processor
//     pReceiveInd[0] = 0;
//     pReceiveNum[0] = Size / ProcNum;
//     for (i = 1; i < ProcNum; i++)
//     {
//         RestRows -= pReceiveNum[i - 1];
//         pReceiveNum[i] = RestRows / (ProcNum - i);
//         pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
//     }
//     // Gather the whole result vector on every processor
//     MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_DOUBLE, pResult,
//                    pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
//     // Free the memory
//     delete[] pReceiveNum;
//     delete[] pReceiveInd;
// }
// // Function for sequential matrix-vector multiplication
// void SerialResultCalculation(double *pMatrix, double *pVector, double *pResult, int Size)
// {
//     int i, j; // Loop variables
//     for (i = 0; i < Size; i++)
//     {
//         pResult[i] = 0;
//         for (j = 0; j < Size; j++)
//             pResult[i] += pMatrix[i * Size + j] * pVector[j];
//     }
// }
// // Function for calculating partial matrix-vector multiplication
// void ParallelResultCalculation(double *pProcRows, double *pVector, double *pProcResult, int Size, int RowNum)
// {
//     int i, j; // Loop variables
//     for (i = 0; i < RowNum; i++)
//     {
//         pProcResult[i] = 0;
//         for (j = 0; j < Size; j++)
//             pProcResult[i] += pProcRows[i * Size + j] * pVector[j];
//     }
// }
// // Function for formatted matrix output
// void PrintMatrix(double *pMatrix, int RowCount, int ColCount)
// {
//     int i, j; // Loop variables
//     for (i = 0; i < RowCount; i++)
//     {
//         for (j = 0; j < ColCount; j++)
//             printf("%7.4f ", pMatrix[i * ColCount + j]);
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
// void TestDistribution(double *pMatrix, double *pVector, double *pProcRows,
//                       int Size, int RowNum)
// {
//     if (ProcRank == 0)
//     {
//         printf("Initial Matrix: \n");
//         PrintMatrix(pMatrix, Size, Size);
//         printf("Initial Vector: \n");
//         PrintVector(pVector, Size);
//     }
//     MPI_Barrier(MPI_COMM_WORLD);
//     for (int i = 0; i < ProcNum; i++)
//     {
//         if (ProcRank == i)
//         {
//             printf("\nProcRank = %d \n", ProcRank);
//             printf(" Matrix Stripe:\n");
//             PrintMatrix(pProcRows, RowNum, Size);
//             printf(" Vector: \n");
//             PrintVector(pVector, Size);
//         }
//         MPI_Barrier(MPI_COMM_WORLD);
//     }
// }
// void TestPartialResults(double *pProcResult, int RowNum)
// {
//     int i; // Loop variables
//     for (i = 0; i < ProcNum; i++)
//     {
//         if (ProcRank == i)
//         {
//             printf("\nProcRank = %d \n Part of result vector: \n", ProcRank);
//             PrintVector(pProcResult, RowNum);
//         }
//         MPI_Barrier(MPI_COMM_WORLD);
//     }
// }
// void TestResult(double *pMatrix, double *pVector, double *pResult,
//                 int Size)
// {
//     // Buffer for storing the result of serial matrix-vector multiplication
//     double *pSerialResult;
//     // Flag, that shows wheather the vectors are identical or not
//     int equal = 0;
//     int i; // Loop variable
//     if (ProcRank == 0)
//     {
//         pSerialResult = new double[Size];
//         SerialResultCalculation(pMatrix, pVector, pSerialResult, Size);
//         for (i = 0; i < Size; i++)
//         {
//             if (pResult[i] != pSerialResult[i])
//                 equal = 1;
//         }
//         if (equal == 1)
//             printf("The results of serial and parallel algorithms "
//                    "are NOT identical. Check your code.");
//         else
//             printf("The results of serial and parallel algorithms "
//                    "are identical.");
//     }
// }
// // Function for computational process termination
// void ProcessTermination(double *pMatrix, double *pVector, double *pResult,
//                         double *pProcRows, double *pProcResult)
// {
//     if (ProcRank == 0)
//         delete[] pMatrix;
//     delete[] pVector;
//     delete[] pResult;
//     delete[] pProcRows;
//     delete[] pProcResult;
// }
// int main(int argc, char *argv[])
// {
//     double *pMatrix; // The first argument - initial matrix
//     double *pVector; // The second argument - initial vector
//     double *pResult; // Result vector for matrix-vector multiplication
//     int Size;        // Sizes of initial matrix and vector
//     double *pProcRows;
//     double *pProcResult;
//     int RowNum;
//     double Start, Finish, DurationParallel, DurationSerial;

//     MPI_Init(&argc, &argv);
//     MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//     MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

//     if (ProcRank == 0)
//     {
//         std::cout << "Number of processes = " << ProcNum << std::endl;
//         // Size = ;
//         Size = atoi(argv[1]);
//         // printf("Matrix size: %d\n", Size);
//     }

//     ProcessInitialization(pMatrix, pVector, pResult, pProcRows, pProcResult, Size, RowNum);

//     Start = MPI_Wtime();
//     auto start = std::chrono::system_clock::now();
//     int64_t scatterTime;

//     DataDistribution(pMatrix, pProcRows, pVector, Size, RowNum, scatterTime);
//     ParallelResultCalculation(pProcRows, pVector, pProcResult, Size, RowNum);
//     ResultReplication(pProcResult, pResult, Size, RowNum);

//     Finish = MPI_Wtime();
//     auto end = std::chrono::system_clock::now();
//     DurationParallel = Finish - Start;

//     // Start = MPI_Wtime();
//     TestResult(pMatrix, pVector, pResult, Size);
//     auto elapsed =
//         std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//     // Finish = MPI_Wtime();

//     // DurationSerial = Finish - Start;
//     if (ProcRank == 0)
//     {
//         printf("\n");
//         printf("Time of parallel execution with scatterv issue = %ld microseconds\n\n", elapsed.count());
//         printf("Time of parallel clean execution = %ld microseconds\n\n", elapsed.count() - scatterTime);

//         // printf("Time of serial execution = %f\n", DurationSerial);
//     }
//     ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcResult);
//     MPI_Finalize();
// }