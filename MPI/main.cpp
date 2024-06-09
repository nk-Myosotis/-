#include <iostream>
#include <mpi.h>
#include <sys/time.h>
#include <pmmintrin.h>
#include <omp.h>

using namespace std;

static const int N = 1000;
static const int thread_count = 4;

float arr[N][N];
float A[N][N];

void init_A(float arr[][N])
{
    for(int i=0;i<N;i++)
    {
        arr[i][i]=1.0;
        for(int j=i+1;j<N;j++)
            arr[i][j]=rand();
    }
    for(int k=0;k<N;k++)
        for(int i=k+1;i<N;i++)
            for(int j=0;j<N;j++)
                arr[i][j]+=A[k][j];
}




void reset_A(float A[][N], float arr[][N])
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = arr[i][j];
}


void print_A(float A[][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}

void ordinary()
{
    reset_A(A, arr);
    timeval t_start;
    timeval t_end;
    gettimeofday(&t_start, NULL);

    for (int k = 0; k < N; k++)
    {
        for (int j = k + 1; j < N; j++)
        {
            A[k][j] = A[k][j] * 1.0 / A[k][k];
        }
        A[k][k] = 1.0;

        for (int i = k + 1; i < N; i++)
        {
            for (int j = k + 1; j < N; j++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
    gettimeofday(&t_end, NULL);
    cout << "ordinary time cost: "
        << 1000 * (t_end.tv_sec - t_start.tv_sec) +
        0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
}



void block_LU(float A[][N], int rank, int num_proc)
{
    int block = N / num_proc;
    int remain = N % num_proc;

    int begin = rank * block;
    int end = rank != num_proc - 1 ? begin + block : begin + block + remain;

    for (int k = 0; k < N; k++)
    {
        //��ǰ�����Լ����̵����񡪡�������ȥ
        if (k >= begin && k < end)
        {
            for (int j = k + 1; j < N; j++)
                A[k][j] = A[k][j] / A[k][k];
            A[k][k] = 1.0;
            ////������Ϣ���򱾽��̺���Ľ��̣�
            //for (int p = rank + 1; p < num_proc; p++)
            //    MPI_Send(&A[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
            //������Ϣ���������������̣�
            for (int p = 0; p < num_proc; p++)
                if (p != rank)
                    MPI_Send(&A[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);

        }
        //��ǰ�в����Լ����̵����񡪡�������Ϣ
        else
        {
            ////������Ϣ������λ���Լ�ǰ��Ľ��̵���Ϣ��
            //int cur_p = k / block;
            //if (cur_p < rank)
            //    MPI_Recv(&A[k], N, MPI_FLOAT, cur_p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //������Ϣ�����������������̵���Ϣ��
            int cur_p = k / block;
             MPI_Recv(&A[k], N, MPI_FLOAT, cur_p, 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

        //��ȥ����
        for (int i = begin; i < end && i < N; i++)
        {
            if (i >= k + 1)
            {
                for (int j = k + 1; j < N; j++)
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                A[i][k] = 0.0;
            }
        }
    }
}


void block_mpi()
{

    timeval t_start;
    timeval t_end;

    int num_proc;//������
    int rank;//ʶ����ý��̵�rank��ֵ��0~size-1

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int block = N / num_proc;
    int remain = N % num_proc;

    //0�Ž��̡������񻮷�
    if (rank == 0)
    {
        reset_A(A, arr);
        gettimeofday(&t_start, NULL);
        //���񻮷�
        for (int i = 1; i < num_proc; i++)
        {
            if (i != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Send(&A[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Send(&A[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }
        }
        block_LU(A, rank, num_proc);
        //������0�Ž����Լ��������������������̴���֮��Ľ��
        for (int i = 1; i < num_proc; i++)
        {
            if (i != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Recv(&A[i * block + j], N, MPI_FLOAT, i, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Recv(&A[i * block + j], N, MPI_FLOAT, i, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        gettimeofday(&t_end, NULL);
        cout << "Block MPI LU time cost: "
            << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
        //print_A(A);
    }

    //��������
    else
    {
        //��0�Ž����Ƚ�������
        if (rank != num_proc - 1)
        {
            for (int j = 0; j < block; j++)
                MPI_Recv(&A[rank * block + j], N, MPI_FLOAT, 0, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            for (int j = 0; j < block + remain; j++)
                MPI_Recv(&A[rank * block + j], N, MPI_FLOAT, 0, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        block_LU(A, rank, num_proc);
        //��0�Ž����������֮�󣬽�������ص�0�Ž���
        if (rank != num_proc - 1)
        {
            for (int j = 0; j < block; j++)
                MPI_Send(&A[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
        else
        {
            for (int j = 0; j < block + remain; j++)
                MPI_Send(&A[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
    }
}
void block_LU_opt(float A[][N], int rank, int num_proc)
{
    __m128 t1, t2, t3;
    int block = N / num_proc;
    int remain = N % num_proc;
    int begin = rank * block;
    int end = rank != num_proc - 1 ? begin + block : begin + block + remain;
    #pragma omp parallel for shared(A) private(t1, t2, t3) schedule(dynamic, 24)
    for (int k = 0; k < N; k++)
    {
        if (k >= begin && k < end)
        {
            float temp1[4] = { A[k][k], A[k][k], A[k][k], A[k][k] };
            t1 = _mm_loadu_ps(temp1);

            for (int j = k + 1; j < N - 3; j += 4)
            {
                t2 = _mm_loadu_ps(A[k] + j);
                t3 = _mm_div_ps(t2, t1);
                _mm_storeu_ps(A[k] + j, t3);
            }
            for (int j = N - N % 4; j < N; j++)
            {
                A[k][j] = A[k][j] / A[k][k];
            }
            A[k][k] = 1.0;
            #pragma omp parallel for schedule(dynamic, 24)
            for (int p = rank + 1; p < num_proc; p++)
                MPI_Send(&A[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
        }
        else
        {
            int cur_p = k / block;
            if (cur_p < rank)
                MPI_Recv(&A[k], N, MPI_FLOAT, cur_p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (int i = begin; i < end && i < N; i++)
        {
            if (i >= k + 1)
            {
                float temp2[4] = { A[i][k], A[i][k], A[i][k], A[i][k] };
                t1 = _mm_loadu_ps(temp2);
                for (int j = k + 1; j <= N - 3; j += 4)
                {
                    t2 = _mm_loadu_ps(A[i] + j);
                    t3 = _mm_loadu_ps(A[k] + j);
                    t3 = _mm_mul_ps(t1, t3);
                    t2 = _mm_sub_ps(t2, t3);
                    _mm_storeu_ps(A[i] + j, t2);
                }
                for (int j = N - N % 4; j < N; j++)
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                A[i][k] = 0;
            }
        }
    }
}

void block_mpi_opt()
{

    timeval t_start;
    timeval t_end;

    int num_proc;//������
    int rank;//ʶ����ý��̵�rank��ֵ��0~size-1

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int block = N / num_proc;
    int remain = N % num_proc;

    //0�Ž��̡������񻮷�
    if (rank == 0)
    {
        reset_A(A, arr);
        gettimeofday(&t_start, NULL);
        //���񻮷�
        for (int i = 1; i < num_proc; i++)
        {
            if (i != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Send(&A[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Send(&A[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            }
        }
        block_LU_opt(A, rank, num_proc);
        //������0�Ž����Լ��������������������̴���֮��Ľ��
        for (int i = 1; i < num_proc; i++)
        {
            if (i != num_proc - 1)
            {
                for (int j = 0; j < block; j++)
                    MPI_Recv(&A[i * block + j], N, MPI_FLOAT, i, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                for (int j = 0; j < block + remain; j++)
                    MPI_Recv(&A[i * block + j], N, MPI_FLOAT, i, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        gettimeofday(&t_end, NULL);
        cout << "Block MPI LU with SSE and OpenMP time cost: "
            << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
    }

    //��������
    else
    {
        //��0�Ž����Ƚ�������
        if (rank != num_proc - 1)
        {
            for (int j = 0; j < block; j++)
                MPI_Recv(&A[rank * block + j], N, MPI_FLOAT, 0, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            for (int j = 0; j < block + remain; j++)
                MPI_Recv(&A[rank * block + j], N, MPI_FLOAT, 0, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        block_LU_opt(A, rank, num_proc);
        //�����������Ž��̷��ؽ��
        if (rank != num_proc - 1)
        {
            for (int j = 0; j < block; j++)
                MPI_Send(&A[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
        else
        {
            for (int j = 0; j < block + remain; j++)
                MPI_Send(&A[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
    }
}
void recycle_LU(float A[][N], int rank, int num_proc)
{
    for (int k = 0; k < N; k++)
    {
        //��ǰ�����Լ����̵����񡪡����г���
        if (int(k % num_proc) == rank)
        {
            for (int j = k + 1; j < N; j++)
                A[k][j] = A[k][j] / A[k][k];
            A[k][k] = 1.0;
            //���������������̷�����Ϣ
            for (int p = 0; p < num_proc; p++)
                if (p != rank)
                    MPI_Send(&A[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
        }
        //��ǰ�в����Լ����̵����񡪡�������Ϣ
        else
        {
            MPI_Recv(&A[k], N, MPI_FLOAT, int(k % num_proc), 2,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (int i = k + 1; i < N; i++)
        {
            if (int(i % num_proc) == rank)
            {
                for (int j = k + 1; j < N; j++)
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                A[i][k] = 0.0;
            }
        }
    }
}
void recycle_mpi()
{
    timeval t_start;
    timeval t_end;

    int num_proc;//������
    int rank;//ʶ����ý��̵�rank��ֵ��0~size-1

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //0�Ž��̡������񻮷�
    if (rank == 0)
    {
        reset_A(A, arr);
        gettimeofday(&t_start, NULL);
        //���񻮷�
        for (int i = 0; i < N; i++)
        {
            int flag = i % num_proc;
            if (flag == rank)
                continue;
            else
                MPI_Send(&A[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
        }
        recycle_LU(A, rank, num_proc);
        //������0�Ž����Լ��������������������̴���֮��Ľ��
        for (int i = 0; i < N; i++)
        {
            int flag = i % num_proc;
            if (flag == rank)
                continue;
            else
                MPI_Recv(&A[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        gettimeofday(&t_end, NULL);
        cout << "Recycle MPI LU time cost: "
            << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
        //print_A(A);
    }
    //��������
    else
    {
        //��0�Ž����Ƚ�������
        for (int i = rank; i < N; i += num_proc)
        {
            MPI_Recv(&A[i], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //ִ������
        recycle_LU(A, rank, num_proc);
        //��0�Ž����������֮�󣬽�������ص�0�Ž���
        for (int i = rank; i < N; i += num_proc)
        {
            MPI_Send(&A[i], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
    }
}
void recycle_LU_opt(float A[][N], int rank, int num_proc)
{
    __m128 t1, t2, t3;
#pragma omp parallel for shared(A) private(t1, t2, t3) schedule(dynamic, 24)
    for (int k = 0; k < N; k++)
    {
        if (int(k % num_proc) == rank)
        {
            float temp1[4] = { A[k][k], A[k][k], A[k][k], A[k][k] };
            t1 = _mm_loadu_ps(temp1);

            for (int j = k + 1; j < N - 3; j += 4)
            {
                t2 = _mm_loadu_ps(A[k] + j);
                t3 = _mm_div_ps(t2, t1);
                _mm_storeu_ps(A[k] + j, t3);
            }
            for (int j = N - N % 4; j < N; j++)
            {
                A[k][j] = A[k][j] / A[k][k];
            }
            A[k][k] = 1.0;
            #pragma omp parallel for schedule(dynamic, 24)
            for (int p = 0; p < num_proc; p++)
                if (p != rank)
                    MPI_Send(&A[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Recv(&A[k], N, MPI_FLOAT, int(k % num_proc), 2,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (int i = k + 1; i < N; i++)
        {
            if (int(i % num_proc) == rank)
            {
                float temp2[4] = { A[i][k], A[i][k], A[i][k], A[i][k] };
                t1 = _mm_loadu_ps(temp2);
                for (int j = k + 1; j <= N - 3; j += 4)
                {
                    t2 = _mm_loadu_ps(A[i] + j);
                    t3 = _mm_loadu_ps(A[k] + j);
                    t3 = _mm_mul_ps(t1, t3);
                    t2 = _mm_sub_ps(t2, t3);
                    _mm_storeu_ps(A[i] + j, t2);
                }
                for (int j = N - N % 4; j < N; j++)
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                A[i][k] = 0;
            }
        }
    }
}

void recycle_mpi_opt()
{
    timeval t_start;
    timeval t_end;

    int num_proc;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        reset_A(A, arr);
        gettimeofday(&t_start, NULL);
        for (int i = 0; i < N; i++)
        {
            int flag = i % num_proc;
            if (flag == rank)
                continue;
            else
                MPI_Send(&A[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
        }
        recycle_LU_opt(A, rank, num_proc);
        for (int i = 0; i < N; i++)
        {
            int flag = i % num_proc;
            if (flag == rank)
                continue;
            else
                MPI_Recv(&A[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        gettimeofday(&t_end, NULL);
        cout << "Recycle MPI LU with SSE and OpenMP time cost: "
            << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;

    }
    else
    {
        for (int i = rank; i < N; i += num_proc)
        {
            MPI_Recv(&A[i], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        recycle_LU_opt(A, rank, num_proc);
        for (int i = rank; i < N; i += num_proc)
        {
            MPI_Send(&A[i], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        }
    }
}

int main()
{
    init_A(arr);

    MPI_Init(NULL, NULL);
    ordinary();
    block_mpi();
    recycle_mpi();
    block_mpi_opt();
    recycle_mpi_opt();

    MPI_Finalize();

}
