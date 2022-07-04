#include<mpi.h>
#include<cmath>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<iostream>
#include<bits/stdc++.h>
#include <climits>
using namespace std;
#include<chrono>
using namespace std::chrono;

vector<int> Dim_Vec; // Vec to Store Matrix Multiplication
string Path = "";   // To Store Optimal Path
char name = 'A';   // Help to Print Optimal Path
int cost = 0;     // To Store Optimal COst

int Total_Matrix = 0;
vector<vector<int>> Matrix_Ranges; // 2D Vector to Store R & C from File fro each matrix

void PrintTree(int i, int j, int n, int *bracket)
{
    // DFS Algorithum to get Optimal Path
	if (i == j){
        Path += name++;
		return;
	}
    Path += "(";
    int Mid = *((bracket+i*n)+j);
    PrintTree(i, Mid, n, bracket);
    PrintTree(Mid+1, j, n, bracket);
    Path += ")";
}

void Optimal_Multiplication()
{
    // Dynamic Programming Code to Get Optimal Path Multipkication
    int n = Dim_Vec.size();
	int Cost[n][n] = {0};
    int bracket[n][n] = {0};
	for(int d=2; d<Dim_Vec.size(); d++)
	{
		for(int i=1; i<Dim_Vec.size()-d+1; i++)
		{
			int j = i+d-1;
			Cost[i][j] = INT_MAX;
			for(int k=i; k<=j-1; k++)
			{
                // Equation to Get Cost 
				int new_cost = Cost[i][k]+Cost[k + 1][j]+Dim_Vec[i-1]*Dim_Vec[k]*Dim_Vec[j];
				if(new_cost < Cost[i][j])
				{
					Cost[i][j] = new_cost;
					bracket[i][j] = k;
				}
			}
		}
	}
    //  To Iterating Bracket tree
    PrintTree(1, n-1, n, (int*)bracket);
    cost = Cost[1][n - 1];
}

void Operation(vector<vector<int>> &A, 
vector<vector<int>> &B, vector<vector<int>> &C, int size, char op)
{
	for(int i=0; i<size; i++)
	{
		for(int j=0; j<size; j++)
		{
            if(op=='S')
            {
                C[i][j] = A[i][j]-B[i][j];
            }else if(op=='A'){
                C[i][j] = A[i][j]+B[i][j];
            }
		}
	}
}

void Strassen_Algorithm(vector<vector<int>> &A, 
vector<vector<int>> &B, vector<vector<int>> &C, int size)
{
	if(size==1){
		C[0][0] = A[0][0]*B[0][0];
		return;
	}else{
		int newSize = size/2;
		vector<int> z(newSize);
		vector<vector<int>> a(newSize, z);
        vector<vector<int>> b(newSize, z);
        vector<vector<int>> c(newSize, z);
        vector<vector<int>> d(newSize, z);
        vector<vector<int>> e(newSize, z);
        vector<vector<int>> f(newSize, z);
        vector<vector<int>> g(newSize, z);
        vector<vector<int>> h(newSize, z);
		
        for(int i=0; i<newSize; i++)
		{
			for(int j=0; j<newSize; j++)
			{
				a[i][j] = A[i][j];
				b[i][j] = A[i][j + newSize];
				c[i][j] = A[i + newSize][j];
				d[i][j] = A[i + newSize][j + newSize];

				e[i][j] = B[i][j];
				f[i][j] = B[i][j + newSize];
				g[i][j] = B[i + newSize][j];
				h[i][j] = B[i + newSize][j + newSize];
			}
		}
        
        vector<vector<int>> c11(newSize, z);
        vector<vector<int>> c12(newSize, z);
        vector<vector<int>> c21(newSize, z);
        vector<vector<int>> c22(newSize, z);
        
        vector<vector<int>> p1(newSize, z);
        vector<vector<int>> p2(newSize, z);
        vector<vector<int>> p3(newSize, z);
        vector<vector<int>> p4(newSize, z);
        vector<vector<int>> p5(newSize, z);
        vector<vector<int>> p6(newSize, z);
        vector<vector<int>> p7(newSize, z);
        
        vector<vector<int>> fResult(newSize, z);
		vector<vector<int>> sResult(newSize, z);
        
        //*p1=a*(f-h)
		Operation(f, h, sResult, newSize, 'S');
		Strassen_Algorithm(a, sResult, p1, newSize);

        //*p2=h*(a+b)
		Operation(a, b, fResult, newSize, 'A');
		Strassen_Algorithm(fResult, h, p2, newSize);

        //*p3=e*(c+d)
		Operation(c, d, fResult, newSize, 'A');
		Strassen_Algorithm(fResult, e, p3, newSize);

        //*p4=d*(g-e)
		Operation(g, e, sResult, newSize, 'S');
		Strassen_Algorithm(d, sResult, p4, newSize);

        //*p5=(a+d)*(e+h)
		Operation(a, d, fResult, newSize, 'A');
		Operation(e, h, sResult, newSize, 'A');
		Strassen_Algorithm(fResult, sResult, p5, newSize);

        //*p6=(b-d)*(g+h)
		Operation(b, d, fResult, newSize, 'S');
		Operation(g, h, sResult, newSize, 'A');
		Strassen_Algorithm(fResult, sResult, p6, newSize);

        //*p7=(a-c)*(e+f)
		Operation(a, c, fResult, newSize, 'S');
		Operation(e, f, sResult, newSize, 'A');
		Strassen_Algorithm(fResult, sResult, p7, newSize);

        // c11=p4+p5+p6-p2
		// c12=p1+p2
		// c21=p3+p4
		// c22=p1-p3+p5-p7

		Operation(p1, p2, c12, newSize, 'A');
		Operation(p3, p4, c21, newSize, 'A'); 

		Operation(p4, p5, fResult, newSize, 'A');
		Operation(fResult, p6, sResult, newSize, 'A');
		Operation(sResult, p2, c11, newSize, 'S');

		Operation(p1, p3, fResult, newSize, 'S');
		Operation(fResult, p5, sResult, newSize, 'A');
		Operation(sResult, p7, c22, newSize, 'S');

		for(int i=0; i<newSize; i++)
		{
			for(int j=0; j<newSize; j++)
			{
				C[i][j] = c11[i][j];
				C[i][j + newSize] = c12[i][j];
				C[i + newSize][j] = c21[i][j];
				C[i + newSize][j + newSize] = c22[i][j];
			}
		}
	}
}

vector<vector<int>> Start_Algo(vector<vector<int>> &A, 
vector<vector<int>> &B, int r1, int c1, int r2, int c2)
{
	int maxSize = max({r1, c1, r2, c2});
	int size = pow(2, int(ceil(log2(maxSize))));

	vector<int> col(size);
	vector<vector<int>> New_A(size, col), New_B(size, col), New_C(size, col);

	for(int i=0; i<r1; i++)
	{
		for(int j=0; j<c1; j++)
		{
			New_A[i][j] = A[i][j];
		}
	}

	for(int i=0; i<r2; i++)
	{
		for(int j=0; j<c2; j++)
		{
			New_B[i][j] = B[i][j];
		}
	}
	Strassen_Algorithm(New_A, New_B, New_C, size);
    return New_C;
}


void Read_File_And_Store(string file_name)
{
    // Reading FIle 
    vector <string> Matrix;
    string line;
    ifstream file(file_name);
    while(getline(file, line))
    {
        // Store All Matrix Information Into Vector 1D
        Matrix.push_back(line);
    }

    // Below Loop Store R & C in 2D Vector for Later Use
    Total_Matrix = Matrix.size();
    for(int i=0; i<Matrix.size(); i++)
    {
        string inter;
        vector <string> tokens;
        stringstream check(Matrix[i]);
        while(getline(check, inter, ' '))
        {
            tokens.push_back(inter);
        }

        int A, B;
        stringstream S1, S2;
        S1<<tokens[0];
        S1>>A;
        S2<<tokens[2];
        S2>>B;

        // Storing Matrix Coordinates
        vector<int> v;
        v.push_back(A);
        v.push_back(B);
        Matrix_Ranges.push_back(v);
        Dim_Vec.push_back(A);
        if(i==Matrix.size()-1){
            Dim_Vec.push_back(B);
        }
    }
    // Getting Optimal Path
    Optimal_Multiplication();
}

int** make_matrix(int R, int C)
{
    // Create 2D Dynamic Matrix and Filling with values
    int **M = new int*[R];
    for(int i=0; i<R; i++)
    {
        M[i] = new int[C];
    }

    for(int i=0; i<R; i++)
    {
        for(int j=0; j<C; j++)
        {
            M[i][j] = 2;
        }
    }
    return M;
}

void Display(int **M, int R, int C)
{
    // 2D pointer Printing Function
    cout<<"Dim : "<<R<<" X "<<C<<endl<<endl;
    cout<<"Matrix"<<endl;
    for(int i=0; i<R; i++)
    {
        for(int j=0; j<C; j++)
        {
            cout<<M[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

int main(int argc, char **argv)
{
    auto start = high_resolution_clock::now();

    // Basic MPI Commands
    int no_of_process, process_Id;
    int no_of_salve, source, dest, rows, offset;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Id);
    MPI_Comm_size(MPI_COMM_WORLD, &no_of_process);

    no_of_salve = no_of_process-1;

    if(process_Id==0)
    {
        MPI_Status status;
        // Reading File
        Read_File_And_Store("file.txt");
        // Printing Optimal Part
        cout<<"Optimal Parenthesization is : "<<Path<<endl;
	    cout<<"Optimal Cost is : "<<cost<<endl<<endl;
        
        //  Creaing 2D Matrix Array and Store Matrix at each index
        int **All_Matrix[Total_Matrix]; 
        for(int i=0; i<Total_Matrix; i++)
        {
            All_Matrix[i] = make_matrix(Matrix_Ranges[i][0], Matrix_Ranges[i][1]);
        }
        cout<<"No. of Matrixes: "<<Total_Matrix<<endl;

        //  Generic Loop to Send TO N Child Process
        for(int m=0; m<Total_Matrix-1; m++)
        {
            // Creating Result Matrix
            int **Result_ptr = make_matrix(Matrix_Ranges[m][0], Matrix_Ranges[m+1][1]);    

            // Matrix That help me to send information to Process        
            int Matrix_A[Matrix_Ranges[m][0]][Matrix_Ranges[m][1]];
            int Matrix_B[Matrix_Ranges[m+1][0]][Matrix_Ranges[m+1][1]];

            // Storing Values in Matrix A
            for(int j=0; j<Matrix_Ranges[m][0]; j++)
            {
                for(int k=0; k<Matrix_Ranges[m][1]; k++)
                {
                    Matrix_A[j][k] = All_Matrix[m][j][k];
                }
            }
            
            // Storing Values in Matrix B
            for(int j=0; j<Matrix_Ranges[m+1][0]; j++)
            {
                for(int k=0; k<Matrix_Ranges[m+1][1]; k++)
                {
                    Matrix_B[j][k] = All_Matrix[m+1][j][k];
                }
            }

            int A_size = Matrix_Ranges[m][0]*Matrix_Ranges[m][1];
            int B_size = Matrix_Ranges[m+1][0]*Matrix_Ranges[m+1][1];

            //  Sending Matrix A Coord
            MPI_Send(&Matrix_Ranges[m][0], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            MPI_Send(&Matrix_Ranges[m][1], 1, MPI_INT, 1, 1, MPI_COMM_WORLD);

            //  Sending Matrix B Coord
            MPI_Send(&Matrix_Ranges[m+1][0], 1, MPI_INT, 1, 2, MPI_COMM_WORLD);
            MPI_Send(&Matrix_Ranges[m+1][1], 1, MPI_INT, 1, 3, MPI_COMM_WORLD);
            //  Sending Matrix A
            MPI_Send(&Matrix_A, A_size, MPI_INT, 1, 4, MPI_COMM_WORLD);
            //  Sending Matrix B 
            MPI_Send(&Matrix_B, B_size, MPI_INT, 1, 5, MPI_COMM_WORLD);

            int Result[Matrix_Ranges[m][0]][Matrix_Ranges[m+1][1]];

            // Receving Result
            int Receive = Matrix_Ranges[m][0]*Matrix_Ranges[m+1][1];
            MPI_Recv(&Result, Receive, MPI_INT, 1, 6, MPI_COMM_WORLD, &status);

            // Storing Result into 2D Matrix for Next Iteration
            int **Output = make_matrix(Matrix_Ranges[m][0], Matrix_Ranges[m+1][1]);
            for(int j=0; j<Matrix_Ranges[m][0]; j++)
            {
                for(int k=0; k<Matrix_Ranges[m+1][1]; k++)
                {
                    Output[j][k] = Result[j][k];
                }
            }

            // Storing 2D Result into 2D ALL_Matrix Array
            All_Matrix[m+1] = Output;
            vector<int> v;
            v.push_back(Matrix_Ranges[m][0]);
            v.push_back(Matrix_Ranges[m+1][1]);
            Matrix_Ranges.at(m+1) = v;

            // Display(Output, Matrix_Ranges[m][0], Matrix_Ranges[m+1][1]);
        }
        // Displaying Final Result
        Display(All_Matrix[Total_Matrix-1], 
        Matrix_Ranges[Total_Matrix-1][0], 
        Matrix_Ranges[Total_Matrix-1][1]);
    }else if(process_Id==1){
        Read_File_And_Store("file.txt");
        MPI_Status status;
        for(int s=0; s<Total_Matrix-1; s++)
        {
            int A_r = 0, A_c = 0;
            MPI_Recv(&A_r, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&A_c, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

            int B_r = 0, B_c = 0;
            MPI_Recv(&B_r, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&B_c, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);

            int A[A_r][A_c];
            int B[B_r][B_c]; 
            
            int A_si = A_r*A_c;
            MPI_Recv(&A, A_si, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);

            int B_si = B_r*B_c;
            MPI_Recv(&B, B_si, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);

            vector<int> A_col(A_c);
            vector<vector<int>> Vec_A(A_r, A_col);

            vector<int> B_col(B_c);
            vector<vector<int>> Vec_B(B_r, B_col);

            for(int i=0; i<A_r; i++)
            {
                for(int j=0; j<A_c; j++)
                {
                    Vec_A[i][j] = A[i][j];
                }
            }

            for(int i=0; i<B_r; i++)
            {
                for(int j=0; j<B_c; j++)
                {
                    Vec_B[i][j] = B[i][j];
                }
            }
            vector<vector<int>> C = Start_Algo(Vec_A, Vec_B, A_r, A_c, B_r, B_c);

            int Result[A_r][B_c];
            for(int i=0; i<A_r; i++)
            {
                for(int j=0; j<B_c; j++)
                {
                    Result[i][j] = C[i][j];
                }
            }

            int Res_Si = A_r*B_c;
            MPI_Send(&Result, Res_Si, MPI_INT, 0, 6, MPI_COMM_WORLD);
        }
    }

    if(process_Id==0){
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout<<"Time Taken: "<<duration.count()<<endl;
    }
    
    MPI_Finalize();
    return 0;
}