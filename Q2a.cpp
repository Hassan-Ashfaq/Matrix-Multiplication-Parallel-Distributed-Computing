#include<mpi.h>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<cstdlib>
#include<climits>
#include<chrono>
using namespace std::chrono;
using namespace std;

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

    Total_Matrix = Matrix.size();
    // Below Loop Store R & C in 2D Vector for Later Use
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

    int no_of_process, process_Id;
    int no_of_salve, source, dest, rows, offset;

    // Basic MPI Commands
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

        //  Generic Loop to Send TO N Child Process
        for(int i=0; i<Total_Matrix-1; i++)
        {
            // Creating Result Matrix
            int **Result_ptr = make_matrix(Matrix_Ranges[i][0], Matrix_Ranges[i+1][1]);            
            int Result[Matrix_Ranges[i][0]][Matrix_Ranges[i+1][1]];

            // Matrix That help me to send information to Process
            int Matrix_A[Matrix_Ranges[i][0]][Matrix_Ranges[i][1]];
            int Matrix_B[Matrix_Ranges[i+1][0]][Matrix_Ranges[i+1][1]];
            
            // Storing Values in Matrix A
            for(int j=0; j<Matrix_Ranges[i][0]; j++)
            {
                for(int k=0; k<Matrix_Ranges[i][1]; k++)
                {
                    Matrix_A[j][k] = All_Matrix[i][j][k];
                }
            }
            
            // Storing Values in Matrix B
            for(int j=0; j<Matrix_Ranges[i+1][0]; j++)
            {
                for(int k=0; k<Matrix_Ranges[i+1][1]; k++)
                {
                    Matrix_B[j][k] = 2; //All_Matrix[i+1][j][k];
                }
            }

            // Spliting Matrix A in N Slave Process
            rows = Matrix_Ranges[i][0]/no_of_salve;
            // cout<<"No.of Rows to Each Salve: "<<rows<<endl;
            
            offset = 0;
            int temp_row = rows;    
            // Sending To Childrens
            int limit = rows;
            // Generic Loop To Send Splitted Matrix into Slave Process
            for(int send=1; send<=no_of_salve; send++)
            {
                // cout<<"Offset: "<<offset<<":: rows: "<<rows<<endl;
                MPI_Send(&Matrix_Ranges[i][0], 1, MPI_INT, send, 1, MPI_COMM_WORLD);
                MPI_Send(&Matrix_Ranges[i][1], 1, MPI_INT, send, 2, MPI_COMM_WORLD);

                MPI_Send(&Matrix_Ranges[i+1][0], 1, MPI_INT, send, 3, MPI_COMM_WORLD);
                MPI_Send(&Matrix_Ranges[i+1][1], 1, MPI_INT, send, 4, MPI_COMM_WORLD);

                int A = Matrix_Ranges[i][1];
                int B = Matrix_Ranges[i+1][0]*Matrix_Ranges[i+1][1];

                MPI_Send(&offset, 1, MPI_INT, send, 5, MPI_COMM_WORLD);
                MPI_Send(&rows, 1, MPI_INT, send, 6, MPI_COMM_WORLD);

                MPI_Send(&Matrix_B, B, MPI_INT, send, 8, MPI_COMM_WORLD);
                MPI_Send(&Matrix_A[offset][0], rows*A, MPI_INT, send, 7, MPI_COMM_WORLD);
                offset += limit; 
                rows += limit;
            }
            // cout<<endl;

            // Receving To Childrens
            offset = 0;
            int N = Matrix_Ranges[i+1][1];
            // Receving Result from Slave Process
            for(int receive=1; receive<=no_of_salve; receive++)
            {
                // cout<<"****************************************"<<endl;
                // cout<<"Rec From: "<<receive<<endl;
                // cout<<"Rows : "<<temp_row<<endl;
                // cout<<"Rec: "<<offset<<" "<<offset+limit<<endl;
                // cout<<"Rec Size : "<<temp_row*N<<endl;
                // cout<<"****************************************"<<endl<<endl;
                MPI_Recv(&Result[offset][0], temp_row*N, MPI_INT, receive, 99, MPI_COMM_WORLD, &status);
                offset += limit;
            }

            // Storing Result into 2D Matrix for Next Iteration
            int **Output = make_matrix(Matrix_Ranges[i][0], Matrix_Ranges[i+1][1]);
            for(int j=0; j<Matrix_Ranges[i][0]; j++)
            {
                for(int k=0; k<Matrix_Ranges[i+1][1]; k++)
                {
                    Output[j][k] = Result[j][k];
                }
            }

            // Storing 2D Result into 2D ALL_Matrix Array
            All_Matrix[i+1] = Output;
            vector<int> v;
            v.push_back(Matrix_Ranges[i][0]);
            v.push_back(Matrix_Ranges[i+1][1]);
            Matrix_Ranges.at(i+1) = v;

            // Display(Output, Matrix_Ranges[i][0], Matrix_Ranges[i+1][1]);
            // break;
        }
        // Displaying Final Result
        Display(All_Matrix[Total_Matrix-1], 
        Matrix_Ranges[Total_Matrix-1][0], 
        Matrix_Ranges[Total_Matrix-1][1]);

    }else if(process_Id>0){
        // Reading File
        Read_File_And_Store("file.txt");
        // Loop To Receive Matrix A chuck and Other Infromation for Muliplication
        for(int p=0; p<Total_Matrix-1; p++)
        {
            MPI_Status status;
            int source = 0;
            int A_r, A_c, B_r, B_c;
            int offset, rows;
            MPI_Recv(&A_r, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&A_c, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);

            MPI_Recv(&B_r, 1, MPI_INT, source, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(&B_c, 1, MPI_INT, source, 4, MPI_COMM_WORLD, &status);
            MPI_Recv(&offset, 1, MPI_INT, source, 5, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, source, 6, MPI_COMM_WORLD, &status);

            int R[A_r][B_c];
            int A[A_r][A_c];
            int B[B_r][B_c];

            int A_S = A_c;
            int B_S = B_r*B_c; 

            MPI_Recv(&B, B_S, MPI_INT, source, 8, MPI_COMM_WORLD, &status);
            MPI_Recv(&A[offset][0], rows*A_S, MPI_INT, source, 7, MPI_COMM_WORLD, &status);
            int new_Row = rows-offset;
            // Multiplying Matrix A Chuck with Whole Matrix B
            for(int i=offset; i<rows; i++) 
            {
                for(int j=0; j<B_c; j++) 
                {
                    R[i][j] = 0;
                    for(int k=0; k<A_c; k++) 
                    {
                        R[i][j] += A[i][k]*B[k][j];
                    }
                }
            }
            // Sending Result Back to Master
            MPI_Send(&R[offset][0], new_Row*B_c, MPI_INT, 0, 99, MPI_COMM_WORLD);
            
            // cout<<"--------------------------------"<<endl;
            // cout<<"Rank : "<<process_Id<<endl;
            // cout<<"Offset : "<<offset<<endl;
            // cout<<"Rows : "<<new_Row<<endl;
            // cout<<"Send : "<<new_Row*B_c<<endl;
            // cout<<"--------------------------------"<<endl<<endl; 
            // break;
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