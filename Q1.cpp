#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<cstdlib>
#include <climits>
using namespace std;
#include<chrono>
using namespace std::chrono;
vector<int> Dim_Vec;

string Path = "";  // To Store the Optimal Matrix Pattern
char name = 'A';  // Help to Print Pattern

void PrintTree(int i, int j, int n, int *bracket)
{
    // Using DFS Algorithum to retrive Optimal pattern
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
//  Implemeted the Dynamic Algothrium To find the Optimal Path
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
                // Formula to Fill the Cost Table
                // Dim_Vec have all matrix dimesionals
				int new_cost = Cost[i][k]+Cost[k + 1][j]+Dim_Vec[i-1]*Dim_Vec[k]*Dim_Vec[j];
				if(new_cost < Cost[i][j])
				{
					Cost[i][j] = new_cost;
					bracket[i][j] = k;
				}
			}
		}
	}
    // Using Depth First Search to print & Get the Optimal Path
    PrintTree(1, n-1, n, (int*)bracket);
    cout<<"Optimal Parenthesization is : "<<Path<<endl;
	cout<<"Optimal Cost is : " << Cost[1][n - 1]<<endl;
}

void Matrix_M(int **R, int **M1, int A_r, int A_c, int **M2, int B_r, int B_c)
{
    // Matrix Multiplication Algorithum
    // Using 2D Pointer to Store thr Result
    for(int i=0; i<A_r; i++) 
    {
        for(int j=0; j<B_c; j++) 
        {
            R[i][j] = 0;
            for(int k=0; k<A_c; k++) 
            {
                R[i][j] += M1[i][k]*M2[k][j];
            }
        }
    }
}

int** make_matrix(int R, int C)
{
    // This Function Create Dynamic Matrix and Fill the Matrix
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
    // Using M 2D Pointer to Print the Matrix Values
    // Row X Columns
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

int main()
{
    auto start = high_resolution_clock::now();
    // Usind 2D vector to Store all Matrix R & C From File
    vector<vector<int>> Matrix_Ranges;
    // To Store Matrix Dims 
    vector <string> Matrix;

    string line;
    ifstream file("file.txt");
    int iterate = 0;
    int **A;
    while(getline(file, line))
    {
        // Pushing Each line into Matrix Vector For later use
        Matrix.push_back(line);
    }

    // Creating 2D Pointer to Store all Matrix from File
    int **All_Matrix[Matrix.size()];

    // Below Loop Get the R & C values and Store into 2D Vector For Later USE
    for(int i=0; i<Matrix.size(); i++)
    {
        string inter;
        vector<string> tokens;
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

        All_Matrix[i] = make_matrix(A, B);
        vector<int> v;
        Dim_Vec.push_back(A);
        v.push_back(A);
        v.push_back(B);

        Matrix_Ranges.push_back(v);
        
        if(i==Matrix.size()-1){
            Dim_Vec.push_back(B);
        }
    }
    // Finding Optimal Multiplication 
    Optimal_Multiplication();

    // Generic Loop To Run For N loop 
    for(int i=0; i<Matrix.size()-1; i++)
    {
        // Making Result Matrix
        int **Result = make_matrix(Matrix_Ranges[i][0], Matrix_Ranges[i+1][1]);
        // SEDN TO  Matrix Multplication Function
        Matrix_M(Result, 
        All_Matrix[i], Matrix_Ranges[i][0], Matrix_Ranges[i][1], 
        All_Matrix[i+1], Matrix_Ranges[i+1][0], Matrix_Ranges[i+1][1]);

        // Store Result into 2D Pointer for Next Iteration
        All_Matrix[i+1] = Result;
        // Updating Coordination Vector
        vector<int> v;
        v.push_back(Matrix_Ranges[i][0]);
        v.push_back(Matrix_Ranges[i+1][1]);
        Matrix_Ranges.at(i+1) = v;
    }

    // Printing Final Result
    Display(All_Matrix[Matrix.size()-1], 
    Matrix_Ranges[Matrix.size()-1][0], 
    Matrix_Ranges[Matrix.size()-1][1]);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    cout<<"Time Taken: "<<duration.count()<<endl;
    return 0;
}