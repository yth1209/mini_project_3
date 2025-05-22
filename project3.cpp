#include<stdio.h>
#include<stdlib.h>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

using namespace std;

void initDoubleSquareArray(int m, bool*** arr){
    *arr = new bool*[m];
    for(int i = 0; i < m; i++){
        (*arr)[i] = new bool[m];
        for(int j = 0; j < m; j++){
            (*arr)[i][j] = false;
        }
    }
}

void freeDoubleSquareArray(int m, bool*** arr){
    for(int i = 0; i < m; i++){
        delete[] (*arr)[i];
    }
    delete[] (*arr);
}

void readBoard(int m, int g, bool*** board){
    string line;
    for(int i = g; i < m+g; i++){
        cin >> line; 
        for(int j = g; j < m+g; j++){
            (*board)[i][j] = (line[j - g] == '#');
        }
    }
}

void printBoard(int m, int g, bool*** board){
    for(int i = g; i < m + g; i++){
        for(int j = g; j < m + g; j++){
            cout << ((*board)[i][j] ? '#' : '.');
        }
        cout << endl;
    }
}

int getNeighCnt(int i, int j, bool*** board) {
    int cntNeigh = 0;

    cntNeigh += (*board)[i][j-1];
    cntNeigh += (*board)[i][j+1];
    
    cntNeigh += (*board)[i-1][j]; 
    cntNeigh += (*board)[i-1][j-1];
    cntNeigh += (*board)[i-1][j+1];
 
    cntNeigh += (*board)[i+1][j];
    cntNeigh += (*board)[i+1][j-1];
    cntNeigh += (*board)[i+1][j+1];
    
    return cntNeigh;
}

void step(int m, int g, bool*** board){
    vector<pair<int,int>> changed_cells;

    for(int i = g; i < m + g; i++){
        for(int j = g; j < m + g; j++){
            int neighCnt = getNeighCnt(i, j, board);

            if((*board)[i][j]){
                if(neighCnt < 2 || neighCnt > 3) changed_cells.push_back({i,j});
            } else {
                if(neighCnt == 3) changed_cells.push_back({i,j});
            }
        }
    }

    for(auto cell: changed_cells){
        (*board)[cell.first][cell.second] = !(*board)[cell.first][cell.second]; 
    }
}

int main(int argc, char* argv[]) {    
    
    int npes, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int d = sqrt(npes);
    if(d*d != npes) {
        cerr << "Number of processes must be a perfect square" << endl;
        return 1;
    }


    int m; // size of board    
    int n; // final generation
    int g; // size of ghost cell
    bool** board; 
    bool** local_board; 

    if(myrank == 0){
        cin >> m >> n >> g; // read m, n, g from txt input
    }
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int board_size = (m+ m % d) + 2 * g; //padding with dimmension
    int l_comp_size = (m + m % d) / d; //local computing size
    int l_board_size = l_comp_size + 2 * g; //local board size

    initDoubleSquareArray(l_board_size, &local_board);
    if(myrank == 0){
        initDoubleSquareArray(board_size, &board);
        readBoard(m, g, &board); // read board from txt input
        
        //SEND the partition of the board to each process
        for(int ri = 0; ri < d; ri++){
            for(int rj = 0; rj < d; rj++){
                int rank = ri * d + rj;
                //SEND (ri*l_comp_size, ri*l_comp_size + l_board_size) * (rj*l_comp_size, rj*l_comp_size + l_board_size)
                for(int i = 0; i < l_board_size; i++){
                    MPI_Send(&board[i + ri*l_comp_size][rj*l_comp_size], l_board_size, MPI_CXX_BOOL, rank, i, MPI_COMM_WORLD);
                }
            }
        }
    }; 

    //RECIEIVE the partition of the board from the root process
    for(int i = 0; i < l_board_size; i++){
        MPI_Status status;
        MPI_Recv(local_board[i], l_board_size, MPI_CXX_BOOL, 0, i, MPI_COMM_WORLD, &status);
    }
    // printBoard(l_comp_size, g, &local_board);


    for(int gen = 0; gen < n; gen++){
        // step(local_m, g, local_board);
        
        //TODO send the ghost cells value to each process
    }

    //TODO gather the board from each process to the root process

    if(myrank == 0){
        // printBoard(m, g, board);
    }

    freeDoubleSquareArray(l_board_size, &local_board);
    freeDoubleSquareArray(board_size, &board);

    MPI_Finalize();
    return 0;
}


    // cout << "myrank: "<< myrank <<  " m: " << m << " n: " << n << " g: " << g << endl;  