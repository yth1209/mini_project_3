#include<stdio.h>
#include<stdlib.h>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <limits>

using namespace std;

void readBoard(int m, int g, vector<vector<bool>>& board){
    string line;
    for(int i = g; i < m+g; i++){
        cin >> line; 
        for(int j = g; j < m+g; j++){
            board[i][j] = (line[j - g] == '#');
        }
    }
}

void printBoard(int m, int g, vector<vector<bool>>& board){
    for(int i = g; i < m + g; i++){
        for(int j = g; j < m + g; j++){
            cout << (board[i][j] ? '#' : '.');
        }
        cout << endl;
    }
    cout << flush;
}

int getNeighCnt(int i, int j, vector<vector<bool>>& board) {
    int m = board.size();

    int cntNeigh = 0;

    cntNeigh += board[i][j-1];
    cntNeigh += board[i][j+1];
    
    cntNeigh += board[i-1][j]; 
    cntNeigh += board[i-1][j-1];
    cntNeigh += board[i-1][j+1];
 
    cntNeigh += board[i+1][j];
    cntNeigh += board[i+1][j-1];
    cntNeigh += board[i+1][j+1];
    
    return cntNeigh;
}

void step(int g, vector<vector<bool>>& board){
    vector<pair<int,int>> changed_cells;

    for(int i = g; i < board.size()-g; i++){
        for(int j = g; j < board.size()-g; j++){
            int neighCnt = getNeighCnt(i, j, board);

            if(board[i][j]){
                if(neighCnt < 2 || neighCnt > 3) changed_cells.push_back({i,j});
            } else {
                if(neighCnt == 3) changed_cells.push_back({i,j});
            }
        }
    }

    for(auto cell: changed_cells){
        board[cell.first][cell.second] = !board[cell.first][cell.second]; 
    }
}

int main(int argc, char* argv[]) {    
    
    int npes, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    int m; // size of board    
    int n; // final generation
    int g; // size of ghost cell
    vector<vector<bool>> board; 
    vector<vector<bool>> local_board; 

    if(myrank == 0){
        cin >> m >> n >> g; // read m, n, g from txt input
    }
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(myrank == 0){
        board = vector<vector<bool>>(m+2*g, vector<bool>(m+2*g, false));
        readBoard(m, g, board); // read board from txt input
    }; 

    //TODO send the partition of the board to each process
    //Partitioned board will be a 2D vector of size (m/root(npes) + 2g)

    for(int gen = 0; gen < n; gen++){
        step(g, board);
        
        //TODO send the ghost cells value to each process
    }

    //TODO gather the board from each process to the root process

    if(myrank == 0){
        printBoard(m, g, board);
    }

    MPI_Finalize();
    return 0;
}


    // cout << "myrank: "<< myrank <<  " m: " << m << " n: " << n << " g: " << g << endl;  