#include<stdio.h>
#include<stdlib.h>
#include <mpi.h>
#include <iostream>
#include <vector>

using namespace std;

void readBoard(vector<vector<bool>>& board){
    int m = board.size();
    for(int i = 0; i < m; i++){
        char line[m];
        cin >> line;
        for(int j = 0; j < m; j++){
            board[i][j] = line[j] == '#';
        }
    }
}

void printBoard(vector<vector<bool>>& board){
    for(int i = 0; i < board.size(); i++){
        for(int j = 0; j < board[i].size(); j++){
            cout << (board[i][j] ? '#' : '.');
        }
        cout << endl;
    }
}

int getNeighCnt(int i, int j, vector<vector<bool>>& board) {
    int m = board.size();

    int cntNeigh = 0;
    if(i != 0) {
        cntNeigh += board[i-1][j]; 
        if(j != 0) cntNeigh += board[i-1][j-1];
        if(j!= m-1) cntNeigh += board[i-1][j+1];
    }
    if (i != m-1) {
        cntNeigh += board[i+1][j];
        if(j != 0) cntNeigh += board[i+1][j-1];
        if(j!= m-1) cntNeigh += board[i+1][j+1];
    }
    if(j != 0) cntNeigh += board[i][j-1];
    if(j != m-1) cntNeigh += board[i][j+1];
    
    return cntNeigh;
}

void step(vector<vector<bool>>& board){
    int m = board.size();

    vector<vector<bool>> newBoard(m, vector<bool>(m));

    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            int neighCnt = getNeighCnt(i, j, board);

            if(board[i][j]){
                if(neighCnt < 2 || neighCnt > 3) newBoard[i][j] = false;
                else newBoard[i][j] = true;
            } else {
                if(neighCnt == 3) newBoard[i][j] = true;
                else newBoard[i][j] = false;
            }
        }
    }
    board = newBoard;
}

int main(int argc, char* argv[]) {    
    
    int npes, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    int m; // size of board    
    int n; // final generation
    int g; // size of ghost cell
    vector<vector<bool>> board; //board

    if(myrank == 0){
        cin >> m >> n >> g; // read m, n, g from txt input
    }
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // cout << "myrank: "<< myrank <<  " m: " << m << " n: " << n << " g: " << g << endl;

    board = vector<vector<bool>> (m, vector<bool>(m));  

    if(myrank == 0){
        readBoard(board); // read board from txt input
    }; 

    //TODO send the partition of the board to each process
    
    for(int gen = 0; gen < n; gen++){
        step(board);
        //TODO send the ghost cells value to each process
    }

    //TODO gather the board from each process to the root process

    if(myrank == 0){
        printBoard(board);
    }

    MPI_Finalize();
    return 0;
}