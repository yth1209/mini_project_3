#include<stdio.h>
#include<stdlib.h>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

using namespace std;

const int PERTITION_TAG = 1;
const int MERGE_TAG = 2;

int npes, myrank;
int d; //number of processes in each dimension

int m; // size of board    
int n; // final generation
int g; // size of ghost cell

int padding_cells; //padding cells for the board to be divisible by d
int board_size; //size of the board after padding
int l_comp_size; //local computing size
int l_board_size; //local board size


inline int idx(int i, int j, int cols) {
    return i * cols + j;
}
inline pair<int,int> rank_ij(int rank, int d) {
    return {rank / d, rank % d};
}

void initSquareDoubleArray(int m, bool*** arr){
    (*arr) = new bool*[m];
    for(int i = 0; i < m; i++){
        (*arr)[i] = new bool[m];
        for(int j = 0; j < m; j++){
            (*arr)[i][j] = false;
        }
    }
}


void freeSquareDoubleArray(int m, bool*** arr){
    for(int i = 0; i < m; i++){
        delete[] (*arr)[i];
    }
    delete[] (*arr);
    *arr = nullptr;
}

void boardToLocalBoard(int m, int g, bool** local_board, int s_i, int s_j, bool*** board){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            (*local_board)[idx(i,j,m)] = (*board)[s_i+i][s_j+j];
        }
    }
}

void localBoardToBoard(int m, int g, bool** local_board, int s_i, int s_j, bool*** board){
    for(int i = g; i < m-g; i++){
        for(int j = g; j < m-g; j++){
            (*board)[s_i+i][s_j+j] = (*local_board)[idx(i,j,m)];
        }
    }
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
    for(int i = g; i < m - g; i++){
        for(int j = g; j < m - g; j++){
            cout << ((*board)[i][j] ? '#' : '.');
        }
        cout << endl;
    }
}

void printBoard(int m, int g, bool** board){
    for(int i = g; i < m - g; i++){
        for(int j = g; j < m - g; j++){
            cout << ((*board)[idx(i,j,m)] ? '#' : '.');
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

void initBoardAndLocalBoard(bool** local_board, bool*** board){
    *local_board = new bool[l_board_size * l_board_size]; //local board by 1D array (for efficent MPI send/recv)

    if(myrank == 0){
        initSquareDoubleArray(board_size, board);
        readBoard(m, g, board); // read board from txt input
        
        //Init the local board of the root process
        boardToLocalBoard(l_board_size, g, local_board, 0, 0, board);

        //SEND the partition of the board to each process
        for(int rank = 1; rank < npes; rank++){
            auto r_ij = rank_ij(rank, d);
            
            bool* send_data = new bool[l_board_size * l_board_size];
            boardToLocalBoard(l_board_size, g, &send_data, r_ij.first*l_comp_size, r_ij.second*l_comp_size, board);
            
            //SEND local board to each process
            MPI_Send(send_data, l_board_size*l_board_size, MPI_CXX_BOOL, rank, PERTITION_TAG, MPI_COMM_WORLD);
            delete[] send_data;
        }
    } else { 
        //RECIEIVE the partition of the board from the root process
        MPI_Status status;
        MPI_Recv(*local_board, l_board_size*l_board_size, MPI_CXX_BOOL, 0, PERTITION_TAG, MPI_COMM_WORLD, &status);
    }
}

void mergeLocalBoardToBoard(bool** local_board, bool*** board){
    //Send the local board result to the root process
    if(myrank == 0){
        localBoardToBoard(l_board_size, g, local_board, 0, 0, board);

        for(int rank = 1; rank < npes; rank++){
            auto r_ij = rank_ij(rank, d);

            MPI_Status merge_status;
            bool* recv_data = new bool[l_board_size * l_board_size];
            MPI_Recv(recv_data, l_board_size*l_board_size, MPI_CXX_BOOL, rank, MERGE_TAG, MPI_COMM_WORLD, &merge_status);
            
            localBoardToBoard(l_board_size, g,  &recv_data, r_ij.first*l_comp_size, r_ij.second*l_comp_size, board);
            
            delete[] recv_data;
        }

        printBoard(board_size-padding_cells, g, board);
        freeSquareDoubleArray(board_size, board);
    } else {
        MPI_Send(*local_board, l_board_size*l_board_size, MPI_CXX_BOOL, 0, MERGE_TAG, MPI_COMM_WORLD);
    }
}

int main(int argc, char* argv[]) {    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    d = sqrt(npes);
    if(d*d != npes) {
        cerr << "Number of processes must be a perfect square" << endl;
        return 1;
    }

    if(myrank == 0){
        cin >> m >> n >> g; // read m, n, g from txt input
    }
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g, 1, MPI_INT, 0, MPI_COMM_WORLD);

    padding_cells = m % d; //padding cells for the board to be divisible by d
    board_size = (m + padding_cells) + 2 * g; //padding with dimmension
    l_comp_size = (m + padding_cells) / d; //local computing size
    l_board_size = l_comp_size + 2 * g; //local board size
    
    bool** board; 
    bool* local_board; 
    
    
    initBoardAndLocalBoard(&local_board, &board);

    for(int gen = 0; gen < n; gen++){
        // step(local_m, g, local_board);
        //TODO send the ghost cells value to each process
    }

    mergeLocalBoardToBoard(&local_board, &board);

    delete[] local_board;

    MPI_Finalize();
    return 0;
}