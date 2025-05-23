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
inline pair<int,int> rank_ij(int rank) {
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

void boardToLocalBoard(int m, bool** local_board, int s_i, int s_j, bool*** board){
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


int getNeighCnt(int i, int j, int board_size, bool** board) {
    int cntNeigh = 0;

    cntNeigh += (*board)[idx(i, j-1, board_size)];
    cntNeigh += (*board)[idx(i, j+1, board_size)];
    cntNeigh += (*board)[idx(i-1, j, board_size)]; 
    cntNeigh += (*board)[idx(i-1, j-1, board_size)];
    cntNeigh += (*board)[idx(i-1, j+1, board_size)];
    cntNeigh += (*board)[idx(i+1, j, board_size)];
    cntNeigh += (*board)[idx(i+1, j-1, board_size)];
    cntNeigh += (*board)[idx(i+1, j+1, board_size)];
    
    return cntNeigh;
}

bool isPaddingCell(int i, int j){
    if(padding_cells == 0) return false; //no padding cells
    auto r_ij = rank_ij(myrank);
    int board_i = i + (r_ij.first * l_comp_size);
    int board_j = j + (r_ij.second * l_comp_size);
    
    if(board_i >= m+g || board_j >= m+g) return true; 
    else return false;
}

void step(bool** local_board){
    vector<pair<int,int>> changed_cells;

    for(int i = g; i < l_comp_size + g; i++){
        for(int j = g; j < l_comp_size + g; j++){
            if(isPaddingCell(i,j)) continue; //skip padding cells
            
            int neighCnt = getNeighCnt(i, j, l_board_size, local_board);

            if((*local_board)[idx(i,j,l_board_size)]){
                if(neighCnt < 2 || neighCnt > 3) changed_cells.push_back({i,j});
            } else {
                if(neighCnt == 3) changed_cells.push_back({i,j});
            }
        }
    }

    for(auto cell: changed_cells){
        (*local_board)[idx(cell.first, cell.second, l_board_size)] = !(*local_board)[idx(cell.first, cell.second, l_board_size)]; 
    }
}

void sendRecvGhostCells(bool** local_board) {
    MPI_Status status;
    MPI_Comm cart_comm;
    int dims[2] = {d, d};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, /*reorder=*/0, &cart_comm);

    int coords[2];
    MPI_Cart_coords(cart_comm, myrank, 2, coords);
    int north, south, west, east;
    MPI_Cart_shift(cart_comm, /*dim=*/0, /*disp=*/1, &north, &south);
    MPI_Cart_shift(cart_comm, /*dim=*/1, /*disp=*/1, &west, &east);

    MPI_Datatype row_type, col_type, corner_type;
    MPI_Type_vector(g, l_comp_size, l_board_size, MPI_CXX_BOOL, &row_type);
    MPI_Type_commit(&row_type);
    MPI_Type_vector(l_comp_size, g, l_board_size, MPI_CXX_BOOL, &col_type);
    MPI_Type_commit(&col_type);
    MPI_Type_vector(g, g, l_board_size, MPI_CXX_BOOL, &corner_type);
    MPI_Type_commit(&corner_type);

    // 상하 방향 교환
    MPI_Sendrecv(
        &(*local_board)[idx(g, g, l_board_size)], 1, row_type, north, 0,
        &(*local_board)[idx(g + l_comp_size, g, l_board_size)], 1, row_type, south, 0,
        cart_comm, &status);
    MPI_Sendrecv(
        &(*local_board)[idx(g + l_comp_size - g, g, l_board_size)], 1, row_type, south, 1,
        &(*local_board)[idx(0, g, l_board_size)], 1, row_type, north, 1,
        cart_comm, &status);

    // 좌우 방향 교환
    MPI_Sendrecv(
        &(*local_board)[idx(g, g, l_board_size)], 1, col_type, west, 2,
        &(*local_board)[idx(g, g + l_comp_size, l_board_size)], 1, col_type, east, 2,
        cart_comm, &status);
    MPI_Sendrecv(
        &(*local_board)[idx(g, g + l_comp_size - g, l_board_size)], 1, col_type, east, 3,
        &(*local_board)[idx(g, 0, l_board_size)], 1, col_type, west, 3,
        cart_comm, &status);

    // 대각선 교환
    int ci = coords[0], cj = coords[1];
    int nw = (ci > 0 && cj > 0) ? (ci - 1) * d + (cj - 1) : MPI_PROC_NULL;
    int ne = (ci > 0 && cj < d - 1) ? (ci - 1) * d + (cj + 1) : MPI_PROC_NULL;
    int sw = (ci < d - 1 && cj > 0) ? (ci + 1) * d + (cj - 1) : MPI_PROC_NULL;
    int se = (ci < d - 1 && cj < d - 1) ? (ci + 1) * d + (cj + 1) : MPI_PROC_NULL;

    MPI_Sendrecv(
        &(*local_board)[idx(g, g, l_board_size)], 1, corner_type, nw, 4,
        &(*local_board)[idx(g + l_comp_size, g + l_comp_size, l_board_size)], 1, corner_type, se, 4,
        cart_comm, &status);
    MPI_Sendrecv(
        &(*local_board)[idx(g, g + l_comp_size - g, l_board_size)], 1, corner_type, ne, 5,
        &(*local_board)[idx(g + l_comp_size, 0, l_board_size)], 1, corner_type, sw, 5,
        cart_comm, &status);

    MPI_Type_free(&row_type);
    MPI_Type_free(&col_type);
    MPI_Type_free(&corner_type);
    MPI_Comm_free(&cart_comm);
}

void initBoardAndLocalBoard(bool** local_board, bool*** board){
    *local_board = new bool[l_board_size * l_board_size]; //local board by 1D array (for efficent MPI send/recv)

    if(myrank == 0){
        initSquareDoubleArray(board_size, board);
        readBoard(m, g, board); // read board from txt input
        
        //Init the local board of the root process
        boardToLocalBoard(l_board_size, local_board, 0, 0, board);

        //SEND the partition of the board to each process
        for(int rank = 1; rank < npes; rank++){
            auto r_ij = rank_ij(rank);
            
            bool* send_data = new bool[l_board_size * l_board_size];
            boardToLocalBoard(l_board_size, &send_data, r_ij.first*l_comp_size, r_ij.second*l_comp_size, board);
            
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
            auto r_ij = rank_ij(rank);

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
        //game of life local computation step
        step(&local_board);
        
        //send the ghost cells value to each process
        sendRecvGhostCells(&local_board);
    }

    mergeLocalBoardToBoard(&local_board, &board);

    delete[] local_board;

    MPI_Finalize();
    return 0;
}