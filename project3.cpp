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
    int m; // size of board    
    int n; // final generation
    int g; // size of ghost cell
    cin >> m >> n >> g; // read m, n, g from txt input

    vector<vector<bool>> board(m, vector<bool>(m)); //board
    readBoard(board);


    for(int gen = 0; gen < n; gen++){
        step(board);
    }

    printBoard(board);
    return 0;
}