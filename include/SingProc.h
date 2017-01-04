#ifndef SINGPROC_H
#define SINGPROC_H
#ifdef USE_MPI
#include <mpi.h>
///Basics class to start a MPI grid
class SingProc{
 private:
  int *Dim;
  int *Period;
  int *Remain;
  int *PCoord;
 public:
  SingProc(int ExtSize,int ExtRank){
    Rank = ExtRank;
    Size = ExtSize;
    NDim = 1;
    Dim = (int *)calloc(NDim,sizeof(int));
    Period = (int *)calloc(NDim,sizeof(int));
    Remain = (int *)calloc(NDim,sizeof(int));
    PCoord = (int *)calloc(NDim,sizeof(int));
    Dim[0] = Size;
    Period[0] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, NDim, Dim, Period, 1, &CommGrid);
    Remain[0] = 1;
    MPI_Cart_sub(CommGrid,Remain,&CommRow);
    Remain[0] = 0;
    MPI_Cart_sub(CommGrid,Remain,&CommCol);
    MPI_Cart_coords(CommGrid,Rank,NDim,PCoord);
    Col = PCoord[0];
    Row = 0;
    //printf("Proc num %d di %d in (%d,%d)\n",Rank,Size,Row,Col);
  }
  SingProc(int ExtSize,int ExtRank,int ExtNRow,int ExtNCol){
    Rank = ExtRank;
    Size = ExtSize;
    NDim = 2;
    Dim = (int *)calloc(NDim,sizeof(int));
    Period = (int *)calloc(NDim,sizeof(int));
    Remain = (int *)calloc(NDim,sizeof(int));
    PCoord = (int *)calloc(NDim,sizeof(int));
    Dim[0] = ExtNRow; Dim[1] = ExtNCol;
    Period[0] = 0;Period[0] = 0;
    //if(NRow*NCol != Size){printf("The number of processor must be a square of a number: n^2\n");exit(0);}
    MPI_Cart_create(MPI_COMM_WORLD, NDim, Dim, Period, 1, &CommGrid);
    Remain[0] = 1;Remain[1] = 0;
    MPI_Cart_sub(CommGrid,Remain,&CommRow);
    Remain[0] = 0;Remain[1] = 1;
    MPI_Cart_sub(CommGrid,Remain,&CommCol);
    MPI_Cart_coords(CommGrid,Rank,NDim,PCoord);
    Row=PCoord[0];
    Col=PCoord[1];
    //    printf("Proc num %d di %d in (%d,%d)\n",Rank,Size,Row,Col);
  }
  int Rank;
  int Size;
  MPI_Comm CommGrid;
  MPI_Comm CommRow;
  MPI_Comm CommCol;
  int Col;
  int Row;
  int NDim;
};
#endif//USE_MPI

#endif//SINGPROC
