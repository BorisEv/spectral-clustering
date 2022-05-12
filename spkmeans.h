#ifndef SPKMEANS_H_
#define SPKMEANS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*converting double 2d array nxm into char* (.4lf format) with comma after each number
 * and \n instead of comma at the end of row. Return the char**/
char* twoDArrToString(double** arr, int n, int m);
/*free memory of 2d array with n number of subarrays */
void arrMemFree(double** arr, int n);
/*initializing 2 darray of doubles of size nxm. Return the 2d array */
double** initArr(int n, int m);
/*multiplying matrix m1 x m2. Return the resulted matrix. It will be of size nxm */
double** matrixMult(double** m1, double** m2, int n, int m);
/*calculate and return transpose matrix of initial matrix nxm*/
double** transposeMat(double** A, int n, int m);
/* converting string of points to 2d double array nxm. The string includes n points
 * with k dimensions. Return the 2d array */
double** strToMatrix(char* str, int n, int k );
/* sorting matrix of eigenvalues and eigenvectors by eigenvalues */
void sortEigenMat(double** arr, int n);

/* return point dimension from the first row in the input file with name filename */
int getPointDimension(char* filename);
/* return number of points from input file */
int getPointsNumber(char* filename);
/* return 2d array of all points from file in which every row defines point */
double** getPointsArray(char* filename, int pointDimension, int pointsNumber);


/* forming and return the weighted adjacency matrix W from X which are in the */
double** wam(double** arrPoints, int pointDimension, int pointsNumber);
/* from given weighted adjacency matrix and it's size calculate and return the diagonal degree matrix */
double** ddg(double** weAdMat, int pointsNumber);
/* from given  weighted adjacency matrix and diagonal degree matrix
 *  calculate and return the the normalized laplacian matrix */
double** lnorm(double** weAdMat, double** diDeMat, int pointsNumber);
/* calculate and return eigenvalues and eigenvectors for symmetric matrix symMat */
double** getJacobiEigenValAndVec(double** symMat, int pointsNumber);

/* find and return the number of clusters k, while eigenVal,
 *  is already sorted array of eigenvalues of lnorm and n is their number */
int findK(double* eigenVal, int n);
/* return U matrix from sorted eigenvalues and eigenvectors. Size of U is n(=pointsNumber)xk */
double** getU(double** eigenValAndVec, int pointsNumber, int k);
/* return T matrix from U matrix. Size of T us n(=pointsNumber)xk */
double** getT(double** U, int pointsNumber, int k);

/* kmeans implementation, gets initial centroids and numbers as string, returns
 * final centroids as string */
char* kmeans(double epsilon, int k, int pointsNumber, int max_iter,char* centroids_str,char* numbers_str);

#endif