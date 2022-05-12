#include "spkmeans.h"

/*converting double 2d array nxm into char* (.4lf format) with comma after each number
 * and \n instead of comma at the end of row. Return the char* */
char* twoDArrToString(double** arr, int n, int m){
	/*initializing needed variables */
	int i;
	int j;
	char* str = (char*)malloc(15*n*m*sizeof(char));
	char* helper =(char*)malloc(15*n*m*sizeof(char));
	if (str == NULL || helper==NULL){
	        printf("An Error Has Occurred");
	        exit(1);
	    }

    /*adding helper+nextNumber to str and then copying str to helper due to compilation
     * warning when calling sprinf on the same variable*/
	for(i = 0; i<n;i++){
		for(j = 0; j<m; j++){
			if(i == 0 && j== 0){
				sprintf(helper, "%.4f,", arr[0][0]);
			}
			else{
				sprintf(str, "%s%.4f", helper, arr[i][j]);
				strcpy(helper,str);
				if(j<m-1){ /* adding comma after each number besides the last in row */
					sprintf(str,"%s,", helper);
				}
				strcpy(helper, str);
			}
		}
		sprintf(str,"%s\n", helper); /* adding \n after each row */
		strcpy(helper, str);
	}
	/* free memory */
	free(helper);
	return str;
}

/* free memory of 2d array with n number of subArrays. Return the 2d array */
void arrMemFree(double** arr, int n){
	int i;
	for(i = 0; i<n;i++){
		free(arr[i]);
	}
	free(arr);
}

/* initializing 2d array of doubles of size nxm. Return the 2d array */
double** initArr(int n, int m){
	int i;
	double** arr = (double**)malloc(n*sizeof(double*));
    if (arr == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
	for(i = 0; i<n; i++){
		arr[i] = (double*)malloc(m*sizeof(double));
	    if (arr[i] == NULL){
	        printf("An Error Has Occurred");
	        exit(1);
	    }
	}
	return arr;
}

/*multiplying matrix m1 x m2. Return the resulted matrix. It will be of size nxm */
double** matrixMult(double** m1, double** m2, int n, int m){
	int i,j,z;
	double** mRes;
	mRes = initArr(n,m);
    /* filling the matrix */
    for(i = 0; i<n; i++){
    	for(j = 0; j<m; j++){
    		mRes[i][j] = 0;
    		for(z = 0; z < n; z++){
    			mRes[i][j] += m1[i][z]*m2[z][j];
    		}
    	}
    }
	return mRes;
}

/*calculate and return transpose matrix of initial matrix nxm*/
double** transposeMat(double** A, int n, int m){
	int i,j;
	double** AT = initArr(m, n);
	for(i = 0; i<m; i++){
		for(j = 0; j<n; j++){
			AT[i][j] = A[j][i];
		}
	}
	return AT;
}

/* converting string of points to 2d double array nxm. The string includes n points
 * with k dimensions. Return the 2d array */
double** strToMatrix(char* str, int n, int k ){
	/* initializing needed variables */
	int i,j;
	double** res = initArr(n,k);
	char* ptr = str-1; /*we use ptr+1 every time to get rid of commas, so initializing it with -1 */
	for(i = 0; i<n;i++){
		for(j=0;j<k;j++){
			res[i][j] = strtod(ptr+1,&ptr);
		}
	}
	return res;
}

/* sorting matrix of eigenvalues and eigenvectors by eigenvalues */
void sortEigenMat(double** arr, int n){
    /* initializing needed variables */
	int i,j,z;
	double temp;
	double* temp2 = (double*)malloc(n*sizeof(double));

	for (i = 0; i < n; ++i)
	{
		for (j = i + 1; j < n; ++j)
		{
			if (arr[0][i] > arr[0][j])
			{
				/*swapping the values*/
				temp = arr[0][i];
				arr[0][i] = arr[0][j];
				arr[0][j] = temp;
				/*swapping the vectors */
				for(z = 0; z<n;z++){
					temp2[z] = arr[z+1][i];
					arr[z+1][i] = arr[z+1][j];
					arr[z+1][j] = temp2[z];
				}
			}
		}
	}
	free(temp2);
}


/* return point dimension from the first row in the input file with name filename */
int getPointDimension(char* filename) {
	/*initializing the needed variables */
	int i;
	char line[7000];
	char comma = ',';
	int counter = 0;
	FILE* file = NULL;
	file = fopen(filename, "r");
	if(file == NULL){
		printf("Invalid Input!");
	 	exit(1);
	    }
	/* counting numbers of commas in the first row +1*/
	fscanf(file,"%s",line);
	counter++;
	for (i = 0; i < (int)strlen(line); i++){
		if(line[i] == comma){
			counter++;
		}
	}
	fclose(file);
    return counter;
}

/* return number of points from input file */
int getPointsNumber(char* filename){
	/* initializing needed variables */
	int pointsNumber = 0;
	char line[7000];
	FILE* file = NULL;
	file = fopen(filename, "r");
	if(file == NULL){
		printf("Invalid Input!");
	 	exit(1);
	    }
	/* counting number of lines in the file, which is number of points */
	while (1){
		fscanf(file,"%s",line);
		if (feof(file))
			break;
		pointsNumber++;
	}
	fclose(file);
	return pointsNumber;
}

/* return 2d array of all points from file in which every row defines point */
double** getPointsArray(char* filename, int pointDimension, int pointsNumber){
	/* initializing needed variables */
	int i,j;
	double** arrPoints;
	char line[7000];
	char *ptr;
	FILE* file = NULL;
	file = fopen(filename, "r");
	if(file == NULL){
		printf("Invalid Input!");
	 	exit(1);
	    }
	arrPoints = initArr(pointsNumber, pointDimension);

	/* for each double* writing a point to it */
	for(i = 0; i<pointsNumber;i++){
		fscanf(file, "%s", line);
	    arrPoints[i][0] = strtod(line, &ptr);
	    /* for each dimension writing number into point */
	    for(j = 1; j< pointDimension; j++){
	    	arrPoints[i][j] = strtod(ptr+1, &ptr);
	    }
	}
	fclose(file);
	return arrPoints;
}



/* forming and return the weighted adjacency matrix W from X which are in the */
double** wam(double** arrPoints, int pointDimension, int pointsNumber){
	/* initializing needed variables */
	int i,j,z;
	double norm; /* norm we will calculate in each iteration */
	double** weAdMat;
	weAdMat = initArr(pointsNumber, pointsNumber);

	/* filling the matrix */
	for(i = 0; i<pointsNumber; i++){
		for(j = 0; j<pointsNumber; j++){
			if(i == j){
				weAdMat[i][j] = 0;
			}
			else{
				/* calculating the norm */
				norm = 0;
				for(z = 0; z<pointDimension; z++){
					norm += pow((arrPoints[i][z]-arrPoints[j][z]),2);
				}
				norm = sqrt(norm);
				weAdMat[i][j] = exp(norm*-0.5);
			}
		}
	}
	return weAdMat;
}


/* from given weighted adjacency matrix and it's size calculate and return the diagonal degree matrix */
double** ddg(double** weAdMat, int pointsNumber){
	/* initializing needed variables */
	int i,j,z;
	double wRowSum; /* sum of i'th row in W matrix, which will be calculated */
	double** diDeMat;
	diDeMat = initArr(pointsNumber, pointsNumber);

	/* filling the matrix */
	for(i = 0; i<pointsNumber; i++){
		for(j = 0; j<pointsNumber; j++){
			if(i != j){
				diDeMat[i][j] = 0;
			}
			else{
				/* calculating the sum of i'th row in weighted adjacency matrix */
				wRowSum = 0;
				for(z = 0; z<pointsNumber; z++){
					wRowSum += weAdMat[i][z];
				}
				diDeMat[i][j] = wRowSum;
			}
		}
	}
	return diDeMat;
}


/* from given  weighted adjacency matrix and diagonal degree matrix
 *  calculate and return the the normalized laplacian matrix */
double** lnorm(double** weAdMat, double** diDeMat, int pointsNumber){
	/*initializing needed variables */
	int i,j;
	double** sqrtD; /*D^-0.5 matrix */
	double** DW; /* multiplication of matrix D by matrix W */
	double** DWD; /* multiplication of matrix DW by matrix D */
	double** normLapMat;
	sqrtD = initArr(pointsNumber, pointsNumber);
	normLapMat = initArr(pointsNumber, pointsNumber);

	/* filling D^-0.5 matrix */
	for(i = 0; i<pointsNumber; i++){
		for(j = 0; j<pointsNumber; j++){
			if(i == j){
				sqrtD[i][j] = 1/sqrt(diDeMat[i][j]);
			}
			else{
				sqrtD[i][j] = diDeMat[i][j];
			}
		}
	}

	DW = matrixMult(sqrtD, weAdMat, pointsNumber, pointsNumber); /* first matrix multiplication */
	DWD = matrixMult(DW, sqrtD, pointsNumber, pointsNumber); /* second matrix multiplication */

	/* filling lnorm matrix */
	for(i = 0; i<pointsNumber; i++){
		for(j = 0; j< pointsNumber; j++){
			if(i == j){
				normLapMat[i][j] = 1-DWD[i][j];
			}
			else{
				normLapMat[i][j] = 0-DWD[i][j];
			}
		}
	}

	/* free memory */
	arrMemFree(DW, pointsNumber);
	arrMemFree(DWD, pointsNumber);
	arrMemFree(sqrtD, pointsNumber);
	return normLapMat;
}


/* calculate and return eigenvalues and eigenvectors for symmetric matrix symMat */
double** getJacobiEigenValAndVec(double** symMat, int pointsNumber){
	/*initializing needed variables */
	int i, j;
	double** eigenValAndVec=  initArr(pointsNumber+1,pointsNumber); /*the first row will be the eigen values and then the V matrix*/
	double** A = initArr(pointsNumber, pointsNumber);
	double** P = initArr(pointsNumber, pointsNumber); /* P matrix */
	double** PT = initArr(pointsNumber, pointsNumber);/* P transpose matrix */
	double** PTA = initArr(pointsNumber, pointsNumber); /* multiplying Ptranspose by A */
	double** V = initArr(pointsNumber, pointsNumber); /* eigenvectors matrix */
	double** Vold; /* to avoid memory leak while multiplying matrices P */
	double offA = 0;
	double offAtag = 0;
	int iter = 0;
	int max_iter = 100;
	double epsilon = pow(10,-5);
	double maxAij; /* pivot element */
	int maxI = 0; /* i index of pivot element */
	int maxJ = 0; /* j index of pivot element */
	double theta;
	char sign;
	double t;
	double c;
	double s;
	/* copying the normLapMat to A, so we dont lose it in iterations */
	for(i = 0; i<pointsNumber; i++){
		for(j = 0; j<pointsNumber; j++){
			A[i][j] = symMat[i][j];
		}
	}
	/*find offA ( the sum of squares of all off-diagonal elements of A */
	for(i =0; i<pointsNumber; i++){
		for(j =0; j<pointsNumber; j++){
			if(i != j){
				offAtag += pow(A[i][j],2);
			}
		}
	}

	/* the main iterations that are counted by iter */
	do{
		offA = offAtag;
		iter++;

		/* find Aij and corresponded i, j */
		maxAij = 0;
		for(i = 0; i<pointsNumber; i++){
			for(j = i+1; j< pointsNumber; j++){
				if(fabs(A[i][j]) >maxAij){
					maxAij = fabs(A[i][j]);
					maxI = i;
					maxJ = j;
				}
			}
		}

		/* calculating c,s */
		theta = (A[maxJ][maxJ]-A[maxI][maxI])/(2*A[maxI][maxJ]);
		sign = (theta<0) ? '-':'+';
		t = 1/(fabs(theta)+sqrt(pow(theta,2)+1));
		if(sign == '-'){
			t = t*-1;
		}
		c = 1/(sqrt(pow(t,2)+1));
		s = t*c;

		/* creating the matrix P */
		arrMemFree(P, pointsNumber);
		P = initArr(pointsNumber, pointsNumber);
		for(i = 0; i<pointsNumber; i++){
			for(j = 0; j<pointsNumber; j++){
				if(i == j){
					if(i==maxI || j == maxJ){
						P[i][j] = c;
					}
					else{
						P[i][j] = 1;
					}
				}
				else{
					if(i == maxI && j == maxJ){
						P[i][j] = s;
					}
					else if(i == maxJ && j== maxI){
						P[i][j] = -1*s;
					}
					else{
						P[i][j] =0;
					}
				}
			}
		}
		/* if it is first iteration - initializing V=P. Else multiplying V=VP */
		if(iter == 1){
			for(i = 0; i<pointsNumber; i++){
				for(j = 0; j<pointsNumber; j++){
					V[i][j] = P[i][j];
				}
			}
		}else{
			Vold = initArr(pointsNumber,pointsNumber);
			/* copying V to Vold */
			for(i = 0; i< pointsNumber; i++){
				for(j =0; j<pointsNumber; j++){
					Vold[i][j] =V[i][j];
				}
			}
			arrMemFree(V, pointsNumber);
			V = matrixMult(Vold, P, pointsNumber, pointsNumber);
			arrMemFree(Vold, pointsNumber);
		}
		/* calculating A' matrix, which we save in A */
		arrMemFree(PT, pointsNumber);
		PT = transposeMat(P, pointsNumber, pointsNumber);
		arrMemFree(PTA, pointsNumber);
		PTA = matrixMult(PT, A, pointsNumber, pointsNumber);
		arrMemFree(A, pointsNumber);
		A = matrixMult(PTA, P,pointsNumber, pointsNumber);
		/* calculating new off(A') */
		offAtag = 0;
		for(i =0; i<pointsNumber; i++){
			for(j =0; j<pointsNumber; j++){
				if(i != j){
					offAtag += pow(A[i][j],2);
				}
			}
		}
	} while(iter<max_iter && (offA-offAtag)>epsilon);

	/* first row of eigenValAndVec are eigenvalues which are the diagonal elements in A */
	for(i = 0; i<pointsNumber; i++){
		eigenValAndVec[0][i] = A[i][i];
	}
	/* eigenvectors are in matrix V, so we write them into eigenValAndVec */
	for(i = 1; i<pointsNumber+1;i++){
		for(j = 0; j<pointsNumber; j++){
			eigenValAndVec[i][j] = V[i-1][j];
		}
	}

	/*free memory*/
	arrMemFree(A, pointsNumber);
	arrMemFree(P, pointsNumber);
	arrMemFree(PT, pointsNumber);
	arrMemFree(PTA, pointsNumber);
	arrMemFree(V, pointsNumber);
	return eigenValAndVec;
}

/* find and return the number of clusters k, while eigenVal,
 *  is already sorted array of eigenvalues of lnorm and n is their number */
int findK(double* eigenVal, int n){
	int i;
	int k;
	double max;
	k = 0;
	max = 0;
	/*getting i of maximum abs(lambda_i-lambda_(i+1)) */
	for (i = 0; i< n/2; i++){
		if(fabs(eigenVal[i]-eigenVal[i+1])>max){
			max = fabs(eigenVal[i]-eigenVal[i+1]);
			k = i+1; /* becauese vector i is in place i+1 in eigenVal array*/
		}
	}
	return k;
}

/* return U matrix from sorted eigenvalues and eigenvectors. Size of U is n(=pointsNumber)xk */
double** getU(double** eigenValAndVec, int pointsNumber, int k){
	/* initializing needed variables */
	int i,j;
	double** U = initArr(pointsNumber,k);
	/* filling the U matrix */
	for(i = 0; i<pointsNumber;i++){
		for(j = 0; j<k; j++){
			U[i][j] = eigenValAndVec[i+1][j];
		}
	}
	return U;
}

/* return T matrix from U matrix. Size of T us n(=pointsNumber)xk */
double** getT(double** U, int pointsNumber, int k){
	/* initializing needed variables */
	int i, j;
	double normUi;
	double** T = initArr(pointsNumber, k);
	/*filling T matrix */
	for(i = 0; i<pointsNumber;i++){
		/* calculating norm for row number i */
		normUi =0;
		for(j = 0; j<k;j++){
			normUi += pow(U[i][j],2);
		}
		normUi = sqrt(normUi);
		/* filling the i row */
		for(j = 0; j<k; j++){
			if(normUi == 0){
				T[i][j] = 0;
			}
			else{
				T[i][j] = U[i][j]/normUi;
			}
		}
	}
	return T;
}

/* kmeans implementation, gets initial centroids and numbers as string, returns
 * final centroids as string */
char* kmeans(double epsilon, int k, int pointsNumber, int max_iter,char* centroids_str,char* numbers_str){
    /* initializing needed variables*/
    int i, j,z;
    /* converting strings to 2d arrays of doubles*/
    double** centroids = strToMatrix(centroids_str,k,k);
    double** numbers = strToMatrix(numbers_str,pointsNumber,k);
    double** clusters = initArr(k,pointsNumber); /* clusters[k*j+i] =0 if no numbers[i] in cluster k */
    char* final_centroids; /* the resulted string */
    /* to check convergence */
    double allCentNorm;
	double prevAllCentNorm;
	double centroidNorm;
	int counter = 0; /* for updating centroids (counter = number of points in cluster */
	/* for finding argmin and jargmin, which is index j of cluster with minimal argument */
	double argmin;
	double tempoarg;
	int jargmin;

    int iter = 0; /*counter of iterations */
    /* the main iterations */
    do{
    	iter++;
    	/* making new clusters */
    	for(i=0;i<k;i++){
    		for(j = 0; j< pointsNumber;j++){
    		    clusters[i][j] = 0;
    		}
    	}

    	/* filling the clusters by finding minimum arg */
    	for(i=0; i<pointsNumber; i++){
    	    /* calculate arg(x_i-mu_j) for the first centroids */
    		argmin = 0;
    		for(z = 0; z<k; z++){
    			argmin += pow((numbers[i][z]-centroids[0][z]),2.0);
    		}
    		jargmin = 0;
    		/* calculating the argmin and jargmin, the index j of centroid with minimal arg */
    		for(j=1;j<k;j++){
    			tempoarg = 0;
    			for(z = 0; z<k; z++){
    			    tempoarg += pow((numbers[i][z]-centroids[j][z]),2.0);
    			}
    			if (tempoarg<argmin){
    				argmin = tempoarg;
    				jargmin = j;
    			}
    		}
    		/* updating cluster */
    		clusters[jargmin][i] = 1;
    	}
    	/* calculate the norm before we change the centroids */
    	prevAllCentNorm = 0;
    	for(j = 0; j<k; j++){
    		centroidNorm = 0;
    		for (z = 0; z<k; z++){
        		centroidNorm += pow(centroids[j][z],2.0);
    		}
    		centroidNorm = sqrt(centroidNorm);
    		prevAllCentNorm += pow(centroidNorm,2);
    	}
    	prevAllCentNorm = sqrt(prevAllCentNorm);
    	/* updating the centroids */
    	for(j=0; j<k;j++){
    		counter = 0;
        	for(z = 0; z<k; z++){
        		centroids[j][z] = 0;
        	}
    		for(i=0; i<pointsNumber; i++){
    			if(clusters[j][i] != 0){
    		    	for( z = 0; z<k; z++){
    		    		centroids[j][z] += numbers[i][z];
    		    	}
        			counter++;
    			}
    		}
    		if(counter == 0){ /* to avoid division by zero, but should never happen */
    		    printf("An Error Has Occurred");
    		    exit(1);
    		}
    		else{
    	    	for(z = 0; z<k; z++){
    	    		centroids[j][z] = centroids[j][z]/counter;
    	    	}
    		}
    	}
    	/* calculating new norm */
    	allCentNorm = 0;
    	for(j=0; j<k; j++){
    		centroidNorm = 0;
    		for ( z = 0; z<k; z++){
        		centroidNorm += pow(centroids[j][z],2.0);
    		}
    		centroidNorm = sqrt(centroidNorm);
    		allCentNorm += pow(centroidNorm,2);
    	}
    	allCentNorm = sqrt(allCentNorm);
    } while(iter<max_iter&&(fabs(allCentNorm-prevAllCentNorm)>epsilon));
    final_centroids = twoDArrToString(centroids,k,k); /* converting final centroids to string */
    /* free memory */
    arrMemFree(centroids,k);
    arrMemFree(numbers, pointsNumber);
    arrMemFree(clusters, k);
    return final_centroids;
}

int main(int argc, char** argv){
	/* initializing needed variables from the arguments*/
	int i,j;
	char* goal = argv[1];
	char* filename = argv[2];
	double** arrPoints;
	/*creating array of all data points */

	int pointDimension = getPointDimension(filename);
	int pointsNumber = getPointsNumber(filename);
	arrPoints = getPointsArray(filename, pointDimension, pointsNumber);
	(void)argc; /*so the compiler won't warn about unused parameter */

	if(strcmp(goal, "wam") == 0){
		/* in case of wam calculating and printing the weighted adjacency matrix */
		double** weAdMat = wam(arrPoints, pointDimension, pointsNumber);
		for(i = 0; i<pointsNumber; i++){
			for(j = 0; j<pointsNumber; j++){
				printf("%.4f", weAdMat[i][j]);
				if(j != pointsNumber-1){
					printf(",");
				}
			}
			printf("\n");
		}
		/* free memory */
		arrMemFree(weAdMat, pointsNumber);
	}

	else if (strcmp(goal, "ddg") == 0){
		/* in case of ddg calculating the weighted adjacency matrix
		 *  and from it calculating and printing the diagonal degree matrix*/
		double** weAdMat = wam(arrPoints, pointDimension, pointsNumber);
		double** diDeMat = ddg(weAdMat, pointsNumber);
		for(i = 0; i<pointsNumber; i++){
			for(j = 0; j<pointsNumber; j++){
				printf("%.4f", diDeMat[i][j]);
				if(j != pointsNumber-1){
					printf(",");
				}
			}
			printf("\n");
		}
		/* free memory */
		arrMemFree(weAdMat, pointsNumber);
		arrMemFree(diDeMat, pointsNumber);
	}
	else if (strcmp(goal, "lnorm") == 0){
		/* in case of lnorm calculating the weighted adjacency matrix,
		 * from it calculating the diagonal degree matrix and from both of them,
		 * calculating and printing the laplacian normalized matrix*/
		double** weAdMat = wam(arrPoints, pointDimension, pointsNumber);
		double** diDeMat = ddg(weAdMat, pointsNumber);
		double** normLapMat = lnorm(weAdMat, diDeMat, pointsNumber);
		for(i = 0; i<pointsNumber; i++){
			for(j = 0; j<pointsNumber; j++){
				printf("%.4f", normLapMat[i][j]);
				if(j != pointsNumber-1){
					printf(",");
				}
			}
			printf("\n");
		}
		/* free memory */
		arrMemFree(weAdMat, pointsNumber);
		arrMemFree(diDeMat, pointsNumber);
		arrMemFree(normLapMat, pointsNumber);
	}
	else if(strcmp(goal, "jacobi") == 0){
		/* in case of jacobi, according to what is written on the forum,
		 * in the we have a symmetric matrix, not date points,
		 * from which we find eigenvectors and eigenvalues and print them */
		double** eigenValAndVec = getJacobiEigenValAndVec(arrPoints, pointsNumber);
		for(i = 0; i<pointsNumber+1; i++){
			for(j = 0; j<pointsNumber; j++){
				if(fabs(eigenValAndVec[i][j])<0.00005){ /* dealing with -0.0000 */
					printf("%.4f", 0.);
				}
				else{
					printf("%.4f", eigenValAndVec[i][j]);
				}
				if(j != pointsNumber-1){
					printf(",");
				}
			}
			printf("\n");
		}
		/* free memory */
		arrMemFree(eigenValAndVec, pointsNumber+1);
	}
	/* free memory */
	arrMemFree(arrPoints, pointsNumber);
	return 0;
}
