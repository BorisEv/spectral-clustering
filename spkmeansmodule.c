#include <Python.h>
#include "spkmeans.c"

/* main function, which we will call in all cases except of kmeans c implementation
 * besides goals "wam", "ddg","lnorm" and "jacobi", it has goal "findT", in this case we
 * calculate the T matrix */
static PyObject* main_capi(PyObject *self, PyObject *args)
{
    /* initializing needed variables, which we get from python */
    char* goal;
    char* filename;
    int k;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "ssi",&goal, &filename, &k)) {
        return NULL; /* If an error has occurred */
    }

	/*creating array of all points */
    double** arrPoints;
	int pointDimension = getPointDimension(filename);
	int pointsNumber = getPointsNumber(filename);
	if(k>pointsNumber){ /* if given k  bigger than N */
		printf("Invalid Input!");
	 	exit(1);
	}
	arrPoints = getPointsArray(filename, pointDimension, pointsNumber);

    char* res = "hi"; /* we wil save here the result of each call and return it as string to python */
	if(strcmp(goal, "wam") == 0){

		/* in case of wam calculating and converting to string the weighted adjacency matrix */
		double** weAdMat = wam(arrPoints, pointDimension, pointsNumber);
        res = twoDArrToString(weAdMat,pointsNumber,pointsNumber);
		/* free memory */
		arrMemFree(weAdMat, pointsNumber);
	}
	else if (strcmp(goal, "ddg") == 0){

		/* in case of ddg calculating the weighted adjacency matrix
		 *  and from it calculating and converting to string the diagonal degree matrix*/
		double** weAdMat = wam(arrPoints, pointDimension, pointsNumber);
		double** diDeMat = ddg(weAdMat, pointsNumber);
        res = twoDArrToString(diDeMat,pointsNumber,pointsNumber);
		/* free memory */
		arrMemFree(weAdMat, pointsNumber);
		arrMemFree(diDeMat, pointsNumber);
	}
	else if (strcmp(goal, "lnorm") == 0){

		/* in case of lnorm calculating the weighted adjacency matrix,
		 * from it calculating the diagonal degree matrix and from both of them,
		 * calculating and converting to string the laplacian normalized matrix*/
		double** weAdMat = wam(arrPoints, pointDimension, pointsNumber);
		double** diDeMat = ddg(weAdMat, pointsNumber);
		double** normLapMat = lnorm(weAdMat, diDeMat, pointsNumber);
        res = twoDArrToString(normLapMat,pointsNumber,pointsNumber);
		/* free memory */
		arrMemFree(weAdMat, pointsNumber);
		arrMemFree(diDeMat, pointsNumber);
		arrMemFree(normLapMat, pointsNumber);
	}
	else if(strcmp(goal, "jacobi") == 0){
		/* in case of jacobi, according to what is written on the forum,
		 * in the input file, and so in the "pointsNumber" array we have
		 * a symmetric matrix, not data points from which we find eigenvectors
		 * and eigenvalues and converting them to string */
		 if(pointsNumber != pointDimension){ /* data includes non symmetric matrix */
		    printf("Invalid Input!");
		    exit(1);
		 }
		double** eigenValAndVec = getJacobiEigenValAndVec(arrPoints, pointsNumber);
		res = twoDArrToString(eigenValAndVec,pointsNumber+1,pointsNumber);
		/* free memory */
		arrMemFree(eigenValAndVec, pointsNumber);
	}
	else if(strcmp(goal, "findT") == 0){
	    /* in this case we want to find T matrix, so we- calculating the weighted adjacency matrix,
		 * from it calculating the diagonal degree matrix, from both of them,
		 * calculating the laplacian normalized matrix, finding the eigenvalues
		 * and eigenvectors, sorting them, finding k if not given, calculating U matrix,
		 * and from it T matrix, which we convert to string */
    	double** U;
		double** T;
		double** weAdMat = wam(arrPoints, pointDimension, pointsNumber);
		double** diDeMat = ddg(weAdMat, pointsNumber);
		double** normLapMat = lnorm(weAdMat, diDeMat, pointsNumber);
		double** eigenValAndVec = getJacobiEigenValAndVec(normLapMat, pointsNumber);
		sortEigenMat(eigenValAndVec, pointsNumber);
		/* if k=0 we need to use the eigengap heuristic */
		if(k==0){
			k = findK(eigenValAndVec[0], pointsNumber);
		}
		if(k==1){ /* k = 1 is not legal */
		    printf("An Error Has Occurred");
    		exit(1);
		}
		U = getU(eigenValAndVec, pointsNumber,k);
		T = getT(U, pointsNumber, k);
		res = twoDArrToString(T,pointsNumber,k);
		arrMemFree(U, pointsNumber);
		arrMemFree(T,pointsNumber);
		arrMemFree(weAdMat, pointsNumber);
		arrMemFree(diDeMat, pointsNumber);
		arrMemFree(normLapMat, pointsNumber);
		arrMemFree(eigenValAndVec, pointsNumber+1);
    }
	/* free memory */
	arrMemFree(arrPoints, pointsNumber);
    return Py_BuildValue("s", res); /*  returning the resulted string we got to python */
}

/* this function will implement kmeans algorithm */
static PyObject* kmeans_capi(PyObject *self, PyObject *args)
{
    /* initializing needed variables, which we get from python */
    int k;
    int pointsNumber;
    int max_iter;
    char* centroids_str;
    char* numbers_str;
    double epsilon;
    char *res;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "diiiss",&epsilon, &k,&pointsNumber,&max_iter, &centroids_str, &numbers_str)) {
        return NULL; /* If an error has occurred */
    }
    /* calculating final centroids */
    res = kmeans(epsilon,k,pointsNumber,max_iter,centroids_str,numbers_str);
    return Py_BuildValue("s", res); /*  returning the resulted string we got to python */
}


static PyMethodDef capiMethods[] = {
    /* the main function */
    {"main",
    (PyCFunction) main_capi,
    METH_VARARGS,
    PyDoc_STR("kmeans algorithm")},
    /* kmeans function */
    {"kmeans",
    (PyCFunction) kmeans_capi,
    METH_VARARGS,
    PyDoc_STR("kmeans algorithm")},
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}