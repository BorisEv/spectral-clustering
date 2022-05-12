import sys
import numpy as np
import spkmeansmodule as spkm
np.random.seed(0)  # setting needed random

# dealing with the arguments
try:
    k = int(sys.argv[1])
except:  # in case of invalid k input
    print("Invalid input!")
    sys.exit()
goals = ["wam","ddg","lnorm","jacobi","spk"]
goal = sys.argv[2]
filename = sys.argv[3]

# checking the arguments and acting accordingly
if (goal not in goals) or (k < 0):
    print("Invalid Input!")
    sys.exit()
else:
    if goal == "spk":
        # in case of spk we calculate T matrix, then implementing k++ algorithm,
        # and call to c kmeans function with initial centroids and Tmatrix
        try:
            # calling for c to get T matrix
            Tmatrix = spkm.main("findT", filename, k)
        except:
            print("An Error Has Occurred")
            sys.exit()
        # converting string to 2d python list of floats "numbers"
        numbers = []
        mat = Tmatrix.splitlines()
        for i in range(len(mat)):
            mat[i] = mat[i].split(",")
        for i in range(len(mat)):
            numbers.append([])
            for j in range(len(mat[i])):
                try:
                    numbers[i].append(float(mat[i][j]))
                except:
                    print("An Error Has Occurred")
                    sys.exit()

        # starting the K++ algorithm
        indexes = []  # saving the indexes
        # initializing needed variables
        N = len(numbers)
        k = len(numbers[0])
        D = [np.inf for i in range(N)]
        P = [0 for i in range(N)]
        initial_centroids = [[] for i in range(k)]

        numbersindexes = [i for i in range(len(numbers))]  # creating indexes because the random function works with 1D
        # initialize the first centroid
        numberindex = np.random.choice(numbersindexes)
        indexes.append(numberindex)  # updating the indexes we output
        initial_centroids[0] = numbers[numberindex]

        # the main iterations
        for i in range(1, k):
            for x in range(N):
                for j in range(i):
                    # calculating the norm
                    norm = 0
                    for z in range(len(numbers[x])):
                        norm += ((numbers[x][z]) - initial_centroids[j][z]) ** 2
                    # finding Dl
                    D[x] = min(D[x], norm)
            sumD = sum(D)

            # calculating the probability distribution
            for x in range(N):
                P[x] = D[x] / sumD
            # choosing next centroid with respect to p and updating indexes we output
            numberindex = np.random.choice(numbersindexes, p=P)
            indexes.append(numberindex)
            initial_centroids[i] = numbers[numberindex]

        # initializing needed variables to kmeans algorithm
        epsilon = 0
        max_iter = 300
        centroids_str = ""
        # converting initial_centroids to string
        for i in initial_centroids:
            for j in i:
                centroids_str += ("%.4f" % j + ",")
            centroids_str = centroids_str[0:len(centroids_str)-1]+ "\n"
        # call kmeans algorithm. centroids_str - initial centroids. Tmatrix - N data numbers
        res = spkm.kmeans(epsilon, k, N, max_iter,centroids_str,Tmatrix)
        # output indexes of initial centroids and then the result
        print(*indexes, sep=",")
        print(res, end="")
    else:
        # if the goal is not spk, call c program with right goal parameter
        res = spkm.main(goal, filename,k)
        # output the result
        print(res, end="")


