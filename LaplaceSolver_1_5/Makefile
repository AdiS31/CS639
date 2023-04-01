CC = icc
CFLAGS = -qopenmp -mkl

main: main.o ConjugateGradients.o Laplacian.o MatVecMultiply.o PointwiseOps.o Reductions.o Utilities.o
	$(CC) $(CFLAGS) -o main main.o ConjugateGradients.o Laplacian.o MatVecMultiply.o PointwiseOps.o Reductions.o Utilities.o
main.o: main.cpp ConjugateGradients.h Timer.h Utilities.h Laplacian.h
	$(CC) $(CFLAGS) -c main.cpp
ConjugateGradients.o: ConjugateGradients.cpp ConjugateGradients.h Laplacian.h PointwiseOps.h Reductions.h Utilities.h Timer.h
	$(CC) $(CFLAGS) -c ConjugateGradients.cpp
Laplacian.o: Laplacian.cpp Laplacian.h CSRMatrixHelper.h MatVecMultiply.h
	$(CC) $(CFLAGS) -c Laplacian.cpp
MatVecMultiply.o: MatVecMultiply.cpp MatVecMultiply.h 
	$(CC) $(CFLAGS) -c MatVecMultiply.cpp
PointwiseOps.o: PointwiseOps.cpp PointwiseOps.h 
	$(CC) $(CFLAGS) -c PointwiseOps.cpp
Reductions.o: Reductions.cpp Reductions.h 
	$(CC) $(CFLAGS) -c Reductions.cpp
Utilities.o: Utilies.cpp Utilities.h 
	$(CC) $(CFLAGS) -c Utilites.cpp