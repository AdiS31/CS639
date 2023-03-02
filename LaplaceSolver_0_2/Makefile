CC = g++
CFLAGS = -fopenmp -Ofast

main: main.o Laplacian.o ConjugateGradients.o PointwiseOps.o Reductions.o Utilities.o
	$(CC) $(CFLAGS) -o main main.o Laplacian.o ConjugateGradients.o PointwiseOps.o Reductions.o Utilities.o
main.o: main.cpp ConjugateGradients.h Timer.h Utilities.h
	$(CC) $(CFLAGS) -c main.cpp
Laplacian.o: Laplacian.cpp Laplacian.h
	$(CC) $(CFLAGS) -c Laplacian.cpp
ConjugateGradients.o: ConjugateGradients.cpp Laplacian.h Parameters.h PointwiseOps.h Reductions.h Utilities.h Timer.h
	$(CC) $(CFLAGS) -c ConjugateGradients.cpp
PointwiseOps.o: PointwiseOps.cpp PointwiseOps.h
	$(CC) $(CFLAGS) -c PointwiseOps.cpp
Reductions.o: Reductions.cpp Reductions.h
	$(CC) $(CFLAGS) -c Reductions.cpp
Utilities.o: Utilities.cpp Utilities.h
	$(CC) $(CFLAGS) -c Utilities.cpp
clean:
	rm -rf *.o