#include "ConjugateGradients.h"
#include "Timer.h"
#include "Utilities.h"

double ComputeLaplacianTime = 0;
double SaxpyTime = 0;
double CopyTime = 0;
double InnerProductTime = 0;
double NormTime = 0;

int main(int argc, char *argv[])
{
    using array_t = float (&) [XDIM][YDIM][ZDIM];

    float *xRaw = new float [XDIM*YDIM*ZDIM];
    float *fRaw = new float [XDIM*YDIM*ZDIM];
    float *pRaw = new float [XDIM*YDIM*ZDIM];
    float *rRaw = new float [XDIM*YDIM*ZDIM];
    float *zRaw = new float [XDIM*YDIM*ZDIM];
    
    array_t x = reinterpret_cast<array_t>(*xRaw);
    array_t f = reinterpret_cast<array_t>(*fRaw);
    array_t p = reinterpret_cast<array_t>(*pRaw);
    array_t r = reinterpret_cast<array_t>(*rRaw);
    array_t z = reinterpret_cast<array_t>(*zRaw);
    
    // Initialization
    {
        Timer timer;
        timer.Start();
        InitializeProblem(x, f);
        timer.Stop("Initialization : ");
    }

    // Call Conjugate Gradients algorithm
	Timer t;
	t.Start();
    ConjugateGradients(x, f, p, r, z);
	t.Stop();
	double total = t.getRuntime();
	//printf("Total ConjugateGradients time: %fms\n", t.getRuntime());
	
	double AvgCL = ComputeLaplacianTime / 256.0;
	double AvgS = SaxpyTime / 256.0;
	double AvgC = CopyTime / 256.0;
	double AvgIP = InnerProductTime / 256.0;
	double AvgN = NormTime / 256.0;
    printf("ComputeLaplacian: Total %fms, Average %fms\nSaxpy: Total %fms, Average %fms\nCopy: Total %fms, Average %fms\nInnerProduct: Total %fms, Average %fms\nNorm: Total %fms, Average %fms\n", ComputeLaplacianTime, AvgCL, SaxpyTime, AvgS, CopyTime, AvgC, InnerProductTime, AvgIP, NormTime, AvgN);
    printf("Total ConjugateGradients time: %fms\n", total);
	printf("Sum of all Kernel time = %fms\n", ComputeLaplacianTime + SaxpyTime + CopyTime + InnerProductTime + NormTime);
	return 0;
}
