#include "Laplacian.h"
#include "Parameters.h"
#include "PointwiseOps.h"
#include "Reductions.h"
#include "Utilities.h"
#include "Timer.h"
#include <iostream>

void ConjugateGradients(
    float (&x)[XDIM][YDIM][ZDIM],
    const float (&f)[XDIM][YDIM][ZDIM],
    float (&p)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const bool writeIterations)
{
	Timer t;
	extern double ComputeLaplacianTime;
	extern double SaxpyTime;
	extern double CopyTime;
	extern double InnerProductTime;
	extern double NormTime;
    // Algorithm : Line 2
	t.Start();
    ComputeLaplacian(x, z);
	t.Stop();
	ComputeLaplacianTime += t.getRuntime();
	
	t.Start();
    Saxpy(z, f, r, -1);
    t.Stop();
	SaxpyTime += t.getRuntime();
	
	float nu = Norm(r);

    // Algorithm : Line 3
    if (nu < nuMax) return;
        
    // Algorithm : Line 4
    t.Start();
	Copy(r, p);
    t.Stop();
	CopyTime += t.getRuntime();
	
	t.Start();
	float rho=InnerProduct(p, r);
    t.Stop();
	InnerProductTime += t.getRuntime();
	
    // Beginning of loop from Line 5
    for(int k=0;;k++)
    {
        std::cout << "Residual norm (nu) after " << k << " iterations = " << nu << std::endl;

        // Algorithm : Line 6
		t.Start();
        ComputeLaplacian(p, z);
		t.Stop();
		ComputeLaplacianTime += t.getRuntime();
		
		t.Start();
        float sigma=InnerProduct(p, z);
		t.Stop();
		InnerProductTime += t.getRuntime();
        // Algorithm : Line 7
        float alpha=rho/sigma;

        // Algorithm : Line 8
		t.Start();
        Saxpy(z, r, r, -alpha);
		t.Stop();
		SaxpyTime += t.getRuntime();
		
		t.Start();
        nu=Norm(r);
		t.Stop();
		NormTime += t.getRuntime();
		
        // Algorithm : Lines 9-12
        if (nu < nuMax || k == kMax) {
			t.Start();
            Saxpy(p, x, x, alpha);
			t.Stop();
			SaxpyTime += t.getRuntime();
			
            std::cout << "Conjugate Gradients terminated after " << k << " iterations; residual norm (nu) = " << nu << std::endl;
            //if (writeIterations) WriteAsImage("x", x, k, 0, 127);
            return;
        }
            
        // Algorithm : Line 13
		t.Start();
        Copy(r, z);
		t.Stop();
		CopyTime += t.getRuntime();
		
		t.Start();
        float rho_new = InnerProduct(z, r);
		t.Stop();
		InnerProductTime += t.getRuntime();
		
        // Algorithm : Line 14
        float beta = rho_new/rho;

        // Algorithm : Line 15
        rho=rho_new;

        // Algorithm : Line 16
		/**
		t.Start();
        Saxpy(p, x, x, alpha);
		t.Stop();
		SaxpyTime += t.getRuntime();
		t.Start();
        Saxpy(p, r, p, beta);
		t.Stop();
		SaxpyTime += t.getRuntime();
		**/
		t.Start();
		Saxpy(p, x, x, alpha, p, r, p, beta);
		t.Stop();
		SaxpyTime += t.getRuntime();
		
        //if (writeIterations) WriteAsImage("x", x, k, 0, 127);
    }
}
