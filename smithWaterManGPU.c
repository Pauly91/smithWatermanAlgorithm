#include <stdio.h>
#include <stdlib.h>


#define idealBlockSize 25

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ int MAX2(int a, int b)
{
	if(a >= b)
		return a;
	else
		return b;
}


__device__ int MAX3(int a, int b, int c)
{
	if(a >= b)
	{
		if(a >= c)
			return a;
		else
			return c;
	}
	else
	{
		if(b >= c)
			return b;
		else
			return c;
	}
}

__device__ int MAXn(int * array, int row, int dnaStrandSize, int length, int Ge)
{
	int max = 0;
	int value;
	int i;
	for (i = 1; i < length; ++i)
	{
		value = (array[row * dnaStrandSize + length - i] - i * Ge);
		if(value > max)
			max = value; 
	}
	return max;
}


__global__ void alignmentFinder(char *dna, char *dnaCompare, int *E_, int *H_,int *F, int *matrixH, int dnaStrandSize, int row)
{

    const int thread_1D_pos = blockIdx.x * blockDim.x + threadIdx.x;
 
	int Gs = 8;
	int Ge = 1;
	int S;

	if(thread_1D_pos  != 0)
	{
		if(dna[thread_1D_pos - 1] == dnaCompare[row - 1])
			S = 5;
		else
			S = -3;
		F[row * (dnaStrandSize + 1) + thread_1D_pos] = MAX2(F[(row - 1) * (dnaStrandSize + 1) + thread_1D_pos], matrixH[(row - 1) * (dnaStrandSize + 1) + thread_1D_pos] - Gs) - Ge;	
		H_[row * (dnaStrandSize + 1) + thread_1D_pos] = MAX3(matrixH[(row - 1) * (dnaStrandSize + 1) + thread_1D_pos - 1] + S, F[row * (dnaStrandSize + 1) + thread_1D_pos],0);
		E_[row * (dnaStrandSize + 1) + thread_1D_pos] = MAXn(H_,row, (dnaStrandSize + 1),thread_1D_pos, Ge);  
		matrixH[row * (dnaStrandSize + 1) + thread_1D_pos] = MAX2(H_[row * (dnaStrandSize + 1) + thread_1D_pos], E_[row * (dnaStrandSize + 1) + thread_1D_pos] - Gs);			
	}


}

int main(int argc, char const *argv[])
{
	int i,j;
	int *matrixH = NULL;
	int *F = NULL;
	int *H_ = NULL;
	int *E_ = NULL;
	
	int dnaStrandSize;
	int dnaTestStrandSize;
	int fscanfReturn;

	int *d_matrixH = NULL;
	int *d_F = NULL;
	int *d_H_ = NULL;
	int *d_E_ = NULL;

	int blockSize;
	

	char *dna = NULL;
	char *dnaCompare = NULL;

	char *d_dna = NULL;
	char *d_dnaCompare = NULL;

    float time;
    cudaEvent_t start, stop;




	FILE *fp = NULL;

	if(!(fp = fopen(argv[1],"r")))
	{
		printf("%s is not open !! \n",argv[1]);
	}

	fscanfReturn = fscanf(fp,"%d",&dnaStrandSize);
	if(fscanfReturn < 0)
		printf("Scanning of Value Failed\n");
	fscanfReturn = fscanf(fp,"%d",&dnaTestStrandSize);
	if(fscanfReturn < 0)
		printf("Scanning of Value Failed\n");

	printf("%d \n",dnaStrandSize);
	printf("%d \n",dnaTestStrandSize);
	
	dna = (char *) malloc(dnaStrandSize * sizeof(char));
	dnaCompare = (char *) malloc(dnaTestStrandSize * sizeof(char));	

	F = (int *) calloc((dnaTestStrandSize + 1) * ( dnaStrandSize + 1) , sizeof(int));
	H_ = (int *) calloc((dnaTestStrandSize + 1) * ( dnaStrandSize + 1) , sizeof(int));
	E_ = (int *) calloc((dnaTestStrandSize + 1) * ( dnaStrandSize + 1) , sizeof(int));
	matrixH = (int *) calloc((dnaTestStrandSize + 1) * ( dnaStrandSize + 1) , sizeof(int));


	for (i = 0; i < dnaStrandSize; ++i)
	{
		fscanfReturn = fscanf(fp," %c ",&dna[i]);
		if(fscanfReturn < 0)
			printf("Scanning of Value Failed\n");
	}

	for (i = 0; i < dnaTestStrandSize; ++i)
	{
		fscanfReturn = fscanf(fp," %c ",&dnaCompare[i]);
		if(fscanfReturn < 0)
			printf("Scanning of Value Failed\n");	
	}
	fclose(fp);

    gpuErrchk( cudaEventCreate(&start) );
    gpuErrchk( cudaEventCreate(&stop) );
    gpuErrchk( cudaEventRecord(start, 0) );

	gpuErrchk(cudaMalloc((void**)&d_E_, (dnaStrandSize + 1)  * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&d_H_ ,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&d_F,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&d_matrixH,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));

	gpuErrchk(cudaMemset(d_E_, 0,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMemset(d_H_, 0,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMemset(d_F, 0,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMemset(d_matrixH, 0,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));


	gpuErrchk(cudaMalloc((void**)&d_dna, dnaStrandSize * sizeof(char)));
	gpuErrchk(cudaMalloc((void**)&d_dnaCompare, dnaTestStrandSize * sizeof(char)));


	gpuErrchk(cudaMemcpy(d_dna, dna, dnaStrandSize * sizeof(char), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_dnaCompare, dnaCompare, dnaTestStrandSize * sizeof(char), cudaMemcpyHostToDevice));

	if((dnaStrandSize + 1) < idealBlockSize)
		blockSize = 1;
	else
		blockSize = idealBlockSize;

	dim3 grid((dnaStrandSize + 1)/blockSize);
	dim3 block(blockSize);

	for (i = 1; i < dnaTestStrandSize + 1; ++i)
	{
		alignmentFinder<<<grid, block>>>(d_dna, d_dnaCompare, d_E_, d_H_, d_F, d_matrixH, dnaStrandSize,i);
	}
	gpuErrchk(cudaMemcpy(matrixH, d_matrixH, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int) , cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(E_, d_E_, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int) , cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(F, d_F, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int) , cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(H_, d_H_, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int) , cudaMemcpyDeviceToHost));


    gpuErrchk( cudaEventRecord(stop, 0) );
    gpuErrchk( cudaEventSynchronize(stop) );
    gpuErrchk( cudaEventElapsedTime(&time, start, stop) );


	if(!(fp = fopen("GPUdata","w")))
	{
		printf("%s is not open !! \n","GPUdata");
		return -1;
	}
	for (i = 0; i < dnaTestStrandSize + 1; ++i)
	{
		for (j = 0; j < dnaStrandSize + 1; ++j)
		{
			fprintf(fp, "%d ",matrixH[i * (dnaStrandSize + 1) + j]);
		}
		fprintf(fp,"\n");
	}
	printf("\n");
	fprintf(fp,"\n");
	fclose(fp);
	printf("\n*** GPU ***\n");
    printf("Time to generate:  %3.1f ms \n", time);
    if((fp = fopen("timingInformation","a+")) == NULL)
    {
        printf("File opening has failed\n");
        return -1;
    }
    fprintf(fp, "%s: DNA Sequence:%d Test DNA Sequence:%d Time:%f\n","GPU",dnaStrandSize,dnaTestStrandSize,time);
    fclose(fp);

	return 0;
}