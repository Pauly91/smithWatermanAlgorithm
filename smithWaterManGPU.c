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

__device__ int MAXn(int * array, int length, int Ge)
{
	int max = 0;
	int value;
	int i;
	for (i = 1; i < length; ++i)
	{
		value = (array[length - i] - i * Ge);
		if(value > max)
			max = value; 
	}
	return max;
}


__global__ void alignmentFinder(char *d_dna, char *d_dnaCompare, int *d_E_, int *d_H_,int *d_F, int *d_matrixH, int dnaTestStrandSize, int dnaStrandSize)
{

// need to work on this 
	int Gs = 8;
	int Ge = 1;
	int S;

	for (int i = 1 , int l = 0; i < dnaTestStrandSize; ++i, ++l) //  startinf from 1 to account for initial zero matrix
	{// l for compare dna count
		for (int j = 1, int k = 0; j < dnaStrandSize; ++j, ++k)
		{
			// k for dna count

			if(dna[k] == dnaCompare[l])
				S = 5;
			else
				S = -3;

			//printf("dna:%c dnaCompare:%c Gs:%d Ge:%d S:%d\n",dna[k],dnaCompare[l],Gs,Ge,S);

			F[i][j] = MAX2(F[i - 1][j], matrixH[i - 1][j] - Gs) - Ge;
			H_[i][j] = MAX3(matrixH[i - 1][j - 1] + S, F[i][j],0);
			E_[i][j] = MAXn(H_[i], j, Ge);
			matrixH[i][j] = MAX2(H_[i][j], E_[i][j] - Gs);
			//printf("-->i:%d j:%d F:%d H_:%d E_:%d matrixH:%d \n",i,j,F[i][j],H_[i][j],E_[i][j],matrixH[i][j]);
		}
	}	
}

int main(int argc, char const *argv[])
{
	int i,j,k,l;
	int **matrixH = NULL;
	int **F = NULL;
	int **H_ = NULL;
	int **E_ = NULL;
	
	int dnaStrandSize;
	int dnaTestStrandSize;
	int fscanfReturn;

	int **d_matrixH = NULL;
	int **d_F = NULL;
	int **d_H_ = NULL;
	int **d_E_ = NULL;
	

	char *dna = NULL;
	char *dnaCompare = NULL;

	char *d_dna = NULL;
	char *d_dnaCompare = NULL;

	int Gs = 0,Ge = 0,gapInit,S;

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

	printf("%d \n",strandSize);
	
	dna = (char *) malloc(dnaStrandSize * sizeof(char));
	dnaCompare = (char *) malloc(dnaTestStrandSize * sizeof(char));	

	F = (int **) calloc((dnaTestStrandSize + 1), sizeof(int *));
	H_ = (int **) calloc((dnaTestStrandSize + 1), sizeof(int *));
	E_ = (int **) calloc((dnaTestStrandSize + 1), sizeof(int *));
	matrixH = (int **) calloc((dnaTestStrandSize + 1 ), sizeof(int *));
	
	for (i = 0; i < dnaTestStrandSize + 1; ++i)
	{
		F[i] = (int*) calloc(dnaStrandSize + 1, sizeof(int));
		H_[i] = (int*) calloc(dnaStrandSize + 1, sizeof(int));
		E_[i] = (int*) calloc(dnaStrandSize + 1, sizeof(int));
		matrixH[i] = (int*) calloc(dnaStrandSize + 1, sizeof(int));

	}


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


	gpuErrchk(cudaMalloc((void**)&d_E_, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&d_H_ ,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&d_F,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&d_matrixH,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));

	gpuErrchk(cudaMemset(d_E_, 0,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMemset(d_H_, 0,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMemset(d_F, 0,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));
	gpuErrchk(cudaMemset(d_matrixH, 0,(dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int)));


	gpuErrchk(cudaMalloc((void**)&d_dna, dnaStrandSize * sizeof(char)));
	gpuErrchk(cudaMalloc((void**)&d_dnaCompare, dnaTestStrandSize * sizeof(char)));


	gpuErrchk(cudaMemcpy(dna, d_dna dnaTestStrandSize * sizeof(char), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dnaCompare, d_dnaCompare dnaTestStrandSize * sizeof(char), cudaMemcpyHostToDevice));

	dim3 grid(dnaTestStrandSize/idealBlockSize, dnaStrandSize/idealBlockSize);
	dim3 block(idealBlockSize, idealBlockSize);

	alignmentFinder<<<grid, block>>>(d_dna, d_dnaCompare, d_E_, d_H_, d_F, d_matrixH, dnaTestStrandSize, dnaStrandSize);

	gpuErrchk(cudaMemcpy(matrixH, d_matrixH, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int) , cudaMemcpyDeviceToHost));


	for (i = 0; i < strandSize + 1; ++i)
	{
		for (j = 0; j < strandSize + 1; ++j)
		{
			printf("%d ",matrixH[i][j]);
		}
		printf("\n");
	}


	return 0;
}