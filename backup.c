#include <stdio.h>
#include <stdlib.h>


#define idealBlockSize 1

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

__device__ int MAXn(int * array, int thread_1D_pos, int dnaStrandSize, int length, int Ge)
{
	int max = 0;
	int value;
	int i;
	for (i = 1; i < length; ++i)
	{
		value = (array[thread_1D_pos * dnaStrandSize + length - i] - i * Ge);
		if(value > max)
			max = value; 
	}
	return max;
}


__global__ void alignmentFinder(char *dna, char *dnaCompare, int *E_, int *H_,int *F, int *matrixH, int dnaTestStrandSize, int dnaStrandSize)
{

// need to work on this


    const int thread_1D_pos = blockIdx.x * blockDim.x + threadIdx.x;
 	
	int i;
	int Gs = 8;
	int Ge = 1;
	int S;
	// if(threadIdx.x == 0)
	// 	printf("width:%d Height:%d \n",blockDim.x * gridDim.x, blockDim.y * gridDim.y);
	
/*

	for (i = 1 , l = 0; i < dnaTestStrandSize; ++i, ++l) //  startinf from 1 to account for initial zero matrix
	{// l for compare dna count
		for (j = 1, k = 0, gapInit = 1; j < dnaStrandSize; ++j, ++k)
*/

   
/*

	Thread 0 - dnaTestSize to colummn
	rest to row
	and then diagonal


*/


	if(thread_1D_pos < dnaStrandSize)
	{
		if(thread_1D_pos  != 0)
		{
			for (i = 1; i < dnaTestStrandSize + 1; ++i) //  startinf from 1 to account for initial zero matrix
			{// l for compare dna count
				F[i * (dnaStrandSize + 1) + thread_1D_pos] = MAX2(F[(i - 1) * (dnaStrandSize + 1) + thread_1D_pos], matrixH[(i - 1) * (dnaStrandSize + 1) + thread_1D_pos] - Gs) - Ge;
			}			
		}		
	}


	if(thread_1D_pos  != 0)
	{
		printf("INDEX:%d\n",thread_1D_pos);
	//  set up threads to not edit column 0

			// for (int j = 0; j < thread_1D_pos; ++j)
			// {
			// 	__syncthreads();
			// }	
		


			if(dna[i - 1] == dnaCompare[thread_1D_pos - 1])
				S = 5;
			else
				S = -3;

		for (i = thread_2D_pos.x,; i < dnaStrandSize ; ++i)
		{
			H_[i * (dnaStrandSize + 1) + i] = MAX3(matrixH[(thread_1D_pos - 1) * (dnaStrandSize + 1) + i -1 ] + S, F[thread_1D_pos * (dnaStrandSize + 1) + i],0);
		}
			
/*

http://stackoverflow.com/questions/23064866/matrix-filling-for-the-smith-waterman-algorithm-in-cuda

paralleize one row at a time !
 ie loop over and call kernel calls


 do it in three kernel codes

 first 

 F
 H_

 second

 E_

 third

 Hi,j


 all the above looping over rows
 

*/




			//printf("dna:%c dnaCompare:%c Gs:%d Ge:%d S:%d\n",dna[k],dnaCompare[l],Gs,Ge,S);
	/*
				F[i][j] = MAX2(F[i - 1][j], matrixH[i - 1][j] - Gs) - Ge;
				H_[i][j] = MAX3(matrixH[i - 1][j - 1] + S, F[i][j],0);
				E_[i][j] = MAXn(H_[i], j, Ge);
				matrixH[i][j] = MAX2(H_[i][j], E_[i][j] - Gs);

	*/


/*

spawn column  width threads

column threads can do  eqa 1
all threads would do equa 2 in a diagonal fashion
all threads for equa 3 in normal fashion
roq fpr equ 3
*/


			F[i * (dnaStrandSize + 1) + thread_2D_pos.x] = MAX2(F[(thread_1D_pos - 1) * (dnaStrandSize + 1) + i], matrixH[(thread_1D_pos - 1) * (dnaStrandSize + 1) + i] - Gs) - Ge;
			 // do this column wise 

			// __syncthreads();
			//printf("-->thread: %d i:%d F:%d F_:%d matrixH:%d \n",thread_1D_pos, i,F[thread_1D_pos * (dnaStrandSize + 1) + i],F[(thread_1D_pos - 1) * (dnaStrandSize + 1) + i],matrixH[(thread_1D_pos - 1) * (dnaStrandSize + 1) + i]);
			// __syncthreads();
			H_[thread_1D_pos * (dnaStrandSize + 1) + i] = MAX3(matrixH[(thread_1D_pos - 1) * (dnaStrandSize + 1) + i -1 ] + S, F[thread_1D_pos * (dnaStrandSize + 1) + i],0);
			//
			// __syncthreads();

			// diagonal ?? 
			E_[thread_1D_pos * (dnaStrandSize + 1) + i] = MAXn(H_,thread_1D_pos, (dnaStrandSize + 1), i, Ge); // check this function for i 
			// 
			// __syncthreads();
			matrixH[thread_1D_pos * (dnaStrandSize + 1) + i] = MAX2(H_[thread_1D_pos * (dnaStrandSize + 1) + i], E_[thread_1D_pos * (dnaStrandSize + 1) + i] - Gs);
			// if(threadIdx.x == 0)
		   // printf("-->thread: %d i:%d F:%d H_:%d E_:%d matrixH:%d \n",thread_1D_pos, i,F[thread_1D_pos * (dnaStrandSize + 1) + i],H_[thread_1D_pos * (dnaStrandSize + 1) + i],E_[thread_1D_pos * (dnaStrandSize + 1) + i],matrixH[thread_1D_pos * (dnaStrandSize + 1) + i]);
			//__syncthreads();
		}
	}			//	
}

int main(int argc, char const *argv[])
{
	int i,j,k,l;
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
	

	char *dna = NULL;
	char *dnaCompare = NULL;

	char *d_dna = NULL;
	char *d_dnaCompare = NULL;



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


	gpuErrchk(cudaMemcpy(d_dna, dna, dnaStrandSize * sizeof(char), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_dnaCompare, dnaCompare, dnaTestStrandSize * sizeof(char), cudaMemcpyHostToDevice));

	// dim3 grid(dnaStrandSize/idealBlockSize, dnaTestStrandSize/idealBlockSize);
	// dim3 block(idealBlockSize, idealBlockSize);

	dim3 grid((dnaStrandSize + 1)/idealBlockSize + (dnaTestStrandSize + 1)/idealBlockSize);
	dim3 block(idealBlockSize);


	alignmentFinder<<<grid, block>>>(d_dna, d_dnaCompare, d_E_, d_H_, d_F, d_matrixH, dnaTestStrandSize, dnaStrandSize);

	gpuErrchk(cudaMemcpy(matrixH, d_matrixH, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int) , cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(E_, d_E_, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int) , cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(F, d_F, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int) , cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(H_, d_H_, (dnaStrandSize + 1) * (dnaTestStrandSize + 1) * sizeof(int) , cudaMemcpyDeviceToHost));


	if(!(fp = fopen("GPUdata","w")))
	{
		printf("%s is not open !! \n","GPUdata");
		return -1;
	}
	for (i = 0; i < dnaTestStrandSize + 1; ++i)
	{
		for (j = 0; j < dnaStrandSize + 1; ++j)
		{
			printf("%d ",matrixH[i * (dnaStrandSize + 1) + j]);
			fprintf(fp, "%d ",matrixH[i * (dnaStrandSize + 1) + j]);
		}
		printf("\n");
		fprintf(fp,"\n");
	}
	printf("\n");
	fprintf(fp,"\n");

	for (i = 0; i < dnaTestStrandSize + 1; ++i)
	{
		for (j = 0; j < dnaStrandSize + 1; ++j)
		{
			printf("%d ",F[i * (dnaStrandSize + 1) + j]);
			fprintf(fp, "%d ",F[i * (dnaStrandSize + 1) + j]);
		}
		printf("\n");
		fprintf(fp,"\n");
	}


	printf("\n");
	fprintf(fp,"\n");

	for (i = 0; i < dnaTestStrandSize + 1; ++i)
	{
		for (j = 0; j < dnaStrandSize + 1; ++j)
		{
			printf("%d ",H_[i * (dnaStrandSize + 1) + j]);
			fprintf(fp, "%d ",H_[i * (dnaStrandSize + 1) + j]);
		}
		printf("\n");
		fprintf(fp,"\n");
	}


	printf("\n");
	fprintf(fp,"\n");

	for (i = 0; i < dnaTestStrandSize + 1; ++i)
	{
		for (j = 0; j < dnaStrandSize + 1; ++j)
		{
			printf("%d ",E_[i * (dnaStrandSize + 1) + j]);
			fprintf(fp, "%d ",E_[i * (dnaStrandSize + 1) + j]);
		}
		printf("\n");
		fprintf(fp,"\n");
	}	
	fclose(fp);

	return 0;
}