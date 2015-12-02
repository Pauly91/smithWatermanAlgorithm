#include <stdio.h>
#include <stdlib.h>


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


int MAX2(int a, int b)
{
	if(a >= b)
		return a;
	else
		return b;
}


int MAX3(int a, int b, int c)
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

int MAXn(int * array, int length, int Ge)
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
	char *dna = NULL;
	char *dnaCompare = NULL;

	int Gs = 8,Ge = 0,S;

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
	printf("%d \n",dnaStrandSize);
	

	fscanfReturn = fscanf(fp,"%d",&dnaTestStrandSize);
	if(fscanfReturn < 0)
		printf("Scanning of Value Failed\n");
	printf("%d \n",dnaTestStrandSize);	

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
    gpuErrchk( cudaEventCreate(&start) );
    gpuErrchk( cudaEventCreate(&stop) );
    gpuErrchk( cudaEventRecord(start, 0) );



	for (i = 1 , l = 0; i < dnaTestStrandSize + 1; ++i, ++l)
	{
		for (j = 1, k = 0; j < dnaStrandSize + 1; ++j, ++k)
		{
			if(dna[k] == dnaCompare[l])
				S = 5;
			else
				S = -3;
			Gs = 8;
			Ge = 1;
			F[i][j] = MAX2(F[i - 1][j], matrixH[i - 1][j] - Gs) - Ge;
			H_[i][j] = MAX3(matrixH[i - 1][j - 1] + S, F[i][j],0);
			E_[i][j] = MAXn(H_[i], j, Ge);
			matrixH[i][j] = MAX2(H_[i][j], E_[i][j] - Gs);
		}
	}


    gpuErrchk( cudaEventRecord(stop, 0) );
    gpuErrchk( cudaEventSynchronize(stop) );
    gpuErrchk( cudaEventElapsedTime(&time, start, stop) );

	if(!(fp = fopen("CPUdata","w")))
	{
		printf("%s is not open !! \n","CPUdata");
		return -1;
	}	
	for (i = 0; i < dnaTestStrandSize + 1; ++i)
	{
		for (j = 0; j < dnaStrandSize + 1; ++j)
		{
			fprintf(fp, "%d ",matrixH[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
	printf("\n");
	fclose(fp);
	printf("\n*** CPU ***\n");
    printf("Time to generate:  %3.1f ms \n", time);
    if((fp = fopen("timingInformation","a+")) == NULL)
    {
        printf("File opening has failed\n");
        return -1;
    }
    fprintf(fp, "%s: DNA Sequence:%d Test DNA Sequence:%d Time:%f\n","CPU",dnaStrandSize,dnaTestStrandSize,time);
    fclose(fp);

	return 0;
}