#include <stdio.h>
#include <stdlib.h>


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

int MAXn(int * array, int strandSize, int Ge)
{
	int max = 0;
	int value;
	int i;
	for (i = 1; i < strandSize; ++i)
	{
		value = (array[strandSize - i] - i * Ge);
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
	int strandSize;
	int fscanfReturn;
	char *dna = NULL;
	char *dnaCompare = NULL;

	int Gs,Ge,gapInit,S;

	FILE *fp = NULL;

	if(!(fp = fopen(argv[1],"r")))
	{
		printf("%s is not open !! \n",argv[1]);
	}

	fscanfReturn = fscanf(fp,"%d",&strandSize);
	if(fscanfReturn < 0)
		printf("Scanning of Value Failed\n");
	printf("%d \n",strandSize);
	
	dna = (char *) malloc(strandSize * sizeof(char));
	dnaCompare = (char *) malloc(strandSize * sizeof(char));	

	F = (int **) calloc((strandSize + 1), sizeof(int *));
	H_ = (int **) calloc((strandSize + 1), sizeof(int *));
	E_ = (int **) calloc((strandSize + 1), sizeof(int *));
	matrixH = (int **) calloc((strandSize + 1 ), sizeof(int *));
	
	for (i = 0; i < strandSize; ++i)
	{
		F[i] = (int*) calloc(strandSize + 1, sizeof(int));
		H_[i] = (int*) calloc(strandSize + 1, sizeof(int));
		E_[i] = (int*) calloc(strandSize + 1, sizeof(int));
		matrixH[i] = (int*) calloc(strandSize + 1, sizeof(int));

	}


	for (i = 0; i < strandSize; ++i)
	{
		fscanfReturn = fscanf(fp," %c ",&dna[i]);
		if(fscanfReturn < 0)
			printf("Scanning of Value Failed\n");
	}

	for (i = 0; i < strandSize; ++i)
	{
		fscanfReturn = fscanf(fp," %c ",&dnaCompare[i]);
		if(fscanfReturn < 0)
			printf("Scanning of Value Failed\n");	
	}
	fclose(fp);



	for (i = 1 , l = 0; i < strandSize; ++i, ++l) //  startinf from 1 to account for initial zero matrix
	{// l for compare dna count
		for (j = 1, k = 0, gapInit = 1; j < strandSize; ++j, ++k)
		{
			// k for dna count
			if(dna[k] == ' ' && dna[j - 1] != dnaCompare[l])// since i and j are starting from 1
			{	
				if(gapInit == 1) // first gap in the entire sequence
				{
					Gs = 8;
					Ge = 0;
					gapInit = 0;
				}
				else
				{
					if(dna[k - 1] != ' ') // gap Starting
					{
						Gs = 8;
						Ge = 0;
					}
					else // extension
					{
						Gs = 0;
						Ge = 1;	
					}
				}
				S = -3;
			}
			else
			{
				if(dna[k] == dnaCompare[l])
					S = 5;
				else
					S = -3;
			}

			F[i][j] = MAX2(F[i - 1][j], matrixH[i - 1][j] - Gs) - Ge;
			H_[i][j] = MAX3(matrixH[i - 1][j - 1] + S, F[i][j],0);
			E_[i][j] = MAXn(H_[i], strandSize, Ge);
			matrixH[i][j] = MAX2(H_[i][j], E_[i][j] - Gs);
		}
	}



	return 0;
}