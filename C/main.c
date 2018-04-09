#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void imprimeMatrizCompleta(double **m, int dim);

void memcpy_d(double *dest, double *src, int _size)
{
	int i;

	for(i = 0; i< _size ;i++)
	{
		dest[i] = src[i];
	}
}

double *substituicaoRegressiva(double **m, int dim)
{
    double *root = (double*)malloc(dim * sizeof(double));
    double sum;
    int i,j,n;

    n = dim - 1;
    root[n] = m[n][dim]/(double)m[n][n];

    for(i = n - 1; i >= 0; i--)
    {
        sum = 0;

        for(j = i + 1; j <= n; j++ )
        {
            sum += m[i][j] * root[j];
        }

        root[i] = (m[i][dim] - sum)/(double)m[i][i];
    }
    return root;
}

void triangularSuperior_p(double **m, int dim)
{
    int i,j,k,l,troca;
    double n;

    double *aux;

	aux = (double*)malloc((dim + 1) * sizeof(double));

    for(i = 0; i < dim; i++)
    {

		troca = -1;
		for(l = i; l < dim; l++)
		{
			if(m[i][i] < fabs(m[l][i]))
			{
				troca = l;
			}
		}

		if(troca != -1)
		{
			memcpy_d(aux, m[troca], dim + 1);
			memcpy_d(m[troca], m[i], dim + 1);
			memcpy_d(m[i], aux, dim + 1);
		}
		imprimeMatrizCompleta(m,dim);

        for(j = i + 1; j < dim; j++)
        {
            n = m[j][i]/(double)m[i][i];

            for(k = 0; k < dim + 1; k++)
            {
                m[j][k] = m[j][k] - n * m[i][k];
            }
        }
    }
}


void triangularSuperior(double **m, int dim)
{
    int i,j,k;
    double n;

    for(i = 0; i < dim; i++)
    {
        for(j = i + 1; j < dim; j++)
        {
            n = m[j][i]/(double)m[i][i];

            for(k = 0; k < dim + 1; k++)
            {
                m[j][k] = m[j][k] - n * m[i][k];
            }
        }
    }
}

double **lerMatrizCompleta(const char *arg, int *dim)
{
    double **m;
    int i,j;

    FILE *arq = fopen(arg,"r");

    if(arq == NULL)
    {
        printf("arquivo nao encontrado\n");
        exit(1);
    }

    fscanf(arq,"%d",dim);

    m = (double**)malloc((*dim) * sizeof(double*));

    for(i = 0; i< *dim; i++)
    {
        m[i] = (double*)malloc((*dim + 1) * sizeof(double));
    }

    for(i = 0; i < *dim; i++)
    {
        for(j = 0; j < *dim + 1; j++)
        {
            fscanf(arq,"%lf", &m[i][j]);
        }
    }
    fclose(arq);
    return m;

}

void imprimeMatrizCompleta(double **m, int dim)
{
    int i,j;

    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim + 1; j++)
        {
            (j == dim)? printf("| %5.2lf\t",m[i][j]): printf("%5.2lf\t",m[i][j]);
        }
        puts("");
    }
    puts("________");
}

void imprimeRaiz(double *r, int dim)
{
    int i;

    puts("\n__roots__\n");

    for(i = 0; i < dim; i++)
    {
        printf("x[%d] = %6.5lf\n",i + 1,r[i]);
    }

    puts("________");
}

int main(int argc, char **argv)
{
    double **matriz;
    double *root;
    int dim;

    matriz = lerMatrizCompleta(argv[1],&dim);
    printf("Matriz de entrada\n");
    imprimeMatrizCompleta(matriz,dim);
    printf("pivotamento\n");
    triangularSuperior_p(matriz,dim);
     printf("Matriz de saida\n");
    imprimeMatrizCompleta(matriz,dim);

    root = substituicaoRegressiva(matriz,dim);
    imprimeRaiz(root,dim);

    return 0;
}
