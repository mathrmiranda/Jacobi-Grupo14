#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ERRO_ACEITAVEL 0.001

int main(int argc, char *argv[])
{
    /**Inicialização das matrizes.**/
    int tam,seed,i,j;
    double *b,*xk,*xk_1, somalinha = 0;

     if ( argc  != 2)
    {
	printf("Argumentos errados. Utilize <ordem da matriz> \n");
	exit(0);
    }

    tam = atoi(argv[1]);
    seed = 2;

    double A[tam][tam];
    b =(double*)malloc(tam*sizeof(double));
    xk =(double*)malloc(tam*sizeof(double));
    xk_1 =(double*)malloc(tam*sizeof(double));


    /**Preenchimento das matrizes A e b **/

    srand(seed);

    //Matriz A
    for(i=0;i<tam;i++)
    {
        for(j=0;j<tam;j++)
        {
            if(i != j)
            {
                A[i][j] = rand() %100;
                somalinha += A[i][j];
            }
        }
        A[i][i] = somalinha + (rand()%10) + 1; //Garante convergencia
        somalinha = 0;
    }

    printf("A = \n");
    for(i=0;i<tam;i++)
    {
        for(j=0;j<tam;j++)
        {
            printf("%.1f ", A[i][j]);
        }
        printf("\n");
    }

    //Matriz b
    printf("\nb = \n");
    for(i=0;i<tam;i++)
    {
        b[i] = rand() % 100;
        printf("%.1f\n", b[i]);
    }
    printf("\n");



    /**METODO DE GAUSS-JACOBI**/

    int k = 0;
    double somaax = 0, precit = 1;
    double precnum = 0, precden = 0;

    //x para a primeira iteraçao (k = 0)
    for(i=0;i<tam;i++)
    {
        xk[i] = b[i]/A[i][i];
    }

    //Começo das iterações
    while(precit > ERRO_ACEITAVEL)
    {
        //Iteração K+1
        for(i=0;i<tam;i++)
        {
            for(j=0;j<tam;j++)
            {
                if(i != j)
                    somaax = somaax + (xk[j]*A[i][j]);
            }
            xk_1[i] = (b[i] - somaax)/A[i][i];

            //Calculando precisao
            if(fabs(xk_1[i]-xk[i]) > precnum)
            {
                precnum = fabs(xk_1[i]-xk[i]);
            }

            if(fabs(xk_1[i]) > precden)
            {
                precden = fabs(xk_1[i]);
            }
            somaax = 0; // reset da soma para proxima linha
        }
        precit = precnum/precden;

        //resetando variaveis para próxima iteração
        precnum = 0;
        precden = 0;
        for(i=0;i<tam;i++)
        {
            xk[i] = xk_1[i];
        }
        k++;
    }

    //Imprimindo resposta
    printf("x (k = %d) = \n", k);
    for(i=0;i<tam;i++)
    {
        printf("%0.3f\n", xk_1[i]);
    }


    /**TESTE DA RESPOSTA**/
    int equa;
    double soma = 0;
    printf("\n");
    printf("qual equacao deseja utilizar para verificar resposta ? (considere comecando em 1)\n");
    scanf("%d", &equa);
    printf("\n");

    equa--;
    for(j=0;j<tam;j++)
    {
        soma += A[equa][j]*xk_1[j];
    }

    printf("resultado pelo metodo Gauss-Jacobi = %.3f\n", soma);
    printf("valor de b%d = %.3f\n", equa+1, b[equa]);

    /**LIBERACAO DE MEMORIA**/
    free(b);
    free(xk);
    free(xk_1);

    return 0;
}
