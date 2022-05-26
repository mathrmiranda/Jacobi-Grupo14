#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define ERRO_ACEITAVEL 0.001

int main(int argc, char *argv[])
{
    /**RECEBENDO ARGUMENTOS**/
    int tam,seed,i,j,n_threads,threads_criadas;
    double *xk_1;

    if ( argc  != 3)
    {
	printf("Argumentos errados. Utilize <ordem da matriz> <numero de threads> \n");
	exit(0);
    }

    tam = atoi(argv[1]);
    n_threads = atoi(argv[2]);
    seed = 2;

    xk_1 =(double*)malloc(tam*sizeof(double));


    /**INICIALIZANDO MATRIZ S[TAM][TAM+2]**/
    double S[tam][tam+2];
    double somalinha = 0;
    srand(seed);

    //Preenchendo matrizes S e xk_1 em paralelo
    #pragma omp parallel for num_threads(n_threads) private(i,j,somalinha) shared(S,tam,xk_1)
        for(i=0;i<tam;i++)
        {
                for(j=0;j<tam;j++)
                {
                    if(i != j)
                    {
                        S[i][j] = rand() %100;
                        somalinha += S[i][j];
                        //printf("thread %d colocou %.1f\n", omp_get_thread_num(), S[i][j]);
                    }
                }
            //printf("thread %d tem somalinha = %.1f\n", omp_get_thread_num(), somalinha);
            S[i][i] = somalinha + (rand()%10) + 1; //Garante convergencia
            S[i][tam+1] = rand() % 100;
            S[i][tam] = S[i][tam+1]/S[i][i];
            xk_1[i] = 0;
            somalinha = 0;
        }

    //Imprimindo S
    printf("S = \n");
    for(i=0;i<tam;i++)
    {
        for(j=0;j<tam+2;j++)
        {
            printf("%.1f ", S[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    /**METODO DE GAUSS-JACOBI**/
    int k = 0;
    double somaax = 0, precit = 1;
    double precnum = 0, precden = 0;

    //Começo das iterações
    while(precit > ERRO_ACEITAVEL)
    {
        //Iteração K+1

        #pragma omp parallel for num_threads(n_threads) private(i,j,somaax) shared(S,tam,precnum,precden,xk_1)
        for(i=0;i<tam;i++)
        {
            //printf("thread %d fara a linha %d\n", omp_get_thread_num(), i);
            for(j=0;j<tam;j++)
            {
                if(i != j){
                    somaax += (S[j][tam]*S[i][j]);

                }
            }
            xk_1[i] = (S[i][tam+1] - somaax)/S[i][i];


            //Calculando precisao
            #pragma omp critical (calculo_precisao)
            {
                if(fabs(xk_1[i]-S[i][tam]) > precnum)
                {
                    precnum = fabs(xk_1[i]-S[i][tam]);
                }
                if(fabs(xk_1[i]) > precden)
                {
                    precden = fabs(xk_1[i]);
                }
            }
            somaax = 0;
        }
        precit = precnum/precden; //precisao
        printf("precisao (k = %d) = %.5f\n",k+1,precit);
        //resetando variaveis para próxima iteração
        precnum = 0;
        precden = 0;

        for(i=0;i<tam;i++)
        {
            S[i][tam] = xk_1[i];
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
        soma += S[equa][j]*xk_1[j];
    }

    printf("resultado pelo metodo Gauss-Jacobi = %.3f\n", soma);
    printf("valor de b%d = %.3f\n", equa+1, S[equa][tam+1]);

    /**LIBERACAO DE MEMORIA**/

    free(xk_1);

    return 0;




















}
