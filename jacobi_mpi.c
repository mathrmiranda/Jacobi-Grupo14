#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define ERRO_ACEITAVEL 0.00001

int main(int argc, char *argv[]){
    int npes, myrank, src, dest, msgtag, ret, tam, seed, i, j,k = 0, flag_parada = 0, n_threads;
    double *xk,*xk_1;
    double somalinha = 0, precnum = 0, precden = 0, precit;
	double *bufrecv;

	/**INICIALIZANDO MPI**/
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	MPI_Get_processor_name(processor_name, &name_len);

	/**RECEBENDO ARGUMENTOS**/
	 if ( argc  != 3)
    {
	printf("Argumentos errados. Utilize <ordem da matriz> <numero de threads>\n");
	exit(0);
    }

    tam = atoi(argv[1]);
    n_threads = atoi(argv[2]);
    seed = 2;
    xk = (double*)malloc((tam)*sizeof(double));

    //printf("Processo %d executando no processador %s\n", myrank, processor_name);

    /**PROCESSO 0**/
    if ( myrank == 0) {
        
        /**INICIALIZANDO E PREENCHENDO AS MATRIZES S**/
        double S[tam][tam+1];

        srand(seed);

        #pragma omp parallel for num_threads(n_threads) private(i,j,somalinha) shared(S,tam)
            for(i=0;i<tam;i++)
            {
                    for(j=0;j<tam;j++)
                    {
                        if(i != j)
                        {
                            S[i][j] = rand() %10;
                            somalinha += S[i][j];
                        }
                    }
                S[i][i] = somalinha + (rand()%10) + 1; //Garante convergencia
                S[i][tam] = rand() % 100;
                somalinha = 0;
            }


            for(i=0;i<tam;i++){
                xk[i] = S[i][tam]/S[i][i]; //Vetor x de resposta da iteracao 0
            }

        /**CALCULANDO QUANTAS LINHAS CADA PROCESSO RECEBERA**/
        int nlinhas, resto, ordem;

        nlinhas = tam/(npes-1);
        resto = tam % (npes-1);

        /**ENVIANDO AS LINHAS DAS MATRIZES S**/
        int linhaidx = 0;
        msgtag = 1;
        ordem = nlinhas;

        for (dest = 1; dest<npes-1; dest++){ //O resto da divisao de linhas ser mandado para o ultimo processo
            for(i=0;i<ordem;i++){
                MPI_Send(&S[linhaidx], tam+1, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
                linhaidx++;
            }
                
        }

        ordem = nlinhas+resto;
        dest = npes-1;
        for(i=0;i<ordem;i++){ //Ultimo processo recebera o resto das linhas que sobraram
                MPI_Send(&S[linhaidx], tam+1, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
                linhaidx++;
        }

        bufrecv = (double*)malloc((nlinhas+resto)*sizeof(double));
        xk_1 = (double*)malloc((tam)*sizeof(double));
        int count = 0;
        
        while(flag_parada == 0){
            
            /**ENVIANDO XK**/
            MPI_Bcast(xk,tam,MPI_DOUBLE,0,MPI_COMM_WORLD);
            
            /**RECEBENDO RESULTADO (xk_1)**/
            
            count = 0;
            for(src=1;src<(npes-1);src++){
                MPI_Recv(bufrecv, nlinhas, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
                for(i=0;i<nlinhas;i++){
                    xk_1[count] = bufrecv[i];
                    count++;
                }
            }
            MPI_Recv(bufrecv, nlinhas+resto, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
            for(i=0;i<nlinhas+resto;i++){
                    xk_1[count] = bufrecv[i];
                    count++;
                }

            //Calculando precisao
            for(i=0;i<tam;i++){
                if(fabs(xk_1[i]-xk[i]) > precnum)
                {
                    precnum = fabs(xk_1[i]-xk[i]);
                }

                if(fabs(xk_1[i]) > precden)
                {
                    precden = fabs(xk_1[i]);
                }         
            }
            precit = precnum/precden;
            printf("erro = %.6f iteracao %d\n", precit,k);

            if(precit > ERRO_ACEITAVEL){
                flag_parada = 0;
            }
                
           
            if(precit < ERRO_ACEITAVEL){
                flag_parada = 1;
            }

            

            for(i=0;i<tam;i++){
                xk[i] = xk_1[i];
            }
            precnum = precden = 0;
            k++;
            MPI_Bcast(&flag_parada,1,MPI_INT,0,MPI_COMM_WORLD);
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
        printf("valor de b%d = %.3f\n", equa+1, S[equa][tam]);
        
	}
	else{

        /**CALCULANDO QUANTAS LINHAS CADA PROCESSO RECEBERA**/
        int nlinhas, resto, ordem;
        double somaax = 0;
        nlinhas = tam/(npes-1);
        resto = tam % (npes-1);

        if(myrank != npes-1)
            ordem = nlinhas;
        else
            ordem = nlinhas+resto;

        /**RECEBER AS LINHAS DA MATRIZ S**/
        src = 0;
        msgtag = 1;
        bufrecv = (double*)malloc((tam+1)*sizeof(double));
        double S[ordem][tam+1];
        for(i=0;i<ordem;i++){
            MPI_Recv(bufrecv, tam+1, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
            for(j=0;j<tam+1;j++){
                S[i][j] = bufrecv[j];
            }
        }
        xk_1 = (double*)malloc(ordem*sizeof(double));
        while(flag_parada == 0){

            MPI_Bcast(xk,tam,MPI_DOUBLE,0,MPI_COMM_WORLD);

            /**METODO DE GAUSS-JACOBI**/
            
            
            #pragma omp parallel for num_threads(n_threads) private(i,j) shared(S,tam,xk_1,xk) reduction(+:somaax)
            for(i=0;i<ordem;i++)
            {
                somaax = 0;
                for(j=0;j<tam;j++)
                {
                    if(j != ((myrank-1)*nlinhas)+i){
                        somaax += (xk[j]*S[i][j]);
                    }
                }
                xk_1[i] = (S[i][tam] - somaax)/S[i][((myrank-1)*nlinhas)+i];
            }

            /**ENVIAR RESULTADO PARA P0**/
            int dest = 0;
            MPI_Send(xk_1, ordem, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);

            MPI_Bcast(&flag_parada,1,MPI_INT,0,MPI_COMM_WORLD);
        }
        
	}

    free(xk);
    free(xk_1);
	free(bufrecv);
	fflush(0);
	ret = MPI_Finalize();
	if (ret == MPI_SUCCESS)
		printf("MPI_Finalize success! \n");
	return(0);

}
