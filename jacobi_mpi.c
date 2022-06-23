#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"


int main(int argc, char *argv[]){
    int npes, myrank, src, dest, msgtag, ret, tam, seed, i, j;
    double *xk;
    double somalinha = 0;
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
	 if ( argc  != 2)
    {
	printf("Argumentos errados. Utilize <ordem da matriz> \n");
	exit(0);
    }

    tam = atoi(argv[1]);
    seed = 2;
    xk = (double*)malloc(tam*sizeof(double));

    /**PROCESSO 0**/
    if ( myrank == 0) {
        printf("Existem %d processos.\n", npes);
        printf("Eu sou o processo %d no processador %s. As matrizes A e b sao as seguintes \n", myrank, processor_name);

        /**INICIALIZANDO E PREENCHENDO AS MATRIZES A E b**/
        double S[tam][tam+1];

        srand(seed);

        for(i=0;i<tam;i++)
        {
                for(j=0;j<tam;j++)
                {
                    if(i != j)
                    {
                        S[i][j] = rand() %100;
                        somalinha += S[i][j];
                    }
                }
            S[i][i] = somalinha + (rand()%10) + 1; //Garante convergencia
            S[i][tam] = rand() % 100;
            somalinha = 0;
        }

        for(i=0;i<tam;i++){
            xk[i] = S[i][tam]/S[i][i];
            printf("%.1f ",xk[i]);
        }

        //Imprimindo S
        printf("\nS = \n");
        for(i=0;i<tam;i++)
        {
            for(j=0;j<tam+1;j++)
            {
                printf("%.1f ", S[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        /**CALCULANDO QUANTAS LINHAS CADA PROCESSO RECEBERA**/
        int nlinhas, resto, ordem;

        nlinhas = tam/(npes-1);
        resto = tam % (npes-1);


        /**ENVIANDO AS LINHAS DAS MATRIZES A E b**/
        int linhaidx = 0;
        msgtag = 1;
        ordem = nlinhas;

		for (dest = 1; dest<npes-1; dest++){ //O resto da divisão de linhas será mandado para o último processo
            for(i=0;i<ordem;i++){
                MPI_Send(&S[linhaidx], tam+1, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
                linhaidx++;
            }
            MPI_Send(xk, tam, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
		}

		ordem = nlinhas+resto;
		dest = npes-1;
        for(i=0;i<ordem;i++){ //Ultimo processo recebera o resto das linhas que sobraram
                MPI_Send(&S[linhaidx], tam+1, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);
                linhaidx++;
        }
        MPI_Send(xk, tam, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);

        /**RECEBENDO RESULTADO (xk_1)**/
        ordem = nlinhas;
        bufrecv = (double*)malloc((tam)*sizeof(double));
        for(src=1;src<npes-1;src++){
            MPI_Recv(bufrecv, ordem, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
            for(i=(src-1)*ordem;i<src*ordem;i++){
                xk[i] = bufrecv[i];
                printf("%.1f ", xk[i]);
            }
        }
        ordem = nlinhas+resto;
        src = npes-1
        MPI_Recv(bufrecv, ordem, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
        for(i=(src-1)*nlinhas;i<src*ordem;i++){
                xk[i] = bufrecv[i];
                printf("%.1f ", xk[i]);
            }

	}
	else{
        /**CALCULANDO QUANTAS LINHAS CADA PROCESSO RECEBERA**/
        int nlinhas, resto, ordem;

        nlinhas = tam/(npes-1);
        resto = tam % (npes-1);

        if(myrank != npes-1)
            ordem = nlinhas;
        else
            ordem = nlinhas+resto;

        /**RECEBER AS LINHAS DAS MATRIZES A E b**/
        src = 0;
        msgtag = 1;
        bufrecv = (double*)malloc((tam+1)*sizeof(double));
        double S[ordem][tam+1];
        for(i=0;i<ordem;i++){
            MPI_Recv(bufrecv, tam+1, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
            for(j=0;j<tam+1;j++){
                S[i][j] = bufrecv[j];
                printf("%.1f ", S[i][j]);
            }
            printf("S[%d](processo %d)\n",i, myrank);
        }

        /**RECEBER xk**/
        MPI_Recv(bufrecv, tam+1, MPI_DOUBLE, src, msgtag, MPI_COMM_WORLD, &status);
        for(i=0;i<tam;i++){
            xk[i] = bufrecv[i];
            printf("%.1f ", xk[i]);
        }

        /**METODO DE GAUSS-JACOBI**/
        int somaax;
        double *xk_1;
        xk_1 = (double*)malloc(ordem*sizeof(double));

        for(i=0;i<ordem;i++)
        {
            for(j=0;j<tam;j++)
            {
                if(i != j){
                    somaax += (xk[j]*S[i][j]);

                }
            }
            xk_1[i] = (S[i][tam] - somaax)/S[i][i];
        }
        printf("xk_1 = ");
        for(i=0;i<ordem;i++)
            printf("%.1f ", xk_1[i]);
        printf("xk_1 (processo %d)\n", myrank);

        /**ENVIAR RESULTADO PARA P0**/
        dest = 0;
        MPI_Send(xk_1, ordem, MPI_DOUBLE, dest, msgtag, MPI_COMM_WORLD);

        free(xk_1);
	}
	free(bufrecv);
    printf("\n\n");
	fflush(0);
	ret = MPI_Finalize();
	if (ret == MPI_SUCCESS)
		printf("MPI_Finalize success! \n");
	return(0);

}
