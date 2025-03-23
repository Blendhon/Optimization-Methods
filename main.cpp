#include "hub_problem.h"

// Definição das variáveis globais
Node nodes[MAX_NODES];
int num_nos;
int num_hubs;
double beta = 1.0, alfa = 0.75, lambda = 1.0;
double mat_distancias[MAX_NODES][MAX_NODES];
double mat_custo[MAX_NODES][MAX_NODES];

int melhor_hub[MAX_HUBS];
double melhor_fo = INFINITY;

// Implementação das funções
double calculate_distance(Node a, Node b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

/*void calculate_distance_matrix() {
    for (int i = 0; i < num_nos; i++) {
        for (int j = 1; j < num_nos; j++) {
            if (i == j) {
                mat_distancias[i][j] = 0;
            }
            else {
                mat_distancias[i][j] = calculate_distance(nodes[i], nodes[j]);
            }
            printf("%.2lf ", mat_distancias[i][j]);
        }
        printf("\n");
    }
}*/

void calculate_distance_matrix() {
    FILE *file = fopen("distancias.txt", "w"); // Abre o arquivo para escrita
    if (!file) {
        perror("Erro ao abrir o arquivo");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < num_nos; i++) {
    	mat_distancias[i][i] = 0;
        for (int j = i + 1; j < num_nos; j++) { // Começa de i + 1 para evitar duplicatas
            double distancia = calculate_distance(nodes[i], nodes[j]);
            mat_distancias[i][j] = distancia;
            mat_distancias[j][i] = distancia; // A matriz é simétrica

            // Escreve no arquivo no formato especificado
            fprintf(file, "Ponto %d - Ponto %d -> %.2lf\n", i, j, distancia);
        }
    }

    fclose(file); // Fecha o arquivo
    printf("Distancias salvas no arquivo 'distancias.txt'.\n");
}

void read_instance(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Erro ao abrir o arquivo");
        exit(EXIT_FAILURE);
    }
    fscanf(file, "%d", &num_nos);
    for (int i = 0; i < num_nos; i++) {
        fscanf(file, "%lf %lf", &nodes[i].x, &nodes[i].y);
    }
    fclose(file);
    calculate_distance_matrix();
}

void initialize_solution(Solution *sol) {
    sol->fo = 0;
    memset(sol->vet_sol, -1, sizeof(sol->vet_sol));  // Inicializa hubs com valores inválidos
    //memset(sol->hubs, -1, sizeof(sol->hubs));  // Inicializa hubs com valores inválidos
    //memset(sol->allocation, -1, sizeof(sol->allocation));  // Inicializa alocações
}

void clone_solution(Solution *original, Solution *clone) {
    memcpy(clone->vet_sol, original->vet_sol, sizeof(original->vet_sol));
    //memcpy(clone->hubs, original->hubs, sizeof(original->hubs));
    //memcpy(clone->allocation, original->allocation, sizeof(original->allocation));
    clone->fo = original->fo;
}

int read_solution(const char *filename, Solution *sol) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Erro ao abrir o arquivo de solução");
        return 0; // Falha ao abrir o arquivo
    }

    // Lê num_nos e p
    if (fscanf(file, "num_nos: %d\tp: %d\n", &num_nos, &num_hubs) != 2) {
        fclose(file);
        return 0; // Formato inválido
    }

    // Lê a FO
    if (fscanf(file, "FO: %lf\n", &sol->fo) != 1) {
        fclose(file);
        return 0; // Formato inválido
    }

    // Lê os hubs
    char line[256];
    fgets(line, sizeof(line), file); // Lê a linha "HUBS: [...]"
    char *token = strtok(line, " [],"); // Divide a linha em tokens
    while (token != NULL) {
        for(int i = 0; i < num_hubs; i++) {
            sol->vet_sol[i] = atoi(token); // Testar
        }
        //sol->hubs[num_hubs++] = atoi(token);
        token = strtok(NULL, " [],");
    }

    // Lê as alocações (OR, H1, H2, DS, CUSTO)
    while (fgets(line, sizeof(line), file)) {
        int ori, h1, h2, ds;
        double custo;
        if (sscanf(line, "%d\t%d\t%d\t%d\t%lf", &ori, &h1, &h2, &ds, &custo) == 5) {
            sol->allocation[ori] = h1; // Assumindo que H1 é o hub alocado para o nó OR
        }
    }

    fclose(file);
    return 1; // Sucesso
}

void heu_cons_ale_gul(Solution *sol, int use_random_seed) {
    //int available_hubs[MAX_NODES];

    for (int i = 0; i < num_nos; i++) {
        sol->vet_sol[i] = i;
        //available_hubs[i] = i;  
        //printf("%d ", available_hubs[i]);
        //printf("%d ", sol->vet_sol[i]);
    }
    //printf("\n");


    // Seleção de hubs
        /* int pos = 3;
        int temp = sol->vet_sol[0];
        sol->vet_sol[0] = sol->vet_sol[pos];
        sol->vet_sol[pos] = temp;

        int pos1 = 5;
        int temp1 = sol->vet_sol[1];
        sol->vet_sol[1] = sol->vet_sol[pos1];
        sol->vet_sol[pos1] = temp1;

        int pos2 = 13;
        int temp2 = sol->vet_sol[2];
        sol->vet_sol[2] = sol->vet_sol[pos2];
        sol->vet_sol[pos2] = temp2;

        int pos3 = 16;
        int temp3 = sol->vet_sol[3];
        sol->vet_sol[3] = sol->vet_sol[pos3];
        sol->vet_sol[pos3] = temp3; */


    for (int i = 0; i < num_hubs; i++) {

        int pos = (use_random_seed) ? rand() % (num_nos - i) : 0;
        int temp = sol->vet_sol[i];
        sol->vet_sol[i] = sol->vet_sol[pos];
        sol->vet_sol[pos] = temp;
    }
        //sol->hubs[i] = available_hubs[pos];
        //available_hubs[pos] = available_hubs[num_nos - i - 1];

    /*for (int i = 0; i < num_nos; i++) {
        //printf("%d ", available_hubs[i]);
        printf("%d ", sol->vet_sol[i]);
    }
    printf("\n");*/
    
    // Alocação
    /* for (int i = num_hubs; i < num_nos; i++) {
        double min_dist = INFINITY;
        for (int j = 0; j < num_hubs; j++) {
            //double dist = mat_distancias[i][sol->hubs[j]];
            double dist = mat_distancias[i][sol->vet_sol[j]];
            if (dist < min_dist) {
                min_dist = dist;
                //sol->allocation[i] = sol->hubs[j];
                sol->vet_sol[i] = sol->vet_sol[j];
                //sol->allocation[i] = sol->vet_sol[j];
            }
        }
    } */
	
	/*{
		sol->vet_sol[0] = 0;
		sol->vet_sol[1] = 2;
		sol->vet_sol[2] = 4;
		sol->vet_sol[3] = 1;
		sol->vet_sol[4] = 3;
	}*/
	
	/*{
		sol->vet_sol[0] = 3;
		sol->vet_sol[1] = 5;
		sol->vet_sol[2] = 13;
		sol->vet_sol[3] = 16;
		sol->vet_sol[4] = 0;
		sol->vet_sol[5] = 1;
		sol->vet_sol[6] = 2;
		sol->vet_sol[7] = 4;
		sol->vet_sol[8] = 6;
		sol->vet_sol[9] = 7;
		sol->vet_sol[10] = 8;
		sol->vet_sol[11] = 9;
		sol->vet_sol[12] = 10;
		sol->vet_sol[13] = 11;
		sol->vet_sol[14] = 12;
		sol->vet_sol[15] = 14;
		sol->vet_sol[16] = 15;
		sol->vet_sol[17] = 17;
		sol->vet_sol[18] = 18;
		sol->vet_sol[19] = 19;
	}*/
	
	/*printf("Hubs escolhidos: ");
    for (int i = 0; i < num_nos; i++) {
        //printf("%d ", available_hubs[i]);
        printf("%d ", sol->vet_sol[i]);
    }
    printf("\n");*/
}

void calculo_fo(Solution& s) {
    s.fo = 0;

    // Inicializar os custos entre nós não-hubs (não-hubs x não-hubs)
    for (int i = num_hubs; i < num_nos; i++) {
        for (int j = i; j < num_nos; j++) {
            mat_custo[s.vet_sol[i]][s.vet_sol[j]] = INT_MAX;
            mat_custo[s.vet_sol[j]][s.vet_sol[i]] = INT_MAX;
        }
    }

    // Inicializar os custos (hubs x hubs) e (hubs x não-hubs)
    for (int k = 0; k < num_hubs; k++) {
        // Hubs
        for (int l = k + 1; l < num_hubs; l++) {
            mat_custo[s.vet_sol[k]][s.vet_sol[l]] = alfa * mat_distancias[s.vet_sol[k]][s.vet_sol[l]];
            mat_custo[s.vet_sol[l]][s.vet_sol[k]] = alfa * mat_distancias[s.vet_sol[l]][s.vet_sol[k]];
        }
        // Não-hubs
        for (int i = num_hubs; i < num_nos; i++) {
            mat_custo[s.vet_sol[k]][s.vet_sol[i]] = lambda * mat_distancias[s.vet_sol[k]][s.vet_sol[i]];
            mat_custo[s.vet_sol[i]][s.vet_sol[k]] = beta * mat_distancias[s.vet_sol[i]][s.vet_sol[k]];
        }
    }

    // Atualizar os custos (hubs x não-hubs) passando por 1 hub intermediário
    for (int k = 0; k < num_hubs; k++) {
        for (int i = num_hubs; i < num_nos; i++) {
            for (int j = 0 + 1; j < num_hubs; j++) {
                mat_custo[s.vet_sol[k]][s.vet_sol[i]] = std::min(mat_custo[s.vet_sol[k]][s.vet_sol[i]],
                (mat_custo[s.vet_sol[k]][s.vet_sol[j]] + mat_custo[s.vet_sol[j]][s.vet_sol[i]]));
                mat_custo[s.vet_sol[i]][s.vet_sol[k]] = std::min(mat_custo[s.vet_sol[i]][s.vet_sol[k]],
                (mat_custo[s.vet_sol[i]][s.vet_sol[j]] + mat_custo[s.vet_sol[j]][s.vet_sol[k]]));
            }
        }
    }

    // Atualizar custos considerando dois hubs intermediários (melhor caminho via hubs)
    for (int i = num_hubs; i < num_nos; i++) {
        for (int j = i; j < num_nos; j++) {
            for (int k = 0; k < num_hubs; k++) {
              mat_custo[s.vet_sol[i]][s.vet_sol[j]] = std::min(mat_custo[s.vet_sol[i]][s.vet_sol[j]],
                  mat_custo[s.vet_sol[i]][s.vet_sol[k]] + mat_custo[s.vet_sol[k]][s.vet_sol[j]]);

              mat_custo[s.vet_sol[j]][s.vet_sol[i]] = std::min(mat_custo[s.vet_sol[j]][s.vet_sol[i]],
                  mat_custo[s.vet_sol[j]][s.vet_sol[k]] + mat_custo[s.vet_sol[k]][s.vet_sol[i]]);
              
                for (int l = 0; l < num_hubs; l++) {
                    mat_custo[s.vet_sol[i]][s.vet_sol[j]] = std::min(mat_custo[s.vet_sol[i]][s.vet_sol[j]],
                         mat_custo[s.vet_sol[i]][s.vet_sol[k]] + mat_custo[s.vet_sol[k]][s.vet_sol[l]] + mat_custo[s.vet_sol[l]][s.vet_sol[j]]);

                    mat_custo[s.vet_sol[j]][s.vet_sol[i]] = std::min(mat_custo[s.vet_sol[j]][s.vet_sol[i]],
                    mat_custo[s.vet_sol[j]][s.vet_sol[l]] + mat_custo[s.vet_sol[l]][s.vet_sol[k]] + mat_custo[s.vet_sol[k]][s.vet_sol[i]]);
                }
            }
        }
    }
    
    // Calcular a função objetivo (FO)
    for (int i = 0; i < num_nos; i++) {
        for (int j = i; j < num_nos; j++) {
            s.fo = std::max(s.fo, mat_custo[s.vet_sol[i]][s.vet_sol[j]]);
            //printf("%.2lf ", mat_custo[s.vet_sol[i]][s.vet_sol[j]]);
        }
        //printf("\n");
    }

    
    if (melhor_fo > s.fo) {
    
    	printf("Hubs escolhidos: ");
    	for (int i = 0; i < num_nos; i++) {
    	    //printf("%d ", available_hubs[i]);
    	    printf("%d ", s.vet_sol[i]);
    	}
    	printf("\n");
    	printf("%.2lf\n", s.fo);
    
    	melhor_fo = s.fo;
	    FILE *file = fopen("mat_custo.txt", "w"); // Abre o arquivo para escrita
	    if (!file) {
	        perror("Erro ao abrir o arquivo");
	        exit(EXIT_FAILURE);
	    }
	    
	    fprintf(file, "HUBS: [");
	    for (int i = 0; i < num_hubs; i++) {
	        fprintf(file, "%d%s", s.vet_sol[i], (i < num_hubs - 1) ? ", " : "");
	    }
	    fprintf(file, "]\n\n");
	    
	    for (int i = 0; i < num_nos; i++) {
	    	for (int j = 0; j < num_nos; j++) {
	    		if (mat_custo[i][j] < 10000)
	    			fprintf(file, "%.2lf\t\t|\t", mat_custo[i][j]);
	    		else
	    			fprintf(file, "%.2lf\t|\t", mat_custo[i][j]);
	    		
			}
	    	fprintf(file, "\n");
		}
	    
	    fprintf(file, "\n%lf", s.fo);
	    fclose(file); // Fecha o arquivo
	}
}

/* double menor_hub_final(int hub, Solution *sol) {
    double menor_hub_final = INFINITY;
    int temp = 0;

    for (int i = 0; i < num_nos; i++) {
        if (mat_distancias[sol->hubs[hub]][i] < menor_hub_final) {
            if (mat_distancias[sol->hubs[hub]][i] != 0) {
                menor_hub_final = mat_distancias[i][sol->hubs[hub]];
                temp = i;
            }
        }
    }

    return menor_hub_final;
} */

/* double menor_hub_hub(int hub, Solution *sol) {
    double menor_hub_hub = INFINITY;
    int temp = 0, hub_temp = 0;

    for (int i = 0; i < num_hubs; i++) {
        if (mat_distancias[sol->hubs[i]][sol->hubs[hub]] < menor_hub_hub) {
            if (mat_distancias[sol->hubs[i]][sol->hubs[hub]] != 0) {
                menor_hub_hub = mat_distancias[i][sol->hubs[hub]];
                temp = sol->hubs[i];
                hub_temp = i;
            }
        }
    }

    return menor_hub_hub * 0.75 + menor_hub_final(hub_temp, sol);
} */

/* double menor_ponto_hub(Solution *sol) {
    double menor_ponto_hub = INFINITY;
    int temp1 = 0, temp2 = 0, temp3;

    for (int i = 0; i < num_nos; i++) {
        for (int j = 0; j < num_hubs; j++) {
            if (mat_distancias[i][sol->hubs[j]] < menor_ponto_hub) {
                if (mat_distancias[i][sol->hubs[j]] != 0) {
                    menor_ponto_hub = mat_distancias[i][sol->hubs[j]];
                    temp1 = i;
                    temp2 = sol->hubs[j];
                    temp3 = j;
                }
            }
        }
    }

    double t = menor_ponto_hub + menor_hub_hub(temp3, sol);
    return t;
} */

void save_solution_details(const char *filename, Solution *sol) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Erro ao salvar a solução");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "n: %d\tp: %d\n", num_nos, num_hubs);
    fprintf(file, "FO: %.2lf\n", sol->fo);
    fprintf(file, "HUBS: [");
    for (int i = 0; i < num_hubs; i++) {
        fprintf(file, "%d%s", sol->vet_sol[i], (i < num_hubs - 1) ? ", " : "");
    }
    fprintf(file, "]\n");
    fprintf(file, "OR\tH1\tH2\tDS\tCUSTO\n");
    for (int i = 0; i < num_nos; i++) {
        for (int j = 0; j < num_nos; j++) {
            int h1 = sol->allocation[i];
            int h2 = sol->allocation[j];
            double tij = beta * mat_distancias[i][h1] + 
                         alfa * mat_distancias[h1][h2] + 
                         lambda * mat_distancias[h2][j];
            fprintf(file, "%d\t%d\t%d\t%d\t%.2lf\n", i, h1, h2, j, tij);
        }
    }
    fclose(file);
}

void display_solution(Solution *sol) {
    printf("\n--- Solucao ---\n");
    printf("n: %d\tp: %d\n", num_nos, num_hubs);
    printf("FO: %.2lf\n", sol->fo);
    printf("HUBS: [");
    for (int i = 0; i < num_hubs; i++) {
        printf("%d%s", sol->vet_sol[i], (i < num_hubs - 1) ? ", " : "");
    }
    printf("]\n");
    printf("OR\tH1\tH2\tDS\tCUSTO\n");
    for (int i = 0; i < num_nos; i++) {
        for (int j = 0; j < 5; j++) {
            int h1 = sol->allocation[i];
            int h2 = sol->allocation[j];
            double tij = beta * mat_distancias[i][h1] + 
                         alfa * mat_distancias[h1][h2] + 
                         lambda * mat_distancias[h2][j];
            printf("%d\t%d\t%d\t%d\t%.2lf\n", i, h1, h2, j, tij);
        }
    }
}

void run_benchmark(const char *filename, int iterations) {
    // Configuração inicial
    read_instance(filename);
    
    // Modo de execução única
    Solution initial_sol;
    initialize_solution(&initial_sol);
    //clock_t start = clock();
    
    for (int i = 0; i < iterations; i++) {
	    heu_cons_ale_gul(&initial_sol, 1);
	 
	    //initial_sol.fo = calculo_fo(initial_sol);
	    calculo_fo(initial_sol);
	}

    //double time_single = (double)(clock() - start) / CLOCKS_PER_SEC;

    // Salvar e exibir solução inicial
    save_solution_details("solucao_inicial.txt", &initial_sol);
    // display_solution(&initial_sol);

    // Modo iterativo
    /* double total_time_heuristic = 0, total_time_fo = 0;
    double maior = 0;
    double temp = 0;

    clock_t start_heu = clock();
    for(int i = 0; i < iterations; i++) {
        Solution temp_sol;
        initialize_solution(&temp_sol);
        heu_cons_ale_gul(&temp_sol, 1);
    }
    total_time_heuristic = (double)(clock() - start_heu) / CLOCKS_PER_SEC;

    clock_t start_fo = clock();
    for (int j = 0; j < iterations; j++) {
        Solution temp_sol;
        initialize_solution(&temp_sol);
        heu_cons_ale_gul(&temp_sol, 1);
        calculo_fo(&temp_sol);
        temp = menor_ponto_hub(&initial_sol);
        if (maior < temp) {
            maior = temp;
        }
    }
    total_time_fo = (double)(clock() - start_fo) / CLOCKS_PER_SEC;

    calculo_fo(&initial_sol);
    printf("\nDescricao do Computador\t?\t\t(seg.)\t\tFO\t\tTempo (seg.)\t\tTempo Sol. Inicial (seg.)\tTempo Calc. FO (seg.)\n");
    printf("R5 4.4GHz - 16GB ram - Windows 11\t97\t\t%.2lf\t%.5lf\t\t\t%.5lf\t\t\t\t%.5lf\n\n", 
           initial_sol.fo,
           time_single,
           total_time_heuristic,
           total_time_fo);
    
    printf("Maior custo encontrado (Funcao Objetivo): %.2lf\n", maior); */

}

int main(int argc, char *argv[]) {
    const char *default_instance = "inst200.txt";
    int default_hub_count = 50;

    const char *instance_file = (argc > 1) ? argv[1] : default_instance;
    num_hubs = (argc > 2) ? atoi(argv[2]) : default_hub_count;

    srand(time(NULL));

    Solution sol;
    if (read_solution("solucao_salva.txt", &sol)) {
        printf("Solução carregada com sucesso!\n");
        display_solution(&sol);
    } else {
        printf("Falha ao carregar a solução.\n");
    }
	
	run_benchmark(instance_file, 1000);
    
    return EXIT_SUCCESS;
}
