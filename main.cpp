#include "hub_problem.h"

// Definição Instancias e Numero de HUBs default
const char *default_instance = "inst50.txt";
int default_hub_count = 10;
// Defina aqui o tempo limite
int time_limit_sec = 290;

// Definição das variáveis globais
Node nodes[MAX_NODES];
int num_nos;
int num_hubs;
double beta = 1.0, alfa = 0.75, lambda = 1.0;
double mat_distancias[MAX_NODES][MAX_NODES];
double mat_custo[MAX_NODES][MAX_NODES];
int melhor_hub[MAX_HUBS];
double melhor_fo = INFINITY;
volatile int stop = 0;

void* time_limit(void* arg) {
    int limit = *(int*)arg;
    sleep(limit);  // Espera pelo tempo limite
    stop = 1;  // Sinaliza para a thread principal parar
    printf("Tempo limite atingido! Encerrando execução...\n");
    pthread_exit(NULL);
}

double calculate_distance(Node a, Node b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

void calculate_distance_matrix() {
	for (int i = 0; i < num_nos; i++) {
		mat_distancias[i][i] = 0;
		for (int j = i + 1; j < num_nos; j++) { // Começa de i + 1 para evitar duplicatas
			double distancia = calculate_distance(nodes[i], nodes[j]);
			mat_distancias[i][j] = distancia;
			mat_distancias[j][i] = distancia; // A matriz é simétrica
		}
	}
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
}

void clone_solution(Solution *original, Solution *clone) {
    memcpy(clone->vet_sol, original->vet_sol, sizeof(original->vet_sol));
    clone->fo = original->fo;
}

int read_solution(const char *filename, Solution *sol) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Erro ao abrir o arquivo de solução");
        return 0;
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
    return 1;
}

void heu_cons_ale_gul(Solution *sol, int use_random_seed) {
    // Inicializa a solução com os índices dos nós
    for (int i = 0; i < num_nos; i++) {
        sol->vet_sol[i] = i;
    }

    // Encontra o nó central
    int center_index = 0;
    int min_max_dist = INT_MAX;

    for (int i = 0; i < num_nos; i++) {
        int max_dist = 0;
        for (int j = 0; j < num_nos; j++) {
            if (mat_distancias[i][j] > max_dist) {
                max_dist = mat_distancias[i][j];
            }
        }
        if (max_dist < min_max_dist) {
            min_max_dist = max_dist;
            center_index = i;
        }
    }

    // Seleciona os hubs
    for (int i = 0; i < num_hubs; i++) {
	    int pos;
	    if (use_random_seed) {
	        // Combinação de estratégia gulosa e aleatória
	        // Considera todos os nós restantes como candidatos
	        int num_candidates = num_nos - i;
	        int farthest_candidate = i;
	        double max_dist = 0.0;
	
	        for (int j = i; j < num_nos; j++) {
	            double dist = mat_distancias[j][center_index];
	            if (dist > max_dist) {
	                max_dist = dist;
	                farthest_candidate = j;
	            }
	        }
	
	        // Seleciona aleatoriamente entre os candidatos mais distantes
	        pos = i + rand() % num_candidates;
	    } else {
	        // Estratégia puramente gulosa (seleciona o mais distante)
	        pos = i;
	        double max_dist = 0.0;
	        for (int j = i; j < num_nos; j++) {
	            double dist = mat_distancias[j][center_index];
	            if (dist > max_dist) {
	                max_dist = dist;
	                pos = j;
	            }
	        }
	    }
	
	    // Troca os elementos
	    int temp = sol->vet_sol[i];
	    sol->vet_sol[i] = sol->vet_sol[pos];
	    sol->vet_sol[pos] = temp;
	}	
}

void calculo_fo(Solution& s, int iterations) {
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
        for (int l = k; l <= num_hubs; l++) {
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
            for (int j = 0; j < num_hubs; j++) {
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
        }
    }

    // Atualizar a melhor FO, se necessário
    if (melhor_fo > s.fo) {
        melhor_fo = s.fo;

		printf("Hubs escolhidos: ");
    	for (int i = 0; i < num_nos; i++) {
    	    printf("%d ", s.vet_sol[i]);
		}
    	printf("\n");
    	printf("N_Int: %d -> FO: %.2lf\n", iterations, s.fo);

        save_solution_details(s);
    }    
}

void save_solution_details(Solution &s) {
    FILE *file = fopen("resultados.txt", "w");
    if (!file) {
        perror("Erro ao salvar a solução");
        exit(EXIT_FAILURE);
    }
    
    // Escreve o cabeçalho
    fprintf(file, "n: %d\tp: %d\n", num_nos, num_hubs);
    fprintf(file, "FO: %.2lf\n", s.fo);
    fprintf(file, "HUBS: [");
    for (int i = 0; i < num_hubs; i++) {
        fprintf(file, "%d%s", s.vet_sol[i], (i < num_hubs - 1) ? ", " : "");
    }
    fprintf(file, "]\n");
    fprintf(file, "OR\tH1\tH2\tDS\tCUSTO\n");
    
    for ( int ori = 0; ori < num_hubs; ori++ ) {
    	for ( int h1 = 0; h1 < num_hubs; h1++ ) {
    		for ( int h2 = 0; h2 < num_hubs; h2++ ) {
    			for ( int des = 0; des < num_hubs; des++ ) {
    				if ( ori == h1 && h1 == h2 && h2 == des);
					else if ( ori == des && h1 == h2 )
						fprintf(file, "%d\t%d\t%d\t%d\t%.2lf\n",
							s.vet_sol[ori], s.vet_sol[h1], s.vet_sol[h2], s.vet_sol[des],
							alfa*mat_distancias[s.vet_sol[ori]][s.vet_sol[h1]]);
				}
			}
		}
	}
    
    for ( int ori = num_hubs; ori < num_nos; ori++ ) {
    	for ( int h1 = 0; h1 < num_hubs; h1++ ) {
    		for ( int h2 = 0; h2 < num_hubs; h2++ ) {
    			for ( int des = num_hubs; des < num_nos; des++ ) {
    				if ( mat_custo[ori][des] + 5 >= beta*mat_distancias[s.vet_sol[ori]][s.vet_sol[h1]] +
							alfa*mat_distancias[s.vet_sol[h1]][s.vet_sol[h2]] +
							lambda*mat_distancias[s.vet_sol[h2]][s.vet_sol[des]] || beta*mat_distancias[s.vet_sol[ori]][s.vet_sol[h1]] +
							alfa*mat_distancias[s.vet_sol[h1]][s.vet_sol[h2]] +
							lambda*mat_distancias[s.vet_sol[h2]][s.vet_sol[des]] < s.fo + 5 )
	    				fprintf(file, "%d\t%d\t%d\t%d\t%.2lf\n",
							s.vet_sol[ori], s.vet_sol[h1], s.vet_sol[h2], s.vet_sol[des],
							beta*mat_distancias[s.vet_sol[ori]][s.vet_sol[h1]] +
							alfa*mat_distancias[s.vet_sol[h1]][s.vet_sol[h2]] +
							lambda*mat_distancias[s.vet_sol[h2]][s.vet_sol[des]]);
					else if ( mat_custo[ori][des] + 5 >= beta*mat_distancias[s.vet_sol[ori]][s.vet_sol[h1]] +
							alfa*mat_distancias[s.vet_sol[h1]][s.vet_sol[h2]] )
	    				fprintf(file, "%d\t%d\t%d\t%d\t%.2lf\n",
							s.vet_sol[ori], s.vet_sol[h1], s.vet_sol[h2], s.vet_sol[des],
							beta*mat_distancias[s.vet_sol[ori]][s.vet_sol[h1]] +
							alfa*mat_distancias[s.vet_sol[h1]][s.vet_sol[h2]]);
					else if ( mat_custo[ori][des] + 5 >= alfa*mat_distancias[s.vet_sol[h1]][s.vet_sol[h2]] +
							lambda*mat_distancias[s.vet_sol[h2]][s.vet_sol[des]] )
	    				fprintf(file, "%d\t%d\t%d\t%d\t%.2lf\n",
							s.vet_sol[ori], s.vet_sol[h1], s.vet_sol[h2], s.vet_sol[des],
							alfa*mat_distancias[s.vet_sol[h1]][s.vet_sol[h2]] +
							lambda*mat_distancias[s.vet_sol[h2]][s.vet_sol[des]]);
				}
			}
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

void* run_benchmark(void* arg) {
    // Configuração inicial
    read_instance(default_instance);
    
    // Modo de execução única
    Solution initial_sol;
    initialize_solution(&initial_sol);
    //clock_t start = clock();
    
    int i = 0;
    while (!stop) {
	    heu_cons_ale_gul(&initial_sol, 1);
	    calculo_fo(initial_sol, i + 1);
	    i++;
	}
	printf("Meta-heurística finalizada.\n");
	// Salvar e exibir solução inicial
    // save_solution_details("solucao_inicial.txt", initial_sol);
    pthread_exit(NULL);

    //double time_single = (double)(clock() - start) / CLOCKS_PER_SEC;

    
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

    const char *instance_file = (argc > 1) ? argv[1] : default_instance;
    num_hubs = (argc > 2) ? atoi(argv[2]) : default_hub_count;

    srand(time(NULL));

    Solution sol;
   /* if (read_solution("solucao_salva.txt", &sol)) {
        printf("Solução carregada com sucesso!\n");
        display_solution(&sol);
    } else {
        printf("Falha ao carregar a solução.\n");
    }*/
	
	pthread_t timer_thread, algorithm_thread;

    // Criação das threads
    pthread_create(&timer_thread, NULL, time_limit, &time_limit_sec);
    pthread_create(&algorithm_thread, NULL, run_benchmark, NULL);

    // Espera pelas threads
    pthread_join(timer_thread, NULL);
    pthread_join(algorithm_thread, NULL);
	  
    return EXIT_SUCCESS;
}
