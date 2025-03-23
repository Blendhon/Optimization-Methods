#include "hub_problem.h"

#include <pthread.h>
#include <unistd.h>

const char *default_instance = "inst100.txt";
int default_hub_count = 50;

volatile int stop = 0;  // Variável compartilhada para sinalizar parada

void* time_limit(void* arg) {
    int limit = *(int*)arg;
    sleep(limit);  // Espera pelo tempo limite
    stop = 1;  // Sinaliza para a thread principal parar
    printf("Tempo limite atingido! Encerrando execução...\n");
    pthread_exit(NULL);
}

// Definição das variáveis globais
Node nodes[MAX_NODES];
int num_nos;
int num_hubs;
double beta = 1.0, alfa = 0.75, lambda = 1.0;
double mat_distancias[MAX_NODES][MAX_NODES];
double mat_custo[MAX_NODES][MAX_NODES];

int melhor_hub[MAX_HUBS];
double melhor_fo = INFINITY;

double calculate_distance(Node a, Node b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

void calculate_distance_matrix() {
/*	FILE *file = fopen("distancias.txt", "w");
	if (!file) {
	    perror("Erro ao abrir o arquivo");
	    exit(EXIT_FAILURE);
	}*/
	
	for (int i = 0; i < num_nos; i++) {
		mat_distancias[i][i] = 0;
		for (int j = i + 1; j < num_nos; j++) { // Começa de i + 1 para evitar duplicatas
			double distancia = calculate_distance(nodes[i], nodes[j]);
			mat_distancias[i][j] = distancia;
			mat_distancias[j][i] = distancia; // A matriz é simétrica
			
			// Escreve no arquivo no formato especificado
			// fprintf(file, "Ponto %d - Ponto %d -> %.2lf\n", i, j, distancia);
		}
	}

	/*fclose(file);
	printf("Distancias salvas no arquivo 'distancias.txt'.\n");*/
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
    for (int i = 0; i < num_nos; i++) {
        sol->vet_sol[i] = i;
    }

    int center_index = 0;
    int min_max_dist = 0x7F800000; // Use um valor grande para inicializar

    for (int i = 0; i < num_nos; i++) {
        int max_dist = 0;
    // Encontra a distância máxima do nó i para qualquer outro nó
        for (int j = 0; j < num_nos; j++) {
            if (mat_distancias[i][j] > max_dist) {
                max_dist = mat_distancias[i][j];
            }
        }
    // Atualiza o centro se a distância máxima for menor
        if (max_dist < min_max_dist) {
            min_max_dist = max_dist;
            center_index = i;
        }
    }

    for (int i = 0; i < num_hubs; i++) {
        int pos;
        if (use_random_seed) {
            // Combinação de estratégia gulosa e aleatória
            // Seleciona aleatoriamente entre os nós mais distantes do centro
            int num_candidates = (num_nos - i) / 2; // Ajuste o número de candidatos conforme necessário
            int farthest_candidate = i;
            double max_dist = 0.0;

        // Encontra os 'num_candidates' nós mais distantes do centro
            for (int j = i; j < num_nos; j++) {
                double dist = mat_distancias[j][center_index]; // Supondo que center_index é o índice do centro
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
                double dist = mat_distancias[j][center_index]; // Supondo que center_index é o índice do centro
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

    /*for (int i = 0; i < num_nos; i++) {
        //printf("%d ", available_hubs[i]);
        printf("%d ", sol->vet_sol[i]);
    }
    printf("\n");*/
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
            //printf("%.2lf ", mat_custo[s.vet_sol[i]][s.vet_sol[j]]);
        }
        //printf("\n");
    }

    //printf("Interação: %d > FO Atual: %.2lf -> Nova FO: %.2lf\n", iterations, melhor_fo, s.fo);
 /*   printf("Hubs escolhidos: ");
    	for (int i = 0; i < num_nos; i++) {
    	    printf("%d ", s.vet_sol[i]);
		}
    	printf("\n");*/
    
    if (melhor_fo >= s.fo) {
    	melhor_fo = s.fo;

    /*	printf("Hubs escolhidos: ");
    	for (int i = 0; i < num_nos; i++) {
    	    //printf("%d ", available_hubs[i]);
    	    printf("%d ", s.vet_sol[i]);
		}
    	printf("\n");*/
    	printf("N_Int: %d -> FO: %.2lf\n", iterations, s.fo);

	/*    FILE *file = fopen("mat_custo.txt", "w"); // Abre o arquivo para escrita
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
	    fclose(file);*/
	}
}

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
    pthread_exit(NULL);

    //double time_single = (double)(clock() - start) / CLOCKS_PER_SEC;

    // Salvar e exibir solução inicial
    // save_solution_details("solucao_inicial.txt", &initial_sol);
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
    /*if (read_solution("solucao_salva.txt", &sol)) {
        printf("Solução carregada com sucesso!\n");
        display_solution(&sol);
    } else {
        printf("Falha ao carregar a solução.\n");
    }*/
	
	pthread_t timer_thread, algorithm_thread;
    int time_limit_sec = 19;  // Defina aqui o tempo limite

    // Criação das threads
    pthread_create(&timer_thread, NULL, time_limit, &time_limit_sec);
    pthread_create(&algorithm_thread, NULL, run_benchmark, NULL);

    // Espera pelas threads
    pthread_join(timer_thread, NULL);
    pthread_join(algorithm_thread, NULL);
	
	/*int i = 0;
	read_instance(instance_file);
    Solution initial_sol;
    initialize_solution(&initial_sol);
	while (melhor_fo > 38000) {
	    i++;
		heu_cons_ale_gul(&initial_sol, 1);
	    calculo_fo(initial_sol, i);
	}*/
    
    return EXIT_SUCCESS;
}
