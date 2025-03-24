#include "hub_problem.h"

// Definição Instancias e Numero de HUBs default
const char *default_instance = "inst200.txt";
int default_hub_count = 3;
// Defina aqui o tempo limite do GRASP
int grasp_time_limit_sec = 60;

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
// Tempo limite do código
int time_limit_sec = 330;

// Função de comparação para candidatos
int compare_candidates(const void *a, const void *b) {
    Candidate *ca = (Candidate *)a;
    Candidate *cb = (Candidate *)b;
    if (ca->dist > cb->dist) return -1;
    if (ca->dist < cb->dist) return 1;
    return 0;
}

// Função auxiliar para trocar elementos
void swap_elements(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Busca local usando VND (Variable Neighborhood Descent)
void local_search(Solution *sol) {
    double current_cost = sol->fo;
    int improved = 1;
    int neighborhood = 1;
    const int max_neighborhood = 3;
    
    while (improved && !stop) {
        improved = 0;
        
        switch (neighborhood) {
            case 1: improved = swap_search(sol); break;
            case 2: improved = insertion_search(sol); break;
            case 3: improved = replacement_search(sol); break;
        }
        
        if (improved) {
            neighborhood = 1;
            current_cost = sol->fo;
        } else if (neighborhood < max_neighborhood) {
            neighborhood++;
        } else {
            break;
        }
    }
}

// Vizinhança Swap (troca um hub com um não-hub)
int swap_search(Solution *sol) {
    double best_cost = sol->fo;
    int best_i = -1, best_j = -1;
    
    for (int i = 0; i < num_hubs; i++) {
        for (int j = num_hubs; j < num_nos; j++) {
            // Faz a troca temporária
            swap_elements(&sol->vet_sol[i], &sol->vet_sol[j]);
            
            // Calcula novo custo
            calculo_fo(*sol); // Passamos 0 pois não queremos imprimir
            
            // Verifica melhora
            if (sol->fo < best_cost) {
                best_cost = sol->fo;
                best_i = i;
                best_j = j;
            }
            
            // Desfaz a troca
            swap_elements(&sol->vet_sol[i], &sol->vet_sol[j]);
            sol->fo = melhor_fo; // Restaura o custo original
        }
    }
    
    // Aplica a melhor troca encontrada
    if (best_i != -1) {
        swap_elements(&sol->vet_sol[best_i], &sol->vet_sol[best_j]);
        sol->fo = best_cost;
        return 1;
    }
    
    return 0;
}

// Vizinhança de Inserção (substitui um hub por um não-hub)
int insertion_search(Solution *sol) {
    double best_cost = sol->fo;
    int best_hub_pos = -1;
    int best_node = -1;
    double current_cost = sol->fo;
    
    for (int i = 0; i < num_hubs; i++) {
        int current_hub = sol->vet_sol[i];
        
        for (int j = num_hubs; j < num_nos; j++) {
            int candidate = sol->vet_sol[j];
            
            // Substitui temporariamente
            sol->vet_sol[i] = candidate;
            
            // Calcula novo custo
            calculo_fo(*sol);
            
            // Verifica melhora
            if (sol->fo < best_cost) {
                best_cost = sol->fo;
                best_hub_pos = i;
                best_node = j;
            }
            
            // Restaura
            sol->vet_sol[i] = current_hub;
            sol->fo = current_cost;
        }
    }
    
    // Aplica a melhor inserção
    if (best_hub_pos != -1) {
        sol->vet_sol[best_hub_pos] = sol->vet_sol[best_node];
        sol->fo = best_cost;
        return 1;
    }
    
    return 0;
}

// Vizinhança de Substituição (troca dois hubs de posição)
int replacement_search(Solution *sol) {
    double best_cost = sol->fo;
    int best_i = -1, best_j = -1;
    double current_cost = sol->fo;
    
    for (int i = 0; i < num_hubs; i++) {
        for (int j = i+1; j < num_hubs; j++) {
            // Faz a troca temporária
            swap_elements(&sol->vet_sol[i], &sol->vet_sol[j]);
            
            // Calcula novo custo
            calculo_fo(*sol);
            
            // Verifica melhora
            if (sol->fo < best_cost) {
                best_cost = sol->fo;
                best_i = i;
                best_j = j;
            }
            
            // Desfaz a troca
            swap_elements(&sol->vet_sol[i], &sol->vet_sol[j]);
            sol->fo = current_cost;
        }
    }
    
    // Aplica a melhor troca
    if (best_i != -1) {
        swap_elements(&sol->vet_sol[best_i], &sol->vet_sol[best_j]);
        sol->fo = best_cost;
        return 1;
    }
    
    return 0;
}

void* time_limit(void* arg) {
    int limit = *(int*)arg;
    sleep(limit);  // Espera pelo tempo limite
    stop = 1;  // Sinaliza para a thread principal parar
    printf("Tempo limite atingido! Encerrando execução...\n");
    exit(1);
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
        return 0;  // Retorna 0 se houver erro ao abrir o arquivo
    }

    // Lê o número de nós e hubs
    if (fscanf(file, "n: %d\tp: %d\n", &num_nos, &num_hubs) != 2) {
        printf("Erro ao ler número de nós e hubs.\n");
        fclose(file);
        return 0;
    }

    // Lê o valor da função objetivo (FO)
    if (fscanf(file, "FO: %lf\n", &sol->fo) != 1) {
        printf("Erro ao ler o valor da função objetivo.\n");
        fclose(file);
        return 0;
    }

    // Lê os hubs
    if (fscanf(file, "HUBS: [") == EOF) {
        printf("Erro ao ler a lista de hubs.\n");
        fclose(file);
        return 0;
    }

    for (int i = 0; i < num_hubs; i++) {
        if (fscanf(file, "%d", &sol->vet_sol[i]) != 1) {
            printf("Erro ao ler os hubs na posição %d.\n", i);
            fclose(file);
            return 0;
        }
        // Ignora vírgula ou fecha colchetes na última iteração
        if (i < num_hubs - 1) fgetc(file);  // Pula a vírgula
    }

    // Pula o fechamento do colchete e a nova linha
    if (fgetc(file) != ']') {
        printf("Erro ao fechar o colchete dos hubs.\n");
        fclose(file);
        return 0;
    }
    fgetc(file);
    fclose(file);
    return 1;
}

void display_solution(const char *filename, Solution *sol) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Erro ao abrir o arquivo de solução");
        return;
    }

    char line[256]; // Buffer para ler linhas do arquivo

    // Lê e exibe linha por linha
    while (fgets(line, sizeof(line), file)) {
        // Remove o caractere de nova linha se existir
        line[strcspn(line, "\n")] = '\0';
        printf("%s\n", line);
    }

    fclose(file);
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

void calculo_fo(Solution& s) {
	
	if (((double)(clock() - start) / CLOCKS_PER_SEC) >= grasp_time_limit_sec - 1) {
        goto label;
	}
	
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
        melhor_tempo = ((double)(clock() - start) / CLOCKS_PER_SEC);

		int i = 0;
		printf("Hubs escolhidos: [");
    	for (i = 0; i < num_hubs - 1; i++) {
    	    printf("%d, ", s.vet_sol[i]);
		}
    	printf("%d]", s.vet_sol[i]);
    	
    	for (int i = 0; i < num_nos; i++) {
			vet_bin[i] = 0;
			for (int j = 0; j < num_hubs; j++)
				if ( i == s.vet_sol[j] )
					vet_bin[i] = 1;
		}
		
    	/*printf(" -> ")
		for (int j = 0; j < num_nos; j++) {
			printf("%d ", vet_bin[j]);
		}*/
		
		printf("\n");
		
    	printf("-> FO: %.2lf\nTempo: %.5lf\n", s.fo, melhor_tempo);
    	
        save_solution_details(s);  	
    }
    
	// Verifica tempo limite
    if (((double)(clock() - start) / CLOCKS_PER_SEC) >= grasp_time_limit_sec - 0.5) {
    	label:
    		
    	double time = ((double)(clock() - start) / CLOCKS_PER_SEC);
	    FILE *file = fopen("tempos.txt", "a");
	    if (!file) {
	        perror("Erro ao abrir tempos.txt");
	        exit(EXIT_FAILURE);
	    }
	    
	    fprintf(file, "inst%d.txt | p = %d -> FO: %.2lf -> Tempo(seg): %.2lf\n", 
	            num_nos, num_hubs, melhor_fo, melhor_tempo);
	    fclose(file);
	    printf("Meta-heurística finalizada em %.5lfseg\n", time);
	    exit(1);
	}        
}

void save_solution_details(Solution &s) {
	double Tij, Tjk = 0.75, Tkl;
	double valor;
	
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
    
	for ( int i = 0; i < num_nos; i++ ) {
		for ( int j = 0; j < num_nos; j++ ) {
			if ( vet_bin[j] == 0 ) {continue;}
			else if ( vet_bin[i] == 0 ) {Tij = 1;}
			else {Tij = 0.75;}
			for ( int k = 0; k < num_nos; k++ ) {
				if (vet_bin[k] == 0) {continue;}
				for ( int l = 0; l < num_nos; l++ ) {
					if (vet_bin[l] == 0) {Tkl = 1;}
					else {Tkl = 0.75;}
					valor = mat_distancias[i][j]*Tij + mat_distancias[j][k]*Tjk + mat_distancias[k][l]*Tkl;
					if ( valor == 0 ) {continue;}
					else if ( mat_custo[i][l] + 1 >= valor && valor >= mat_custo[i][l] - 1 )
						fprintf(file, "%d\t%d\t%d\t%d\t%.2lf\n", i, j, k, l, valor);
				}
			}
		}
	}
	
    fclose(file);
}

void* run_benchmark(void* arg) {
    read_instance(default_instance);
    Solution initial_sol;
    initialize_solution(&initial_sol);
    
    srand(time(NULL));
    
    start = clock();
    while (!stop) {
        // Fase de construção
        heu_cons_ale_gul(&initial_sol, 1);
        
        // Fase de busca local
        local_search(&initial_sol);
               
        // Avaliação da solução
        calculo_fo(initial_sol);
    }
    
    printf("Tempo limite atingido.\n");
    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    const char *instance_file = (argc > 1) ? argv[1] : default_instance;
    num_hubs = (argc > 2) ? atoi(argv[2]) : default_hub_count;

    Solution sol;
    
	/*if (read_solution("resultados.txt", &sol)) {
        printf("Solução carregada com sucesso!\n");
    } else {
        printf("Falha ao carregar a solução.\n");
    }*/
    
    // display_solution("resultados.txt", &sol);
	
	pthread_t timer_thread, algorithm_thread;

    // Criação das threads
    pthread_create(&timer_thread, NULL, time_limit, &time_limit_sec);
    pthread_create(&algorithm_thread, NULL, run_benchmark, NULL);

    // Espera pelas threads
    pthread_join(timer_thread, NULL);
    pthread_join(algorithm_thread, NULL);
	  
    return EXIT_SUCCESS;
}
