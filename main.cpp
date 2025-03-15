#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_NODES 200
#define MAX_HUBS 50
#define INFINITY 1e9

typedef struct {
    double x, y;
} Node;

typedef struct {
    int hubs[MAX_HUBS];
    int hub_count;
    double objective_function;
    int allocation[MAX_NODES];
    clock_t start_time;
    clock_t end_time;
} Solution;

Node nodes[MAX_NODES];
int n, p;
double beta = 1.0, alpha = 0.75, lambda = 1.0;
double distance_matrix[MAX_NODES][MAX_NODES];
double vetor_hub[MAX_HUBS];
double vetor_nao_hub[MAX_NODES];

// Funções básicas
double calculate_distance(Node a, Node b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

void calculate_distance_matrix() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                distance_matrix[i][j] = 0;
            }
            else {
                distance_matrix[i][j] = calculate_distance(nodes[i], nodes[j]);
            }
        }
    }
}

// Leitura de instância com parâmetros
void read_instance(const char *filename, int hub_count) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Erro ao abrir o arquivo");
        exit(EXIT_FAILURE);
    }
    fscanf(file, "%d", &n);
    for (int i = 0; i < n; i++) {
        fscanf(file, "%lf %lf", &nodes[i].x, &nodes[i].y);
    }
    fclose(file);
    p = hub_count;
    calculate_distance_matrix();
}

// Funções de solução
void initialize_solution(Solution *sol) {
    sol->hub_count = p;
    sol->objective_function = 0;
    memset(sol->hubs, -1, sizeof(sol->hubs));  // Inicializa hubs com valores inválidos
    memset(sol->allocation, -1, sizeof(sol->allocation));  // Inicializa alocações
}

void clone_solution(Solution *original, Solution *clone) {
    memcpy(clone->hubs, original->hubs, sizeof(original->hubs));
    memcpy(clone->allocation, original->allocation, sizeof(original->allocation));
    clone->objective_function = original->objective_function;
}

// Heurística construtiva (sem semente aleatória para inst200)
void heu_cons_ale_gul(Solution *sol, int use_random_seed) {
    int available_hubs[MAX_NODES];
    for (int i = 0; i < n; i++) available_hubs[i] = i;

    // Seleção de hubs
    for (int i = 0; i < p; i++) {
        int pos = (use_random_seed) ? rand() % (n - i) : 0;
        sol->hubs[i] = available_hubs[pos];
        available_hubs[pos] = available_hubs[n - i - 1];
    }

    // Alocação
    for (int i = 0; i < n; i++) {
        double min_dist = INFINITY;
        for (int j = 0; j < p; j++) {
            double dist = distance_matrix[i][sol->hubs[j]];
            if (dist < min_dist) {
                min_dist = dist;
                sol->allocation[i] = sol->hubs[j];
            }
        }
    }
}

double menor_hub_final(int hub, Solution *sol) {
    double menor_hub_final = INFINITY;
    int temp = 0;

    for (int i = 0; i < n; i++) {
        if (distance_matrix[sol->hubs[hub]][i] < menor_hub_final) {
            if (distance_matrix[sol->hubs[hub]][i] != 0) {
                menor_hub_final = distance_matrix[i][sol->hubs[hub]];
                temp = i;
            }
        }
    }

    return menor_hub_final;
}

double menor_hub_hub(int hub, Solution *sol) {
    double menor_hub_hub = INFINITY;
    int temp = 0, hub_temp = 0;

    for (int i = 0; i < sol->hub_count; i++) {
        if (distance_matrix[sol->hubs[i]][sol->hubs[hub]] < menor_hub_hub) {
            if (distance_matrix[sol->hubs[i]][sol->hubs[hub]] != 0) {
                menor_hub_hub = distance_matrix[i][sol->hubs[hub]];
                temp = sol->hubs[i];
                hub_temp = i;
            }
        }
    }

    return menor_hub_hub * 0.75 + menor_hub_final(hub_temp, sol);
}

double menor_ponto_hub(Solution *sol) {
    double menor_ponto_hub = INFINITY;
    int temp1 = 0, temp2 = 0, temp3;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < sol->hub_count; j++) {
            if (distance_matrix[i][sol->hubs[j]] < menor_ponto_hub) {
                if (distance_matrix[i][sol->hubs[j]] != 0) {
                    menor_ponto_hub = distance_matrix[i][sol->hubs[j]];
                    temp1 = i;
                    temp2 = sol->hubs[j];
                    temp3 = j;
                }
            }
        }
    }

    double t = menor_ponto_hub + menor_hub_hub(temp3, sol);
    return t;
}

// Cálculo da FO com temporização
double compute_objective(Solution *sol) {
    double max_cost = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int k = sol->allocation[i];
            int l = sol->allocation[j];
            double cost = beta * distance_matrix[i][k] + 
                         alpha * distance_matrix[k][l] + 
                         lambda * distance_matrix[l][j];
            if (cost > max_cost) max_cost = cost;
        }
    }
    sol->objective_function = max_cost;
    return max_cost;
}

// Salvar solução em arquivo com detalhes
void save_solution_details(const char *filename, Solution *sol) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Erro ao salvar a solução");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "n: %d\tp: %d\n", n, p);
    fprintf(file, "FO: %.2lf\n", sol->objective_function);
    fprintf(file, "HUBS: [");
    for (int i = 0; i < p; i++) {
        fprintf(file, "%d%s", sol->hubs[i], (i < p - 1) ? ", " : "");
    }
    fprintf(file, "]\n");
    fprintf(file, "OR\tH1\tH2\tDS\tCUSTO\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int h1 = sol->allocation[i];
            int h2 = sol->allocation[j];
            double tij = beta * distance_matrix[i][h1] + 
                         alpha * distance_matrix[h1][h2] + 
                         lambda * distance_matrix[h2][j];
            fprintf(file, "%d\t%d\t%d\t%d\t%.2lf\n", i, h1, h2, j, tij);
        }
    }
    fclose(file);
}

// Exibir solução na tela
void display_solution(Solution *sol) {
    printf("\n--- Solucao ---\n");
    printf("n: %d\tp: %d\n", n, p);
    printf("FO: %.2lf\n", sol->objective_function);
    printf("HUBS: [");
    for (int i = 0; i < p; i++) {
        printf("%d%s", sol->hubs[i], (i < p - 1) ? ", " : "");
    }
    printf("]\n");
    printf("OR\tH1\tH2\tDS\tCUSTO\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 5; j++) {
            int h1 = sol->allocation[i];
            int h2 = sol->allocation[j];
            double tij = beta * distance_matrix[i][h1] + 
                         alpha * distance_matrix[h1][h2] + 
                         lambda * distance_matrix[h2][j];
            printf("%d\t%d\t%d\t%d\t%.2lf\n", i, h1, h2, j, tij);
        }
    }
}

// Rotina de testes e métricas
void run_benchmark(const char *filename, int hub_count, int iterations) {
    // Configuração inicial
    read_instance(filename, hub_count);
    
    // Modo de execução única
    Solution initial_sol;
    initialize_solution(&initial_sol);
    clock_t start = clock();
    heu_cons_ale_gul(&initial_sol, 1);
    initial_sol.objective_function = compute_objective(&initial_sol);
    double time_single = (double)(clock() - start) / CLOCKS_PER_SEC;

    // Salvar e exibir solução inicial
    save_solution_details("solucao_inicial.txt", &initial_sol);
    display_solution(&initial_sol);

    // Modo iterativo
    double total_time_heuristic = 0, total_time_fo = 0;
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
        for(int i = 0; i < 1000 ; i++) {
            heu_cons_ale_gul(&temp_sol, 1);
        }
        compute_objective(&temp_sol);
        temp = menor_ponto_hub(&initial_sol);
        if (maior < temp) {
            maior = temp;
        }
    }
    total_time_fo = (double)(clock() - start_fo) / CLOCKS_PER_SEC;

    compute_objective(&initial_sol);
    printf("\nDescricao do Computador\t?\t\t(seg.)\t\tFO\t\tTempo (seg.)\t\tTempo Sol. Inicial (seg.)\tTempo Calc. FO (seg.)\n");
    printf("R5 4.4GHz - 16GB ram - Windows 11\t97\t\t%.2lf\t%.5lf\t\t\t%.5lf\t\t\t\t%.5lf\n\n", 
           initial_sol.objective_function,
           time_single,
           total_time_heuristic,
           total_time_fo);
    
    printf("Maior custo encontrado (Funcao Objetivo): %.2lf\n", maior);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Uso: %s <arquivo_instancia> <num_hubs>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Definir a semente do gerador de números aleatórios
    srand(time(NULL));

    run_benchmark(argv[1], atoi(argv[2]), 1000);
    return EXIT_SUCCESS;
}
