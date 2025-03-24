#ifndef HUB_PROBLEM_H
#define HUB_PROBLEM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h> // Para DBL_MAX
#include <algorithm>
#include <pthread.h>
#include <unistd.h>
#include <limits.h>

// Constantes
#define MAX_NODES 200
#define MAX_HUBS 50
#define INFINITY DBL_MAX

// Estruturas
typedef struct node {
    double x, y;
} Node;

typedef struct solution {
    int hubs[MAX_HUBS];
    int vet_sol[MAX_NODES];
    double fo;
    int allocation[MAX_NODES];
} Solution;

// Variáveis globais
extern int time_limit_sec;
extern Node nodes[MAX_NODES];
extern int num_nos;
extern int num_hubs;
extern double beta, alfa, lambda;
extern double mat_distancias[MAX_NODES][MAX_NODES];
extern double mat_custo[MAX_NODES][MAX_NODES];
extern int melhor_hub[MAX_HUBS];
extern double melhor_fo;
aint vet_bin[MAX_NODES];

// Protótipos das funções
double calculate_distance(Node a, Node b);
void calculate_distance_matrix();
void read_instance(const char *filename);
void initialize_solution(Solution *sol);
void clone_solution(Solution *original, Solution *clone);
int read_solution(const char *filename, Solution *sol);
void heu_cons_ale_gul(Solution *sol, int use_random_seed);
void calculo_fo(Solution &sol, int iterations);
void save_solution_details(Solution &sol);
void display_solution(Solution *sol);
void* run_benchmark(void* arg);
void* time_limit(void* arg);

#endif // HUB_PROBLEM_H
