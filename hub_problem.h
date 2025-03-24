#ifndef HUB_PROBLEM_H
#define HUB_PROBLEM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>
#include <limits.h>

// Constantes
#define MAX_NODES 200
#define MAX_HUBS 50

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

// Estrutura para candidatos (usada na busca local)
typedef struct candidate{
    int index;
    double dist;
} Candidate;

// Protótipos das funções de busca local
void local_search(Solution *sol);
int swap_search(Solution *sol);
int insertion_search(Solution *sol);
int replacement_search(Solution *sol);

// Funções auxiliares
int compare_candidates(const void *a, const void *b);
void swap_elements(int *a, int *b);

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
int vet_bin[MAX_NODES];
clock_t start;
double melhor_tempo;

// Protótipos das funções
double calculate_distance(Node a, Node b);
void calculate_distance_matrix();
void read_instance(const char *filename);
void initialize_solution(Solution *sol);
void clone_solution(Solution *original, Solution *clone);
int read_solution(const char *filename, Solution *sol);
void heu_cons_ale_gul(Solution *sol, int use_random_seed);
void calculo_fo(Solution &sol);
void save_solution_details(Solution &sol);
void display_solution(const char *filename, Solution *sol);
void* run_benchmark(void* arg);
void* time_limit(void* arg);

#endif // HUB_PROBLEM_H
