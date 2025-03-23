#ifndef HUB_PROBLEM_H
#define HUB_PROBLEM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h> // Para DBL_MAX
#include <algorithm>

// Constantes
#define MAX_NODES 200
#define MAX_HUBS 50
#define INFINITY DBL_MAX

// Estruturas
typedef struct {
    double x, y;
} Node;

typedef struct {
    int hubs[MAX_HUBS];
    int vet_sol[MAX_NODES];
    double fo;
    int allocation[MAX_NODES];
} Solution;

// Variáveis globais
extern Node nodes[MAX_NODES];
extern int n;
extern double beta, alpha, lambda;
extern double distance_matrix[MAX_NODES][MAX_NODES];
//extern double vetor_hub[MAX_HUBS];
//extern double vetor_nao_hub[MAX_NODES];

// Protótipos das funções
double calculate_distance(Node a, Node b);
void calculate_distance_matrix();
void read_instance(const char *filename);
void initialize_solution(Solution *sol);
void clone_solution(Solution *original, Solution *clone);
int read_solution(const char *filename, Solution *sol);
void heu_cons_ale_gul(Solution *sol, int use_random_seed);
double calculo_fo(Solution *sol, int iterations);
void save_solution_details(const char *filename, Solution *sol);
void display_solution(Solution *sol);
void* run_benchmark(void* arg);
void calcular_fo_g(Solution& s);

#endif // HUB_PROBLEM_H
