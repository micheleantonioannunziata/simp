#include "../mymat/mymat.h"
#include <stdio.h>

typedef struct sp *simp;

// cattura le info dal file
simp catch_problem(FILE*);

// esegue il test di ottimalità, restuitisce 1 se viene passato
// altrimenti restituisce 0 e assegna a entering_index 
// l'indice della variabile che deve entrare in base
int optimality_test(simp s, int *entering_index);

// esegue il test di illimitatezza rispetto alla variabile x_index che deve entrare in base;
// se tutti i rappoti minimi sono non positivi allora c'è un ottimo illimitato
// altrimenti imposta exiting_index alla variabile che deve uscire 
// (ossia quella a cui è associato il rapporto minimo) 
// e restituisce il valore di x_index
float unboundedness_test(simp s, const int index, int *exiting_index);

// calcola il coefficiente di costo ridotto z_index - c_index
// z_index = vettore trasposto dei coefficienti di costo di base 
//              x matrice inversa di base
//              x colonna numero index della matrice dei coefficienti tecnologici
// c_index = coefficienti di costo associato alla variabile x_index 
float get_reduced_cost_coeff_by_index(simp s, const int index);

// aggiorna la base facendo entrare x_entering_index ed uscie x_exiting_index
// aggiorna la matrice di base, il vettore degli indici in base,
// il vettore degli indici fuori base, il vettore dei coefficienti di costo in base,
// il vettore delle variabili in base ed il valore ottimo
void update_base(simp *s, const int entering_index, const int exiting_index);

// algoritmo del simplesso
void solve_simplex(simp s);

// metodo delle due fasi
void solve_problem(simp);

// prima fase del metodo delle due fasi
// costruisce e risolve il problema artificiale
// restituisce il vettore degli indici delle variabili in base
// della soluzione ottima di base
int* first_phase(simp);

// restituisce 1 se la matrice identità è presente nella matrice
// dei coefficienti tecnologici, altrimenti 0
int is_identity_in_tech_coeffs(simp);

// restituisce il vettore degli indici delle colonne dell'identità
// presente nella matrice dei coefficienti tecnologici
int* find_identity_columns_in_tech_coeffs(simp);

// restuituisce una copia
simp get_copy_simp(simp);

// restituisce il numero di variabili
const int get_num_variables(simp);

// restituisce il numero di vincoli
const int get_num_constraints(simp);

// restituisce il vettore dei termini noti
matrix get_known_terms(simp);

// restituisce la matrice dei coefficienti tecnoligici
matrix get_tech_coeffs(simp);

// restituisce il vettore dei coefficienti di costo
matrix get_cost_coeffs(simp);

// restituisce il vettore dei coefficienti di costo in base
matrix get_cost_coeffs_base(simp);

// restituisce la matrice di base
matrix get_base(simp);

// restituisce il vettore delle variabili in base
matrix get_x_base_value(simp);

// restituisce il numero di righe della matrice diìei coefficienti tecnologici
const int get_rows_tech_coeffs(simp);

// restituisce il numero di colonne della matrice diìei coefficienti tecnologici
const int get_columns_tech_coeffs(simp);

// restituisce gli indici delle variabili in base
int* get_base_indices(simp);

// restituisce gli indici delle colonne dell'identità presenti 
// nella matrice dei coefficienti tecnologici
int* get_identity_columns_in_tech_coeffs(simp);

// restituisce il numero di colonne dell'identità presenti 
// nella matrice dei coefficienti tecnologici
int get_num_identity_columns_in_tech_coeffs(simp);

// restituisce l'i-esimo indice delle variabili in base
const int get_base_index(simp s, const int i);

// restituisce l'i-esimo indice delle variabili fuori base
const int get_out_base_index(simp s, const int i);

// imposta l'i-esimo indice delle variabili in base
void set_base_index(simp s, const int i, const int value);

// restituisce gli indici delle variabili fuori base
int* get_out_base_indices(simp);

// restituisce il valore ottimo
float get_z_value(simp);

// imposta l'i-esimo indice delle variabili fuori base
void set_out_base_index(simp s, const int i, const int value);

// verifica se x_index è in base
int is_in_base(simp s, const int index);

// printer
void print_matrix_by_info(matrix, char*);
void print_general_info(simp);
void print_num_variables(simp);
void print_num_constraints(simp);
void print_known_terms(simp);
void print_cost_coeffs(simp);
void print_cost_coeffs_base(simp);
void print_tech_coeffs(simp);
void print_base_indices(simp);
void print_out_base_indices(simp);
void print_base(simp);
void print_x_base_value(simp);
void print_z_value(simp);
void print_solution_info(simp);
void print_identity_columns_in_tech_coeffs(simp);

void verify_base_indices(simp);

// costruisce la matrice di base a partire dagli indici di base
matrix build_base_by_indices(simp);

// costruisce il vettore degli indici fuori base a partire dagli indici di base
int* build_out_base_indices(simp);

// restituisce il valore ottimo
// vettore trasposto dei coefficienti di costo in base 
//      x vettore delle variabili in base
float build_z_value_by_base(simp);

// restituisce il vettore delle variabili in base
// matrice inversa di base x vettore dei termini noti
matrix build_x_base_value(simp);

// restituisce il vettore dei coefficienti di costo in base
matrix build_cost_coeffs_base(simp);

// dealloca tutto ciò che concerne la base 
// al netto del vettore degli indici di base
void clean_base_related_except_base_indices(simp);

// dealloca tutto ciò che concerne la base
void clean_base_related(simp);

// deallocazione generica
void clean_ptr(void**);

// deallocazione totale
void clean_all(simp);