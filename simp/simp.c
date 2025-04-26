#include "./simp.h"
#include "../parser/parser.h"
#include "../mymat/mymat.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

struct sp {
  int num_variables, num_constraints; // n: numero variabili, m: numero vincoli
  matrix cost_coeffs; // vattore dei coefficienti di costo (n, 1)
  matrix tech_coeffs; // matrice dei coefficienti tecnologici (m, n)
  matrix known_terms; // vettori dei termini noti (m, 1)
  matrix cost_coeffs_base; // vettore dei coefficienti di costo in base (m, 1)

  int *base_indices, *out_base_indices; // indici delle variabili in base (m), indici delle variabili fuori base (n - m)
  matrix base; // matrice di base (m, m)

  matrix x_base_value; // vettore delle variabili in base
  float z_value; // valore della funzione obiettivo
};

simp catch_problem(FILE *f) {
  if (!f) {
    perror("cannot use a null file pointer");
    exit(EXIT_FAILURE);
  }

  simp s;

  // alloca e controlla
  if (!(s = calloc(1, sizeof(struct sp)))) {
    perror("bad allocation");
    exit(EXIT_FAILURE);
  }

  // verifica funzione obiettivo
  verify_goal_function(f);

  // cattura info dal file utilizzando il parser
  s->num_variables = parse_num_variables(f);
  s->tech_coeffs = parse_tech_coeffs(f);
  s->num_constraints = get_rows_tech_coeffs(s);
  s->cost_coeffs = parse_cost_coeffs(f);
  s->known_terms = parse_known_terms(f);

  // contralla che il vettore dei termini noti non abbia componenti negative
  if (has_negative_component(get_known_terms(s))) {
    fprintf(stderr, "error: known terms vector has a negative component");
    print_matrix(get_known_terms(s));
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  // parsing base
  s->base_indices = parse_base_indices(f);

  // verifica congruenza indici di base
  verify_base_indices(s);

  // costruisci matrice di base
  s->base = build_base_by_indices(s);

  // verifica invertibilità
  if (determinant(get_base(s)) == .0) {
    fprintf(stderr, "error: base not invertible");
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  // costruisci indici fuori base
  s->out_base_indices = build_out_base_indices(s);

  // costruisci coefficienti di costo in base
  s->cost_coeffs_base = build_cost_coeffs_base(s);

  // calcola valore funzione obiettivo
  s->z_value = build_z_value_by_base(s);

  // calcola soluzione ottima di base
  s->x_base_value = build_x_base_value(s);

  return s;
}

int optimality_test(simp s, int *entering_index) {
  // num_variabili - num_vincoli variabili fuori base
  int *out_base_indices = get_out_base_indices(s),
    num_out_base = get_num_variables(s) - get_num_constraints(s),
    best_index = -1, i, var_index;

  float max_reduced_cost = -1e9, reduced_cost;

  // per ogni variabile fuori base
  for (i = 0; i < num_out_base; i++) {
      var_index = out_base_indices[i];

      // calcola il coefficiente di costo ridotto z_j - c_j
      reduced_cost = get_reduced_cost_coeff_by_index(s, var_index);

      // preserva il coefficiente di costo ridotto massimo
      if (reduced_cost > max_reduced_cost) {
          max_reduced_cost = reduced_cost;
          best_index = var_index;
      }
  }

  // se sono tutti non positivi, test passato
  if (max_reduced_cost <= 0)  return 1;

  // altrimenti entra in base la variabile con coefficiente
  // di costo ridotto massimo
  if (entering_index != NULL)
    *entering_index = best_index;

  return 0;
}

float unboundedness_test(simp s, const int index, int *exiting_index) {
  // verifics indici
  if (index <= 0 || index > get_num_variables(s)) {
      fprintf(stderr, "error: index %d out of bound in problem composed by %d variables",
              index, get_num_variables(s));
      clean_all(s);
      exit(EXIT_FAILURE);
  }

  // verifica index in base
  if (is_in_base(s, index)) {
      fprintf(stderr, "error: cannot calculate reduced cost coeff for index %d in base", index);
      print_base(s);
      clean_all(s);
      exit(EXIT_FAILURE);
  }

  // b_overline = inversa_matrice_base x vettore_termini_noti
  // y_index = inversa_matrice_base x colonna numero index matrice coefficienti tecnologici
  matrix base_inv = inverse(get_base(s)),
    b_overline = get_x_base_value(s),
    y_index = row_by_column_multiplication(base_inv, get_column(get_tech_coeffs(s), index));

  int *base_indices = get_base_indices(s), min_index = -1, i;
  const int num_constraints = get_num_constraints(s);

  // min_ratio = min {b_overline[i] / y_index[i]: y_index[i] > 0}
  float min_ratio = 1e9, y_value, ratio;

  for (i = 1; i <= num_constraints; i++) {
    y_value = get_matrix_element(y_index, i, 1);
    printf("\ny%d%d: %.2f\n", i, index, y_value);
    if (y_value > 0) {
        ratio = get_matrix_element(b_overline, i, 1) / y_value;
        if (ratio < min_ratio) {
            min_ratio = ratio;
            min_index = base_indices[i - 1];
        }
    }
  }

  free(base_inv);

  // ottimo illimitato, tutti i rapporti minimi sono <= 0
  if (min_index == -1) {
      fprintf(stderr, "unbounded optimum\n");
      print_matrix(y_index);
      clean_all(s);
      exit(EXIT_FAILURE);
  }

  free(y_index);

  // esce dalla base la variabile a cui è associata il rapporto minimo
  if (exiting_index) *exiting_index = min_index;

  // la variabile x_index entra con valore min_ratio
  return min_ratio;
}

// coefficiente di costo ridotto z_index - c_index
float get_reduced_cost_coeff_by_index(simp s, const int index) {
  if (index <= 0 || index > get_num_variables(s)) {
    fprintf(stderr, "error: index %d out of bound in problem composed by %d variables",
      index, get_num_variables(s));
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  if (is_in_base(s, index)) {
    fprintf(stderr, "error: cannot calculate reduced cost coeff for index %d in base", index);
    print_base(s);
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  // cbt: vettore trasposto dei coefficienti di costo di base
  // base_inv: matrice inversa di base
  // tech_coeffs_column: colonna numero index della matrice dei coefficienti tecnologici
  // c_index: coefficienti di costo della variabile x_index
  // z_index: cbt x base_inv x tech_coeffs_column
  float z_index, c_index = get_matrix_element(get_cost_coeffs(s), index, 1), res;
  matrix cbt = transpose(get_cost_coeffs_base(s)),
    base_inv = inverse(get_base(s)),
    tech_coeffs_column = get_column(get_tech_coeffs(s), index);

  matrix tmp = row_by_column_multiplication(cbt, base_inv);
  matrix res_matrix = row_by_column_multiplication(tmp, tech_coeffs_column);

  z_index = get_matrix_element(res_matrix, 1, 1);

  res = z_index - c_index;
  printf("\nz%d - c%d: %.2f\n", index, index, res);

  free(tech_coeffs_column);
  free(cbt);
  free(base_inv);
  free(tmp);
  free(res_matrix);

  return res;
}

// algoritmo del simplesso
void solve_simplex(simp s) {
  int entering_index, exiting_index, iteration = 0;
  float min_ratio;

  while (1) {
    printf("\n\n---------------------------------------------------------------\n\n");
    printf("\niteration num. %d\n", ++iteration);

    // stampa info generali
    print_general_info(s);

    // test di ottimalità
    if (optimality_test(s, &entering_index)) {
        printf("current base solution is optimal\n");
        break;
    }

    // test di illimitatezza
    min_ratio = unboundedness_test(s, entering_index, &exiting_index);

    printf("current base solution isn't optimal, x%d exists and x%d joins the base with value %.2f\n",
      exiting_index, entering_index, min_ratio);

    printf("updating base ...\n");

    // aggiorna base
    update_base(&s, entering_index, exiting_index);
  }
}

// getters
const int get_num_variables(simp s) { return s->num_variables; }
const int get_num_constraints(simp s) { return s->num_constraints; }
matrix get_known_terms(simp s) { return s->known_terms; }
matrix get_tech_coeffs(simp s) { return s->tech_coeffs; }
matrix get_cost_coeffs(simp s) { return s->cost_coeffs; }
const int get_rows_tech_coeffs(simp s) { return get_rows(s->tech_coeffs); }
const int get_columns_tech_coeffs(simp s) { return get_columns(s->tech_coeffs); }
int* get_base_indices(simp s) { return s->base_indices; }
int* get_out_base_indices(simp s) { return s->out_base_indices; }
matrix get_base(simp s) { return s->base; }
matrix get_cost_coeffs_base(simp s) { return s->cost_coeffs_base; }
matrix get_x_base_value(simp s) { return s->x_base_value; }
float get_z_value(simp s) { return s->z_value; }

const int get_base_index(simp s, const int i) {
  if (i <= 0 || i > get_num_constraints(s)) {
    fprintf(stderr, "error getting base index %d of a problem with %d constraints", i, get_num_constraints(s));
    clean_all(s);
    exit(EXIT_FAILURE);
  }
  return s->base_indices[i - 1];
}

void set_base_index(simp s, const int i, const int value) {
  if (i <= 0 || i > get_num_constraints(s)) {
    fprintf(stderr, "error getting base index %d of a problem with %d constraints", i, get_num_constraints(s));
    clean_all(s);
    exit(EXIT_FAILURE);
  }
  s->base_indices[i - 1] = value;
}

const int get_out_base_index(simp s, const int i) {
  int num_variables = get_num_variables(s), num_constraints = get_num_constraints(s);
  if (i <= 0 || i > num_variables + num_constraints) {
    fprintf(stderr, "error getting out base index %d of a problem with %d variables and %d constraints",
      i, num_variables, num_constraints);
    clean_all(s);
    exit(EXIT_FAILURE);
  }
  return s->out_base_indices[i - 1];
}

void print_general_info(simp s) {
  print_num_variables(s);
  print_num_constraints(s);
  print_cost_coeffs(s);
  print_tech_coeffs(s);
  print_known_terms(s);
  print_base_indices(s);
  print_out_base_indices(s);
  print_base(s);
  print_cost_coeffs_base(s);
  print_solution_info(s);
}

void print_num_variables(simp s) {
  printf("\nnum variables: %d (from x1 to x%d)\n", get_num_variables(s), get_num_variables(s));
}

void print_num_constraints(simp s) {
  printf("\nnum constraints: %d\n", get_num_constraints(s));
}

void print_matrix_by_info(matrix m, char *s) {
  printf("\n%s (%d x %d):", s, get_rows(m), get_columns(m));
  print_matrix(m);
}

void print_known_terms(simp s) { print_matrix_by_info(get_known_terms(s), "known terms"); }
void print_cost_coeffs(simp s) { print_matrix_by_info(get_cost_coeffs(s), "cost coeffs"); }
void print_cost_coeffs_base(simp s) { print_matrix_by_info(get_cost_coeffs_base(s), "cost coeffs in base"); }
void print_tech_coeffs(simp s) { print_matrix_by_info(get_tech_coeffs(s), "tech coeffs"); }
void print_base(simp s) { print_matrix_by_info(get_base(s), "base matrix"); }
void print_x_base_value(simp s) { print_matrix_by_info(get_x_base_value(s), "x base value"); }
void print_z_value(simp s) { printf("\nz min value: %.1f\n", get_z_value(s)); }
void print_solution_info(simp s) { print_x_base_value(s); print_z_value(s); }

void print_base_indices(simp s) {
  matrix m = alloc_matrix(1, get_num_constraints(s));

  int j;
  for (j = 1; j <= get_num_constraints(s); j++)
    set_matrix_element(m, 1, j, get_base_index(s, j));

  print_matrix_by_info(m, "base indices");
  destroy_matrix(m);
}

void print_out_base_indices(simp s) {
  int num_columns = get_num_variables(s) - get_num_constraints(s);
  matrix m = alloc_matrix(1, num_columns);

  int j;
  for (j = 1; j <= num_columns; j++)
    set_matrix_element(m, 1, j, get_out_base_index(s, j));

  print_matrix_by_info(m, "out base indices");
  destroy_matrix(m);
}

void verify_base_indices(simp s) {
  int num_variables = get_num_variables(s),
    num_constraints = get_num_constraints(s), i, value;

  for (i = 1; i <= num_constraints; i++) {
    value = get_base_index(s, i);

    if (value <= 0 || value > num_variables) {
      fprintf(stderr, "error: base index in pos %d has value %d in a problem composed by variables from x1 to x%d",
        i, value, num_variables);
      clean_all(s);
      exit(EXIT_FAILURE);
    }
  }
}

matrix build_base_by_indices(simp s) {
  return submatrix_by_columns(get_tech_coeffs(s), get_base_indices(s), get_num_constraints(s));
}

int is_in_base(simp s, const int index) {
  int *base_indices = get_base_indices(s),
    num_constraints = get_num_constraints(s),
    num_variables = get_num_variables(s), i;

  if (index <= 0 || index > num_variables) {
    fprintf(stderr, "error evaluating index %d is in base for a problem composed by %d variables",
      index, num_variables);
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < num_constraints; i++)
    if (base_indices[i] == index) return 1;

  return 0;
}

// aggiorna la base facendo entrare x_entering_index ed uscire x_exiting_index
void update_base(simp *s, const int entering_index, const int exiting_index) {
  const int num_constraints = get_num_constraints(*s), num_variables = get_num_variables(*s);

  // verifica indice
  if (entering_index <= 0 || entering_index > num_variables) {
    fprintf(stderr, "error: cannot enter in base index %d out of bound in a problem composed by %d variables\n",
            entering_index, num_variables);
    clean_all(*s);
    exit(EXIT_FAILURE);
  }

  // verifica indice
  if (exiting_index <= 0 || exiting_index > num_variables) {
    fprintf(stderr, "error: cannot remove index %d from base, out of bound in a problem composed by %d variables\n",
            exiting_index, num_variables);
    clean_all(*s);
    exit(EXIT_FAILURE);
  }

  // verifica indice entrante fuori base
  if (is_in_base(*s, entering_index)) {
    fprintf(stderr, "error: index %d is already in base\n", entering_index);
    print_base(*s);
    clean_all(*s);
    exit(EXIT_FAILURE);
  }

  // verifica indice uscente in base
  if (!is_in_base(*s, exiting_index)) {
    fprintf(stderr, "error: index %d is not in base\n", exiting_index);
    print_base(*s);
    clean_all(*s);
    exit(EXIT_FAILURE);
  }

  int i;

  // modifica vettore indici di base
  for (i = 1; i <= num_constraints; i++)
    if (get_base_index(*s, i) == exiting_index) {
      set_base_index(*s, i, entering_index);
      break;
    }

  // pulisci
  clean_base_related_except_base_indices(*s);

  // ricostruisci tutto
  (*s)->base = build_base_by_indices(*s);
  (*s)->out_base_indices = build_out_base_indices(*s);
  (*s)->cost_coeffs_base = build_cost_coeffs_base(*s);
  (*s)->x_base_value = build_x_base_value(*s);
  (*s)->z_value = build_z_value_by_base(*s);
}


int* build_out_base_indices(simp s) {
  int *base_indices = get_base_indices(s),
    num_constraints = get_num_constraints(s),
    num_variables = get_num_variables(s), *out_base;

  // alloca e controlla
  if (!(out_base = calloc(num_variables - num_constraints, sizeof(int)))) {
    perror("bad allocation\n");
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  int i, k = 0;

  // per ogni variabile,
  // se non è in base si aggiunge nel vettore degli indici fuori base
  for (i = 1; i <= num_variables; i++)
    if (!is_in_base(s, i))  out_base[k++] = i;

  return out_base;
}

matrix build_cost_coeffs_base(simp s) {
  int *base_indices = get_base_indices(s),
    num_constraints = get_num_constraints(s), i, base_var;
  matrix cost_coeffs = get_cost_coeffs(s);
  matrix res = alloc_matrix(num_constraints, 1);
  float coeff;

  // si aggiunge a res solo i coefficienti di costo delle variabili
  // i cui indici sono presenti nel vettore degli indici di base
  for (i = 1; i <= num_constraints; i++) {
    base_var = get_base_index(s, i);
    if (base_var <= 0 || base_var > get_rows(cost_coeffs)) {
      fprintf(stderr, "error: invalid base variable index x%d\n", base_var);
      clean_all(s);
      exit(EXIT_FAILURE);
    }
    coeff = get_matrix_element(cost_coeffs, base_var, 1);
    set_matrix_element(res, i, 1, coeff);
  }

  return res;
}

float build_z_value_by_base(simp s) {
  matrix cbt = transpose(get_cost_coeffs_base(s)),
    base_inv = inverse(get_base(s)), b = get_known_terms(s);

  matrix tmp = row_by_column_multiplication(cbt, base_inv);
  matrix res_matrix = row_by_column_multiplication(tmp, b);

  float res = get_matrix_element(res_matrix, 1, 1);

  free(cbt); free(base_inv); free(tmp); free(res_matrix);

  return res;
}

// vettore delle variabili in base
matrix build_x_base_value(simp s) {
  matrix base_inv = inverse(get_base(s));
  matrix res = row_by_column_multiplication(base_inv, get_known_terms(s));

  free(base_inv);

  // controlla ammissibilità
  if (has_negative_component(res)) {
    fprintf(stderr, "base solution ineligible\n");
    print_base_indices(s);
    print_matrix_by_info(res, "x base value");
    free(res);
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  return res;
}

// deallocazione generica
void clean_ptr(void **ptr) {
  if (*ptr) { free(*ptr); *ptr = NULL; }
}

void clean_base_related_except_base_indices(simp s) {
  if (!s) return;

  clean_ptr((void**) &s->base);
  clean_ptr((void**) &s->out_base_indices);
  clean_ptr((void**) &s->cost_coeffs_base);
  clean_ptr((void**) &s->x_base_value);
}

void clean_base_related(simp s) {
  if (!s) return;

  clean_base_related_except_base_indices(s);
  clean_ptr((void**) &s->base_indices);
}

void clean_all(simp s) {
  if (!s) return;

  clean_ptr((void**) &s->cost_coeffs);
  clean_ptr((void**) &s->tech_coeffs);
  clean_ptr((void**) &s->known_terms);
  clean_base_related(s);

  free(s);
}
