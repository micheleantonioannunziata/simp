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

  // vettore degli indici delle colonne dell'identità presenti nella matrice dei coefficienti tecnologici
  int *identity_columns_in_tech_coeffs, num_identity_columns_in_tech_coeffs;

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
    fclose(f);
    exit(EXIT_FAILURE);
  }

  // verifica funzione obiettivo
  verify_goal_function(f);

  // cattura info dal file utilizzando il parser
  s->num_variables = parse_num_variables(f);
  s->tech_coeffs = parse_tech_coeffs(f);

  s->identity_columns_in_tech_coeffs = find_identity_columns_in_tech_coeffs(s);

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

  s->base_indices = NULL;
  s->out_base_indices = NULL;
  s->base = NULL;
  s->cost_coeffs_base = NULL;
  s->x_base_value = NULL;

  return s;
}

int optimality_test(simp s, int *entering_index) {
  // num_variabili_fuori_base = num_variabili - num_vincoli
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
      exit(EXIT_FAILURE);
  }

  // verifica index in base
  if (is_in_base(s, index)) {
      fprintf(stderr, "error: cannot calculate reduced cost coeff for index %d in base", index);
      print_base(s);
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

// risoluzione problema con metodo delle due fasi
void solve_problem(simp s) {
  int i, *identity_columns_in_tech_coeffs = get_identity_columns_in_tech_coeffs(s);

  const int num_identity_columns_in_tech_coeffs = get_num_identity_columns_in_tech_coeffs(s),
    num_constraints = get_num_constraints(s), num_variables = get_num_variables(s);

  matrix original_tech_coeffs = get_copy_matrix(get_tech_coeffs(s));

  // se nella matrice dei coefficienti tecnologici c'è l'identità,
  // parti da quella come base
  if (num_identity_columns_in_tech_coeffs == num_constraints) {
    printf("\nin tech coeffs index there is already identity matrix (%d x %d), so let's start from it\n",
      num_constraints, num_constraints);

    // indici di base
    for (i = 1; i <= num_constraints; i++)
      set_base_index(s, i, identity_columns_in_tech_coeffs[i - 1]);

    print_base_indices(s);
  }
  else {
    printf("\nin tech coeffs index there is NOT identity matrix (%d x %d), so let's do two phases method\n",
      num_constraints, num_constraints);
    // metodo delle due fasi
    int *base_indices = first_phase(s), second_phase = 1;

    // verifica se c'è una variabile artificiale in base
    for (i = 0; i < num_constraints; i++)
      if (base_indices[i] > num_variables) {
        printf("\nin the optimal base solution of the artificial problem there is at least one artifial variable (x%d),\nso the original problem is ineligible\n",
          base_indices[i]);
        return;
      }

    // seconda fase
    if (!(s->base_indices = (int*) calloc(num_constraints, sizeof(int)))) {
      perror("bad allocation\n");
      free(base_indices);
      clean_all(s);
      exit(EXIT_FAILURE);
    }

    for (i = 1; i <= num_constraints; i++)
      set_base_index(s, i, base_indices[i - 1]);

    clean_ptr((void**) &s->tech_coeffs);
    s->tech_coeffs = original_tech_coeffs;

    printf("\nstarting second phase\n");
  }

  // ricostruisci tutto
  clean_base_related_except_base_indices(s);
  s->base = build_base_by_indices(s);
  s->out_base_indices = build_out_base_indices(s);
  s->cost_coeffs_base = build_cost_coeffs_base(s);
  s->x_base_value = build_x_base_value(s);
  s->z_value = build_z_value_by_base(s);

  solve_simplex(s);

  print_x_base_value(s);
  print_z_value(s);
}

int* first_phase(simp s) {
  simp art = get_copy_simp(s);

  const int num_constraints = get_num_constraints(art),
    num_variable_original_problem = get_num_variables(art),
    num_identity_columns_in_tech_coeffs = get_num_identity_columns_in_tech_coeffs(art);

  int i, k, *identity_columns_in_tech_coeffs = get_identity_columns_in_tech_coeffs(art);

  // al problema artificiale si aggiungono variabili per ottenere
  // la matrice identità
  const int num_variable_artificial_problem = num_variable_original_problem +
     (num_constraints - num_identity_columns_in_tech_coeffs);
  art->num_variables = num_variable_artificial_problem;

  // allocazione coefficienti di costo
  clean_ptr((void**) &art->tech_coeffs);
  art->cost_coeffs = alloc_matrix(num_variable_artificial_problem, 1);

  // coefficienti di costo delle variabili del problema originale nulli
  for (i = 1; i <= num_variable_artificial_problem; i++)
    if (i <= num_variable_original_problem)
      set_matrix_element(get_cost_coeffs(art), i, 1, 0);
    // coefficienti di costo delle variabili artificiali pari ad 1
    else  set_matrix_element(get_cost_coeffs(art), i, 1, 1);

  // allocazione matrice dei coefficienti tecnologici
  clean_ptr((void**) &art->tech_coeffs);

  // array per tenere traccia di quali colonne dell'identità sono già presenti
  int* identity_positions;
  if (!(identity_positions = (int*) calloc(num_constraints + 1, sizeof(int)))) {
    fprintf(stderr, "bad allocation\n");
    free(identity_columns_in_tech_coeffs);
    clean_all(s);
    clean_all(art);
    exit(EXIT_FAILURE);
  }

  // marca posizioni delle colonne dell'identità già presenti
  matrix current_column;
  int position;

  for (i = 0; i < num_identity_columns_in_tech_coeffs; i++) {
    current_column = get_column(get_tech_coeffs(s), identity_columns_in_tech_coeffs[i]);
    position = verify_identity_column(current_column);
    identity_positions[position] = 1;  // segna posizione come già presente
    destroy_matrix(current_column);
  }

  // ora identity_positions[i] == 1 indica che la i-esima colonna
  // dell'identità è già presente

  // crea la matrice dei coefficienti tecnologici per il problema artificiale
  // prima copia la matrice originale
  matrix artificial_tech_coeffs = get_copy_matrix(get_tech_coeffs(s));

  // identifica le colonne dell'identità che mancano
  // e crea le relative colonne artificiali
  matrix identity_matrix = identity(num_constraints), id_current_column, temp;
  matrix missing_columns = NULL;
  int num_missing_columns = 0;

  // per ogni possibile colonna dell'identità
  for (i = 1; i <= num_constraints; i++)
    // se non presente
    if (!identity_positions[i]) {
      id_current_column = get_column(identity_matrix, i);

      if (missing_columns == NULL)  missing_columns = id_current_column;
      else {
        matrix temp = horizontal_concatenation(missing_columns, id_current_column);
        destroy_matrix(id_current_column);
        destroy_matrix(missing_columns);
        missing_columns = temp;
      }

      printf("identity column %d added to tech coeffs matrix\n", i);
      num_missing_columns++;
    }

  // concatena la matrice dei coefficienti tecnologici con le colonne mancanti dell'identità
  if (missing_columns != NULL) {
    art->tech_coeffs = horizontal_concatenation(artificial_tech_coeffs, missing_columns);
    destroy_matrix(artificial_tech_coeffs);
    destroy_matrix(missing_columns);
  } else  art->tech_coeffs = artificial_tech_coeffs;

  destroy_matrix(identity_matrix);
  free(identity_positions);

  // indici di base - indici variabili artificiali
  k = 1;

  if (!(art->base_indices = (int*) calloc(num_constraints, sizeof(int)))) {
    perror("bad allocation\n");
    clean_all(art);
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  // indici delle colonne dell'identità iniziali
  for (i = 1; i <= num_identity_columns_in_tech_coeffs; i++)
    set_base_index(art, k++, identity_columns_in_tech_coeffs[i - 1]);

  // più le variabili artificiali
  for (i = num_variable_original_problem + 1; i <= num_variable_artificial_problem; i++)
    set_base_index(art, k++, i);

  clean_base_related_except_base_indices(art);
  art->base = build_base_by_indices(art);
  art->out_base_indices = build_out_base_indices(art);

  art->cost_coeffs_base = alloc_matrix(num_constraints, 1);
  // coefficienti di costo in base pari a 1
  initizialize_matrix(get_cost_coeffs_base(art), 1);

  art->x_base_value = build_x_base_value(art);
  art->z_value = build_z_value_by_base(art);

  printf("\nstarting first phase\n");
  print_general_info(art);

  solve_simplex(art);

  int *base_optimal_indices_copy;

  if (!(base_optimal_indices_copy =  calloc(num_constraints, sizeof(int)))) {
    perror("bad allocation\n");
    clean_all(art);
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  for (i = 1; i <= num_constraints; i++)
    base_optimal_indices_copy[i - 1] = get_base_index(art, i);

  clean_all(art);

  return base_optimal_indices_copy;
}

// algoritmo del simplesso
void solve_simplex(simp s) {
  int entering_index, exiting_index, iteration = 0;
  float min_ratio;

  while (1) {
    printf("\n\n---------------------------------------------------------------\n\n");
    printf("\niteration num. %d\n", ++iteration);

    print_solution_info(s);

    print_base_indices(s);

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

// restituisce il vettore degli indici delle colonne dell'identità
// presente nella matrice dei coefficienti tecnologici
int* find_identity_columns_in_tech_coeffs(simp s) {
  return find_identity_columns(get_tech_coeffs(s), &s->num_identity_columns_in_tech_coeffs);
}

// getters
simp get_copy_simp(simp s) {
  simp res;

  int i, num_constraints = get_num_constraints(s),
    num_variables = get_num_variables(s);

  if (!(res = calloc(1, sizeof(struct sp)))) {
    perror("bad allocation\n");
    clean_all(s);
    exit(EXIT_FAILURE);
  }

  res->num_variables = num_variables;
  res->num_constraints = num_constraints;

  if (s->cost_coeffs)   res->cost_coeffs = get_copy_matrix(get_cost_coeffs(s));
  if (s->tech_coeffs)   res->tech_coeffs = get_copy_matrix(get_tech_coeffs(s));
  if (s->base)          res->base = get_copy_matrix(get_base(s));
  if (s->x_base_value)  res->x_base_value = get_copy_matrix(get_x_base_value(s));
  if (s->known_terms)   res->known_terms = get_copy_matrix(get_known_terms(s));
  res->z_value = get_z_value(s);

  if (s->base_indices) {
    if (!(res->base_indices = (int*) calloc(num_constraints, sizeof(int)))) {
      perror("bad allocation\n");
      clean_all(s);
      exit(EXIT_FAILURE);
    }
    for (i = 1; i <= num_constraints; i++)
      set_base_index(res, i, get_base_index(s, i));
  }

  if (s->out_base_indices) {
    if (!(res->out_base_indices = (int*) calloc(num_variables - num_constraints, sizeof(int)))) {
      perror("bad allocation\n");
      clean_all(s);
      exit(EXIT_FAILURE);
    }

    for (i = 1; i <= num_variables - num_constraints; i++)
      set_out_base_index(res, i, get_out_base_index(s, i));
  }

  if (s->base)              res->base = get_base(s);
  if (s->cost_coeffs_base)  res->cost_coeffs_base = get_cost_coeffs(res);

  res->num_identity_columns_in_tech_coeffs = get_num_identity_columns_in_tech_coeffs(s);

  if (s->identity_columns_in_tech_coeffs) {
    if (!(res->identity_columns_in_tech_coeffs = (int*) calloc(get_num_identity_columns_in_tech_coeffs(res), sizeof(int)))) {
      perror("bad allocation\n");
      clean_all(s);
      exit(EXIT_FAILURE);
    }

    for (i = 0; i < get_num_identity_columns_in_tech_coeffs(res); i++)
    res->identity_columns_in_tech_coeffs[i] = s->identity_columns_in_tech_coeffs[i];
  }

  return res;
}

const int get_num_variables(simp s) { return s->num_variables; }
const int get_num_constraints(simp s) { return s->num_constraints; }
matrix get_known_terms(simp s) { return s->known_terms; }
matrix get_tech_coeffs(simp s) { return s->tech_coeffs; }
matrix get_cost_coeffs(simp s) { return s->cost_coeffs; }
const int get_rows_tech_coeffs(simp s) { return get_rows(s->tech_coeffs); }
const int get_columns_tech_coeffs(simp s) { return get_columns(s->tech_coeffs); }
int* get_base_indices(simp s) { return s->base_indices; }
int* get_out_base_indices(simp s) { return s->out_base_indices; }
int* get_identity_columns_in_tech_coeffs(simp s) { return s->identity_columns_in_tech_coeffs; }
int get_num_identity_columns_in_tech_coeffs(simp s) { return s->num_identity_columns_in_tech_coeffs; }
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
    fprintf(stderr, "error setting base index %d to value %d of a problem with %d constraints",
      i, value, get_num_constraints(s));
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

void set_out_base_index(simp s, const int i, const int value) {
  int num_variables = get_num_variables(s), num_constraints = get_num_constraints(s);
  if (i <= 0 || i > num_variables - num_constraints) {
    fprintf(stderr, "error setting out base index %d to value %d of a problem with %d variables and %d constraints",
      i, value, num_variables, num_constraints);
    clean_all(s);
    exit(EXIT_FAILURE);
  }
  s->out_base_indices[i - 1] = value;
}

void print_general_info(simp s) {
  printf("\n\n---------------------------------------------------------------\n\n");
  printf("general info of the problem");
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
  printf("\n\n---------------------------------------------------------------\n\n");
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

void print_identity_columns_in_tech_coeffs(simp s) {
  int num_cols = get_num_identity_columns_in_tech_coeffs(s), j;
  matrix m = alloc_matrix(1, num_cols);


  for (j = 1; j <= num_cols; j++)
    set_matrix_element(m, 1, j, s->identity_columns_in_tech_coeffs[j - 1]);

  print_matrix_by_info(m, "identity columns in tech coeffs");
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
          destroy_matrix(res);
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

  destroy_matrix(cbt);
  destroy_matrix(base_inv);
  destroy_matrix(tmp);
  destroy_matrix(res_matrix);

  return res;
}

// vettore delle variabili in base
matrix build_x_base_value(simp s) {
  matrix base_inv = inverse(get_base(s));
  matrix res = row_by_column_multiplication(base_inv, get_known_terms(s));

  destroy_matrix(base_inv);

  if (has_negative_component(res)) {
    fprintf(stderr, "base solution ineligible\n");
    print_base_indices(s);
    clean_all(s);
    destroy_matrix(base_inv);
    destroy_matrix(res);
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
  clean_ptr((void**) &s->identity_columns_in_tech_coeffs);
  clean_base_related(s);

  free(s);
}
