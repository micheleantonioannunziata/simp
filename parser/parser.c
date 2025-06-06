#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../mymat/mymat.h"
#include "parser.h"

#define MAX 256

// restituisce il numero di variabili
int parse_num_variables(FILE *f) {
  rewind(f); // riavvolgi
  char line[MAX], *token;
  int max_index = 0, var_index;

  // per ogni riga letta
  while (fgets(line, sizeof(line), f)) {
    // leggi spazio per spazio
    token = strtok(line, " ");
    while (token != NULL) {
      // leggi indice
      if (sscanf(token, "%*[^x]x%d", &var_index) == 1) {
        // aggiorna indice massimo
        if (var_index > max_index)
            max_index = var_index;
      }
      token = strtok(NULL, " ");
    }
  }

  return max_index;
}

// restituisce il vettore colonna dei coefficienti di costo
matrix parse_cost_coeffs(FILE* f) {
  rewind(f); // riavvolgi
  char line[MAX];

  // leggi la prima riga
  fgets(line, sizeof(line), f);

  int num_variables = parse_num_variables(f), index;
  float coeff;

  matrix res = alloc_matrix(num_variables, 1);

  // scandisci la funzione obiettivo spazio per spazio
  char *token = strtok(line, " ");
  while (token) {
    // non considera min
    if (strcmp(token, "min") != 0) {

      // lettura coefficiente razionale
      if (sscanf(token, "%fx%d", &coeff, &index) == 2)
        set_matrix_element(res, index, 1, coeff);

      // lettura -1x
      else if (sscanf(token, "-x%d", &index) == 1)
        set_matrix_element(res, index, 1, -1.0);

      // lettura +1x
      else if (sscanf(token, "+x%d", &index) == 1 || sscanf(token, "x%d", &index) == 1)
        set_matrix_element(res, index, 1, 1.0);

      else {
        fprintf(stderr, "error: not valid token '%s' in goal function\n", token);
        exit(EXIT_FAILURE);
      }
    }
    token = strtok(NULL, " ");
  }

  return res;
}

// verifica congruenza funziona obiettivo
void verify_goal_function(FILE *f) {
  rewind(f);
  char line[MAX];
  fgets(line, sizeof(line), f);

  // controlla che la stringa inizi con "min"
  if (strncmp(line, "min", 3)) {
    fprintf(stderr, "error: cost function must start with 'min' to be in standard form\n");
    exit(EXIT_FAILURE);
  }

  int var_seen[MAX] = {0}, index;
  float coeff;
  char *token = strtok(line, " ");
  token = strtok(NULL, " ");

  while (token) {
    // verifica coefficienti
    if (sscanf(token, "%fx%d", &coeff, &index) == 2);
    else if (sscanf(token, "+x%d", &index) == 1)  coeff = 1.0;
    else if (sscanf(token, "-x%d", &index) == 1)  coeff = -1.0;
    else if (sscanf(token, "x%d", &index) == 1)   coeff = 1.0;
    else {
        fprintf(stderr, "error: invalid term in goal function: '%s'\n", token);
        exit(EXIT_FAILURE);
    }

    // verifica indici
    if (index <= 0 || index >= MAX) {
        fprintf(stderr, "error: variable x%d out of bound (1 - %d)\n", index, MAX);
        exit(EXIT_FAILURE);
    }

    // gestione duplicati
    if (var_seen[index]) {
        fprintf(stderr, "error: duplicate variable x%d in goal function\n", index);
        exit(EXIT_FAILURE);
    }

    var_seen[index] = 1;
    token = strtok(NULL, " ");
  }
}

// restituisce la matrice dei coefficienti tecnologici
matrix parse_tech_coeffs(FILE *f) {
  rewind(f);
  char line[MAX];
  fgets(line, sizeof(line), f); // salta prima riga (funzione obiettivo)
  long pos;
  int num_constraints = parse_num_contraints(f, &pos);

  int num_vars = parse_num_variables(f);
  matrix res = alloc_matrix(num_constraints, num_vars);
  fseek(f, pos, SEEK_SET); // ritorna alla posizione iniziale (dopo la funzione obiettivo)

  int row = 1;
  float coeff;
  int index;
  char *token;

  // per ogni riga tranne quella della base
  while (fgets(line, sizeof(line), f)) {
    token = strtok(line, " ");

    // scandisci spazio per spazio
    while (token != NULL) {
      // lettura coefficiente razionale
      if (sscanf(token, "%fx%d", &coeff, &index) == 2)  set_matrix_element(res, row, index, coeff);
      // lettura coefficiente +1
      else if (sscanf(token, "+x%d", &index) == 1 || sscanf(token, "x%d", &index) == 1)
        set_matrix_element(res, row, index, 1);
      // lettura coefficiente -1
      else if (sscanf(token, "-x%d", &index) == 1)
        set_matrix_element(res, row, index, -1.0);
      // interrompi lettura riga ad '='
      else if (strcmp(token, "=") == 0)
        break;
      else {
        fprintf(stderr, "error: unknown token '%s' in constraint %d\n", token, row);
        exit(EXIT_FAILURE);
      }

      token = strtok(NULL, " ");
    }

    row++;
  }

  return res;
}

const int parse_num_contraints(FILE *f, long *pos) {
  rewind(f);
  char line[MAX];
  fgets(line, sizeof(line), f); // salta funzione obiettivo

  int num_constraints = 0;
  *pos = ftell(f); // posizione attuale

  // per ogni riga
  while (fgets(line, sizeof(line), f))  num_constraints++;

  return num_constraints;
}

// termini noti
matrix parse_known_terms(FILE *f) {
  rewind(f);
  char line[MAX];

  // 1 per ogni vincolo
  long pos;
  int num_constraints = parse_num_contraints(f, &pos);

  matrix res = alloc_matrix(num_constraints, 1);
  fseek(f, pos, SEEK_SET);

  float value;
  int row = 1;
  while (fgets(line, sizeof(line), f)) {
      if (sscanf(strrchr(line, '=') + 1, "%f", &value) == 1)
          set_matrix_element(res, row, 1, value);
      row++;
  }

  return res;
}
