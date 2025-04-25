#include "mymat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct mt {
  int rows, columns; // num_righe e num_colonne
  float *data;
};

// allocazione matrice
matrix alloc_matrix(const int rows, const int columns) {
  matrix m;

  // alloca spazio e controlla
  if (!(m = calloc(1, sizeof(struct mt)))) {
    perror("bad allocation\n");
    exit(EXIT_FAILURE);
  }

  // assegnamento righe e colonne
  m->rows = rows;
  m->columns = columns;

  // alloca spazio per i valori e controlla
  if (!(m->data = calloc(get_rows(m) * get_columns(m), sizeof(float)))) {
    perror("bad allocation\n");
    exit(EXIT_FAILURE);
  }

  // inizializza a zero
  initizialize_matrix(m, 0);

  return m;
}

// num_righe
const int get_rows(matrix m) {  return m->rows; }
// num_colonne
const int get_columns(matrix m) {  return m->columns; }

// j_esima colonna
matrix get_column(matrix m, const int column) {
  // verifica indice
  if (column <= 0 || column > get_columns(m)) {
    fprintf(stderr, "error during request column %d in a (%d x %d) matrix\n",
      column, get_rows(m), get_columns(m));
    exit(EXIT_FAILURE);
  }

  int i;
  matrix res = alloc_matrix(get_rows(m), 1); // m righe, 1 colonna

  // cattura elementi della colonna column della matrice m e imposta in res
  for (i = 1; i <= get_rows(m); i++)
      set_matrix_element(res, i, 1, get_matrix_element(m, i, column));

  return res;
}

// i-esima riga
matrix get_row(matrix m, const int row) {
  // verifica indice
  if (row <= 0 || row > get_rows(m)) {
    fprintf(stderr, "error: requested row %d in a (%d x %d) matrix\n",
      row, get_rows(m), get_columns(m));
    exit(EXIT_FAILURE);
  }

  int j;
  matrix res = alloc_matrix(1, get_columns(m)); // 1 riga, n colonne

  // cattura elementi della riga row della matrice m e imposta in res
  for (j = 1; j <= get_columns(m); j++)
      set_matrix_element(res, 1, j, get_matrix_element(m, row, j));

  return res;
}

// restituire il valore dell' elemento m[i][j] utilizzando la formula i * num_colonne + j
float get_matrix_element(matrix m, const int rows_index, const int columns_index) {
  // verifica indici
  if (rows_index <= 0 || columns_index <= 0 || rows_index > get_rows(m) || columns_index > get_columns(m)) {
    fprintf(stderr, "error during getting element in pos [%d][%d] of a (%d x %d) matrix",
        rows_index, columns_index, get_rows(m), get_columns(m));
    exit(EXIT_FAILURE);
  }

  // -1 per astrarre l'indicizzazione a partire da 0
  return m->data[(rows_index - 1) * get_columns(m) + (columns_index - 1)];
}

// impostare il valore dell'elemento m[i][j] utilizzando la formula i * num_colonne + j
void set_matrix_element(matrix m, const int rows_index, const int columns_index, const float value) {
  if (rows_index <= 0 || columns_index <= 0 || rows_index > get_rows(m) || columns_index > get_columns(m)) {
    fprintf(stderr, "error during setting element in pos [%d][%d] to %.1f of a (%d x %d) matrix",
        rows_index, columns_index, value, get_rows(m), get_columns(m));
    exit(EXIT_FAILURE);
  }

  // -1 per astrarre l'indicizzazione a partire da 0
  m->data[(rows_index - 1) * get_columns(m) + (columns_index - 1)] = value;
}

// inizializza tutti i componenti a value
void initizialize_matrix(matrix m, const float value) {
  int i, j;

  // per ogni riga
  for (i = 1; i <= get_rows(m); i++)
    // per ogni colonna
    for (j = 1; j <= get_columns(m); j++)
      // imposta il valore a value
      set_matrix_element(m, i, j, value);
}

// stampa
void print_matrix(matrix m) {
  int i, j;
  float value;

  // stampa indici colonne
  printf("\npos");
  for (i = 1; i <= get_columns(m); i++) printf("\t[%d]", i);
  printf("\n");

  // per ogni riga
  for (i = 1; i <= get_rows(m); i++) {
    // stampa indice riga
    printf("[%d]", i);

    // per ogni colonna
    for (j = 1; j <= get_columns(m); j++) {
      // cattura valore
      value = get_matrix_element(m, i, j);

      // se non ha parte frazionaria stampa solo la parte intera
      if (value == (int) value)  printf("\t %d", (int) value);
      else  printf("\t %.1f", value);
    }

    printf("\n");
  }
}

// leggi da stdin
void scan_matrix(matrix m) {
	int i, j;
  float value;

	for (i = 1; i <= get_rows(m); i++)
		for (j = 1; j <= get_columns(m); j++) {
			printf("insert element in pos. [%d][%d]\n>>> ", i, j);
			scanf("%f", &value);
      set_matrix_element(m, i, j, value);
		}
}

// verifica se m è una matrice quadrata
int is_quadratic(matrix m) {
  return get_columns(m) == get_rows(m);
}

// calcola il determinante
float determinant(matrix m) {

  // se m non quadrata, allora non esiste il determinante
  if (!is_quadratic(m)) {
    fprintf(stderr, "error: matrix must be square to compute determinant\n");
    exit(EXIT_FAILURE);
  }

  int dim = get_rows(m);

  // passo base
  // se m ha solo un componente allora il determinante è il componente stesso
  if (dim == 1)
    return get_matrix_element(m, 1, 1);

  // passo base
  // se m è 2 x 2 si applica il metodo di sarrus
  // (prodotto dei coefficienti in diagonale meno
  // il prodotto dei restanti due coefficienti)
  if (dim == 2)
    return get_matrix_element(m, 1, 1) * get_matrix_element(m, 2, 2)
         - get_matrix_element(m, 1, 2) * get_matrix_element(m, 2, 1);

  // sviluppo di laplace
  // considerati tutti i coefficienti di una certa riga o di una certa colonna,
  // il determinante è uguale alla somma dei prodotti di un dato coefficiente
  // per il corrispondente complemento algebrico

  // si definisce complemento algebrico dell'elemento m[i][j] il prodotto tra il
  // determinante della matrice complementare rispetto alla riga i e alla
  // colonna j, e (-1)^{i+j}
  float det = 0, sign;
  int j;
  matrix comp;

  // considera la prima colonna
  for (j = 1; j <= dim; j++) {
    // complementare
    comp = complementary(m, 1, j);

    // calcola segno
    sign = ((1 + j) % 2 == 0) ? 1.0f : -1.0f;

    // chiamata ricorsiva
    det += sign * get_matrix_element(m, 1, j) * determinant(comp);

    // dealloca matrice complementare
    destroy_matrix(comp);
  }

  return det;
}

// la matrice complementare di m rispetto alla riga i e alla colonna j
// si ottiene escludendo la riga i e la colonna j di m
matrix complementary(matrix m, const int row_c, const int col_c) {
	int i, j, h, k, rows = get_rows(m), cols = get_columns(m);

  // verifica indici
  if (row_c <= 0 || row_c > rows || col_c <= 0 || col_c > cols) {
    fprintf(stderr, "error: complementary indices (%d, %d) out of bounds for matrix (%d x %d)\n",
            row_c, col_c, rows, cols);
    exit(EXIT_FAILURE);
  }

  matrix res = alloc_matrix(rows - 1, cols -1); // m - 1 x n - 1

  // per ogni riga - h: indice riga per la matrice res
	for (i = h = 1; i <= rows; i++) {
    // escludi row_c
		if (i != row_c) {
      // per ogni colonna - k: indice colonna per la matrice res
			for (j = k = 1; j <= cols; j++)
				// escludi col_c
        if (j != col_c)
          set_matrix_element(res, h, k++, get_matrix_element(m, i, j));
			h++;
		}
	}

	return res;
}

matrix scalar_matrix_multiplication(matrix m, const float scalar) {
  matrix res = alloc_matrix(get_rows(m), get_columns(m));
  int i, j;

  // moltiplica ogni elemento di m per scalar
  for (i = 1; i <= get_rows(res); i++)
    for (j = 1; j <= get_columns(res); j++)
      set_matrix_element(res, i, j, scalar * get_matrix_element(m, i, j));

  return res;
}

// costruisci l'identità - matrice diagonale scalare in cui
// tutti gli elementi della diagonale sono pari a 1
matrix identity(const int dim) {
  int i, j;
  matrix identity = alloc_matrix(dim, dim); // inizializza a 0

  // imposta 1 sulla diagonale
  for (i = 1; i <= dim; i++)
      for (j = 1; j <= dim; j++)
          if (i == j) set_matrix_element(identity, i, j, 1);

  return identity;
}

// verifica se la matrice passata è l'identità
const int is_identity(matrix m) {
  if (!is_quadratic(m))
    return 0;

  matrix id = identity(get_rows(m));
  const int res = matrix_equality(m, id);
  destroy_matrix(id);

  return res;
}

// restituisce la matrice trasposta che sii ottiene
// scambiando le righe con le colonne
matrix transpose(matrix m) {
  int i, j;
  int rows = get_rows(m), cols = get_columns(m);

  matrix trp = alloc_matrix(cols, rows);

  for (i = 1; i <= rows; i++)
    for (j = 1; j <= cols; j++)
      set_matrix_element(trp, j, i, get_matrix_element(m, i, j));

  return trp;
}

// effettua la somma componente per componente tra matrici se sign = 1,
// se sign = 0 effettua la differenza componente per componente
matrix matrix_addition_helper(matrix m1, matrix m2, const int sign) {
  if (sign < 0 || sign > 1) {
    fprintf(stderr, "error: sign %d for matrix_addition_helper not valid", sign);
    exit(EXIT_FAILURE);
  }

  // verifica compatibilità
  if (get_rows(m1) != get_rows(m2) || get_columns(m1) != get_columns(m2)) {
    fprintf(stderr, "error during matrix addition between (%d x %d) and (%d x %d) matrices",
        get_rows(m1), get_columns(m1), get_rows(m2), get_columns(m2));
    exit(EXIT_FAILURE);
  }

  int i, j, rows = get_rows(m1), columns = get_columns(m1);
  float value1, value2;
  matrix res = alloc_matrix(rows, columns);

  for (i = 1; i <= rows; i++)
    for (j = 1; j <= columns; j++) {
      value1 = get_matrix_element(m1, i, j);
      value2 = get_matrix_element(m2, i, j);

      // se sign = 0 allora cambia segno a value2
      if (!sign) value2 *= -1;

      set_matrix_element(res, i, j, value1 + value2);
    }

  return res;
}

matrix matrix_addition(matrix m1, matrix m2) {
  return matrix_addition_helper(m1, m2, 1);
}

matrix matrix_subtraction(matrix m1, matrix m2) {
  return matrix_addition_helper(m1, m2, 0);
}

// moltiplicazione righe per colonne
// m1 (m, n), m2 (n x t) ====> m1 x m2 (m, t)
matrix row_by_column_multiplication(matrix m1, matrix m2) {
  // al fine di poter calcolare il prodotto righe per colonne tra due matrici,
  // il numero delle righe di una deve essere uguale al numero di colonne dell'altra
  if (get_columns(m1) != get_rows(m2)) {
    fprintf(stderr, "incompatible dimensions for matrix multiplication (%d x %d) * (%d x %d)\n",
            get_rows(m1), get_columns(m1), get_rows(m2), get_columns(m2));
    exit(EXIT_FAILURE);
  }

  int i, j, k;
  int rows_m1 = get_rows(m1); // m
  int cols_m1 = get_columns(m1); // n
  int cols_m2 = get_columns(m2); // t
  float sum = .0;

  matrix result = alloc_matrix(rows_m1, cols_m2); // m x t

  for (i = 1; i <= rows_m1; i++)
    for (j = 1; j <= cols_m2; j++) {
      sum = .0;
      for (k = 1; k <= cols_m1; k++)
        sum += get_matrix_element(m1, i, k) * get_matrix_element(m2, k, j);
      set_matrix_element(result, i, j, sum);
    }

  return result;
}

// uguaglianza componente per componente tra matrici
const int matrix_equality(matrix m1, matrix m2) {
  if (get_rows(m1) != get_rows(m2) || get_columns(m1) != get_columns(m2))
    return 0;

  int i, j;

  for (i = 1; i <= get_rows(m1); i++)
    for (j = 1; j <= get_columns(m1); j++)
      if (get_matrix_element(m1, i, j) != get_matrix_element(m1, i, j))
        return 0;

  return 1;
}

// calcolo dell'inversa
// l'inversa di una matrice quadrata esiste solo se il suo determinante è non nullo.
// se esiste è definita come: (1 / determinante) ⋅x matrice trasposta dei cofattori.
// la matrice dei cofattori è quella costituita da tutti i complementi algebrici
matrix inverse(matrix m) {
  // verifica matrice quadrata
  if (!is_quadratic(m)) {
    fprintf(stderr, "error: matrix must be square to compute the inverse.\n");
    exit(EXIT_FAILURE);
  }

  int dim = get_rows(m);
  float det = determinant(m);

  // determinante diverso da 0
  if (fabs(det) < 1e-6 * dim) {
    fprintf(stderr, "error: matrix is (almost) singular, cannot compute inverse.\n");
    exit(EXIT_FAILURE);
  }

  matrix cof = alloc_matrix(dim, dim); // matrice dei cofattori
  matrix comp, cof_trp, inv;
  float comp_det, sign;
  int i, j;

  // per ogni elemento [i, j], si calcola il cof[i][j] = (-1)^{i+j} x det(comp[i][j])
  for (i = 1; i <= dim; i++) {
    for (j = 1; j <= dim; j++) {
      // matrice complementare escludendo riga i e colonna j
      comp = complementary(m, i, j);
      // determinante complementare
      comp_det = determinant(comp);

      // calcolo segno
      sign = ((i + j) % 2 == 0) ? 1.0f : -1.0f;

      // imposta elemento [i][j] dei cofattori
      set_matrix_element(cof, i, j, sign * comp_det);

      // deallocazione
      destroy_matrix(comp);
    }
  }

  // trasposta dei cofattori
  cof_trp = transpose(cof);
  destroy_matrix(cof);

  // moltiplicazione scalare tra la trasposta dei cofattori e (1 / determinante)
  inv = scalar_matrix_multiplication(cof_trp, 1.0f / det);
  destroy_matrix(cof_trp);

  return inv;
}

// verifica se m ha una componente negativa
const int has_negative_component(matrix m) {
  int i, j;

  for (i = 1; i <= get_rows(m); i++)
    for (j = 1; j <= get_columns(m); j++)
      if (get_matrix_element(m, i, j) < 0)
        return 1;

  return 0;
}

// estrae una sottomatrice dalle colonne specificate dell'input m
// la matrice risultante ha lo stesso numero di righe di m e num_columns colonne
// le colonne vengono specificate tramite un array di indici (base 1)
matrix submatrix_by_columns(matrix m, int *columns, const int num_columns) {
  int total_cols = get_columns(m), total_rows = get_rows(m), i, j;

  // verifica indici
  if (num_columns <= 0 || num_columns > total_cols) {
      fprintf(stderr, "error: trying to extract %d columns from a (%d x %d) matrix\n",
              num_columns, total_rows, total_cols);
      exit(EXIT_FAILURE);
  }

  // verifica indici dell'array columns ed eventuali duplicati
  for (i = 0; i < num_columns; i++) {
      if (columns[i] <= 0 || columns[i] > total_cols) {
          fprintf(stderr, "error: invalid column index %d in a (%d x %d) matrix\n",
                  columns[i], total_rows, total_cols);
          exit(EXIT_FAILURE);
      }

      for (j = i + 1; j < num_columns; j++)
          if (columns[i] == columns[j]) {
              fprintf(stderr, "error: duplicate column index %d found\n", columns[i]);
              exit(EXIT_FAILURE);
          }
  }

  matrix res = alloc_matrix(total_rows, num_columns); // m righe, num_columns colonne
  float value;

  // copia ogni elemento delle colonne selezionate nella nuova matrice
  for (i = 1; i <= total_rows; i++)
      for (j = 0; j < num_columns; j++) {
          value = get_matrix_element(m, i, columns[j]);
          set_matrix_element(res, i, j + 1, value);
      }

  return res;
}

// concatenazione orizzontale tra matrici
// m1 (m, n), m2 (m, t) ===> m1!m2 (m, n + t)
matrix horizontal_concatenation(matrix m1, matrix m2) {
  int rows = get_rows(m1), cols1 = get_columns(m1), cols2 = get_columns(m2);

  // verifica compatibilità
  if (rows != get_rows(m2)) {
      fprintf(stderr, "error: cannot horizontally concatenate (%d x %d) and (%d x %d) matrices\n",
              rows, cols1, get_rows(m2), cols2);
      exit(EXIT_FAILURE);
  }

  matrix res = alloc_matrix(rows, cols1 + cols2);

  int i, j;
  float value;

  for (i = 1; i <= rows; i++)
    for (j = 1; j <= cols1 + cols2; j++) {
      value = (j <= cols1)
        ? get_matrix_element(m1, i, j)
        : get_matrix_element(m2, i, j - cols1);
      set_matrix_element(res, i, j, value);
    }

  return res;
}

// deallocazione
void destroy_matrix(matrix m) {
  free(m->data);
  free(m);
}
