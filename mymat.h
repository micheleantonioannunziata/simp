// tipo matrice: indicizzata da 1 a m per le righe e da 1 a n per le colonne
typedef struct mt *matrix;

// alloca una matrice m x n (m righe, n colonne)
matrix alloc_matrix(const int m, const int n);

// restituisce il numero di righe della matrice
const int get_rows(matrix m);

// restituisce il numero di colonne della matrice
const int get_columns(matrix m);

// restituisce la i-esima riga, 1 <= i <= rows
matrix get_row(matrix m, const int i);

// restituisce la j-esima colonna, 1 <= j <= columns
matrix get_column(matrix m, const int j);

// restituisce il valore in posizione [i, j], 1 <= i <= rows, 1 <= j <= columns
float get_matrix_element(matrix m, const int i, const int j);

// imposta il valore in posizione [i, j], 1 <= i <= rows, 1 <= j <= columns
void set_matrix_element(matrix m, const int i, const int j, const float value);

// inizializza tutti gli elementi della matrice al valore value
void initizialize_matrix(matrix m, const float value);

// stampa gli elementi della matrice su stdout
void print_matrix(matrix m);

// legge gli elementi della matrice da stdin
void scan_matrix(matrix m);

// verifica se la matrice è quadrata (stesso numero di righe e colonne)
int is_quadratic(matrix m);

// calcola il determinante della matrice
float determinant(matrix m);

// restituisce la matrice complementare rispetto all'elemento [i, j], 1 <= i <= rows, 1 <= j <= columns
matrix complementary(matrix m, const int i, const int j);

// moltiplica tutti gli elementi della matrice per uno scalare
matrix scalar_matrix_multiplication(matrix m, const float scalare);

// genera la matrice identità di dimensione dim
matrix identity(const int dim);

// verifica se la matrice è l'identità
const int is_identity(matrix m);

// restituisce la trasposta della matrice
matrix transpose(matrix m);

// somma due matrici componente per componente
matrix matrix_addition(matrix a, matrix b);

// sottrae due matrici componente per componente
matrix matrix_subtraction(matrix a, matrix b);

// moltiplica due matrici (righe per colonne)
matrix row_by_column_multiplication(matrix a, matrix b);

// verifica l'uguaglianza tra due matrici
const int matrix_equality(matrix a, matrix b);

// calcola la matrice inversa
matrix inverse(matrix m);

// restituisce una sottomatrice basata su un array di colonne
matrix submatrix_by_columns(matrix m, int *columns, const int count);

// verifica se esiste almeno un elemento negativo nella matrice
const int has_negative_component(matrix m);

// concatena due matrici orizzontalmente
matrix horizontal_concatenation(matrix a, matrix b);

// libera la memoria occupata dalla matrice
void destroy_matrix(matrix m);
