// effettua una serie di controlli sulla funzione obiettivo parsata
// verifica che non ci siano variabili duplicate, ad esempio -x1 + 4x1
// controlla il 'min' iniziale
// controlla la validità dei coefficienti - il segno deve esser non separato alla variabile
// controlla la validità dei indici delle variabili
void verify_goal_function(FILE*);

// restituisce, leggendo dal file, il numero totale di variabili del problema
int parse_num_variables(FILE*);

// restituisce, leggendo dal file il vettore colonna dei coefficienti di costo
matrix parse_cost_coeffs(FILE*);

// restituisce, leggendo dal file, la matrice dei coefficienti tecnologici
matrix parse_tech_coeffs(FILE*);

// restituisce, leggendo dal file, il vettore colonna dei termini noti
matrix parse_known_terms(FILE*);

// restituisce il numero di vincoli e imposta pos alla posizione iniziale del puntatore a file (dopo la funzione obiettivo)
const int parse_num_contraints(FILE *f, long *pos) ;