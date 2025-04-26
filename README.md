# simp

Il progetto **simp** è una banale implementazione in linguaggio C del **metodo del simplesso**. Tale algoritmo consente di risolvere all'ottimo problemi di programmazione matematica lineare purché siano in forma standard di minimo e purché sia data in input una base ammissibile. 

La soluzione proposta si basa su una serie di ipotesi fondamentali: la matrice dei coefficienti tecnologici deve essere a rango pieno, ovvero il numero di vincoli deve essere strettamente inferiore al numero di variabili; la funzione obiettivo deve essere convessa; inoltre, il problema deve essere fornito in uno specifico formato descritto di seguito.

## Moduli

### 1. **mymat**
Il modulo **mymat** è responsabile della gestione delle matrici, un aspetto fondamentale del metodo del simplesso. Esso contiene le funzioni per:
- Allocare, deallocare e manipolare matrici.
- Eseguire operazioni matematiche elementari su matrici, come la somma componente per componente, la moltiplicazione righe per colonne la trasposizione e l'inversione.
- Determinare la matrice complementare rispetto una certa riga ed una certa colonna, calcolare il determinante di una matrice.

### 2. **parser**
Il modulo **parser** si occupa della lettura e dell'analisi del file di input che contiene il problema di programmazione lineare. Le funzioni principali del modulo sono:
- Lettura del file `prob.txt`.
- Estrazione e salvataggio dei dati riguardanti la funzione obiettivo, coefficienti di costo, coefficienti tecnologici, termini noti, variabili, indici di base.
- Verifica della correttezza del formato del file di input.

### 3. **simp**
Il modulo **simp** implementa l'algoritmo del simplesso vero e proprio. Le sue funzioni includono:
- Costruzione della base iniziale.
- Esecuzione delle iterazioni del metodo del simplesso.
- Determinazione della soluzione ottima o la diagnosi di un problema non risolvibile.

Il modulo **simp** interagisce con **mymat** per le operazioni sulle matrici e con **parser** per leggere e caricare il problema da risolvere.

## Formato file input

Il programma richiede un file di input chiamato `prob.txt`, che contiene il problema di programmazione lineare da risolvere. Il file deve seguire il formato descritto di seguito.

### Funzione obiettivo

La prima riga del file definisce la funzione obiettivo. Essa deve essere nel formato `min c1x1 + c2x2 + ... + cn*xn`.
Dove `c1, c2, ..., cn` sono i coefficienti di costo delle variabili `x1, x2, ..., xn`.

### Vincoli

Ogni riga successiva definisce una restrizione del problema nella forma `a1x1 + a2x2 + ... + an*xn = b`.
Dove `a1, a2, ..., an` sono i coefficienti tecnologici delle variabili, e `b` è il termine noto del vincolo.

### Variabili in base

L'ultima riga del file contiene l'elenco degli indici delle variabili in base, in formato `B = {i1, i2, ..., im}`

### Esempio di `prob.txt`

Un esempio di file `prob.txt` potrebbe essere il seguente:

```
min 3x1 -2x2
 -x1 +x2 +x3 = 4
 x1 +x4 = 3
 x1 +3x2 -x5 = 3
B = {2, 3, 4}
```

In questo esempio:
- La funzione obiettivo è `3x1 - 2x2`.
- Ci sono tre vincoli:
  - `-x1 +x2 +x3 = 4`
  - `x1 +x4 = 3`
  - `x1 +3x2 -x5 = 3`
- Ci sono cinque variabili: `x1, x2, x3, x4, x5`.
- Il vettore dei termini noti è `[4, 3, 3]`.
- La matrice dei coefficienti tecnologici è
  ```
  -1  1  1   0  0 
   1  0  0   1  0 
   1  3  0   0 -1 
  ```
- Le variabili in base sono **x2, x3, e x4**.

### Dettagli formato

- Le variabili sono indicate con `x1`, `x2`, ecc., e devono essere numerate consecutivamente a partire da 1.
- Le espressioni della funzione obiettivo e dei vincoli devono essere espresse in forma lineare.
- Gli operatori `+` e `-` devono essere preceduti da uno spazio.
- Le variabili in base sono identificate dalla lista di indici fornita nell'ultima riga.

## Utilizzo

1. **Compilazione progetto**:
   Per compilare il progetto, basta utilizzare il comando `make` nella directory del progetto:
   ```bash
   make
   ```

2. **Preparazione file input**: Crea il file `prob.txt` nel formato descritto sopra.

3. **Esecuzione programma**: Una volta compilato il progetto e preparato il file di input, è possibile eseguire il programma con il comando `./simp`.
