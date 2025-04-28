#include <stdio.h>
#include "./simp/simp.h"
#include "./mymat/mymat.h"

int main (void) {
  FILE *f = fopen("prob.txt", "r");
  simp s = catch_problem(f);

  solve_problem(s);

  fclose(f);
  clean_all(s);
}
