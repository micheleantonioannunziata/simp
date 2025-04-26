simp.exe: clean mymat.o parser.o simp.o main.o
	gcc mymat.o parser.o simp.o main.o -o simp.exe

main.o: main.c
	gcc -c main.c

simp.o: ./simp/simp.c
	gcc -c ./simp/simp.c -o simp.o

parser.o: ./parser/parser.c
	gcc -c ./parser/parser.c -o parser.o

mymat.o: ./mymat/mymat.c
	gcc -c ./mymat/mymat.c -o mymat.o

clean:
	rm -f *.o *.exe