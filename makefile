obj_dir = obj

simp.exe: $(obj_dir) mymat.o parser.o simp.o main.o
	gcc $(obj_dir)/mymat.o $(obj_dir)/parser.o $(obj_dir)/simp.o $(obj_dir)/main.o -o simp.exe

main.o: main.c
	gcc -c main.c -o $(obj_dir)/main.o

simp.o: ./simp/simp.c
	gcc -c ./simp/simp.c -o $(obj_dir)/simp.o

parser.o: ./parser/parser.c
	gcc -c ./parser/parser.c -o $(obj_dir)/parser.o

mymat.o: ./mymat/mymat.c
	gcc -c ./mymat/mymat.c -o $(obj_dir)/mymat.o

clean:
	rm -rf $(obj_dir) *.exe

# crea obj/ se non esiste
$(obj_dir):
	mkdir -p $(obj_dir)

# ogni .o dipende da obj/
$(obj_dir)/%.o: %.c | $(obj_dir)
	gcc -c $< -o $@
