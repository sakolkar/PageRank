main:
	@mpicc -g -Wall -lm -o main src/main.c src/Lab4_IO.c

datatrim:
	@gcc -o datatrim src/datatrim.c

serialtester:
	@gcc -lm -o serialtester src/serialtester.c src/Lab4_IO.c

clean:
	@rm -f main data_input data_output serialtester datatrim
