main:
	@gcc -g -Wall -o main src/main.c src/Lab4_IO.c -lm

datatrim:
	@gcc -o datatrim src/datatrim.c

serialtester:
	@gcc -o serialtester src/serialtester.c src/Lab4_IO.c -lm

clean:
	@rm -f main data_input data_output serialtester datatrim
