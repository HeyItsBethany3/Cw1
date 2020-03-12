CC = g++ # compiler
.PHONY: all compile clean
all: compile run

# Cpp file names
file1 = LSQ
file2 = Driver


compile: ${file1}.cpp ${file2}.cpp
	${CC} -g -c -o ${file1}.o ${file1}.cpp
	${CC} -g -c -o ${file2}.o ${file2}.cpp
	${CC} -g -o Execute ${file1}.o ${file2}.o

run: Execute
	./Execute

clean:
	rm *.o Execute
