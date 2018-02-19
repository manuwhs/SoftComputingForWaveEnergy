ELM: ./bin/ELM.o ./bin/gsl_func.o ./bin/main.o ./bin/general_func.o ./bin/GA.o ./bin/SA.o ./bin/FP.o ./bin/Graph.o ./bin/AIS.o ./bin/process_param.o ./bin/output_func.o ./bin/aux_func.o ./bin/threads.o ./bin/ES.o ./bin/bitvector_func.o
	gcc -L/usr/lib -o "ELM"  ./bin/FP.o ./bin/ELM.o ./bin/gsl_func.o ./bin/threads.o ./bin/main.o ./bin/Graph.o ./bin/SA.o ./bin/AIS.o  ./bin/output_func.o ./bin/general_func.o ./bin/ES.o ./bin/bitvector_func.o ./bin/GA.o ./bin/process_param.o ./bin/aux_func.o -lgsl -lgslcblas -lm -lpthread 
	rm -f ./codes/*.o
	rm -f ./headers/*.gch

./bin/main.o: main.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -c -lpthread main.c  -o ./bin/main.o

./bin/ELM.o: ./codes/ELM.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/ELM.c -o ./bin/ELM.o
	
./bin/gsl_func.o: ./codes/gsl_func.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/gsl_func.c -o ./bin/gsl_func.o
	
./bin/general_func.o: ./codes/general_func.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/general_func.c -o ./bin/general_func.o
	
./bin/aux_func.o: ./codes/aux_func.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/aux_func.c -o ./bin/aux_func.o
	
./bin/threads.o: ./codes/threads.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/threads.c -o ./bin/threads.o
	
./bin/output_func.o: ./codes/output_func.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/output_func.c -o ./bin/output_func.o
	
./bin/AIS.o: ./codes/AIS.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/AIS.c -o ./bin/AIS.o

./bin/GA.o: ./codes/GA.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/GA.c  -o ./bin/GA.o
	
./bin/SA.o: ./codes/SA.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/SA.c -o ./bin/SA.o
	
./bin/Graph.o: ./codes/Graph.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/Graph.c -o ./bin/Graph.o
	
./bin/ES.o: ./codes/ES.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/ES.c -o ./bin/ES.o
	
./bin/FP.o: ./codes/FP.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/FP.c -o ./bin/FP.o

./bin/bitvector_func.o: ./codes/bitvector_func.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/bitvector_func.c -o ./bin/bitvector_func.o
		
./bin/process_param.o: ./codes/process_param.c ELMheader.h
	gcc -I/usr/include/gsl -O0 -g3 -Wall -c ./codes/process_param.c -o ./bin/process_param.o

.PHONY: clean
clean:
	rm -f bin/*.o
	rm -f *.txt
	rm -R Results/RMSE/*

