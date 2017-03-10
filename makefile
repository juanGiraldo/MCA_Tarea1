energias.pdf : plotE.py energia.txt
	python plotE.py

energia.txt : a.out
	(time ./a.out 1) #> tiempo.txt 2>&1
	(time ./a.out 2) #>> tiempo.txt 2>&1
	(time ./a.out 4) #>> tiempo.txt 2>&1

a.out : Fermi_Pasta.c
	gcc -fopenmp Fermi_Pasta.c -lm

clean :
	rm -rf energias.pdf a.out tiempo.txt
