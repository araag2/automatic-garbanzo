/*--------------------------------------------------------------------
| Sistemas Operativos 2017/2018 1o Semestre 1a
| Artur Guimaraes, Nº86389  Diogo Fernandes, Nº86408
---------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <pthread.h>

#include "matrix2d.h"

/*--------------------------------------------------------------------
| Function: Macros
| Abstraction of some concepts for the programmer/ making the code look
| better.
---------------------------------------------------------------------*/

#define size_formula    	  ((getColumns)/getNThreads) // Formula for each slice, we could've joined the two of them
#define size_formula_main 	  (N/trab)                   // but it would cost us some abstraction.
#define getColumns 			  args->Columns
#define getNThreads 		  args->numThreads
#define getNIterations 		  args->numIterations
#define getLine               args->myLine
#define getMaxD               args->maxD
#define getMatrix             matrix
#define getMatrix_Aux         matrix_aux
#define ResetCounter(i)       i=0

    
/*--------------------------------------------------------------------
| Function: thread_struct
| The struct that saves the threads arguments
|
| Columns-> Number of columns in the matrix
|
| numThreads-> The number of threads that exist. It's used so the thread
| has a notion wether it's the last one or not. Also serves the purpose
| of calculating the size_formula, which discovers the size of each slice
| of matrix.
|
| numIterations-> The number of iterations each thread is running the heat
| simulation. As this number grows, so does the stability of the system
| and the median.
|
| myLine-> This field contains a variable of type int which defines the starting line
| for the thread's slice of simulation. It's easier to assign a line to each
| thread. Other implementations would calculate this through a specific thread
| id and total Colum number math (admiting the square matrix conditions).
| 
| maxD-> maxD defines the minimum heat difference that can occur in a given
| iteration. If this condition is not met, the program should stop calculating 
| heat and exit every thread after the given iteration finishes.
---------------------------------------------------------------------*/

typedef struct {
  int Columns;
  int numThreads;
  int numIterations;
  int myLine;
  double maxD;
} args_mythread_t;

/*--------------------------------------------------------------------
|||||||||||||||||||||||||GLOBAL VARIABLES|||||||||||||||||||||||||||||
---------------------------------------------------------------------*/

DoubleMatrix2D *matrix;				//Global matrix
DoubleMatrix2D *matrix_aux;			//Global matrix
pthread_mutex_t mutex;				//Mutex for the critical region
pthread_cond_t cond;				//Condition variable for all threads except the last
int wthreads;						//How many threads did the math in this iteration.
int ithreads;						//You can only stop waiting if all threads have completed the math for this iteration. Incremented by the last thread
int ref_maxD;     					//If =0, condition not met, if =1 condition is met, if =-1, program should shutdown orderly.
/*--------------------------------------------------------------------
| Function: MatrixInitialization
| Sets the starting values of the edges on the matrix. Matrix needs to be
| created first. 
| tSup-> 	First line
| tLeft->	First column
| tInf->    Last line
| tRight->	Last column
---------------------------------------------------------------------*/

void MatrixInitiation(DoubleMatrix2D *matrix,int N,double tSup,double tInf,double tLeft,double tRight){
  dm2dSetLineTo (matrix, 0, tSup);
  dm2dSetLineTo (matrix, N+1, tInf);
  dm2dSetColumnTo (matrix, 0, tLeft);
  dm2dSetColumnTo (matrix, N+1, tRight);
}

/*--------------------------------------------------------------------
| Function: parse_integer_or_exit
| Parses an integer using sscanf or exits the program without returning 
| the value
---------------------------------------------------------------------*/

int parse_integer_or_exit(char const *str, char const *name)
{
  int value;
  if(sscanf(str, "%d", &value) != 1) {
    fprintf(stderr, "\nError in argument \"%s\".\n\n", name);
    exit(1);
  }
  return value;
}

/*--------------------------------------------------------------------
| Function: parse_double_or_exit
| Parses a double using sscanf or exits the program without returning 
| the value
---------------------------------------------------------------------*/

double parse_double_or_exit(char const *str, char const *name)
{
  double value;
  if(sscanf(str,"%lf",&value) != 1) {
    fprintf(stderr,"\nError in argument\"%s\".\n\n",name);
    exit(1);
  }
  return value;
}

/*--------------------------------------------------------------------
| Function: parse_multiple_or_exit
| Parses a integer that must be a multiple of N using sscanf or exits 
| the program without returning the value
---------------------------------------------------------------------*/

int parse_multiple_or_exit(char const *str, char const *name, int N)
{
  int value;
  if(sscanf(str, "%d", &value) != 1 && N%value != 0) {
    fprintf(stderr, "\nError in argument \"%s\".\n\n", name);
    exit(1);
  }
  return value;
}

/*--------------------------------------------------------------------
| Function: Thread_simul
|
| General Info: This is the function where all the works is done. Each thread runs this function
| at the same time (give or take) and does the math based on our heat formulas for it's slice. All
| results are written at the same time on the current matrix_aux, since we read the results from the
| last iteration from the current matrix. Both matrix_aux and matrix are global, being a shared memory
| area for all threads.
| 
| Synchronizing: Since threads read the previous iteration values from one thread and write the values
| for the current iteration on the other one, it is possible for all threads to execute the math for their
| slice at the same time. The only condition for this is to guarantee that all threads are in the same iteration
| at a given time (or if there differences, that each thread waits for the rest).
|
| Critical Section: Our critical section is guarded by a mutex to assure only one thread accesses it at any given
| time. This section has 2 main mechanisms to assure everything runs smoothly: 
|
|		wthreads -> A  globalvariable of type intcalled wthreads, which counts the number of threads that have reached that 
|		area and compares it to the total number of threads, to identify which thread is the last onein any given iteration. 
|		If the thread is not the last one, then she must wait for all others to complete. It's reset at the end of every 
|		iteration by the last thread to reach the if condition after the locked mutex.
|
|		ithreads -> The threads wait in a while loop, which checks a global int variable called ithreads. 
| 		Serves to compare it's value to the current iteration to see if the thread should be waiting (to 
|		solve eventual spurious wakeups) or if it should move on. This variable is incremented by the last
|		of the threads in any given iteration, after it switched the matrix's pointers around (so they swap
|		the matrix's in which they read and write), and just before the broadcast to wake all other threads up,
|		to minimize any potential problems.
|		
| Conditional variable: We use a pthread_cond_t type variable for the threads to wait upon paired with a pthread_mutex_t,
| making sure all threads wait for others to complete the iteration and that the last thread wakes up all others threads.		
---------------------------------------------------------------------*/

 void *Thread_simul(void *a) {

  args_mythread_t *args   = (args_mythread_t *) a;

  DoubleMatrix2D *temp;

  int iter,i,j;
  double value,dif;

  for(iter=0; iter < getNIterations; iter++){
     for (i=getLine; i<getLine + size_formula; i++){                                       
      for (j = 1; j<getColumns+1; j++) {                                     
        value = (dm2dGetEntry(getMatrix,i-1,j) + dm2dGetEntry(getMatrix,i+1,j) +
        dm2dGetEntry(getMatrix,i,j-1) + dm2dGetEntry(getMatrix,i,j+1)) / 4.0;
        dm2dSetEntry(getMatrix_Aux,i,j,value);
        if(!(ref_maxD)){
          dif=value-dm2dGetEntry(getMatrix,i,j);
          if(dif>=getMaxD)
            ref_maxD=1;
        }
      }
     }                                                                     


    if(pthread_mutex_lock(&mutex)!=0){
      fprintf(stderr,"\nError blocking mutex\n");
      return NULL;												// Opting to use return NULL instead of pthread_exit() because of a slight valgrind
    }															// problem associated with it.


    if(wthreads++<getNThreads-1){
      while(ithreads!=(iter+1))									// While ithreads is different than iter+1 symbolizes a wait condition for the last Thread to wake the waiting thread up
      	pthread_cond_wait(&cond,&mutex);
 	}

    else{	
       ResetCounter(wthreads);   				 				// Reset the thread counter which they use to find out if they're the last thread that was to run that given iteration

      if(!(ref_maxD)){
        ref_maxD=-1;

      }
      else                  
        ResetCounter(ref_maxD);			 						// Reseting the variable which checks if the maxD condition was met
      

      temp = getMatrix_Aux;				 						// Switching the pointers around
      getMatrix_Aux = getMatrix;                                                 
      getMatrix = temp;
      
      ithreads++;												// Increases the counter which let's the threads getout of the waiting period. Used
      pthread_cond_broadcast(&cond);							// to prevent possible spurious wakeup.	
    }

    if(pthread_mutex_unlock(&mutex) != 0) {						// Smart unlock with an error message
      fprintf(stderr, "\nError unblocking mutex\n");
      return NULL;
    }

    if(ref_maxD==-1)								    		// If the thread makes plus two itterations after stopping dynamicly it exits
        return NULL;
  }

  return NULL;
 }

/*--------------------------------------------------------------------
| Function: Master_thread
|
| General Info:
| This function server as the master thread of the function, being responsible 
| for initializing and terminating each thread, as well as sending the initial
| data of the problem. It also serves to initialize the matrixes, mutexes, condition
| variables and destroying/freeing them. Has appropriate error messages for each
| fail condition.
---------------------------------------------------------------------*/
void Master_thread(int N, double tLeft, double tSup, double tRight, double tInf, int Iterations, int trab, double maxD, char* fichS, int periodoS){

  args_mythread_t *slave_args;
  pthread_t       *slaves;
  int i;

  FILE* file = fopen(fichS, "r");
  if (file != NULL){
    matrix = readMatrix2dFromFile(file, N+2, N+2);
    fclose(file);
  }

  if(matrix == NULL){
    matrix = dm2dNew(N+2, N+2);
    MatrixInitiation(matrix, N, tSup, tInf, tLeft, tRight);
  }
  
  matrix_aux = dm2dNew(N+2, N+2);
  dm2dCopy(matrix_aux, matrix);

  if (matrix == NULL || matrix_aux == NULL) {
    fprintf(stderr, "\nError: Could not alocate the memory to generate the matrix or the matrix_aux\n\n");
    return;
  }

  dm2dCopy(matrix_aux, matrix);
 
  slave_args = (args_mythread_t*)malloc(trab*sizeof(args_mythread_t));              // Allocates the necessary space for everything
  slaves     = (pthread_t*)malloc(trab*sizeof(pthread_t));
  
  if(pthread_mutex_init(&mutex, NULL) != 0) {                                       // Making sure the mutex initializes well
        fprintf(stderr, "\nError initializing mutex\n");
        return;
  }

  if(pthread_cond_init(&cond, NULL) != 0) {                                         // Making sure the conditional variable initializes well
        fprintf(stderr, "\nError initializing the condition\n");
        return;
  }

  for (i=0; i<trab; i++){                                                         // Creates the each thread and associates it's values
      slave_args[i].Columns = N;
      slave_args[i].numThreads = trab;
      slave_args[i].numIterations = Iterations;
      slave_args[i].myLine = 1+i*size_formula_main;
      slave_args[i].maxD=maxD;
      pthread_create(&slaves[i], NULL, Thread_simul, &slave_args[i]);   
  }

  for(i=0; i<trab; i++){                                                        // Wait for each thread to terminate
    pthread_join(slaves[i],NULL);
  }


  if(pthread_cond_destroy(&cond) != 0){                                        // Making sure we destroy the conditional variable well
        fprintf(stderr, "\nError destroying variable\n");
        return;
  }

  if(pthread_mutex_destroy(&mutex) != 0){                                      // Making sure we destroy the mutex well 
        fprintf(stderr, "\nError destroying mutex\n");
        return;
  }

  dm2dPrint(matrix);

  file = fopen(fichS, "w");
  writeMatrixToFile(matrix, file, N+2);
  fclose(file);

  file = fopen(fichS, "r");
  matrix = readMatrix2dFromFile(file, N+2, N+2);
  fclose(file);


  dm2dPrint(matrix);
  free(slave_args);
  free(slaves);
  dm2dFree(matrix);                                                              // Free's everything after the threads have been terminated
  dm2dFree(matrix_aux); 
  return;
}

/*--------------------------------------------------------------------
| Function: main
|
| Nothing of big importance to note, some parsers and conditions to make
| sure the inputs actually make sense in the context of our work. Sends
| all the work to a different function for simplicity sake and to tidy
| things up a bit.
|
---------------------------------------------------------------------*/

int main (int argc, char** argv) {


  if(argc != 11) {																  //Just checking the number of arguments
    fprintf(stderr, "\nInvalid number of arguments while parsing\n");
    fprintf(stderr, "Arguments: heatSim N tLeft tSup tRight tInf Iterations trab maxD\n\n");
    return 1;
  }

  /* argv[0] = program name */
  int N = parse_integer_or_exit(argv[1], "N");
  double tLeft = parse_double_or_exit(argv[2], "tLeft");
  double tSup = parse_double_or_exit(argv[3], "tSup");
  double tRight = parse_double_or_exit(argv[4], "tRight");
  double tInf = parse_double_or_exit(argv[5], "tInf");
  int Iterations = parse_integer_or_exit(argv[6], "Iterations");                    // Parses everything we need to parse
  int trab = parse_multiple_or_exit(argv[7], "trab",N);
  double maxD = parse_double_or_exit(argv[8], "maxD");
  char* fichS = argv[9];
  int periodoS = parse_integer_or_exit(argv[10], "periodoS");

  

  fprintf(stderr, "\nArguments:\n"                                                 
	" N=%d tLeft=%.1f tSup=%.1f tRight=%.1f tInf=%.1f Iterations=%d trab=%d maxD=%.4f fichS=%s periodoS=%d \n",
	N, tLeft, tSup, tRight, tInf, Iterations, trab, maxD, fichS, periodoS);

  if(N < 1 || tLeft < 0 || tSup < 0 || tRight < 0 || tInf < 0 || Iterations < 1 || N%trab !=0 || maxD<0) { // Error conditions and error messages
    fprintf(stderr, "\nError: Invalid arguments\n"
	" Just a quick reminder that N >= 1, temperatures >= 0.0, Iterations >= 1, N must be a multiple of trab , maxD>=0.0\n\n");
    return 1;
  }

  Master_thread(N,tLeft,tSup,tRight,tInf,Iterations,trab,maxD, fichS, periodoS);                 // Start's the master thread with everything it needs

  return 0;
}
