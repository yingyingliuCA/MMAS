//Author: Ying Ying Liu
//Last Modified: 20160221

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#define bool int
#define true 1
#define false 0

/* constants for a random number generator, for details see numerical recipes in C */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
long int seed = 12345678;
int n; //#vertices
int m; //#edges
int numants;
int niteration;
double alpha;
double beta;
double rho;
double epsilon;
int *input_graph;

char graphFile[50];
char graphType[1];

/****************** Graph input *******************************/
//int sol_ref[12] = {1, 0, 0, 0, 6, 7, 7, 5, 7, 10, 9, 9};
//int membership_ref[n];
//int membership_ref[12] = {1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3};

/********** Parameters used by CEC paper *********************/
//double alpha = 1;
//double beta = 2;
//double rho = 0.8;
//double epsilon = 0.005;

/************ Varaibles used by CEC paper ********************/

//Variables for heuristic information calculation
double *eta; // heuristic information for formula 2
double *cor; // pearson correlation
double *stand_dev; // standard deviations for n rows
double *avg_deg; // averages for n rows
int *deg; // degrees for n vertices

//Variables for ACO
double *pheromone;
int *sol_ants;
int *membership_ants;
int *global_best_sol;
int *global_best_mem;
int *local_best_sol;
double *mod_ants;
double global_best_mod;
double local_best_mod;
double tau_min;
double tau_max;

/************** function declaration *************************/

/** Miscellaneous **/
double ran01( long *idum );
void printDoubleArray(double* array, int size);
void printAdjMatrix();
void useSampleGraph();
void generate_rand_graph(int v);
int findLineLen(char* fin);
void readDimacsGraphFromFile(char* dimacsGraphFile);
void readRmatGraphFromFile(char *rmatGraphFile);

/** Functions for Heurisitc Information **/
void calculate_stand_dev();
void calculate_cor();
void calculate_eta();

/** Functions for finding membership vector and modularity calculation **/
//int insert(int i, int num_mem, int q[][n], int count, int visited[], int mem[]);
//void bfs(int vec[], int mem[]);
int insert(int i, int num_mem, int *q, int count, int *visited, int *mem);
void bfs(int *vec, int *mem);
double calculate_mod(int mem[]);

/** ACO functions **/
void initialize_pheromone();
void clear_sol();
void solution_construction(int ant);
void update_iteration_best();
void update_global_best();
void pheromone_evaparation();
void update_pheromone();

/************** Program *************************/

double ran01( long *idum )
/*
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}

void printDoubleArray(double* array, int size)
{
  int i;
  for(i=0; i<size; i++)
  {
    printf("%.3f ",array[i]);
  }
  printf("\n");
}

void printAdjMatrix()
{
  int i;
  printf("Adjacency Matrix:\n");
  for (i = 0; i < n*n; i++)
  {
    if (i%n == 0)
    {
      printf("\n");
    }
    printf("%d ", input_graph[i]);
  }
  printf("\n");
}

void useSampleGraph()
{
  n = 12;
  m = 21;
  //n = 15;
  input_graph = (int *)malloc(n*n*sizeof(int));
  deg = (int *)malloc(n*sizeof(int));
  avg_deg = (double *)malloc(n*sizeof(double));
  int known_graph[144]= 
  {0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0,
  0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0};

  /*int known_graph[225]= {1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						   1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						   1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						   1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						   1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						   0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
						   0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
						   0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
						   0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
						   0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
						   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
						   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
						   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
						   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
						   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1};*/

  int i, j, tmp;

  for (i=0; i<n; i++)
  {
  	tmp = 0;
  	for (j=0; j<n; j++)
  	{
  		input_graph[i*n+j] = known_graph[i*n+j];
  		tmp += known_graph[i*n+j];
  	}
  	deg[i] = tmp;
  	avg_deg[i] = tmp*1.0/n;

  }
}

void generate_rand_graph(int v)
{
  int i, j, k;
  n = v;
  m = 0;
  input_graph = (int *)malloc(n*n*sizeof(int));
  deg = (int *)malloc(n*sizeof(int));
  avg_deg = (double *)malloc(n*sizeof(double));

  for(i=0; i<n; i++){
    for(j=0; j<=i; j++){
      //k = rand()%100+1;
      k=(int)(ran01(&seed)*100);
      if(k<=50)
        { input_graph[i*n+j] = 0;}
      else
        { 
        	if (i!=j)
        	{
        		input_graph[i*n+j] = 1;
          		input_graph[j*n+i] = 1;
          		deg[i*n+j]++;
          		deg[j*n+i]++;   
          		m++;      		
        	}
        }
    }
  }
  for(i=0; i<n; i++)
  {
    avg_deg[i] = deg[i]*1.0/n;
  }

  printf("generate random graph ok\n");
}

int findLineLen(char* fin)
{
  FILE *fp;
  int c, lnlen=0, maxlen=0, count=0;
  fp = fopen(fin, "r");
  if (fp == NULL)
  {
      fprintf(stderr, "Unable to open file:%s\n",fin);
      exit(-1);
  } 

  c = fgetc(fp);
  while(c!= EOF)
  {
    if(c!='\n')
    {
      lnlen++;
    }
    else
    {
      if(maxlen<lnlen)
      {
        maxlen = lnlen;
      }
      //printf("\nnode %d: line length: %d\n",count, lnlen);
      count++;
      lnlen = 0;
    }
    c = fgetc(fp);
  }
  //printf("file %s max line length = %d\n",fin, maxlen);
  fclose(fp);
  return maxlen;    
}

//Function to read benchmark graphs from 10th Dimacs Graph clustering and partitioning competition 
void readDimacsGraphFromFile(char* dimacsGraphFile)
{
  FILE *fp;
  char *line, *currstr, *nextstr; 
  int c, lnlen=0; 
  int nVertices, nEdges, fmt, nUniqueEdges=0, i, j, index1, index2, vertextID=0, edgeID, firstline=1, minusOne=1;
  int *xadj, *adjncy;

  lnlen = findLineLen(dimacsGraphFile);
  printf("%s line len:%d\n",dimacsGraphFile,lnlen);
    
  fp = fopen(dimacsGraphFile,"r");
  if (fp == NULL) {
    fprintf(stderr, "Unable to open dimacsGraphFile:%s\n",dimacsGraphFile);
    exit(-1);
  }

  line = malloc(lnlen*sizeof(char));

  while (fgets(line, lnlen, fp) != NULL) {
    if (*line == '%') continue; //comments
    else
    { 
      if (firstline == 1)
      {
        sscanf(line,"%d %d %d", &nVertices, &nEdges, &fmt); //All the test graphs so far have 0 weights. 
        if (nVertices <=0 || nEdges < 0)
        {
            fprintf(stderr, "nVertices or nEdges in %s should be greater than 0\n",dimacsGraphFile);
            exit(-1);         
        }
        printf("graph %s, %d vertices, %d edges\n",dimacsGraphFile,nVertices,nEdges);
              //initialize adjacency matrix
        n = nVertices;
        m = nEdges; 
        input_graph = (int *)malloc(n*n*sizeof(int));
      	deg = (int *)malloc(n*sizeof(int));
      	avg_deg = (double *)malloc(n*sizeof(double));
        firstline = 0;
      }
      else
      {
        //printf("vertextID=%d:%s", vertextID, line);
        char *data = line;
        int offset;
        while(sscanf(data, " %d%n", &edgeID,&offset) == 1)
        {
          data += offset;
          if (input_graph != NULL && vertextID < nVertices)
          {
             //i = (start-1) * n + (end-1);
             i = vertextID * n + edgeID-1;
             //printf("vertextID=%d, edgeID-1 = %d, i=%d\n",vertextID, edgeID-1, i);
             //There may be repeated edges in the rmat, so don't have to count them again
             if(input_graph[i] == 0)
             {
               //symmetric adjacency matrix
               j = (edgeID-1) * n + vertextID;
               //printf("j=%d\n",j);
               input_graph[i] = 1;
               input_graph[j] = 1;
	             deg[vertextID] ++;
           	   deg[edgeID-1] ++;              
               nUniqueEdges ++;
             }
          }
          //printf("vertextID=%d, edgeID=%d, offset=%d\n",vertextID,edgeID,offset);
        }
        if (!(*line == ' ' || *line == '\n')) // Without this condition, the graph file may not be read correctly due to empty lines
          vertextID ++;  
      }

    }  
  }

  fclose(fp);
  if (input_graph == NULL)
  {
    fprintf(stderr, "Error in assigning adjacency matrix from dimacsGraphFile:%s\n", dimacsGraphFile);
    exit(-1);
  }
  else
  {
    for(i=0; i<n; i++)
    {
      avg_deg[i] = deg[i]*1.0/n;
    }
    //m = nUniqueEdges; 
    printf("Read graph file: %s\n",dimacsGraphFile);
    /*printf("nUniqueEdges: %d\n",nUniqueEdges);
    printAdjMatrix(g, n);
    printf("\ndegree:\n");
    printIntArray(deg, n);
    printf("\naverage degree:\n");
    printDoubleArray(avg_deg, n);*/
  }

}

void readRmatGraphFromFile(char *rmatGraphFile)
{
  /* read parameters from file */
  FILE *fp;
  char line[128], var[2], var2[2];
  int start, end, nVertices, nEdges, nUniqueEdges=0, i, j, index1, index2, firstline=1, minusOne=1;

  fp = fopen(rmatGraphFile,"r");
  if (fp == NULL) {
    fprintf(stderr, "Unable to open rmatGraphFile:%s\n",rmatGraphFile);
    exit(-1);
  }
  while (fgets(line, sizeof (line), fp) != NULL) {
    if (*line == 'c') continue;  /* comment */
    else if (*line == 'p')
    {
      sscanf(line,"%s %s %d %d",var, var2, &nVertices, &nEdges);
      //printf("var=%s, var2=%s, nVertices = %d, nEdges = %d\n", var, var2, nVertices, nEdges);

      //initialize adjacency matrix
      n = nVertices;
      printf("n=%d\n",n);
      input_graph = (int *)malloc(n*n*sizeof(int));
      deg = (int *)malloc(n*sizeof(int));
      avg_deg = (double *)malloc(n*sizeof(double));
      //printAdjMatrix(g, n);
    }
    else if (*line == 'a')
    {
      sscanf(line,"%s%d %d",var, &start, &end);
      //RMAT output file added 1 to the actual graph vertice numbers to start with 1 instead of 0
      //Check whether vertex number start with 0 or 1
      if (firstline == 1)
      {
        if (start == 0)
        {  minusOne = 0;
        }
        firstline = 0;
      }
      if (minusOne == 1)
      {
        start = start - 1;
        end = end - 1;
      }
      //printf("%d -> %d\n", start, end);
      if (input_graph != NULL)
      {
         i = start * n + end;
         //There may be repeated edges in the rmat, so don't have to count them again
         if(input_graph[i] == 0)
         {
           //symmetric adjacency matrix
           j = end * n + start;
           input_graph[i] = 1;
           input_graph[j] = 1;
           nUniqueEdges ++;
           deg[start] ++;
           deg[end] ++;
         }
      }

    }
  }
  if (input_graph == NULL)
  {
    fprintf(stderr, "Error in assigning adjacency matrix from rmatGraphFile:%s\n", rmatGraphFile);
    exit(-1);
  }
  else
  {
    for(i=0; i<n; i++)
    {
      avg_deg[i] = deg[i]*1.0/n;
    }
    m = nUniqueEdges; 
    printf("Read graph file: %s\n",rmatGraphFile);
    /*printf("nUniqueEdges: %d\n",nUniqueEdges);
    printAdjMatrix(g, n);
    printf("\ndegree:\n");
    printIntArray(deg, n);
    printf("\naverage degree:\n");
    printDoubleArray(avg_deg, n);*/
  }

}


void calculate_stand_dev()
{
	int i, j;
	double temp;
	//printf("\nstand dev: \n");
	for(i=0; i< n; i++)
	{
		stand_dev[i] = 0; // clear the element
		temp = stand_dev[i];
		for(j=0; j<n; j++)
		{
			temp = temp + pow((input_graph[i*n+j] - avg_deg[i]), 2.0);
		}
		temp = sqrt(temp / n);
		stand_dev[i] =  temp;
		//printf("%.2f ", stand_dev[i]);
	}
}

void calculate_cor()
{
	int i, j, k;
	double temp;
	calculate_stand_dev();
	//printf("\ncorrelation:\n");
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			cor[i*n+j] = 0;
			temp = 0;
			for(k=0; k<n; k++)
			{
				temp = temp + (input_graph[i*n+k] - avg_deg[i])*(input_graph[j*n+k] - avg_deg[j]);
			}
			temp = temp / (n * stand_dev[i] * stand_dev[j]);
			cor[i*n+j] = temp;
			//printf("%.2f ", cor[i][j]);
		}
		//printf("\n");
	}
}

void calculate_eta()
{
	int i, j;
	calculate_cor();
	//printf("\nheuristic information:\n");
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			eta[i*n+j] = 1 / (1 + exp( -1.0 * cor[i*n+j]));
			//printf("%.2f ", eta[i][j]);
		}
		//printf("\n");
	}
}

//int insert(int i, int num_mem, int q[][n], int count, int visited[], int mem[])
int insert(int i, int num_mem, int *q, int count, int *visited, int *mem)
{
	visited[i] = 1;
	mem[i] = num_mem;
	//q[num_mem-1][count] = i;
	q[(num_mem-1)*n+count] = i;
	count++;
	//printf("mem[%d]=%d\n", i, num_mem);
	return count; // return number of elements in the same membership
}

//Breadth first search to find the membership vector according to the solution vector.
//void bfs(int vec[], int mem[])
void bfs(int *vec, int *mem)
{
	int i, j, k, temp;
	int count = 0, num_mem = 0;
	/*
	int visited[n], q[n][n];
	for(i = 0; i < n; i ++)
	{
		visited[i] = 0;
		for(j = 0; j < n; j ++)
		{
			q[i][j] = -1;
		}
	}*/
	int *visited, *q; 
	visited = (int *)malloc(n*sizeof(int)); 
	q = (int *)malloc(n*n*sizeof(int));

	for(i = 0; i < n; i ++)
	{
		count = 0;
		if(!visited[i])
		{
			num_mem++;
			//printf("i=%d, vec[%d]=%d\n",i,i,vec[i]);
			//For the first unvisited vertex, add both the vertex and its locus into the queue
			count = insert(i, num_mem, q, count, visited, mem);
			//printf("inserted %d\n",i);
			count = insert(vec[i], num_mem, q, count, visited, mem);
			//printf("inserted %d\n",vec[i]);
			for(j = 0; j < n; j ++)
			{
				//match all vertices and all elements in the queue
				k = 0;
				while(!visited[j] && k < count)
				{
					//temp = q[num_mem-1][k];
					temp = q[(num_mem-1)*n+k];
					//since the solution is a 1d array instead of 2d matrix, we need to be careful
					//to look at all possiblities of connected elements
					if( j == temp || vec[j] == temp || j == vec[temp] || vec[j] == vec[temp])
					{
						count = insert(j, num_mem, q, count, visited, mem);
					}
					else
					{
						k++;
					}

				}
			}
		}
	}
	
	/*printf("membership found:\n");
	for(i = 0; i < n; i++)
	{
		printf("%d ", mem[i]);
	}*/
}

double calculate_mod(int mem[])
{
	double mod, temp;
	int i, j;
	mod = 0;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(mem[i] == mem[j])
			{
				//mod = mod + (input_graph[i*n+j] - deg[i]*deg[j]*1.0/(2*n));
				mod = mod + (input_graph[i*n+j] - deg[i]*deg[j]*1.0/(2*m));
			}
			else
			{
				mod = mod + 0;
			}
		}
	}
	//mod = mod / (2*n);
	mod = mod / (2*m);
	//printf("modularity: %.3f\n", mod);
	return mod;
}

void initialize_pheromone()
{
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			pheromone[i*n+j] = 1;
		}
	}
	printf("\npheromone initialization ok\n");
}

void clear_sol()
{
	int i, j;
	for(i=0; i< numants; i++)
	{
		for (j=0; j< n; j++)
		{
			sol_ants[i*n+j] = -1;
			membership_ants[i*n+j] = -1;
		}
	}
	//printf("clear solution done\n");
}

void solution_construction(int ant)
{
  int i, j, r, s, v, sCount, vCount, index;
	int pos[n];
	int start_v = ant % n;
	double pk[n];
	double umax, pksum, usum, rnd;

	v = start_v;
	s = -1;
	sCount = 0;
	//printf("\nant %d:\n",ant);
	// choose the solution edge
	while(sCount < n)
	{
		//roulette wheel selection
		pksum = 0;
		vCount = 0;
		for (i = 0; i < n; i ++)
		{
			if (input_graph[v*n+i])
			{
				pk[vCount] = pow(pheromone[v*n+i], alpha) * pow(eta[v*n+i], beta);
				pos[vCount] = i;
				pksum += pk[vCount];
				vCount ++;
			}
		}
		//printf("v %d probability: \n", v);
		for (i = 0; i < vCount; i ++)
		{
			pk[i] = pk[i] / pksum;
			//printf("pos %d: prob %.3f\n", pos[i], pk[i]);
		}

		while(s == -1)
		{
			index = -1;
			usum = 0;
			rnd = ran01(&seed);
			//printf("rnd=%.3f\n", rnd);
			while (usum < rnd){
				index ++;
				usum += pk[index];
				//printf("index = %d, usum = %.4f, rnd = %.4f\n", index, usum, rnd);
			}
			s = pos[index];
			//printf("s=%d\n", s);
		}

		sol_ants[ant*n+v] = s;
		//printf("ant=%d,n=%d,v=%d,sol_ants[%d]=%d\n",ant,n,v,ant*n+v,s);
		v = (v + 1) % n;
		s = -1;
		sCount ++;
	} // end while

	/*printf("\nant %d solution:\n", ant);
	for(i=0; i< n; i++)
	{
		printf("\nant=%d, n=%d, i=%d:",ant,n,i);
		printf("sol:%d \n", sol_ants[ant*n+i]);
	}*/
	//printf("End of solution construction\n");

}

void update_iteration_best()
{
	int i, index=0;
	double temp=0;

	for (i = 0; i < numants; i++)
	{
		//printf("finding membership for ant %d: sol_ants[%d]=%d, membership_ants[%d]=%d\n", i, i*n, sol_ants[i*n], i*n, membership_ants[i*n]);
		bfs(&sol_ants[i*n], &membership_ants[i*n]);
		mod_ants[i] =  calculate_mod(&membership_ants[i*n]);
		if (mod_ants[i] > temp)
		{
			temp = mod_ants[i];
			index = i;
		}
	}
	local_best_mod = temp;
	for(i = 0; i < n; i ++)
	{
		local_best_sol[i] = sol_ants[index*n+i];
	}
	//printf("iteration best mod: %.3f\n",local_best_mod); 

}

void update_global_best()
{
	int i;
	if (local_best_mod > global_best_mod)
	{
		global_best_mod = local_best_mod;
		tau_max = global_best_mod / (1 - rho);
		tau_min = epsilon * tau_max;
		for (i = 0; i < n; i ++)
		{
			global_best_sol[i] = local_best_sol[i];
		}
		printf("global best mod: %.3f\n", global_best_mod);
	}

	//printf("local_best_mod: %.3f, global_best_mod:%.3f\n", local_best_mod, global_best_mod);
}

void pheromone_evaparation()
{
	int i, j;
	for(i = 0; i < n; i ++)
	{
		for (j = 0; j < n; j ++)
		{
			pheromone[i*n+j] = rho * pheromone[i*n+j];
			if (pheromone[i*n+j] < tau_min)
			{
				pheromone[i*n+j] = tau_min;
			}
			// unnecessary, but add it in case
			if (pheromone[i*n+j] > tau_max)
			{
				pheromone[i*n+j] = tau_max;
			}

		}
	}
}

void update_pheromone()
{
	int i, j;
	pheromone_evaparation();
	for(i = 0; i < n; i ++)
	{
		j = local_best_sol[i];
		pheromone[i*n+j] += local_best_mod;
		if (pheromone[i*n+j] > tau_max)
		{
			pheromone[i*n+j] = tau_max;
		}
		// unnecessary, but add it in case
		if (pheromone[i*n+j] < tau_min)
		{
			pheromone[i*n+j] = tau_min;
		}
	}

}

void usage() {

	fprintf(stderr, "aco_clustering [-options]\n");
        fprintf(stderr, "\t-c ###  graph file to use\n");
        fprintf(stderr, "\t-t ###  2 graph types to use: r (for rmat) or m (for metis)\n");
        fprintf(stderr, "\t-a ###  number of ants\n");
        fprintf(stderr, "\t-i ###  number of iterations\n");
        fprintf(stderr, "\t-l ###  alpha\n");
        fprintf(stderr, "\t-b ###  beta\n");
        fprintf(stderr, "\t-r ###  rho\n");
        fprintf(stderr, "\t-e ###  epsilon\n");
        fprintf(stderr, "For example: ./aco_clustering -c karate.gr -t r -a 30 -i 40 -l 1 -b 2 -r 0.8 -e 0.005\n");
	exit(-1);
}


void parseUserInput(int argc, char** argv) {
	int c;
  //printf("argc = %d\n",argc);
    int numantsSpecified = 0;
  	int niterationSpecified = 0;
    int alphaSpecified = 0;
    int betaSpecified = 0;
    int rhoSpecified = 0;
    int epsilonSpecified = 0;  

	if (argc == 1) {
		usage();
	} else if (argc <= 21) {

		while((c = getopt(argc, argv, "c:t:a:i:l:b:r:e:")) != -1) {

			switch (c) {

		  case 'c':
			strcpy(graphFile, "testGraphs/");
			strcat(graphFile, optarg);
			break;

		  case 't':
			strcpy(graphType, optarg);
			if ( !(strcmp(graphType, "r")==0 || strcmp(graphType, "m")==0))
				usage();
			break;

      case 'a':
			numants = atol(optarg);
			if(numants <= 0)
			  usage();
			else
			  numantsSpecified = 1;
			break;

      case 'i':
			niteration = atol(optarg);
			if(niteration <= 0)
			  usage();
			else
			  niterationSpecified = 1;
			break;

      case 'l':
      alphaSpecified = 1;  
      alpha = (double)atof(optarg);
      break;     

      case 'b':
      betaSpecified = 1;  
      beta = (double)atof(optarg);
      break; 

      case 'r':
      rhoSpecified = 1;  
      rho = (double)atof(optarg);
      break;   

      case 'e':
      epsilonSpecified = 1;  
      epsilon = (double)atof(optarg);
      break;                 

		default:
			usage();
			}

      if(numantsSpecified == 0)
        numants = 30;
      if(niterationSpecified == 0)
        niteration = 100;
      if(alphaSpecified == 0)
        alpha = 1; 
      if(betaSpecified == 0)
        beta = 2; 
      if(rhoSpecified == 0)
        rho = 0.8;
      if(epsilonSpecified == 0)
        epsilon = 0.005;
		}

	} else {
		fprintf(stderr, "Invalid input arguments\n");
		usage();
	}

  printf("parameters:\ngraphFile:%s,numants:%d,niteration:%d,alpha:%.3f,beta:%.3f,rho:%.3f,epsilon:%.3f\n",
    graphFile, numants, niteration,alpha,beta,rho,epsilon);

}


int main(int argc, char *argv[]){

	clock_t start = clock(), diff;

 	//working index variables
  int i, j, k;

 	//Generate graph and calculate modularity for sample solution
	parseUserInput(argc, argv);
	if (strcmp(graphType, "r")==0)
		readRmatGraphFromFile(graphFile);
	if (strcmp(graphType, "m")==0)
		readDimacsGraphFromFile(graphFile);
	//generate_rand_graph();
	//useSampleGraph();

	//printAdjMatrix();

	//Calculate heurisitic information for solution construction

	stand_dev = (double *)malloc(n*sizeof(double));
	cor = (double *)malloc(n*n*sizeof(double));
	eta = (double *)malloc(n*n*sizeof(double));
	pheromone = (double *)malloc(n*n*sizeof(double));
	sol_ants = (int *)malloc(numants*n*sizeof(int));
	membership_ants = (int *)malloc(numants*n*sizeof(int));
	mod_ants = (double *)malloc(numants*sizeof(double));
	global_best_sol = (int *)malloc(n*sizeof(int));
	local_best_sol = (int *)malloc(n*sizeof(int));
	global_best_mem = (int *)malloc(n*sizeof(int));

	calculate_eta();
	//printf("eta:\n");
	//printDoubleArray(eta, n*n);

	//Initialization
	initialize_pheromone();
	global_best_mod = 0;
	local_best_mod = 0;
	printf("\nniteration:%d,numants=%d\n",niteration,numants);
	//Start ACO
	for(i = 0; i < niteration; i++)
	{
		clear_sol();
		for(j = 0; j < numants; j ++)
		{
			solution_construction(j);
		}
		//printf("Iteration %d\n",i);
		update_iteration_best();
		update_global_best();
		update_pheromone();
	} //End ACO

	printf("\nACO result: \nmodularity: %.3f", global_best_mod);
	/*printf("\nglobal best solution:\n");
	for(i = 0; i < n; i++)
	{
		printf("%d ", global_best_sol[i]);
	}*/
	bfs(global_best_sol, global_best_mem);
  k = 0; //Find the largest cluster id
	printf("\nglobal best membership:\n");
	for(i = 0; i < n; i++)
	{
		printf("%d,", global_best_mem[i]);
    if(k<global_best_mem[i])
      k = global_best_mem[i]; 
	}
  k++;//Number of clusters

	diff = clock() - start;
	int msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("\nTime taken: %d seconds %d milliseconds", msec/1000, msec%1000);

  //write results to file
	  FILE *fp;
	  fp = fopen("acoexp.csv","a");
	  if(fp==NULL)
	  {
	    fp = fopen("acoexp.csv","w+");
	    fprintf(fp, "graphFile,numV,numE,numants,niteration,msec,mod,k,alpha,beta,rho,epsilon\n");
	  }
	  	fprintf(fp, "%s,%d,%d,%d,%d,%d,%.5f,%d,%.5f,%.5f,%.5f,%.5f\n",
	    graphFile,n,m,numants,niteration,msec,global_best_mod,k,alpha,beta,rho,epsilon);
	  fclose(fp);


}// end main
