/* PageRank */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* allocate one object of given type */
#define NEW(type) ((type*)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type (TAB)*/
#define NEW_A(num,type) ((type*)calloc((size_t)(num),(size_t)sizeof(type)))

typedef unsigned int u_int;
typedef double Real;

/* vector definition */
typedef struct
{
  u_int dim;
  Real *ve;
} VEC;

/* matrix definition */
typedef struct
{
  u_int m, n;
  Real **me;
} MAT;

/* sparse matrix definition */
typedef struct graphnode
{
  u_int col;
  Real val;
  struct graphnode *next;
} NODE;

typedef struct
{
  u_int m, n;
  NODE *rows;
} SMAT;

/* v_get -- gets a VEC of dimension 'dim'
   Precondition: size >= 0
   Postcondition: initialized to zero */
VEC *v_get(u_int size)
{
  VEC *v;
  
  if ( (v = NEW(VEC)) == (VEC *)NULL )
  {
    fprintf(stderr, "v_get memory error");
    exit(-1);
  }
  
  v->dim = size;
  if ( (v->ve = NEW_A(size,Real)) == (Real *)NULL )
  {
    free(v);
    fprintf(stderr, "v_get memory error");
    exit(-1);
  }
  
  return (v);
}

/* v_free -- returns VEC & associated memory back to memory heap */
int v_free(VEC *vec)
{
  if ( vec == (VEC *)NULL )
    return (-1);
  
  if ( vec->ve == (Real *)NULL ) 
  {
    free(vec);
  }
  else
  {
    free(vec->ve);
    free(vec);
  }
  
  return (0);
}

SMAT *sm_get(u_int m, u_int n)
{
  SMAT *G;

  if ( (G = NEW(SMAT)) == (SMAT *)NULL )
  {
    fprintf(stderr, "sm_get memory error");
    exit(-1);
  }
  
  G->m = m ; G->n = n;

  if ((G->rows = NEW_A(m,NODE)) == (NODE *)NULL )
  {
    free(G);
    fprintf(stderr, "sm_get memory error");
    exit(-1);
  }

  for (u_int i=0 ; i<G->m ; i++)
    (G->rows[i]).val = -1;

  return (G);
}

int sm_free(SMAT *G)
{
  if ( G == (SMAT *)NULL )
    return (-1);
  
  if ( G->rows == (NODE *)NULL ) 
  {
    free(G);
  }
  else
  {
    NODE *n0;
    NODE *n1;
    for (u_int i=0 ; i<G->m ; i++)
    {
      n0 = &(G->rows[i]);
      if (n0->val < 0.0) break; /* empty line */
      n0 = n0->next;
      while (n0->val >= 0.0)
      {
        n1 = n0->next;
        free(n0);
        n0 = n1;
      }
      free(n0);
    }
    free(G->rows);
    free(G);
  }
  
  return (0);
}

NODE *sm_add(NODE *n0, u_int c, Real v)
{
  NODE *n1;
  n0->col = c;
  n0->val = v;
  if ( (n1 = NEW(NODE)) == (NODE *)NULL )
  {
    fprintf(stderr, "sm_add memory error");
    exit(-1);
  }
  n1->val = -1;
  n0->next = n1;
  return (n1);
}

/* m_get -- gets an mxn matrix by dynamic memory allocation 
   Precondition: m>=0 && n>=0
   Postcondition: initialized to zero */
MAT *m_get(u_int m, u_int n)
{
  MAT *g;
  
  if ( (g = NEW(MAT)) == (MAT *)NULL )
  {
    fprintf(stderr, "m_get memory error");
    exit(-1);
  }
  
  g->m = m ; g->n = n;

  if ((g->me = NEW_A(m,Real*)) == (Real **)NULL )
  {
    free(g);
    fprintf(stderr, "m_get memory error");
    exit(-1);
  }
  
  for ( int i = 0; i < m; i++ )
    if ( (g->me[i] = NEW_A(n,Real)) == (Real *)NULL )
    {
      fprintf(stderr, "m_get memory error");
      exit(-1);
    }
  
  return (g);
}

/* m_free -- returns MAT & associated memory back to memory heap */
int m_free(MAT *mat)
{
  if ( mat == (MAT *)NULL )
    return (-1);
  
  for ( int i = 0; i < mat->m; i++ )
    if ( mat->me[i] != (Real *)NULL ) free(mat->me[i]);

  if ( mat->me != (Real **)NULL ) free(mat->me);
  
  free(mat);
  
  return (0);
}

/* m_input -- file input of matrix */
MAT *m_input(FILE *fp)
{
  MAT *g;
  u_int m,n,val;
  
  /* get dimension */
  if ( fscanf(fp," Matrix: %u by %u",&m,&n) < 2 )
  {
    fprintf(stderr, "m_input error reading dimensions");
    exit(-1);
  }
  
  /* allocate memory */
  g = m_get(m,n);
  
  /* get entries */
  for ( u_int i = 0 ; i < m; i++ )
  {
    if ( fscanf(fp," row %u:",&val) < 1 )
    {
      fprintf(stderr, "m_input error reading line %u", i);
      exit(-1);
    }
    for ( u_int j = 0; j < n; j++ )
      if ( fscanf(fp,"%lf",&g->me[i][j]) < 1 )
      {
        fprintf(stderr, "m_input error reading line %u col %u", i, j);
        exit(-1);
      }
  }
  
  return (g);
}

/* sm_input -- file input of sparse matrix */
SMAT *sm_input(FILE *fp)
{
  SMAT *g;
  u_int m,n,row;
  Real col;
  NODE *n0;
  
  /* get dimension */
  if ( fscanf(fp," SparseMatrix: %u by %u",&m,&n) < 2 )
  {
    fprintf(stderr, "sm_input error reading dimensions");
    exit(-1);
  }
  
  g = sm_get(m,n);
  
  /* get entries */
  for ( u_int i = 0 ; i < m; i++ )
  {
    if ( fscanf(fp," row %u:",&row) < 1 )
    {
      fprintf(stderr, "sm_input error reading line %u", i);
      exit(-1);
    }
    n0 = &(g->rows[i]);
    for (;;)
    {
     if ( fscanf(fp,"%lf",&col) < 1 )
      {
        fprintf(stderr, "sm_input error reading line %u col x", i);
        exit(-1);
      }
      if (col < 0.0) break;
      n0 = sm_add(n0, (u_int)col, 1.0);
    }
  }
  
  return (g);
}

static char *format = "%1.5g ";

void sm_output(FILE *fp, SMAT *G)
{
  NODE *n0;

  fprintf(fp, "SparseMatrix: %d by %d\n", G->m, G->n);
  for ( u_int i = 0 ; i < G->m ; i++ )
  {
    fprintf(fp, "row %u: ", i); 
    n0 = &(G->rows[i]);
    while (n0->val >= 0.0)
    {
      fprintf(fp,format,(Real)n0->col);
      n0 = n0->next;
    }
    fprintf(fp,"-1\n");
  }
}

/* m_output -- file output of matrix 
   Precondition: Memory already allocated for the matrix */
void m_output(FILE *fp, MAT *g)
{
   fprintf(fp, "Matrix: %d by %d\n", g->m, g->n);
   for ( u_int i = 0 ; i < g->m ; i++ )
   {
     fprintf(fp,"row %u: ", i);
     for ( u_int j = 0 ; j < g->n ; j++ )
     {
       fprintf(fp,format,g->me[i][j]);
     }
     putc('\n',fp);
   }
}

/* v_output -- file output of vector */
void v_output(FILE *fp, VEC *v)
{
  fprintf(fp, "Vector: %d\n", v->dim);
  for (u_int i=0 ; i<v->dim ; i++) fprintf(fp,format,v->ve[i]);
  putc('\n',fp);
}

/* m_cp -- copy matrix M in OUT
   Precondition: memory is already allocated for M and OUT
   Precondition: sizes of M and OUT must match*/
void m_cp(MAT *M, MAT *OUT)
{
  for ( u_int i = 0; i < M->m; i++ )
    memmove(&(OUT->me[i][0]), &(M->me[i][0]), (M->n)*sizeof(Real));
}

/* v_cp -- copy vector v in out
   Precondition: memory is already allocated for v and out*/
void v_cp(VEC *v, VEC *out)
{
  memmove(&(out->ve[0]), &(v->ve[0]), (v->dim)*sizeof(Real));
}

/* v_zero -- zero a vector
   Precondition: memory is already allocated for v */
void v_zero(VEC *v)
{
  memset(v->ve,0,v->dim*sizeof(Real));
}

void m_zero(MAT *g)
{
	for(u_int i = 0; i < g->m; ++i) 
	{
		for(u_int j = 0; j < g->n; ++j) 
		{
				g->me[i][j] = 0;
		}
	}
}

MAT* mToH(const MAT *g)
{
	MAT *new = m_get(g->m, g->n);
	for(u_int i = 0; i < g->m; ++i) 
	{
		u_int countNotNull = 0;
		for(u_int j = 0; j < g->n; ++j) 
		{
			if (g->me[i][j] != 0) countNotNull++; 
			new->me[i][j] = g->me[i][j];
		}
		if (countNotNull != 0) {
			for(u_int j = 0; j < g->n; ++j) 
			{
				new->me[i][j] /= countNotNull ; 
			}
		}
		
	}
	
	return new;
}

VEC* multiplyMToV (MAT *g, VEC *v)
{
	VEC *res = v_get(v->dim);
	v_zero(res);
	for (u_int j = 0; j < g->n; ++j)
	{
		for (u_int i = 0; i < g->m; ++i)
		{
			res->ve[j] += v->ve[i]*g->me[i][j];
		}
	}
	return res;	
}

MAT* multiplyMToM (MAT *m1, MAT* m2)
{
	MAT* res = m_get(m1->m, m1->n);
	
	for(u_int i = 0; i < res->m; ++i) 
	{
		for(u_int j = 0; j < m2->n; ++j) 
		{
			for(u_int l = 0; l < m2->n; ++l) 
			{
			res->me[i][j] += m1->me[i][l]*m2->me[l][j]; 
			}
		}
	}
	return res;
}

MAT* mToPowerN (MAT *g, u_int n)
{
	MAT* res = m_get(g->m, g->n);
	for (u_int i = 0; i < g->m; ++i)
	{
		for (u_int j = 0; j < g->n; ++j)
		{
			if(i==j)
			{
				res->me[i][j]=1;
			}
			else
			{
				
			res->me[i][j]=0;
			}
		}
	}
	
	for (u_int i = 0; i < n; ++i)
	{
		MAT* temp=res;
		res = multiplyMToM(res, g);
		m_free(temp);
	}
	return res;
}
void removeAbs(MAT *g)
{
	for (u_int i = 0; i < g->m; ++i)
	{
		u_int j ;
		for (j= 0; j < g->n && g->me[i][j]==0; ++j)
		{
			
		}
		if(j==g->n)
		{
			for (j= 0; j < g->n ; ++j)
			{
				g->me[i][j]=1.0;
			}
		}
	}
}

VEC * returnA(MAT *g)
{
	VEC *a=v_get(g->m);
	for (u_int i = 0; i < g->m; ++i)
	{
		u_int j ;
		for (j= 0; j < g->n && g->me[i][j]==0; ++j)
		{
			
		}
		if(j==g->n)
		{
			a->ve[i]=1;
		}
	}
	return a;
}

void multiplyMToConst (MAT* g, Real a) 
{
	for (u_int i = 0; i < g->m; ++i)
	{
		for (u_int j = 0; j < g->n; ++j)
		{
			g->me[i][j] *= a;
		} 
	}
}

MAT* matriceTelep (u_int n) 
{
	MAT* res = m_get(n,n);
	for (u_int i = 0; i < n; ++i)
	{
		for (u_int j= 0; j < n; ++j)
		{
			res->me[i][j] = 1.0/n;
		}
	}
	return res;
}

void main()
{
	Real alfa = 0.85;
  VEC *v;
  v = v_get(3);
  v_output(stdout,v);
  printf("size: %d\n", v->dim);
  v_free(v);

  MAT *G;
  G = m_get(3,3);
  printf("nb-line: %d ; nb-col: %d\n", G->m, G->n);
  m_free(G);

  FILE *fp;
 
  fp = fopen("dataset/g.dat", "r");
  G = m_input(fp);
  fclose(fp);
  m_output(stdout,G);
  
  v=v_get(G->m);
 /* typedef struct
{
  u_int dim;
  Real *ve;
} VEC;*/

	for(u_int i=0;i<v->dim;i++)
	{
		v->ve[i]=1.0/v->dim;
	}
	
	v_output(stdout, v);
	
	// noeuds absorbants
	VEC * a = returnA(G);
	//removeAbs(G);
	
	//stocausticite
	MAT *H = mToH(G);
	
	//ergodicite
	
	MAT * telep = matriceTelep(H->m);
	VEC * res=v_get(H->m);
	v_cp(v,res);
	
	for (u_int n = 0; n < 1000; ++n)
	{
		 v_output(stdout, res);
		for (u_int i = 0; i < H->m; ++i)
		{
			
			
			Real totalTemp=0;
			Real totalTempA=0;
			for (u_int j = 0; j < H->n; ++j)
			{
				totalTemp += res->ve[j] * H->me[i][j];
				totalTempA += res->ve[j] * a->ve[j];
			}
			res->ve[i]=alfa*totalTemp+(((totalTempA)*alfa)+1-alfa)/H->m;
			
		}
		
	}
	
	/*for (u_int i = 0; i < H->m; ++i)
	{
		for (u_int j = 0; j < H->n; ++j)
		{
			
			H->me[i][j] = H->me[i][j] * (alfa) +  (1-alfa)/(H->m);
		}
	}*/
	//m_output(stdout, H);
	
	//convergence
	
	//H = mToPowerN(H, 2000);
	//m_output(stdout, H);
	//VEC *resultat = multiplyMToV(H, v);
	//v_output(stdout, resultat);
  
    v_output(stdout, res);

  
  //MAT *Hcarre = multiplyMToM(H, H);
	//m_output(stdout, Hcarre);
    
  //MAT *HPow = mToPowerN(H, 200);

	//m_output(stdout, HPow);
  
  //v_free(resultat);
  //v_free(v);
  //m_free(H);
  //m_free(G);
  //m_free(HPow);


  SMAT *SG;
  fp = fopen("dataset/genetic.dat","r");
  SG = sm_input(fp);
  sm_output(stdout, SG);
  
  VEC * aBis=v_get(SG->m);
  for(u_int i=0;i<SG->m;i++)
  {
	  NODE * cur=&(SG->rows[i]);
	  u_int nb=0;
	  while(cur->col!=-1)
	  {
		  cur=cur->next;
		  nb++;
		}
		if(nb!=0)
		{
			  while(cur->col!=-1)
			{
				cur->val=1.0/nb;
			cur=cur->next;
			
			}
		}
		else
		{
			//TODO add à a
			aBis->ve[i]=1;
		}
	}
  	for (u_int n = 0; n < 1000; ++n)
	{
		
		for (u_int i = 0; i < H->m; ++i)
		{
		}
	}
  /*typedef struct graphnode
{
  u_int col;
  Real val;
  struct graphnode *next;
} NODE;

typedef struct
{
  u_int m, n;
  NODE *rows;
} SMAT;
*/
  
  
  fclose(fp);
  
  
  
  
  
  
  fp = fopen("dataset/test.dat","w");
  sm_output(fp,SG);
  sm_free(SG);
  
  
  
  
  

  exit(0);
}
