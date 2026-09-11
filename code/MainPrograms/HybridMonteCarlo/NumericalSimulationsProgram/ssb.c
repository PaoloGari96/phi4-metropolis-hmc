
/* 
 *   File ssb.c
 */ 

#define CONTROL
#include "phi4.h"
#include "metropolis.h"
#include "math.h"
#include "observables.h"

/*  
 *  data structures to store all the parameters of the algorithm
 */

static metro_params_t metro_params;
static hmc_params_t hmc_params;
static act_params_t act_params;
static leapfrog_params_t leapfrog_params;
static int seed;


static int get_val(FILE* fp, char *str, char* fmt,  void* val)
{
   char c[128];

   if(1!=fscanf(fp,"%s",c))
   {
      fprintf(stderr,"Error reading input file at %s\n",str);
      exit(1);
   }

   if(strcmp(str,c)!=0)
   {
      fprintf(stderr,"Error reading input file expected %s found %s\n",str,c);
      exit(1);
   }

   if(1!=fscanf(fp,fmt,val))
   {
      fprintf(stderr,"Error reading input file at %s\n",str);
      fprintf(stderr,"Cannot read value format %s\n",fmt);
      exit(1);
   }

   return 0;
}


static int read_input(char *input)
{
   FILE* fp;

   fp=fopen(input,"r");

   if (fp==NULL) 
   {
      fprintf(stderr, "Cannot open input file %s \n",input);
      exit(1);
   }

   get_val(fp, "kappa",       "%lf",&act_params.kappa    );
   get_val(fp, "lambda",      "%lf",&act_params.lambda   );
   get_val(fp, "ntherm",      "%i" ,&metro_params.ntherm );
   get_val(fp, "nsweep",      "%i" ,&metro_params.nsweep );
   get_val(fp, "ntraj",       "%i" ,&hmc_params.ntraj );
   get_val(fp, "N",           "%i" ,&leapfrog_params.N );
   get_val(fp, "delta",       "%lf",&metro_params.delta  );
   get_val(fp, "seed",        "%i" ,&seed                );
   get_val(fp, "naccu",       "%i" ,&metro_params.naccu  );


   return 0;
}


int main(int argc, char* argv[])
{
   int i,j;
   double tzero=1;
   double m,sm,qm;
   double h;
   int N;
   double dt;

   if (argc != 2) 
   {
	   fprintf(stderr,"Number of arguments not correct\n");
	   fprintf(stderr,"Usage: %s <infile> \n",argv[0]);
	   exit(1);
   }

   /* get the parameters from the input file */
   read_input(argv[1]);

   /* initialize random number generator */
   rlxd_init(2,seed);

   /* initialize the nearest neighbor field */
   hopping(hop);

   /* initialize phi field */
   ranlxd(phi,V);


   N = leapfrog_params.N;
   dt=tzero/N;

   for(i=0;i<metro_params.ntherm;i++){
      hmc(&act_params,&metro_params,dt,N);
   }


   for(i=0;i<hmc_params.ntraj;i++){
      m=0;
      sm=0;
      qm=0;

      h=0;
      for(j=0;j<metro_params.naccu;j++){
         m+=Mabs();
         sm+=Mtwo();
         qm+=Mfour();    

         h+=hmc(&act_params,&metro_params,dt,N);
      }
      printf("%d %e %e %e %e \n", i, m/metro_params.naccu, sm/metro_params.naccu, qm/metro_params.naccu, h/metro_params.naccu);
   }

 
   return 0;
}

