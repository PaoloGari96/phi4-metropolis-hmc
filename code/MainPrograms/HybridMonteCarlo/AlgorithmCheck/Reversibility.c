#define CONTROL
#include "phi4.h"
#include "metropolis.h"
#include "math.h"
#include "observables.h"


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
   int i;
   double phi_old[V],mom_old[V];
   double h_old;
   double deltamom, deltaphi;
   double tzero=1;
   int N=10;
   double dt;

   deltamom=deltaphi=0;

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

   /* initialize mom field */
   ranlxd(mom,V);

   N = leapfrog_params.N; 
   dt = tzero/N;


   for(i=0;i<V;i++){
      phi_old[i]=phi[i];
      mom_old[i]=mom[i];
   }
   
   h_old=hamiltonian(&act_params);
   Leapfrog(&act_params,dt,N);

   for(i=0;i<V;i++){
      mom[i]=-mom[i];
   }

   Leapfrog(&act_params,dt,N);

   for(i=0;i<V;i++){
      mom[i]=-mom[i];
   }

   for(i=0;i<V;i++){
      deltaphi+=fabs(phi_old[i]-phi[i]);
      deltamom+=fabs(mom_old[i]-mom[i]);
   }

   deltaphi/=V;
   deltamom/=V;

   printf("%e %e %e\n", deltaphi, deltamom, hamiltonian(&act_params)-h_old);
 
   return 0;
}

