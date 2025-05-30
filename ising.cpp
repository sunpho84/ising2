#include <random>
#include <vector>
#include <cstdio>
#include <chrono>
#include <omp.h>

using namespace std;

/** For reference, critical beta is .4407228 */
const double Beta=0.29;
size_t L;
size_t N;
const int seed=124634;
 int nConfs;

struct timer
{
  static timer timerCost;
  
  static auto now()
  {
    return chrono::high_resolution_clock::now();
  }
  
  size_t tot=0;
  
  size_t n=0;
  
  chrono::high_resolution_clock::time_point from;
  
  void start()
  {
    n++;
    from=now();
  }
  
  void stop()
  {
    tot+=chrono::duration_cast<chrono::nanoseconds>(now()-from).count();
  }
  
  double get(bool sub=true)
  {
    double res=tot;
    if(sub)
      res-=timerCost.tot*(double)n/timerCost.n;
    
    return res/1e9;
  }
};

timer timer::timerCost;

/** Computes the energy relative a given site */
int computeEnSite(const vector<int>& conf,
		  const size_t& iSite)
{
  int en=0;
  
  /** Coverts to coordinate */
  const int y=iSite/L;
  const int x=iSite%L;
  
  /** Computes the neighbor in the two directions */
  const int neighSiteXUp=y*L+(x+1)%L;
  const int neighSiteYUp=((y+1)%L)*L+x;
  const int neighSiteXDw=y*L+(x+L-1)%L;
  const int neighSiteYDw=((y+L-1)%L)*L+x;
  
  en-=conf[neighSiteXUp]*conf[iSite]+conf[neighSiteYUp]*conf[iSite];
  en-=conf[neighSiteXDw]*conf[iSite]+conf[neighSiteYDw]*conf[iSite];
  
  return en;
}

/** Computes the energy */
int computeEn(const vector<int>& conf)
{
  int en=0;
  
#pragma omp parallel for reduction(+:en)
  for(size_t iSite=0;iSite<N;iSite++)
    {
      /** Coverts to coordinate */
      const int y=iSite/L;
      const int x=iSite%L;
      
      /** Computes the neighbor in the two directions */
      const int neighSiteX=y*L+(x+1)%L;
      const int neighSiteY=((y+1)%L)*L+x;
      
      en-=conf[neighSiteX]*conf[iSite]+conf[neighSiteY]*conf[iSite];
    }
  
  return en;
}

/** Computes the magnetization */
double computeMagnetization(const vector<int>& conf)
{
  int mag=0;
  
#pragma omp parallel for reduction(+:mag)
  for(size_t iSite=0;iSite<N;iSite++)
    mag+=conf[iSite];
  
  return (double)mag/N;
}

int main(int narg,char** arg)
{
  if(narg<3)
    {
      fprintf(stderr,"Use %s L nConfs\n",arg[0]);
      exit(1);
    }
  
  L=atoi(arg[1]);
  nConfs=atoi(arg[2]);
  N=L*L;
  
  timer totalTime;
  
  for(size_t i=0;i<100000;i++)
    {
      timer::timerCost.start();
      timer::timerCost.stop();
    }
  
 printf("NThreads: %d\n",omp_get_max_threads());
  
#ifdef PLOT
  /** Open the plot */
  FILE* gp=popen("gnuplot","w");
  fprintf(gp,"unset key\n");
  fprintf(gp,"set style fill solid\n");
#endif
  /** Open the measurement file*/
  FILE* measFile=fopen("meas.txt","w");
  
  /** Configuration */
  vector<int> conf(N);
  
  /** Random number generator */
  // mt19937_64 gen(seed);
  
  /** Local random number generator */
  std::vector<mt19937_64> locGen(N);
  for(size_t iSite=0;iSite<N;iSite++)
    locGen[iSite].seed(seed+iSite);
  
  /** Creates the first configuration */
  for(size_t i=0;i<N;i++)
    conf[i]=binomial_distribution<int>(1,0.5)(locGen[i])*2-1;
  
  totalTime.start();
  
  /** Produce nConfs */
  for(int iConf=0;iConf<nConfs;iConf++)
    {
      /** Updates all sites */
      for(int par=0;par<2;par++)
	{
#pragma omp parallel for
	  for(size_t iSite=0;iSite<N;iSite++)
	    {
	      const int y=iSite/L;
	      const int x=iSite%L;
	      const int pSite=(x+y)%2;

	      if(par==pSite)
		{
		  /** Store the initial state of the site */
		  const int oldVal=conf[iSite];
		  
		  /** Computes energy of the old conf, relative to the site*/
		  const int enBefore=computeEnSite(conf,iSite);
		  
		  /** Draw a fair binomial distribitution, 0 or 1 with probability 50%, and assigns to the site */
		  conf[iSite]=binomial_distribution(1,0.5)(locGen[iSite])*2-1;
		  
		  /** Computes energy of the new conf, rleative to the site */
		  const int enAfter=computeEnSite(conf,iSite);
		  
		  /** Computes energy difference */
		  const int eDiff=enAfter-enBefore;
		  
		  /** Computes the acceptance probability */
		  const double pAcc=std::min(1.0,exp(-Beta*eDiff));
		  
		  /** Draw acceptance from the binomial distribution with probability pAcc */
		  const int acc=binomial_distribution<int>(1,pAcc)(locGen[iSite]);
		  
		  /** If not accepted, bring back the original site*/
		  if(not acc)
		    conf[iSite]=oldVal;
		}
	    }
	}
      
#ifdef PLOT
      /** Plot the configuration */
      fprintf(gp,"plot '-' w boxxyerror\n");
      for(size_t site=0;site<N;site++)
	if(conf[site]==-1)
	  fprintf(gp,"%lg %lg 0.5 0.5\n",site%L+0.5,int(site/L)+0.5);
      fprintf(gp,"e\n");
      fflush(gp);
#endif
      
      /** Computes the magnetization */
      const double mag=computeMagnetization(conf);
      
      /** Computes the energy */
      const double ene=computeEn(conf);
      
      /** Print the measurement */
      fprintf(measFile,"%lg %lg\n",mag,ene);
      fflush(measFile);
      
      /** Report progress */
      printf("Progress: %d/%d\n",iConf,nConfs);
    }
  
#ifdef PLOT
  pclose(gp);
#endif
  
  fclose(measFile);
  
  totalTime.stop();
  
  printf("Duration: %lg s\n",totalTime.get());
  
  return 0;
}
