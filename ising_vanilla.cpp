#include <random>
#include <vector>
#include <cstdio>

using namespace std;

const double Beta=0.9;//.4407228;
const size_t L=200;
const size_t N=L*L;
const int seed=124634;
const int nConfs=1000;

/// Computes the energy
int computeEn(const vector<int>& conf)
{
  int en=0;
  
  for(size_t iSite=0;iSite<N;iSite++)
    {
      int y=iSite/L;
      int x=iSite%L;
      
      int neighSiteX=y*L+(x+1)%L;
      int neighSiteY=((y+1)%L)*L+x;
      
      en-=conf[neighSiteX]*conf[iSite]+conf[neighSiteY]*conf[iSite];
    }
  
  return en;
}

/// Computes the magnetization
double computeMagnetization(const vector<int>& conf)
{
  int mag=0;
  for(size_t iSite=0;iSite<N;iSite++)
    mag+=conf[iSite];
  
  return (double)mag/N;
}

int main()
{
  /** Open the plot */
  FILE* gp=popen("gnuplot","w");
  fprintf(gp,"unset key\n");
  fprintf(gp,"set style fill solid\n");
  
  /** Open the measurement file*/
  FILE* measFile=fopen("meas.txt","w");
  
  /** Configuration */
  vector<int> conf(N);
  
  /** Random number generator */
  mt19937_64 gen(seed);
  
  /** Creates the distribution */
  for(size_t i=0;i<N;i++)
    conf[i]=binomial_distribution<int>(1,0.5)(gen)*2-1;
  
  /** Produce nConfs */
  for(int iConf=0;iConf<nConfs;iConf++)
    {
      /** Updates all sites */
      for(size_t iSite=0;iSite<N;iSite++)
	{
	  /** Computes energy before the change */
	  const int enBefore=computeEn(conf);
	  
	  /** Draw a fair binomial distribitution, 0 or 1 with probability 50% */
	  binomial_distribution siteDistr(1,0.5);
	  if(siteDistr(gen)==0)
	    conf[iSite]=-1;
	  else
	    conf[iSite]=+1;
	  
	  /** Configuration */
	  vector<int> newConf=conf;
	  
	  /** Computes energy after the change */
	  const int enAfter=computeEn(newConf);
	  
	  /** Computes energy difference */
	  const int eDiff=enAfter-enBefore;
	  
	  /** Computes the acceptance probability */
	  const double pAcc=std::min(1.0,exp(-Beta*eDiff));
	  
	  /** Draw acceptance from the binomali distribution with probability pAcc */
	  const int acc=binomial_distribution<int>(1,pAcc)(gen);
	  
	  /** If not accepted, bring back the original site*/
	  if(acc)
	    newConf=conf;
	}
      
      /** Plot the configuration */
      fprintf(gp,"plot '-' w boxxyerror\n");
      for(size_t site=0;site<N;site++)
	if(conf[site]==-1)
	  fprintf(gp,"%lg %lg 0.5 0.5\n",site%L+0.5,int(site/L)+0.5);
      fprintf(gp,"e\n");
      fflush(gp);
      
      /** Computes the magnetization */
      const double mag=computeMagnetization(conf);
      
      /** Computes the energy */
      const double ene=computeEn(conf);
      
      /** Print the measurement */
      fprintf(measFile,"%lg %lg\n",mag,ene);
    }
  
  pclose(gp);
  fclose(measFile);
  
  return 0;
}
