#include <random>
#include <unistd.h>
#include <vector>
#include <iostream>

using namespace std;

int computeEn(vector<int>& conf,int L,int N)
{
  int en=0;
  for(int iSite=0;iSite<N;iSite++)
    {
      int y=iSite/L;
      int x=iSite%L;
      
      int neighSiteX=y*L+(x+1)%L;
      int neighSiteY=((y+1)%L)*L+x;
      
      en-=conf[neighSiteX]*conf[iSite]+conf[neighSiteY]*conf[iSite];
    }
  
  return en;
}

double computeMagnetization(vector<int>& conf,int L,int N)
{
  int mag=0;
  for(int iSite=0;iSite<N;iSite++)
    mag+=conf[iSite];
  
  return (double)mag/N;
}

int main()
{
  FILE* gp=popen("gnuplot","w");
  fprintf(gp,"unset key\n");
  fprintf(gp,"set style fill solid\n");
  
  double beta=0.48;
  int L=100;
  int N=L*L;
  int seed=23534634;
  
  vector<int> conf(N);
  
  mt19937 gen(seed);
  
  for(int& c : conf)
    c=+1;
  
  int nConfs=100;
  
  /** Produce nConfs */
  for(int iConf=0;iConf<nConfs;iConf++)
    {
      /** Update each esite*/
      for(int iSite=0;iSite<N;iSite++)
	{
	  // cout<<"Looping on site "<<iSite<<endl; 
	  
	  int backupSiteState=conf[iSite];
	  
	  // cout<<"Before: "<<conf[iSite]<<endl;
	  int enBefore=computeEn(conf,L,N);
	  // cout<<"enBefore: "<<enBefore<<endl;
	  binomial_distribution siteDistr(1,0.5);
	  if(siteDistr(gen)==0)
	    conf[iSite]=-1;
	  else
	    conf[iSite]=+1;
	  
	  // cout<<"After: "<<conf[iSite]<<endl;
	  int enAfter=computeEn(conf,L,N);
	  // cout<<"enAfter: "<<enAfter<<endl;
	  
	  int eDiff=enAfter-enBefore;
	  // cout<<"eDiff: "<<eDiff<<endl;
	  
	  if(eDiff<=0)
	    // cout<<"Accepted as energy is decreasing"<<endl
	    ;
	  else
	    {
	      double pAcc=exp(-beta*eDiff);
	      
	      // cout<<"Pacc: "<<pAcc<<endl;
	      binomial_distribution<int> distrAcc(1,pAcc);
	      
	      int acc=distrAcc(gen);
	      // cout<<"acc: "<<acc<<endl;
	      
	      if(acc==0)
		{
		  conf[iSite]=backupSiteState;
		  // cout<<"Not accepted"<<endl;
		}
	      // else
	      // 	cout<<"Accepted"<<endl;
	    }
	}
      
      fprintf(gp,"plot '-' w boxxyerror\n");
      for(int site=0;site<N;site++)
	if(conf[site]==-1)
	  fprintf(gp,"%lg %lg 0.5 0.5\n",site%L+0.5,int(site/L)+0.5);
      fflush(gp);
      
      double mag=computeMagnetization(conf,L,N);
      cout<<"Mag "<<mag<<endl;
      
      //sleep(1);
    }
  
  pclose(gp);
  
  return 0;
}
