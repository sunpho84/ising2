#include <random>
#include <unistd.h>
#include <vector>
#include <iostream>

using namespace std;

#define PLOT

#include <chrono>

size_t totalTime=0;
size_t computeEnTime=0;

auto now()
{
  return chrono::high_resolution_clock::now();
}

size_t timeFrom(chrono::high_resolution_clock::time_point from)
{
  return chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now()-from).count();
}

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

int computeSiteEn(vector<int>& conf,int L,int N,int iSite)
{
  int en=0;
      int y=iSite/L;
      int x=iSite%L;
      
      int neighSiteX=y*L+(x+1)%L;
      int neighSiteX2=y*L+(x+L-1)%L;
      int neighSiteY=((y+1)%L)*L+x;
      int neighSiteY2=((y+L-1)%L)*L+x;
      
      en-=conf[neighSiteX]*conf[iSite]+conf[neighSiteY]*conf[iSite];
      en-=conf[neighSiteX2]*conf[iSite]+conf[neighSiteY2]*conf[iSite];
  
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
#ifdef PLOT
  FILE* gp=popen("gnuplot","w");
  fprintf(gp,"unset key\n");
  fprintf(gp,"set style fill solid\n");
#endif
  
  double beta=1.4407228;
  int L=20;
  int N=L*L;
  int seed=124634;
  
  vector<int> conf(N);
  
  mt19937 gen(seed);
  
  for(int& c : conf)
    c=binomial_distribution<int>(1,0.5)(gen)*2-1;
  
  int nConfs=10;
  
  auto beginProgTime=now();
  
  /** Produce nConfs */
  for(int iConf=0;iConf<nConfs;iConf++)
    {
      /** Update each esite*/
      for(int iSite=0;iSite<N;iSite++)
	{
	  // cout<<"Looping on site "<<iSite<<endl; 
	  
	  int backupSiteState=conf[iSite];
	  
	  // cout<<"Before: "<<conf[iSite]<<endl;
	  int enBefore=computeSiteEn(conf,L,N,iSite);
	  // cout<<"enBefore: "<<enBefore<<endl;
	  binomial_distribution siteDistr(1,0.5);
	  if(siteDistr(gen)==0)
	    conf[iSite]=-1;
	  else
	    conf[iSite]=+1;
	  
	  // cout<<"After: "<<conf[iSite]<<endl;
	  auto beginEnMeas=now();
	  int enAfter=computeSiteEn(conf,L,N,iSite);
	  // cout<<"enAfter: "<<enAfter<<endl;
	  computeEnTime+=timeFrom(beginEnMeas);
	  
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

#ifdef PLOT
      fprintf(gp,"plot '-' w boxxyerror\n");
      for(int site=0;site<N;site++)
	if(conf[site]==-1)
	  fprintf(gp,"%lg %lg 0.5 0.5\n",site%L+0.5,int(site/L)+0.5);
      fprintf(gp,"e\n");
      fflush(gp);
#endif
      
      double mag=computeMagnetization(conf,L,N);
      cout<<"Mag "<<mag<<endl;
      
      //sleep(1);
    }
  
  totalTime+=timeFrom(beginProgTime);
  
  cout<<"TotalTime: "<<totalTime/1e9<<" s"<<endl;
  cout<<"ComputeEnergyTime: "<<computeEnTime/1e9<<" s"<<endl;
  
#ifdef PLOT
  pclose(gp);
#endif
  
  return 0;
}
