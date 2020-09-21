#include <vector>
#include <random>
#include <stdlib.h>
#include <algorithm>
#include "agents.h"
#include "Parameters.h"

using namespace std;

void initialize();
void simulationStep(int step,int infectSpread,int numInfectRepeat,float oxyHeal,int numRecurInj, int numABX);
void cellStep(int infectSpread,int numInfectRepeat,float oxyHeal,int numRecurInj);
void recurrentInjury(int step,int numRecurInj);
void giveABX(int step, int *numABX);

extern void injure_infectionFRD(int inj_number);
extern void updateSystemOxy(int step);
extern void evaporate();
extern void diffuse();
extern void recur_injury();
extern void applyAntibiotics();

vector<EC> ecArray;
vector<int> ecIndexes;
vector<pmn> pmnArray;
vector<mono> monoArray;
vector<TH0> TH0array;
vector<TH1> TH1array;
vector<TH2> TH2array;
vector<pmn_marrow> pmn_marrowArray;
vector<mono_marrow> mono_marrowArray;
vector<TH0_germ> TH0_germArray;
vector<TH1_germ> TH1_germArray;
vector<TH2_germ> TH2_germArray;
mt19937 generator;
uniform_int_distribution<int> distribution10k(0,9999);
uniform_int_distribution<int> distribution1000(0,999);
uniform_int_distribution<int> distribution100(0,99);
uniform_int_distribution<int> distribution50(0,49);
uniform_int_distribution<int> distribution12(0,11);
uniform_int_distribution<int> distribution10(0,9);
uniform_int_distribution<int> distribution9(0,8);
uniform_int_distribution<int> distribution8(0,7);
uniform_int_distribution<int> distribution5(0,4);
uniform_int_distribution<int> distribution3(0,2);
uniform_int_distribution<int> distribution2(0,1);

float d_oxy[numTimeSteps],d_infect[numTimeSteps],d_TNF[numTimeSteps],d_sTNFr[numTimeSteps],
d_IL10[numTimeSteps],d_IL6[numTimeSteps],d_GCSF[numTimeSteps],d_pTH1[numTimeSteps],
d_pTH2[numTimeSteps],d_IFNg[numTimeSteps],d_PAF[numTimeSteps],d_IL1[numTimeSteps],
d_IL4[numTimeSteps],d_IL8[numTimeSteps],d_IL12[numTimeSteps],d_IL1ra[numTimeSteps],
d_sIL1r[numTimeSteps];
float system_oxy,oxyDeficit,totalInfection,total_TNF,total_sTNFr,total_IL10,
total_IL6,total_GCSF,total_proTH1,total_proTH2,total_IFNg,total_PAF,
total_IL1,total_IL4,total_IL8,total_IL12,total_sIL1r,total_IL1ra;
float PAFmult,TNFmult,sTNFrmult,IL1mult,sIL1rmult,IL1ramult,IFNgmult,IL4mult,
IL8mult,IL10mult,IL12mult,GCSFmult;
int d_pmn[numTimeSteps],d_mono[numTimeSteps],d_TH1[numTimeSteps],d_TH2[numTimeSteps];
int cellGrid[101][101];

extern const int cellCapacity,xDim,yDim,injuryStep,parameterInput,numTimeSteps;
extern const float antibioticMultiplier;

extern "C" int mainSimulation(float oxyHeal, int infectSpread,
	int numRecurInj, int numInfectRepeat, int numCytokines, float* internalParameterization){

		int i,step,iend,dOut,count6hr,antibiotic1,antibiotic2,istep,x1,x2,inj_number;
		//  int* numABX;
		int numABX;
		int INJ[2]={1,40};
		//cout<<"inj_number="<<inj_number<<"oxy="<<oxyHeal<<"IS="<<infectSpread<<"\n";
		int sum[2]={0,0};
		int flag;

		for(x1=0;x1<2;x1++){
			inj_number=INJ[x1];
			for(x2=0;x2<25;x2++){
				flag=0;
				generator.seed(x2);
				initialize();
				PAFmult=internalParameterization[0];
				TNFmult=internalParameterization[1];
				IL1mult=internalParameterization[2];
				IFNgmult=internalParameterization[3];
				IL4mult=internalParameterization[4];
				IL8mult=internalParameterization[5];
				IL10mult=internalParameterization[6];
				IL12mult=internalParameterization[7];
				GCSFmult=internalParameterization[8];
				dOut=0;
				step=0;
				istep=0;
				ecIndexes.clear();
				for(i=0;i<xDim*yDim;i++){
					ecIndexes.push_back(i);
				}
				// numAntibiotics=0;
				// numABX=&numAntibiotics;
				count6hr=0;
				antibiotic1=0;
				antibiotic2=0;
				numABX=0;
				injure_infectionFRD(inj_number);
				for(i=0;i<numTimeSteps;i++){
					if(step==0){}
					step++;
					istep++;
					if(step==injuryStep){step=1;}
					antibiotic1++;
					antibiotic2++;

					updateSystemOxy(istep);
					if(oxyDeficit>8161||(oxyDeficit<5&&i>0)){
						if((writeEverything==1)||(dailyOutput==1)){
							for(iend=dOut;iend<numTimeSteps;iend++){
								d_oxy[iend]=-1.0;
								d_infect[iend]=-1.0;
								d_TNF[iend]=-1.0;
								d_sTNFr[iend]=-1.0;
								d_IL10[iend]=-1.0;
								d_GCSF[iend]=-1.0;
								d_pTH1[iend]=-1.0;
								d_pTH2[iend]=-1.0;
								d_pmn[iend]=-1;
								d_mono[iend]=-1;
								d_TH1[iend]=-1;
								d_TH2[iend]=-1;
								d_IFNg[iend]=-1.0;
								d_PAF[iend]=-1.0;
								d_IL1[iend]=-1.0;
								d_IL4[iend]=-1.0;
								d_IL8[iend]=-1.0;
								d_IL12[iend]=-1.0;
								d_sIL1r[iend]=-1.0;
								d_IL1ra[iend]=-1.0;

							}
						}
						if(oxyDeficit>=8161){
//							return 0;
								break;
						}
						else{
							sum[x1]++;
//							cout<<"TEST5"<<"\n";
							flag=1;
							break;
						}
					}

					if(((dailyOutput==1)&&(count6hr==50))||(writeEverything==1)){
						d_oxy[dOut]=oxyDeficit;
						d_infect[dOut]=totalInfection;
						d_TNF[dOut]=total_TNF;
						d_sTNFr[dOut]=total_sTNFr;
						d_IL10[dOut]=total_IL10;
						d_GCSF[dOut]=total_GCSF;
						d_pTH1[dOut]=total_proTH1;
						d_pTH2[dOut]=total_proTH2;
						d_pmn[dOut]=pmnArray.size();
						d_mono[dOut]=monoArray.size();
						d_TH1[dOut]=TH1array.size();
						d_TH2[dOut]=TH2array.size();
						d_IFNg[dOut]=total_IFNg;
						d_PAF[dOut]=total_PAF;
						d_IL1[dOut]=total_IL1;
						d_IL4[dOut]=total_IL4;
						d_IL8[dOut]=total_IL8;
						d_IL12[dOut]=total_IL12;
						d_sIL1r[dOut]=total_sIL1r;
						d_IL1ra[dOut]=total_IL1ra;
						dOut++;
					}
					simulationStep(i,infectSpread,numInfectRepeat,oxyHeal,numRecurInj,numABX);

				}
				if(oxyDeficit<8161&&flag==0){
//					cout<<"TEST10"<<"\n";
					sum[x1]++;
				}
			}
		}
//		cout<<"sum1="<<sum[0]<<" sum2="<<sum[1]<<"\n";
		if(sum[0]==sum[1]){
			return 0;
		}
		else{
			return 1;
		}
	}

void simulationStep(int step,int infectSpread,int numInfectRepeat,float oxyHeal,int numRecurInj, int numABX){
  cellStep(infectSpread,numInfectRepeat,oxyHeal,numRecurInj);
  evaporate();
  recurrentInjury(step,numRecurInj);
  diffuse();
  if(antibioticMultiplier>0){
    giveABX(step, &numABX);
  }
}

void cellStep(int infectSpread,int numInfectRepeat,float oxyHeal,int numRecurInj){
  int length,j;
  length=TH0array.size();
  if(length>0){
    shuffle(TH0array.begin(),TH0array.end(),generator);}
  j=0;
  while(j<length){
    TH0array[j].TH0function(j);
    j++;
    length=TH0array.size();}
  length=ecArray.size();
  shuffle(ecIndexes.begin(),ecIndexes.end(),generator);
  j=0;
  while(j<length){
    ecArray[ecIndexes[j]].inj_function(infectSpread,numInfectRepeat);
    ecArray[ecIndexes[j]].ECfunction(oxyHeal);
    j++;
    length=ecArray.size();}
  length=pmnArray.size();
  if(length>0){
    shuffle(pmnArray.begin(),pmnArray.end(),generator);}
  j=0;
  while(j<length){
    pmnArray[j].pmn_function(j);
    j++;
    length=pmnArray.size();}
  length=monoArray.size();
    if(length>0){
      shuffle(monoArray.begin(),monoArray.end(),generator);}
  j=0;
  while(j<length){
    monoArray[j].mono_function(j);
    j++;
    length=monoArray.size();}
  length=TH1array.size();
    if(length>0){
      shuffle(TH1array.begin(),TH1array.end(),generator);}
  j=0;
  while(j<length){
    TH1array[j].TH1function(j);
    j++;
    length=TH1array.size();}
  length=TH2array.size();
    if(length>0){
      shuffle(TH2array.begin(),TH2array.end(),generator);}
  j=0;
  while(j<length){
    TH2array[j].TH2function(j);
    j++;
    length=TH2array.size();}
  length=pmn_marrowArray.size();
    if(length>0){
      shuffle(pmn_marrowArray.begin(),pmn_marrowArray.end(),generator);}
  j=0;
  while(j<length){
    pmn_marrowArray[j].pmn_marrow_function();
    j++;
    length=pmn_marrowArray.size();}

  length=mono_marrowArray.size();
    if(length>0){
      shuffle(mono_marrowArray.begin(),mono_marrowArray.end(),generator);}
  j=0;
  while(j<length){
    mono_marrowArray[j].mono_marrow_function();
    j++;
    length=mono_marrowArray.size();}

  length=TH1_germArray.size();
    if(length>0){
      shuffle(TH1_germArray.begin(),TH1_germArray.end(),generator);}
  j=0;
  while(j<length){
    TH1_germArray[j].TH1_germ_function();
    j++;
    length=TH1_germArray.size();}

  length=TH2_germArray.size();
    if(length>0){
      shuffle(TH2_germArray.begin(),TH2_germArray.end(),generator);}
  j=0;
  while(j<length){
    TH2_germArray[j].TH2_germ_function();
    j++;
    length=TH2_germArray.size();}

  length=TH0_germArray.size();
    if(length>0){
      shuffle(TH0_germArray.begin(),TH0_germArray.end(),generator);}
  j=0;
  while(j<length){
    TH0_germArray[j].TH0_germ_function();
    j++;
    length=TH0_germArray.size();
  }
}

void recurrentInjury(int step, int numRecurInj){
  int i;
  if(step==injuryStep-1){
    for(i=1;i<=numRecurInj;i++){
      recur_injury();
    }
  }
}

void giveABX(int step, int *numABX){
  if((step%injuryStep==102)&&(*numABX<1100)){
    applyAntibiotics();
    *numABX++;
  }
}
