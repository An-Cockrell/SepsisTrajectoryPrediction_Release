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
void updateTrajectoryOutput(float allSignals[][numTimeSteps], int q);

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

float allSignals[20][numTimeSteps];
float allSignalsReturn[20*numTimeSteps];
float system_oxy,oxyDeficit,totalInfection,total_TNF,total_sTNFr,total_IL10,
total_IL6,total_GCSF,total_proTH1,total_proTH2,total_IFNg,total_PAF,
total_IL1,total_IL4,total_IL8,total_IL12,total_sIL1r,total_IL1ra;
float PAFmult,TNFmult,sTNFrmult,IL1mult,sIL1rmult,IL1ramult,IFNgmult,IL4mult,
IL8mult,IL10mult,IL12mult,GCSFmult;
int cellGrid[101][101];

extern const int cellCapacity,xDim,yDim,injuryStep,parameterInput,numTimeSteps;
extern const float antibioticMultiplier;

extern "C" float* mainSimulation(float oxyHeal, int infectSpread,
				 int numRecurInj, int numInfectRepeat, int inj_number, int seed, int numCytokines, int seed2, int seed3, int seed4, int seed5,int seed6, int seed7, int seed8, int seed9, int seed10, float* internalParameterization){

		int i,step,iend,jend,dOut,count6hr,antibiotic1,antibiotic2,istep,x1,x2,k,j;
		int numABX;
		int flag;
		flag=0;
		generator.seed(seed);
//		cout<<"Test10\n";
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
		count6hr=0;
		antibiotic1=0;
		antibiotic2=0;
		numABX=0;
		injure_infectionFRD(inj_number);
		updateSystemOxy(istep);
		for(i=0;i<numTimeSteps;i++){
			step++;
			istep++;
			if(step==injuryStep){step=1;}
			if(istep==1650 && seed2!=-1){generator.seed(seed2);}
			if(istep==1850 && seed3!=-1){generator.seed(seed3);}
			if(istep==1950 && seed4!=-1){generator.seed(seed4);}
			if(istep==2050 && seed5!=-1){generator.seed(seed5);}
			if(istep==2150 && seed6!=-1){generator.seed(seed6);}
			if(istep==2250 && seed7!=-1){generator.seed(seed7);}
			if(istep==2350 && seed8!=-1){generator.seed(seed8);}
			if(istep==2450 && seed9!=-1){generator.seed(seed9);}
			if(istep==2550 && seed10!=-1){generator.seed(seed10);}
			antibiotic1++;
			antibiotic2++;
			updateTrajectoryOutput(allSignals,i);
			simulationStep(i,infectSpread,numInfectRepeat,oxyHeal,numRecurInj,numABX);
			updateSystemOxy(istep);
			if(oxyDeficit>8161||(oxyDeficit<5&&i>0)){
				if((writeEverything==1)||(dailyOutput==1)){
					for(iend=i+1;iend<numTimeSteps;iend++){
						for(jend=0;jend<20;jend++){
								allSignals[jend][iend]=-1.0;
						}
					}
					break;
				}
				if(oxyDeficit>=8161){
						break;
				}
			}
		}
		k=0;
		for(j=0;j<20;j++){
			for(i=0;i<numTimeSteps;i++){
				allSignalsReturn[k]=allSignals[j][i];
				k++;
			}
		}
		return allSignalsReturn;
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

// void updateTrajectoryOutput(float allSignals[], int q){
// 	allSignals[0*numTimeSteps+q]=oxyDeficit;
// 	allSignals[1*numTimeSteps+q]=totalInfection;
// 	allSignals[2*numTimeSteps+q]=total_TNF;
// 	allSignals[3*numTimeSteps+q]=total_sTNFr;
// 	allSignals[4*numTimeSteps+q]=total_IL10;
// 	allSignals[5*numTimeSteps+q]=total_GCSF;
// 	allSignals[6*numTimeSteps+q]=total_proTH1;
// 	allSignals[7*numTimeSteps+q]=total_proTH2;
// 	allSignals[8*numTimeSteps+q]=pmnArray.size();
// 	allSignals[9*numTimeSteps+q]=monoArray.size();
// 	allSignals[10*numTimeSteps+q]=TH1array.size();
// 	allSignals[11*numTimeSteps+q]=TH2array.size();
// 	allSignals[12*numTimeSteps+q]=total_IFNg;
// 	allSignals[13*numTimeSteps+q]=total_PAF;
// 	allSignals[14*numTimeSteps+q]=total_IL1;
// 	allSignals[15*numTimeSteps+q]=total_IL4;
// 	allSignals[16*numTimeSteps+q]=total_IL8;
// 	allSignals[17*numTimeSteps+q]=total_IL12;
// 	allSignals[18*numTimeSteps+q]=total_sIL1r;
// 	allSignals[19*numTimeSteps+q]=total_IL1ra;
// }

void updateTrajectoryOutput(float allSignals[][numTimeSteps], int q){
	allSignals[0][q]=oxyDeficit;
	allSignals[1][q]=totalInfection;
	allSignals[2][q]=total_TNF;
	allSignals[3][q]=total_sTNFr;
	allSignals[4][q]=total_IL10;
	allSignals[5][q]=total_GCSF;
	allSignals[6][q]=total_proTH1;
	allSignals[7][q]=total_proTH2;
	allSignals[8][q]=pmnArray.size();
	allSignals[9][q]=monoArray.size();
	allSignals[10][q]=TH1array.size();
	allSignals[11][q]=TH2array.size();
	allSignals[12][q]=total_IFNg;
	allSignals[13][q]=total_PAF;
	allSignals[14][q]=total_IL1;
	allSignals[15][q]=total_IL4;
	allSignals[16][q]=total_IL8;
	allSignals[17][q]=total_IL12;
	allSignals[18][q]=total_sIL1r;
	allSignals[19][q]=total_IL1ra;
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
