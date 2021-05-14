// EC3,EC5,CA1 are implemented.
// in area1.cpp, input from CA3 is induced. (2016.3.29)
// random connection between pop. in CA1
// interarea connection induced.
// introduce input from CA3:  here, CA3 input is exchanged at Tach (correspodning to place cells representing before the junction and those after that.)
//  *4.cpp: last update 2016.6.5
//  intorudce EC3 dynaics 
//  *5.cpp: last update 2016.6.6
// introduce calcium current in EC5 neurons  (see Saravanan 2015)
//  *6.cpp: last update 2016.6.17
// EC3,5 neuron: SC -> pyr
// VIP introduced
//  *7.cpp: last update 2017.4.7
// VIP -> SOM is introduced.
// irrelevant place information is introduced.
// last update 2017.6.16
// VIP neurons are introduced as Poisson neurons
// PV receives inputs from CA3 not EC3
// delay time is changed  (in previous ver, delay is dependent on hstp. this ver, I fixed it)
// update 2017.6.23:  ver 9
// introduce connections between ext neurons in pre-layer to inh neurons in post layer.
// introduce calcium dynamics and CAN channel
// update 2017.7.6
// change Ach schedule.
// update 2017.7.14
// introduce the right arm
// update 2017.7.15
// VIP oscillaiton is introduced
// check point  gahp in EC3 is modified 0.6 -> 0.3 in this ver.
// check point  G3_CA1_EC5 in this ver.
// update 2017.10.30
// induce decoder 

#include <iostream>
#include <vector>
#include <string>
#include <deque>
#include <set>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include "mt.c"

#define VER_PROG 9
//#define readmode 0

using namespace std;
double const pi = 3.14159265;
double const e = 2.71828182;
double const hstp = 0.02;
double const dbl_eps=1.0e-10;

int osc_flg=0;   // osc_flg==1:  make oscillation input flat


// for longer ver.
double const T=7000.0;
double const Tach_strt =1500;
double const Tach_end_str = 3500;  // strong Ach application is terminated at this time
double const Tach_inc = 5000;
double const Tach_end = 7000;  
double       Tec3off_strt =15000;
double const Tec3off_end  = Tec3off_strt+15000;
double const Tca3_C  =  500;  //300
double const Tca3_L  =  1500;  //800
double const Tca3_o  =  3500; //2300
double const Tca3_C1  = 5000; //3300
double const Tca3_end  = 7000;
double const Tcnct_strt=200;   // in beginning, inter-layer connection is quenched.




double geta=0;
double Gpv_ach = 0.015;
double Icnst_PV=0.0;


//check point
double const tmout=50;





// ECIII
int const NE_EC3 = 240;
int const NI_EC3 = 160;  //check point
int const N_EC3  = NE_EC3 + NI_EC3;
int const NE_EC3_L =120;
int const NE_EC3_R =240;

// cAB: A:post ensemble B:pre ensemble
// GAB: A:post ensemble B:pre ensemble
double const GEE_EC3 = 0.01; 
double const GEI_EC3 = 2.5; //check point 1.6->1.5->2.0->2.5
double const GIE_EC3=  0.5; // check point 0.2 -> 0.25 ->0.5
double const GII_EC3 = 0.01;

double const cEE_EC3 = 0.025;
double const cEI_EC3 = 0.5;
double const cIE_EC3 = 0.3;
double const cII_EC3 = 0.3;

double const Icnst_EC3= 0; 

double const rec3     = 35; //[Hz]
double const Gn_EC3_E = 0.02; 
double const Gn_EC3_I = 0.015; 
double const Gn_EC3_IE = 0.005;  
double const Gn_EC3_II = 0.02; 



// ECV
int const NE_EC5 = 400;
int const NI_EC5 = 120;
int const N_EC5  = NE_EC5 + NI_EC5;
int const NE_EC5_L =200;  
int const NE_EC5_R =400;  

double GEE_EC5 = 0.0;  // an argument of this script
double GEI_EC5 = 1.0;
double GIE_EC5 = 0.5;
double const GII_EC5 = 0.01;
double const cEE_EC5 = 0.005;
double const cIE_EC5 = 0.3;
double const cEI_EC5 = 0.5; 
double const cII_EC5 = 0.3; 

double const Icnst_EC5= 0.0; 

double const rec5     = 30;  
double const Gn_EC5_E = 0.005;  //check point 2017.11.26  0.0055
double const Gn_EC5_I = 0.005;
double const Gn_EC5_IE = 0.01;
double const Gn_EC5_II = 0.015;  




// CA1
int const NE_CA1 = 480;  
double const rIE = 0.5;
int const NPV    = NE_CA1*rIE;
double const rTE = 0.5;
int const NOLM   = NE_CA1*rTE;
int const NVIP   = NE_CA1*rTE;

int const N_CA1  = NE_CA1+NPV+NOLM;
int const NE_CA1_L=120;
int const NE_CA1_R=240;
int const NE_CA1_other=360;

double const GEE_CA1 = 0;   // Gpost,pre
double const GEI_CA1 = 2.2;
double const GET_CA1 = 0.01;
double const GIE_CA1 = 0.2;
double const GII_CA1 = 0.3; 
double const GIT_CA1 = 0.5;
double const GTE_CA1 = 0.5;
double const GTI_CA1 = 0.2;
double const GTT_CA1 = 0.0;

double const cEE_CA1 = 0.2;   // Gpost,pre
double const cEI_CA1 = 0.5; 
double const cET_CA1 = 0.1;
double const cIE_CA1 = 0.5;
double const cII_CA1 = 0.5;
double const cIT_CA1 = 0.5;
double const cTE_CA1 = 0.3;
double const cTI_CA1 = 0.5;
double const cTT_CA1 = 0.0;

double const Golm_ach = 0.02;
double const Icnst_E=-0.1;
double const Icnst_T=0.2; 

// noise into CA1
double const rca1   =  35;  //[Hz] //check point 30 -> 35
double const Gn_CA1_TE =  0.08;  
double const Gn_CA1_TI =  0.01;
double const Gn_CA1_IE =  0.02; 
double const Gn_CA1_II =  0.005;
double const Gn_CA1_EE =  0.008;
double const Gn_CA1_EI =  0.015; //check point 

double const rvip   =  20;  //[Hz]



// CA3
int const NE_CA3  = 480; 
int const NE_CA3_L= 120;
int const NE_CA3_R= 240;
int const NE_CA3_other=360;
int const N_CA3 = NE_CA3; // in this program, we have 3 clusters of ca3 neurons corresponding to each location

int const N_ALN_EC3=0;
int const N_ALN_EC5=N_EC3+5;
int const N_ALN_CA1=N_ALN_EC5+N_EC5+5;
int const N_all    =N_ALN_CA1+N_CA1+NVIP+5;

double const rca3        =3;  //[Hz]
double const GE_CA3_CA1 = 1.0; 
double const GI_CA3_CA1 = 0.05; 
double const cEE_CA3_CA1 = 0.3;// check G3_CA3_CA1 in HH_gennet.h. here, we introduce non-learning network
double const cIE_CA3_CA1 = 0.5;




// external poisson neurons
double const Gec2     = 0.4;  //[uS] 
double const Gms_ec5   =  0.1; // from MS to ECV
double const Gms    =  3.7;



// from ECIII to CA1
double const GEE_EC3_CA1 = 0.001; 
double const cEE_EC3_CA1 = 0.01; //check point gemma_dec   in G3_EC3_CA1(),  we modified prob of connections between unrelated groups. check it
//ver=="no_EC3toPV", these variables are not used
double const GIE_EC3_CA1 = 0.05; 
double const cIE_EC3_CA1 = 0.3;
// effect by OLM
double const Asup        = 0.1;  
double const Tsup        = 5.0;  // [ms] 
// // nonlinear interaction in CA1
double const G_Caspk = 0.01; 
double const tdcy_ca=50.0;
double const Thst_ca=tdcy_ca*2; // a time we keep information of synaptic conductance.int const rca_st
int const    Nhst_ca=int(Thst_ca/hstp)+1;
double const Tca_spk_overlay=100.0;
double const Tamp        =10; // [ms] amplification range of CA3 and CA1 interactions 
double const Tampec3     =15; // [ms] amplification range of EC3 and CA1 interactions 
double const Tplat       =20; // [ms] duration of plat by EC3


// from ECV to ECIII
double const GEE_EC5_EC3 = 0.01; 
double const cEE_EC5_EC3 = 0.0125; 
#define GIE_EC5_EC3 0
#define cIE_EC5_EC3 0 


// from CA1 to ECV
double GEE_CA1_EC5 = 0.0025;    // argument
double const cEE_CA1_EC5 = 0.005;
#define GIE_CA1_EC5 0.01 
#define cIE_CA1_EC5 0.3 



// Calcium dynamics
double gCAN_def = 0.0;  // [mS/cm^2] 
double const Cca_th    =50;//check
double const Ca_L=0.0003;
double const Ca_H=0.004; //check point 0.0045 -> 0.004
double const rhigh=1.5; //check point 1.0 -> 1.5 -> 1.2
double const rlow=0.2;



// for decoder
int Flg_gamma = 0;  // 0: t>6000 with gamma, 1: t>6500 with gamma, 2: nogamma

int const NE_Dec   =100;
int const NI_Dec   =100;
int const N_Dec    =NE_Dec+NI_Dec;
int const NE_Dec_L =50;
int const NE_Dec_R =100;
int const NI_Dec_L =150;
int const NI_Dec_R =200;

double const GEE_Dec = 0.0; // in spike transmission function, there is no transmission between excitatory neurons by definition 
double const GEI_Dec = 1.5; //check point 3.0->2.5->2.0
double const GIE_Dec=  0.0; 
double const GII_Dec = 0.05;

double const cEE_Dec = 0.0;
double const cEI_Dec = 0.8;
double const cIE_Dec = 0.0;
double const cII_Dec = 0.6;

double const GEE_EC3_Dec=0.1;//0.15->0.1 
double const GEE_CA1_Dec=0.4;
double const GIE_EC3_Dec=1.0;//0.6->1.0->1.2 ->1.0
double const GIE_CA1_Dec=0.4;//0.2 -> 0.4
double const cEE_EC3_Dec=0.03; 
double const cEE_CA1_Dec=0.2; 
double const cIE_EC3_Dec=0.4; 
double const cIE_CA1_Dec=0.8;

double const Gn_Dec_E = 0.005; 
double const Gn_Dec_I = 0.03;
double Gn_Dec_IE = 0.05;  
double const Gn_Dec_II = 0.006; 
double const rdec      = 35; //[Hz] 
double Icnst_Dec =-0.2;
double amp_dec=1.0;

double Aosc_gamma=30; //check point
double Aosc_gamma_EC3=20; //check point
double mE_EC3_gamma  =0.02;
double mE_CA1_gamma  =0.005;
int    N_gamma_EC3  = 80;
int    N_gamma_CA1  = 80;
double Gn_EC3_gamma =0.013;
double Gn_CA1_gamma =0.01;






// common parameters
double const Thst_syn=100.0; // a time we keep information of synaptic conductance.   check point5 50.0 -> 100 to match with Thst_ca
int const Nhst=int(Thst_syn/hstp)+1;


const int Nns = 40;
const int FLG_noise=1;
double Vex = 2.0; //mV
double rex = 100.0; //Hz
double rex_e=0.0; //Hz only to excitatory neurons
double rex_i =0.0; //Hz

double Vosc= 2.0;
const int Nosc = 10;
double Aosc=5; 

double Aca3=30; //[Hz] //check point 30->20->30
double freq_ext=10; //[Hz]
double freq_gamma_h=80; //[Hz]

double osc_scl=1.0;

// synapse parameters [*10ms]
double const tdcy_E=5.3;
double const tris_E=0.05;
double const tdcy_I=9.1;
double const tris_I=0.07;
double const tdcy_T=22.0;
double const tris_T=2.0;

//check point in gamma_dec
double const tdcy_E_dec=4.0;
double const tdcy_I_dec=4.0;



double VE = 0.0;
double VI = -80.0;
double Vr = -30.0;
double Vth = 0.0;
double Vmax = 20.0;

//time constants 
double const tmE = 20.0;
double const tmI = 10.0;
//check point
//double ts = 2.0;
double ts = 2.0;
double tgaba = 8.0;
double tnmda = 100.0;
double tref = 5.0;

//synaptic delay
// within area
double dEmax_wti = 2.0;  // [ms]
double dEmin_wti = 0.0;
double dImax_wti = 2.0;
double dImin_wti = 0.0;

// between area
double dEmax_btw = 15.0;  // [ms]
double dEmin_btw = 10.0;  //
double dImax_btw = 15.0;  // this value is not used
double dImin_btw = 5.0;  // this value is not used

int    const N_plat      =(dEmax_btw+Tplat)/(double)hstp;



vector<double> dvec;
vector<int> ivec;
deque<double> ddeque;
deque<int> ideque;
typedef vector< vector<double> > DBLMAT;

double dice(){
  //	return rand()/(RAND_MAX + 1.0);
  return genrand_real2();
}


typedef struct str_var_EC3_tmp
{
  vector<double> v,m,n,h,ahp,rf,rs;
  vector<double> spts; //previous spike time
  vector< deque<double> > gE,gI,gtest;
} STR_VAR_EC3;

typedef struct str_var_ca1_tmp
{
  //caspk  vector<double> v,m,n,h,a,b,r,rf,rs,ahp,mca,hca;
  vector<double> v,m,n,h,a,b,r,rf,rs,ahp,t_Ca_spk;
  vector<double> spts,olmspts,olmspts1,olmspts2; //previous spike time
  vector< deque<double> > ca3spts,ec3spts,gE,gI,gT;
  vector< deque<int> > t_plat,if_plat;
  vector<int>   if_spk;
} STR_VAR_CA1;


typedef struct str_var_EC5_tmp
{
  vector<double> v,m,n,h,ahp,Cal,KM,Napm,Naph,KAa,KAb,Ca,mCAN,KC,hrate_can;
  vector<double> spts; //previous spike time
  vector< deque<double> > gE,gI;
} STR_VAR_EC5;

typedef struct str_var_Dec_tmp
{
  vector<double> v,m,n,h,ahp;
  vector<double> spts; //previous spike time
  vector< deque<int> >    recCA1_spts,recEC3_spts;
  vector< deque<double> > gE,gI,gE_CA1;
} STR_VAR_Dec;


typedef struct str_para_tmp
{
  double I;
  double gCAN;
} STR_PARA;


#define PVinact 1
#include"HH_memdyn.h"
#include"HH_gennet.h"





double ngn(){
	double u = dice(); double v = dice();
	return sqrt(-2.0*log(u))*cos(2.0*pi*v);
}

double v_to_g(double v, double gm){
	return gm*v;
}	



double calc_nrnd(void){
  static int sw=0;
  static double t,u;

  if(sw==0){
    sw=1;
    t=sqrt(-2*log(1-dice())); u = 2*pi*dice();
    return t*cos(u);
  }
  else{
    sw=0;
    return t*sin(u);
  }

}



void init_var_EC3(STR_VAR_EC3* pstr_var_tmp){
  for(int i = 0; i < N_EC3; i++){
    pstr_var_tmp->v.push_back(VI+dice()*(Vr-VI));
    pstr_var_tmp->n.push_back(0.2*dice());
    pstr_var_tmp->h.push_back(1.0-0.2*dice());
    pstr_var_tmp->m.push_back(0.2*dice());
    pstr_var_tmp->ahp.push_back(dice());
    pstr_var_tmp->rf.push_back(dice());
    pstr_var_tmp->rs.push_back(dice());
    pstr_var_tmp->spts.push_back(-1000.0);
  }
  // input current into CA1 neurons: dual exponential function (cf principles of comput. neurosci.)
  for(int i=0; i<N_EC3;i++){
    pstr_var_tmp->gE.push_back(ddeque);
    pstr_var_tmp->gI.push_back(ddeque);
    for( int j=0; j<Nhst;j++){
      pstr_var_tmp->gE[i].push_back(0.0);
      pstr_var_tmp->gI[i].push_back(0.0);
    }
    //check point
    pstr_var_tmp->gtest.push_back(ddeque);
    for( int j=0; j<Nhst;j++)
      pstr_var_tmp->gtest[i].push_back(0.0);;
  }
  
  return;
}


void init_var_CA1(STR_VAR_CA1* pstr_var_tmp){
  // caspk  double dtmp;
  for(int i = 0; i < N_CA1; i++){
    pstr_var_tmp->v.push_back(VI+dice()*(Vr-VI));
    pstr_var_tmp->n.push_back(0.2*dice());
    pstr_var_tmp->h.push_back(1.0-0.2*dice());
    pstr_var_tmp->m.push_back(0.2*dice());
    pstr_var_tmp->a.push_back(dice());
    pstr_var_tmp->b.push_back(dice());
    pstr_var_tmp->r.push_back(dice());
    pstr_var_tmp->rf.push_back(dice());
    pstr_var_tmp->rs.push_back(dice());
    pstr_var_tmp->ahp.push_back(0.0);
    /* caspk
       dtmp=0.0;
    pstr_var_tmp->mca.push_back(dtmp);
    pstr_var_tmp->hca.push_back(1-dtmp);
    */
    pstr_var_tmp->t_Ca_spk.push_back(-1000.0);
    pstr_var_tmp->spts.push_back(-1000.0);
    pstr_var_tmp->ca3spts.push_back(ddeque);
    for(int j=0;j<3;j++) pstr_var_tmp->ca3spts[i].push_back(-1000.0-j*100.0);
    pstr_var_tmp->ec3spts.push_back(ddeque);
    for(int j=0;j<3;j++) pstr_var_tmp->ec3spts[i].push_back(-1000.0-j*100.0);
    pstr_var_tmp->olmspts.push_back(-1000.0-500);
    pstr_var_tmp->olmspts1.push_back(-1000.0-500);
    pstr_var_tmp->olmspts2.push_back(-1000.0-500);
    pstr_var_tmp->t_plat.push_back(ideque);
    pstr_var_tmp->if_plat.push_back(ideque);
    for(int j=0;j<N_plat;j++){
      pstr_var_tmp->t_plat[i].push_back(0);
      pstr_var_tmp->if_plat[i].push_back(0);
    }
    pstr_var_tmp->if_spk.push_back(0);
  }
  // input current into CA1 neurons: dual exponential function (cf principles of comput. neurosci.)
  for(int i=0; i<N_CA1;i++){
    pstr_var_tmp->gE.push_back(ddeque);
    pstr_var_tmp->gI.push_back(ddeque);
    pstr_var_tmp->gT.push_back(ddeque);
    for( int j=0; j<Nhst;j++){
      pstr_var_tmp->gE[i].push_back(0.0);
      pstr_var_tmp->gI[i].push_back(0.0);
      pstr_var_tmp->gT[i].push_back(0.0);
    }
  }

  return;
}


void init_var_EC5(STR_VAR_EC5* pstr_var_tmp){
  for(int i = 0; i < N_EC5; i++){
    pstr_var_tmp->v.push_back(VI+dice()*(Vr-VI));
    pstr_var_tmp->n.push_back(0.2*dice());
    pstr_var_tmp->h.push_back(1.0-0.2*dice());
    pstr_var_tmp->m.push_back(0.2*dice());
    pstr_var_tmp->ahp.push_back(0.2); // check point (if set these values between 0 and 1, spike timings are quite similar)  
    pstr_var_tmp->Cal.push_back(0.0);
    pstr_var_tmp->KM.push_back(0.0); // check point (if set these values between 0 and 1, spike timings are quite similar)  
    pstr_var_tmp->Napm.push_back(0.0); // check point (if set these values between 0 and 1, spike timings are quite similar)  
    pstr_var_tmp->Naph.push_back(0.0); // check point (if set these values between 0 and 1, spike timings are quite similar)  
    pstr_var_tmp->KAa.push_back(0.0); // check point (if set these values between 0 and 1, spike timings are quite similar)  
    pstr_var_tmp->KAb.push_back(0.0); // check point (if set these values between 0 and 1, spike timings are quite similar)  
    pstr_var_tmp->Ca.push_back(0.0);
    pstr_var_tmp->mCAN.push_back(0.0);
    pstr_var_tmp->KC.push_back(0.0);

    pstr_var_tmp->spts.push_back(-1000.0);
    pstr_var_tmp->hrate_can.push_back(1.0);
  }
  // input current into CA1 neurons: dual exponential function (cf principles of comput. neurosci.)
  for(int i=0; i<N_EC5;i++){
    pstr_var_tmp->gE.push_back(ddeque);
    pstr_var_tmp->gI.push_back(ddeque);
    for( int j=0; j<Nhst;j++){
      pstr_var_tmp->gE[i].push_back(0.0);
      pstr_var_tmp->gI[i].push_back(0.0);
    }
  }
  
  return;
}

// for decoder
void init_var_Dec(STR_VAR_Dec* pstr_var_tmp){
  for(int i = 0; i < N_Dec; i++){
    pstr_var_tmp->v.push_back(VI+dice()*(Vr-VI));
    pstr_var_tmp->n.push_back(0.2*dice());
    pstr_var_tmp->h.push_back(1.0-0.2*dice());
    pstr_var_tmp->m.push_back(0.2*dice());
    pstr_var_tmp->ahp.push_back(0.2);
    pstr_var_tmp->spts.push_back(-1000.0);
  }
  for(int i=0; i<N_Dec;i++){
    pstr_var_tmp->gE.push_back(ddeque);
    pstr_var_tmp->gE_CA1.push_back(ddeque);
    pstr_var_tmp->gI.push_back(ddeque);
    for( int j=0; j<Nhst;j++){
      pstr_var_tmp->gE[i].push_back(0.0);
      pstr_var_tmp->gE_CA1[i].push_back(0.0);
      pstr_var_tmp->gI[i].push_back(0.0);
    }
  }
  return;
}



void calc(int ik,int it,double c_Ach_aft, double Ach_test){
  double dtmp;
  int itmp;
  double t;
  double lfp;
  string ver="no_EC3toPV";
  double c_Ach=0.; // consentration of Ach
  
  //  srand(ik*31);
  init_genrand(ik*31);

  //  vector< vector<double> > 
  DBLMAT G_EC3   = calc_G3("EC3");   // G[i][j]: i:post j:pre
  DBLMAT G_EC5   = calc_G3("EC5");
  DBLMAT G_CA1   = calc_G_CA1("CA1");

  DBLMAT G_EC3_CA1=calc_G3_EC3_CA1("EC3","CA1",ver);
  DBLMAT G_CA3_CA1=calc_G3_CA3_CA1("CA3","CA1",ver);
  DBLMAT G_CA1_EC5=calc_G3_CA1_EC5("CA1","EC5");
  DBLMAT G_EC5_EC3=calc_G3_EC5_EC3("EC5","EC3");
  
  //  vector< vector<double> > 
  vector< vector<int> >    d_EC3    = calc_d(N_EC3,NE_EC3);  //d[i][j]: i:post j:pre
  vector< vector<int> >    d_EC5    = calc_d(N_EC5,NE_EC5);
  vector< vector<int> >    d_CA1    = calc_d(N_CA1,NE_CA1);

  vector< vector<int> >    d_EC3_CA1    = calc_d_btw(N_EC3,NE_EC3,N_CA1,NE_CA1);
  vector< vector<int> >    d_CA1_EC5    = calc_d_btw(N_CA1,NE_CA1,N_EC5,NE_EC5);
  vector< vector<int> >    d_EC5_EC3    = calc_d_btw(N_EC5,NE_EC5,N_EC3,NE_EC3);

  // for decoder
  DBLMAT G_Dec    =calc_Dec("Dec","Dec");
  DBLMAT G_EC3_Dec=calc_Dec("EC3","Dec");
  DBLMAT G_CA1_Dec=calc_Dec("CA1","Dec");
  vector< vector<int> >    d_Dec        = calc_d_dec(N_Dec,N_Dec,'w');  // w: within a network, b: between networks
  vector< vector<int> >    d_EC3_Dec    = calc_d_dec(N_EC3,N_Dec,'b');
  vector< vector<int> >    d_CA1_Dec    = calc_d_dec(N_CA1,N_Dec,'b');

  for(int i=0;i<N_Dec;i++){
    for(int j=0;j<d_EC3_Dec[i].size();j++)d_EC3_Dec[i][j]+=(int)(12.5/hstp); //check point dec
  }


  vector< vector<int> > Goutidx_EC3;  // Goutidx[i][j]: i:pre j:post
  for(int i = 0; i < N_EC3; i++){
    Goutidx_EC3.push_back(ivec);
    for(int j = 0; j < N_EC3; j++){
      if( G_EC3[j][i] > 0.0 ) Goutidx_EC3[i].push_back(j);
    }
  }

  vector< vector<int> > Goutidx_EC5;
  for(int i = 0; i < N_EC5; i++){
    Goutidx_EC5.push_back(ivec);
    for(int j = 0; j < N_EC5; j++){
      if( G_EC5[j][i] > 0.0 ) Goutidx_EC5[i].push_back(j);
    }
  }
  
  vector< vector<int> > Goutidx_CA1;
  for(int i = 0; i < N_CA1; i++){
    Goutidx_CA1.push_back(ivec);
    for(int j = 0; j < N_CA1; j++){
      if( G_CA1[j][i] > 0.0 ) Goutidx_CA1[i].push_back(j);
    }
  }


  /*
  for(int i = NE_CA1+NPV; i < N_CA1; i++){
    cout<<i<< "| ";
    for(int jidx = 0; jidx < Goutidx_CA1[i].size(); jidx++){
      int j = Goutidx_CA1[i][jidx];
      if(j<NE_CA1)  cout<< j << " " ;
    }
    cout<<endl;
  }
  exit(0);
  */

  vector< vector<int> > Goutidx_EC3_CA1;
  for(int i = 0; i < N_EC3; i++){
    Goutidx_EC3_CA1.push_back(ivec);
    for(int j = 0; j < N_CA1; j++){
      if( G_EC3_CA1[j][i] > 0.0 ) Goutidx_EC3_CA1[i].push_back(j);
    }
  }

  vector< vector<int> > Goutidx_CA1_EC5;
  for(int i = 0; i < N_CA1; i++){
    Goutidx_CA1_EC5.push_back(ivec);
    for(int j = 0; j < N_EC5; j++){
      if( G_CA1_EC5[j][i] > 0.0 ) Goutidx_CA1_EC5[i].push_back(j);
    }
  }
    
  
  vector< vector<int> > Goutidx_CA3_CA1;
  for(int i = 0; i < N_CA3; i++){
    Goutidx_CA3_CA1.push_back(ivec);
    for(int j = 0; j < NE_CA1; j++){
      if( G_CA3_CA1[j][i] > 0.0 )  Goutidx_CA3_CA1[i].push_back(j);
    }
    for(int j = NE_CA1; j < NE_CA1+NPV; j++){
      if( G_CA3_CA1[j][i] > 0.0 )  Goutidx_CA3_CA1[i].push_back(j);
    }
  }

  vector< vector<int> > Goutidx_EC5_EC3;
  for(int i = 0; i < N_EC5; i++){
    Goutidx_EC5_EC3.push_back(ivec);
    for(int j = 0; j < N_EC3; j++){
      if( G_EC5_EC3[j][i] > 0.0 ) Goutidx_EC5_EC3[i].push_back(j);
    }
  }



  // setting for decoder
  vector< vector<int> > Goutidx_Dec;
  for(int i = 0; i < N_Dec; i++){
    Goutidx_Dec.push_back(ivec);
    for(int j = 0; j < N_Dec; j++){
      if( G_Dec[j][i] > 0.0 ) Goutidx_Dec[i].push_back(j);
    }
  }

  vector< vector<int> > Goutidx_EC3_Dec;
  for(int i = 0; i < N_EC3; i++){
    Goutidx_EC3_Dec.push_back(ivec);
    for(int j = 0; j < N_Dec; j++){
      if( G_EC3_Dec[j][i] > 0.0 ) Goutidx_EC3_Dec[i].push_back(j);
    }
  }

  vector< vector<int> > Goutidx_CA1_Dec;
  for(int i = 0; i < N_CA1; i++){
    Goutidx_CA1_Dec.push_back(ivec);
    for(int j = 0; j < N_Dec; j++){
      if( G_CA1_Dec[j][i] > 0.0 ) Goutidx_CA1_Dec[i].push_back(j);
    }
  }
  /**************/


  /**********   read a spike data   ***************/
#ifdef readmode
  ostringstream ossd; 

  if(Flg_gamma==0)
    ossd <<  "sample_rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdec4long"<<Ach_test<<Aosc_gamma<<Aosc_gamma_EC3<<"0.10.0020.018080gamma12"; 
  else if(Flg_gamma==1)
    ossd <<  "sample_rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdec4long"<<Ach_test<<Aosc_gamma<<Aosc_gamma_EC3<<"0.10.0020.018080gamma1"; 
  else if(Flg_gamma==2)
    ossd <<  "sample_rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdec4long"<<Ach_test<<Aosc_gamma<<Aosc_gamma_EC3<<"0.10.0020.018080gamma3"; 
  else if(Flg_gamma==3)
    ossd <<  "sample_rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdec4long"<<Ach_test<<Aosc_gamma<<Aosc_gamma_EC3<<"0.10.0020.018080nogamma"; 



  string fstr=ossd.str();   ifstream ifs(fstr.c_str()); 

  if(ifs.fail()) {
    cerr << "File do not exist."<<fstr.c_str()<<"\n" ;
    exit(0);
  }


  string str;

  DBLMAT spikeEC3;  // output filter
  DBLMAT spikeCA1;  // output filter
  for(int i=0;i<NE_EC3;i++)   spikeEC3.push_back(dvec);
  for(int i=0;i<NE_CA1;i++)   spikeCA1.push_back(dvec);

  while(true){
    getline(ifs, str);
    sscanf(str.data(), "%lg %d ", &dtmp, &itmp);
    if(itmp>=N_ALN_EC3 && itmp<N_ALN_EC3+NE_EC3)spikeEC3[itmp-N_ALN_EC3].push_back(dtmp);
    if(itmp>=N_ALN_CA1 && itmp<N_ALN_CA1+NE_CA1)spikeCA1[itmp-N_ALN_CA1].push_back(dtmp);
    if(ifs.eof()) break;
  }
    

  /**********************/


  // check point 
  ostringstream osst;
  ostringstream ossf;

  if(Flg_gamma==0){
    osst <<  "rast_8080gamma12"<<ik<<it<<Icnst_Dec<<Gn_Dec_IE<<".dat";
    ossf <<  "dyn_8080gamma12"<<ik<<it<<Icnst_Dec<<Gn_Dec_IE<<".dat";
  }
  else if(Flg_gamma==1){
    osst <<  "rast_8080gamma13"<<ik<<it<<Icnst_Dec<<Gn_Dec_IE<<".dat";
    ossf <<  "dyn_8080gamma13"<<ik<<it<<Icnst_Dec<<Gn_Dec_IE<<".dat";
  }
  else if(Flg_gamma==2){
    osst <<  "rast_8080gamma3"<<ik<<it<<Icnst_Dec<<Gn_Dec_IE<<".dat";
    ossf <<  "dyn_8080gamma3"<<ik<<it<<Icnst_Dec<<Gn_Dec_IE<<".dat";
  }
  else if(Flg_gamma==3){
    osst <<  "rast_8080nogamma"<<ik<<it<<Icnst_Dec<<Gn_Dec_IE<<".dat";
    ossf <<  "dyn_8080nogamma"<<ik<<it<<Icnst_Dec<<Gn_Dec_IE<<".dat";
  }


#else
  ostringstream osst;
  ostringstream ossf;
  if(osc_flg==1){
    osst <<  "rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<osc_scl<<".datdeclong_notheta";
    ossf <<  "dyn_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<osc_scl<<".datdeclong_notheta";
  }
 else if(Tach_end_str>Tec3off_strt){
   osst <<  "rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_block"<<Ach_test; 
   ossf <<  "dyn_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_block"<<Ach_test; 
 }
 else if(Tach_inc>Tec3off_strt){
   osst <<  "rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_block3"<<Ach_test; 
   ossf <<  "dyn_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_block3"<<Ach_test; 
 }
 else if( T > Tec3off_strt){
   osst <<  "rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_block5"<<Ach_test<<"0.1"; 
   ossf <<  "dyn_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_block5"<<Ach_test<<"0.1"; 
  }
  else if(Icnst_PV<-1.0){
    osst <<  "rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_PVinact";
    ossf <<  "dyn_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_PVinact";
  }
  else if(geta>0.0){
    osst <<  "rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_disord";
    ossf <<  "dyn_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdeclong_disord";
  }
  else{
    osst <<  "rast_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdec4long"<<Ach_test<<Aosc_gamma<<Aosc_gamma_EC3<<Icnst_Dec<<Gn_Dec_IE<<"0.01"<<N_gamma_EC3<<N_gamma_CA1<<"gamma15_1"; 
    ossf <<  "dyn_EC5EC3" <<GEE_EC5_EC3 <<"_CA1EC5"<<GEE_CA1_EC5<<"_EC3_CA1"<<GEE_EC3_CA1<<"_Gpv"<<Gpv_ach<<"_cCA1EC5"<<cEE_CA1_EC5<<"_EC5"<<GEE_EC5<<"_gCAN"<<gCAN_def<<c_Ach_aft<<"GEI_EC5"<<GEI_EC5<<"GIE_EC5"<<GIE_EC5<<ik<<"_"<<it<<Gn_EC3_I<<".datdec4long"<<Ach_test<<Aosc_gamma<<Aosc_gamma_EC3<<Icnst_Dec<<Gn_Dec_IE<<"0.01"<<N_gamma_EC3<<N_gamma_CA1 <<"gamma15_1"; 
  }

#endif

  string fstrt = osst.str(); ofstream ofst;   ofst.open( fstrt.c_str() );
  fstrt = ossf.str(); ofstream ofsf; ofsf.open( fstrt.c_str() );

  STR_VAR_EC3 str_var_EC3_tmp;
  STR_VAR_CA1 str_var_CA1_tmp;
  STR_VAR_EC5 str_var_EC5_tmp;
  STR_PARA   str_paras_tmp;
  
  STR_VAR_EC3* pstr_var_EC3;
  STR_VAR_CA1* pstr_var_CA1;
  STR_VAR_EC5* pstr_var_EC5;
  STR_PARA*    pstr_paras;
  pstr_var_EC3 = &str_var_EC3_tmp;
  pstr_var_CA1 = &str_var_CA1_tmp;
  pstr_var_EC5 = &str_var_EC5_tmp;
  pstr_paras   = &str_paras_tmp;
  

  // for decoder
  STR_VAR_Dec  str_var_Dec_tmp;
  STR_VAR_Dec* pstr_var_Dec;
  pstr_var_Dec = &str_var_Dec_tmp;
  
  for(int i=0;i<N_Dec;i++){
    pstr_var_Dec->recCA1_spts.push_back(ideque);
    itmp=0;
    for(int j=0;j<d_CA1_Dec[i].size();j++) itmp= itmp>d_CA1_Dec[i][j] ? itmp:d_CA1_Dec[i][j];
    for(int j=0;j<itmp;j++)    pstr_var_Dec->recCA1_spts[i].push_back(0);
    
    pstr_var_Dec->recEC3_spts.push_back(ideque);
    itmp=0;
    for(int j=0;j<d_EC3_Dec[i].size();j++) itmp= itmp>d_EC3_Dec[i][j] ? itmp:d_EC3_Dec[i][j];
    for(int j=0;j<itmp;j++)                pstr_var_Dec->recEC3_spts[i].push_back(0);
  }
  
  
  
  
  //  srand(itmp);
  itmp=ik*(it+1)*313+ik+it;
  init_genrand(itmp);


  init_var_EC3(pstr_var_EC3);
  init_var_CA1(pstr_var_CA1);
  init_var_EC5(pstr_var_EC5); 
  init_var_Dec(pstr_var_Dec); 


  vector<double> vgEI;

  double curr_dyn_E[Nhst],curr_dyn_I[Nhst],curr_dyn_T[Nhst];
  double curr_dyn_E_dec[Nhst],curr_dyn_I_dec[Nhst];
  for(int i=0; i<Nhst;i++){
    curr_dyn_E[i]=(1/(tris_E-tdcy_E))*(exp(-i*hstp/tris_E)-exp(-i*hstp/tdcy_E));
    curr_dyn_I[i]=(1/(tris_I-tdcy_I))*(exp(-i*hstp/tris_I)-exp(-i*hstp/tdcy_I));
    curr_dyn_T[i]=(1/(tris_T-tdcy_T))*(exp(-i*hstp/tris_T)-exp(-i*hstp/tdcy_T));
    curr_dyn_E_dec[i]=(1/(tris_E-tdcy_E_dec))*(exp(-i*hstp/tris_E)-exp(-i*hstp/tdcy_E_dec));
    curr_dyn_I_dec[i]=(1/(tris_I-tdcy_I_dec))*(exp(-i*hstp/tris_I)-exp(-i*hstp/tdcy_I_dec));
  }
  double Ca_g[Nhst_ca]; //check point time scale 100 -> 50
  for(int i=0;i<5.0/hstp;i++)           Ca_g[i]=i/(5.0/hstp);
  for(int i=5.0/hstp;i<(35.0/hstp);i++) Ca_g[i]=1-(i-5.0/hstp)/(60.0/hstp);
  for(int i=(35.0/hstp);i<Nhst_ca;i++)  Ca_g[i]=0.5*exp(-(i-(35.0/hstp))/(30.0/hstp));
  


  for(t = 0.0; t < T+hstp; t += hstp){

    // cholinergic control
    if( t>= Tach_strt && t<Tach_end_str)    c_Ach=1.0-0.8*exp(-(t-Tach_strt)/200.0);
    else if( t>=Tach_end_str && t<Tach_inc) c_Ach=c_Ach_aft+(1.0-c_Ach_aft)*exp(-(t-Tach_end_str)/200.0);
    else if( t>=Tach_inc && t<Tach_end)     c_Ach=Ach_test-(Ach_test-c_Ach_aft)*exp(-(t-Tach_inc)/200.0);//check point
    else                  c_Ach=0.2; //check point
      

#ifndef readmode
    // VIP input driven by Ach
    double Gosc;
    Gosc=rvip*(sin(2*pi*(freq_ext*t/1000.0+0.5-0.1))+0.3); 
    // S=3.03, scale keeping S is 0.4827
    if(osc_flg==1)
      Gosc=rvip*(1+0.3)*osc_scl;
      //Gosc=rvip*0.4827;
    
    if(Gosc<0.0) Gosc=0.0;

    for(int j=0;j<NVIP;j++){
      // inh into PV
      for(int i = NE_CA1; i < NE_CA1+NPV; i++){
	if( dice() < Gosc/(1000.0/hstp)){
	  if(i==NE_CA1)ofst << t << " " << j+N_all-NVIP << endl; 
	  for( int k=0; k<Nhst; k++)
	    pstr_var_CA1->gI[i][k] += c_Ach*Gpv_ach*curr_dyn_I[k];
	}
      }
      // inh into OLM
      for(int i = NE_CA1+NPV; i < N_CA1; i++){
	if( dice() < Gosc/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_CA1->gI[i][k] += c_Ach*Golm_ach*curr_dyn_I[k];
	}
      }
    }
    
    
    if(c_Ach>0.0)
      pstr_paras->gCAN=c_Ach*gCAN_def;
    else
      pstr_paras->gCAN=0;



    for(int i = 0; i < NE_EC3; i++){ 
      for(int j=0;j<Nns;j++){
	if( dice() < rec3/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_EC3->gE[i][k] += Gn_EC3_E*curr_dyn_E[k];
	}
      }
      for(int j=0;j<Nns;j++){
	if( dice() < rec3/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_EC3->gI[i][k] += Gn_EC3_I*curr_dyn_I[k];
	}
      }
    }

    for(int i = NE_EC3; i < N_EC3; i++){ 
      for(int j=0;j<Nns;j++){
	if( dice() < rec3/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_EC3->gE[i][k] += Gn_EC3_IE*curr_dyn_E[k];
	}
      }
      for(int j=0;j<Nns;j++){
	if( dice() < rec3/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_EC3->gI[i][k] += Gn_EC3_II*curr_dyn_I[k];
	}
      }
    }

	
    for(int i = 0; i < NE_EC5; i++){ 
      for(int j=0;j<Nns;j++){
	if( dice() < rec5/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_EC5->gE[i][k] += Gn_EC5_E*curr_dyn_E[k];
	}
      }
      for(int j=0;j<Nns;j++){
	if( dice() < rec5/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_EC5->gI[i][k] += Gn_EC5_I*curr_dyn_I[k];
	}
      }
    }

    for(int i = NE_EC5; i < N_EC5; i++){ 
      for(int j=0;j<Nns;j++){
	if( dice() < rec5/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++){
	    pstr_var_EC5->gE[i][k] += Gn_EC5_IE*curr_dyn_E[k];
	  }
	}
      }
      for(int j=0;j<Nns;j++){	
	if( dice() < rec5/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++){
	    pstr_var_EC5->gI[i][k] += Gn_EC5_II*curr_dyn_I[k];
	  }
	}
      }
    }

    
    for(int i = 0; i < NE_CA1; i++){ 
      for(int j=0;j<Nns;j++){
	if( dice() < rca1/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_CA1->gE[i][k] += Gn_CA1_EE*curr_dyn_E[k];
	}
      }
      for(int j=0;j<Nns;j++){
	if( dice() < rca1/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_CA1->gI[i][k] += Gn_CA1_EI*curr_dyn_I[k];
	}
      }
    }
    for(int i = NE_CA1; i < NE_CA1+NPV; i++){ 
      for(int j=0;j<Nns;j++){
	if( dice() < rca1/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_CA1->gE[i][k] += Gn_CA1_IE*curr_dyn_E[k];
	}
      }
      for(int j=0;j<Nns;j++){
	if( dice() < rca1/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_CA1->gI[i][k] += Gn_CA1_II*curr_dyn_I[k];
	}
      }
    }
    for(int i = NE_CA1+NPV; i < N_CA1; i++){ 
      for(int j=0;j<Nns;j++){
	if( dice() < rca1/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_CA1->gE[i][k] += Gn_CA1_TE*curr_dyn_E[k];
	}
      }
      for(int j=0;j<Nns;j++){
	if( dice() < rca1/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_CA1->gI[i][k] += Gn_CA1_TI*curr_dyn_I[k];
	}
      }
    }
#endif

    // noise to decoder
    for(int i = 0; i < NE_Dec; i++){ 
      for(int j=0;j<Nns;j++){
	if( dice() < rdec/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_Dec->gE[i][k] += Gn_Dec_E*curr_dyn_E_dec[k];
	}
      }
      for(int j=0;j<Nns;j++){
	if( dice() < rdec/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_Dec->gI[i][k] += Gn_Dec_I*curr_dyn_I_dec[k];
	}
      }
    }

    for(int i = NE_Dec; i < N_Dec; i++){ 
      for(int j=0;j<Nns;j++){
	if( dice() < rdec/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_Dec->gE[i][k] += Gn_Dec_IE*curr_dyn_E_dec[k];
	}
      }
      for(int j=0;j<Nns;j++){
	if( dice() < rdec/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)
	    pstr_var_Dec->gI[i][k] += Gn_Dec_II*curr_dyn_I_dec[k];
	}
      }
    }


    



#ifndef  readmode     

    // current from external noisy stimuli.		  //check point
    if(t>6500){
      Gosc=sin(2*pi*(freq_gamma_h*t/1000.0));
      if(Gosc>=0.0){
	for(int i = 0; i < NE_EC3; i++){ 
	  pstr_var_EC3->gE[i][1]+=Gn_EC3_gamma*Gosc;
	  pstr_var_EC3->gtest[i][1]+=Gn_EC3_gamma*Gosc;
	}
      }
      if(Gosc<0.0){
	for(int i = 0; i < NE_EC3; i++){ 
	  pstr_var_EC3->gI[i][1]+=Gn_EC3_gamma*(-10*Gosc);
	}
      }
      if(Gosc>=0.0){
	for(int i = 0; i < NE_CA1; i++) 
	  pstr_var_CA1->gE[i][1]+=Gn_CA1_gamma*Gosc;
      }
      if(Gosc<0.0){
	for(int i = 0; i < NE_CA1; i++) 
	  pstr_var_CA1->gI[i][1]+=Gn_CA1_gamma*(-10*Gosc);
      }
    }



    // oscillatory input  phase=0: bottom in LFP in CA1-SP
    Gosc=Aosc*(sin(2*pi*(freq_ext*t/1000.0+0.6-0.1))+2.0);  // this is presumably from MS gabaergic which is common in OLM and EC5 //check point 0.5 -> 0.6(r48) -> 0.5 (r47) offset 1.0 ->1.5 1.3(r47) 2.0(r49)
    if(osc_flg==1)
      Gosc=Aosc*(1+1)*osc_scl;
    //      Gosc=Aosc*1.0;
    //S=2*pi, the balance scale is 1
    if(Gosc<0.0) Gosc=0.0;

    // inh into OLM
    for(int j=0;j<Nosc;j++){
      for(int i = NE_CA1+NPV; i < N_CA1; i++){
	if( dice() < Gosc/(1000.0/hstp)){
	  if(i==NE_CA1+NPV)ofst << t << " " << j+N_all << endl; 
	  for( int k=0; k<Nhst; k++)  pstr_var_CA1->gI[i][k] += Gms*curr_dyn_I[k];
	}
      }
    }

    // inh into Inh in EC5
    for(int j=0;j<Nosc;j++){
      for(int i = NE_EC5; i < N_EC5; i++){
	if( dice() < Gosc/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)  pstr_var_EC5->gI[i][k] += Gms_ec5*curr_dyn_I[k];
	}
      }
    }

    // osc to PV (Borhegyi,JNS,2004)
    Gosc=Aosc*(sin(2*pi*(freq_ext*t/1000.0+0.1))+0.8);  // this is presumably from MS gabaergic which is common in OLM and EC5 //check point  off set 1.0 -> 0.8
    if(osc_flg==1)
      Gosc=Aosc*(1+0.8)*osc_scl;
    for(int j=0;j<Nosc;j++){
      for(int i = NE_CA1; i < NE_CA1+NPV; i++){
	if( dice() < Gosc/(1000.0/hstp)){
	  for( int k=0; k<Nhst; k++)  pstr_var_CA1->gI[i][k] += Gms*curr_dyn_I[k];
	}
      }
    }


    //osc input from EC2 to inh in EC3
    double tmpoffset=0.6; //offset 0.5 -> 0.6
    Gosc=Aosc*(sin(2*pi*(freq_ext*t/1000.0+0.3-0.1+geta))+tmpoffset);  
    if(osc_flg==1)
      Gosc=Aosc*(1+tmpoffset)*osc_scl;  //check point
    //      Gosc=Aosc*0.6775;  //check point
    //S=4.257, the scale = 0.6775
    if(Gosc<0.0) Gosc=0.0;
    
    for(int j=0;j<Nosc;j++){  
      for(int i = NE_EC3; i < N_EC3; i++){
	if( dice() < Gosc/(1000.0/hstp)){ //check point
	  for( int k=0; k<Nhst; k++)
	    pstr_var_EC3->gE[i][k] += Gec2*curr_dyn_E[k];
	}
      }
    }



    // CA3 input
    if(t>=Tca3_C && t<Tca3_end){
      Gosc=Aca3*(sin(2*pi*(freq_ext*t/1000.0+0.2))+1.0);  //check point support 1.3 -> 0.5->1.0, shift 0.15 -> 0.25->0.2
      if(osc_flg==1)
	Gosc=Aca3*(1+0.3)*osc_scl;
      //Gosc=Aca3*0.4827;
      // the scale 0.4827
      if(Gosc<0.0) Gosc=0.0;
      int spk_flg=0;

      for(int i = 0; i < N_CA3; i++){
	spk_flg=0;
	if( ( ((t>Tca3_C && t<=Tca3_L) || (t>Tca3_C1 && t<Tca3_end)) && (NE_CA3_R<=i && i <NE_CA3_other)  ) ||  // center position and neurons in center cluster
	    (  (t>Tca3_L && t<=Tca3_o)                               &&  NE_CA3_L> i                      ) ||  // left position and neurons in left cluster
	    (  (t>Tca3_o && t<=Tca3_C1)                              && (NE_CA3_other<=i && i<NE_CA3  )   )     // other 
	    ){
	  if( dice() < (Gosc+rca3)/(1000.0/hstp) )  spk_flg=1;
	}
	else{
	  if( dice() < rca3/(1000.0/hstp) ) spk_flg=1;
	}
	
	if(spk_flg==1){
	  ofst << t << " " << i+N_all+Nosc*2 << endl;
	  for(int jidx = 0; jidx < Goutidx_CA3_CA1[i].size(); jidx++){
	    int j = Goutidx_CA3_CA1[i][jidx];
	    for( int k=0; k<Nhst; k++)  pstr_var_CA1->gE[j][k] += G_CA3_CA1[j][i]*curr_dyn_E[k];
	    pstr_var_CA1->ca3spts[j].pop_front();  pstr_var_CA1->ca3spts[j].push_back(t);
	  }
	}
      }
    }

    // current through recurrent connections
    for(int i = 0; i < N_EC3; i++){
      if( pstr_var_EC3->v[i] > Vth  && t - pstr_var_EC3->spts[i] > tref ){
	ofst << t << " " << i+N_ALN_EC3 << endl; 
	pstr_var_EC3->spts[i] = t;

	for(int jidx = 0; jidx < Goutidx_EC3[i].size(); jidx++){
	  int j = Goutidx_EC3[i][jidx];
	  if( i < NE_EC3 ){
	    for( int k=0; k<Nhst-d_EC3[j][i]; k++)
	      pstr_var_EC3->gE[j][k+d_EC3[j][i]] += G_EC3[j][i]*curr_dyn_E[k];
	  }
	  else if(i<N_EC3){
	    for( int k=0; k<Nhst-d_EC3[j][i]; k++)
	      pstr_var_EC3->gI[j][k+d_EC3[j][i]] += G_EC3[j][i]*curr_dyn_I[k];
	  }
	}

	if(t<Tec3off_strt || t>Tec3off_end){   // temporaly inhibit to EC3 to observe its effect
	  for(int j=0;j<N_CA1;j++)pstr_var_CA1->if_spk[j]=0;
	  for(int jidx = 0; jidx < Goutidx_EC3_CA1[i].size(); jidx++){
	    int j = Goutidx_EC3_CA1[i][jidx];
	    pstr_var_CA1->ec3spts[j].pop_front();  pstr_var_CA1->ec3spts[j].push_back(t+d_EC3_CA1[j][i]*hstp);
	    pstr_var_CA1->if_spk[j]=1;
	    double Aamp_olm= (t-pstr_var_CA1->olmspts[j]<Tsup) ? Asup:1; 
	    for( int k=0; k<Nhst-d_EC3_CA1[j][i]; k++)
	      pstr_var_CA1->gE[j][k+d_EC3_CA1[j][i]] += Aamp_olm*G_EC3_CA1[j][i]*curr_dyn_E[k];
	  }
	}

	// for decoder
	if(i < NE_EC3){
	  for(int jidx = 0; jidx < Goutidx_EC3_Dec[i].size(); jidx++){
	    int j = Goutidx_EC3_Dec[i][jidx];
	    for(int k=0; k<(int)(2.0/hstp);k++) pstr_var_Dec->recEC3_spts[j][d_EC3_Dec[j][i]-k]+=1;
	    for(int k=0; k<Nhst-d_EC3_Dec[j][i]; k++)
	      pstr_var_Dec->gE[j][k+d_EC3_Dec[j][i]] += G_EC3_Dec[j][i]*curr_dyn_E_dec[k];
	  }
	}


      }
    }


    // output from EC5
    for(int i = 0; i < N_EC5; i++){
      if( pstr_var_EC5->v[i] > Vth  && t - pstr_var_EC5->spts[i] > tref ){
	ofst << t << " " << i+N_ALN_EC5 << endl; 
	pstr_var_EC5->spts[i] = t;
	
	
	for(int jidx = 0; jidx < Goutidx_EC5[i].size(); jidx++){
	  int j = Goutidx_EC5[i][jidx];

	  if( i < NE_EC5 ){
	    for( int k=0; k<Nhst-d_EC5[j][i]; k++)
	      pstr_var_EC5->gE[j][k+d_EC5[j][i]] += G_EC5[j][i]*curr_dyn_E[k];
	  }
	  else if(i<N_EC5){
	    for( int k=0; k<Nhst-d_EC5[j][i]; k++)
	      pstr_var_EC5->gI[j][k+d_EC5[j][i]] += G_EC5[j][i]*curr_dyn_I[k];
	  }
	}
	if(t>Tcnct_strt){
	  for(int jidx = 0; jidx < Goutidx_EC5_EC3[i].size(); jidx++){
	    int j = Goutidx_EC5_EC3[i][jidx];
	    for( int k=0; k<Nhst-d_EC5_EC3[j][i]; k++)
	      pstr_var_EC3->gE[j][k+d_EC5_EC3[j][i]] += G_EC5_EC3[j][i]*curr_dyn_E[k];
	  }
	}
      }
    }

    // non-linear integration of EC3 and CA3 inputs.
    for(int i = 0; i < NE_CA1; i++){
      pstr_var_CA1->t_plat[i].pop_front();      pstr_var_CA1->t_plat[i].push_back(0);
      pstr_var_CA1->if_plat[i].pop_front();     pstr_var_CA1->if_plat[i].push_back(0);
      if(pstr_var_CA1->if_spk[i]==1 && pstr_var_CA1->ec3spts[i][2]-pstr_var_CA1->ec3spts[i][0]<Tampec3){  //detect burst (3 spikes within Tampec3) from ECIII to CA1
	int tmp_nstp=int((pstr_var_CA1->ec3spts[i][2]-t)/(double)hstp);
	pstr_var_CA1->t_plat[i][tmp_nstp]=1;
      }
      if(pstr_var_CA1->t_plat[i][0]==1 && t-pstr_var_CA1->olmspts2[i]>Tampec3){  // check plataeu (olmspts2 can be later than t, but this is quite small so here neglected. )
	int tmp_nplat=Tplat/(double)hstp;
	for(int j=0;j<tmp_nplat;j++)pstr_var_CA1->if_plat[i][j]=1;
      }
      if(t-pstr_var_CA1->ca3spts[i][0]<Tamp && pstr_var_CA1->if_plat[i][0]==1){
	if(t-pstr_var_CA1->t_Ca_spk[i]>Tca_spk_overlay){
	  for( int k=0; k<Nhst; k++) pstr_var_CA1->gE[i][k] += G_Caspk*Ca_g[k];
	  pstr_var_CA1->t_Ca_spk[i]=t;
	}
	ofst << t << " " << i+5000 << endl; 
      }
    }



    // output CA1
    for(int i = 0; i < N_CA1; i++){
      if( pstr_var_CA1->v[i] > Vth  && t - pstr_var_CA1->spts[i] > tref ){
	ofst << t << " " << i+N_ALN_CA1 << endl; 
	pstr_var_CA1->spts[i] = t;
	
	for(int jidx = 0; jidx < Goutidx_CA1[i].size(); jidx++){
	  int j = Goutidx_CA1[i][jidx];

	  if( i < NE_CA1 ){
	    for( int k=0; k<Nhst-d_CA1[j][i]; k++)
	      pstr_var_CA1->gE[j][k+d_CA1[j][i]] += G_CA1[j][i]*curr_dyn_E[k];
	  }
	  else if(i<NE_CA1+NPV){
	    for( int k=0; k<Nhst-d_CA1[j][i]; k++)
	      pstr_var_CA1->gI[j][k+d_CA1[j][i]] += G_CA1[j][i]*curr_dyn_I[k];
	  }
	  else{  
	    if(j<NE_CA1){//projections from OLM to Pyr inhibit only inputs from EC3. 
	      pstr_var_CA1->olmspts2[j]=pstr_var_CA1->olmspts1[j];
	      pstr_var_CA1->olmspts1[j]=pstr_var_CA1->olmspts[j];
	      pstr_var_CA1->olmspts[j]=t+d_CA1[j][i]*hstp;
	    }
	    else{
	      for( int k=0; k<Nhst-d_CA1[j][i]; k++)
		pstr_var_CA1->gT[j][k+d_CA1[j][i]] += G_CA1[j][i]*curr_dyn_T[k];
	    }
	  }
	}
	if( i < NE_CA1 ){	 
	  for(int jidx = 0; jidx < Goutidx_CA1_EC5[i].size(); jidx++){
	    int j = Goutidx_CA1_EC5[i][jidx];
	    for( int k=0; k<Nhst-d_CA1_EC5[j][i]; k++)    pstr_var_EC5->gE[j][k+d_CA1_EC5[j][i]] += G_CA1_EC5[j][i]*curr_dyn_E[k];
	  }
	}

	// for decoder
	if(i < NE_CA1){
	  for(int jidx = 0; jidx < Goutidx_CA1_Dec[i].size(); jidx++){
	    int j = Goutidx_CA1_Dec[i][jidx];
	    for(int k=0; k<(int)(2.0/hstp);k++)      pstr_var_Dec->recCA1_spts[j][d_CA1_Dec[j][i]-k]+=1;
	    for(int k=0; k<Nhst-d_CA1_Dec[j][i]; k++)
	      pstr_var_Dec->gE[j][k+d_CA1_Dec[j][i]] += G_CA1_Dec[j][i]*curr_dyn_E_dec[k];
	  }
	}

      }
    }
    
    
    if( ((int)floor(t/hstp))%1000 == 0 ){
      double nearestspiketime = 0.0;
      for(int i = 0; i < N_CA1; i++){
	if( pstr_var_CA1->spts[i] > nearestspiketime ) nearestspiketime = pstr_var_CA1->spts[i];
      }
      if( t-nearestspiketime > 1000.0 ) break;
    }
#else   // in readmode, imaginary spikes to Dec neurons
    for(int i =0;i<NE_CA1;i++){
      for(int l=0;l<spikeCA1[i].size();l++){
	if(fabs(spikeCA1[i][l]-t)<0.01){
	  ofst << t << " " << i+3000+N_Dec+NE_EC3+100 << endl; 
	  for(int jidx = 0; jidx < Goutidx_CA1_Dec[i].size(); jidx++){
	    int j = Goutidx_CA1_Dec[i][jidx];
	    for(int k=0; k<(int)(1.0/hstp);k++)    pstr_var_Dec->recCA1_spts[j][d_CA1_Dec[j][i]-k]+=1;
	    for(int k=0; k<Nhst-d_CA1_Dec[j][i]; k++){
	      pstr_var_Dec->gE[j][k+d_CA1_Dec[j][i]] += G_CA1_Dec[j][i]*curr_dyn_E_dec[k];
	      pstr_var_Dec->gE_CA1[j][k+d_CA1_Dec[j][i]] += G_CA1_Dec[j][i]*curr_dyn_E_dec[k];
	    }
	  }
	  break;
	}
      }
    }

    for(int i =0;i<NE_EC3;i++){
      for(int l=0;l<spikeEC3[i].size();l++){
	if(fabs(spikeEC3[i][l]-t)<0.01){
	  ofst << t << " " << i+3000+N_Dec << endl; 
	  for(int jidx = 0; jidx < Goutidx_EC3_Dec[i].size(); jidx++){
	    int j = Goutidx_EC3_Dec[i][jidx];
	    for(int k=0; k<(int)(1.0/hstp);k++)	      pstr_var_Dec->recEC3_spts[j][d_EC3_Dec[j][i]-k]+=1;
	    for(int k=0; k<Nhst-d_EC3_Dec[j][i]; k++)
	      pstr_var_Dec->gE[j][k+d_EC3_Dec[j][i]] += G_EC3_Dec[j][i]*curr_dyn_E_dec[k];
	  }
	  break;
	}
      }
    }

#endif



    // for decoder
    for(int i = 0; i < NE_Dec; i++){
      if(pstr_var_Dec->recEC3_spts[i][0]*pstr_var_Dec->recCA1_spts[i][0]>0.0){
	for(int k=0; k<Nhst; k++){
	  pstr_var_Dec->gE[i][k] += 5.0*((GEE_CA1_Dec/(NE_CA1_L*cEE_CA1_Dec))+(GEE_EC3_Dec/(NE_EC3_L*cEE_EC3_Dec)))*curr_dyn_E_dec[k];
	}

	for(int k=0; k<(int)(1.0/hstp);k++)  pstr_var_Dec->recEC3_spts[i][k]=0;
	for(int k=0; k<(int)(1.0/hstp);k++)  pstr_var_Dec->recCA1_spts[i][k]=0;
      }
      
      pstr_var_Dec->recEC3_spts[i].pop_front();pstr_var_Dec->recEC3_spts[i].push_back(0);
      pstr_var_Dec->recCA1_spts[i].pop_front();pstr_var_Dec->recCA1_spts[i].push_back(0);
    }
    
    for(int i = 0; i < N_Dec; i++){
      if( pstr_var_Dec->v[i] > Vth  && t - pstr_var_Dec->spts[i] > tref ){
	ofst << t << " " << i+3000 << endl; 
	pstr_var_Dec->spts[i] = t;
	/*  //check point
	for(int jidx = 0; jidx < Goutidx_Dec[i].size(); jidx++){
	  int j = Goutidx_Dec[i][jidx];
	  if( i < NE_Dec ){
	  for( int k=0; k<Nhst-d_Dec[j][i]; k++)
	      pstr_var_Dec->gE[j][k+d_Dec[j][i]] += G_Dec[j][i]*curr_dyn_E[k];
	      }
	*/
	if(i>=NE_Dec){
	  for(int jidx = 0; jidx < Goutidx_Dec[i].size(); jidx++){
	    int j = Goutidx_Dec[i][jidx];
	    for( int k=0; k<Nhst-d_Dec[j][i]; k++)
	      pstr_var_Dec->gI[j][k+d_Dec[j][i]] += G_Dec[j][i]*curr_dyn_I_dec[k];
	  }
	}
      }
    }




#ifndef readmode
    /*************************************************/

    /**********             LFP         **************/

    /*************************************************/

    //calc. LFP as sum of gE,gI of OLM,gI of CA1_E 
    lfp=0;
    for(int i = 0; i < NE_CA1; i++)
      lfp +=pstr_var_CA1->gI[i][0];
    for(int i = NE_CA1+NPV; i < N_CA1; i++){
      lfp+=pstr_var_CA1->gI[i][0];
      lfp-=pstr_var_CA1->gE[i][0];
    }
    

    /*************************************************/
    
    /******           1 time step ahead        *******/

    /*************************************************/

    // 1 step ahead (neurons in EC3)
    for(int i = 0; i < N_EC3; i++){
      pstr_var_EC3->gE[i].pop_front();    pstr_var_EC3->gE[i].push_back(0.0);
      pstr_var_EC3->gI[i].pop_front();    pstr_var_EC3->gI[i].push_back(0.0);			
      pstr_var_EC3->gtest[i].pop_front();    pstr_var_EC3->gtest[i].push_back(0.0);
      double I=-1*(pstr_var_EC3->gE[i][0]*(pstr_var_EC3->v[i]-VE)
		   +pstr_var_EC3->gI[i][0]*(pstr_var_EC3->v[i]-VI));

      if(i < NE_EC3){
	int n_DOF=7;
	double x[n_DOF];
	I+=Icnst_EC3;
	pstr_paras->I=I;
	
	x[0]=pstr_var_EC3->v[i];
	x[1]=pstr_var_EC3->n[i];
	x[2]=pstr_var_EC3->h[i];
	x[3]=pstr_var_EC3->m[i];
	x[4]=pstr_var_EC3->ahp[i];
	x[5]=pstr_var_EC3->rf[i];
	x[6]=pstr_var_EC3->rs[i];
	RG4(x,n_DOF,&t,pstr_paras,difffunc_EC3_E_pyr);
	pstr_var_EC3->v[i]=x[0];
	pstr_var_EC3->n[i]=x[1];
	pstr_var_EC3->h[i]=x[2];
	pstr_var_EC3->m[i]=x[3];
	pstr_var_EC3->ahp[i]=x[4];	
	pstr_var_EC3->rf[i]=x[5];
	pstr_var_EC3->rs[i]=x[6];
      }
      else{
	int n_DOF=3;
	double x[n_DOF];
	pstr_paras->I=I;
	
	x[0]=pstr_var_EC3->v[i];
	x[1]=pstr_var_EC3->n[i];
	x[2]=pstr_var_EC3->h[i];
	RG4(x,n_DOF,&t,pstr_paras,difffunc_CA1_PV);
	pstr_var_EC3->v[i]=x[0];
	pstr_var_EC3->n[i]=x[1];
	pstr_var_EC3->h[i]=x[2];
      }
    }


    // 1 step ahead (neurons in EC5)
    
    for(int i = 0; i < N_EC5; i++){
      pstr_var_EC5->gE[i].pop_front();    pstr_var_EC5->gE[i].push_back(0.0);
      pstr_var_EC5->gI[i].pop_front();    pstr_var_EC5->gI[i].push_back(0.0);			

      
      double I=-1*(pstr_var_EC5->gE[i][0]*(pstr_var_EC5->v[i]-VE)
		   +pstr_var_EC5->gI[i][0]*(pstr_var_EC5->v[i]-VI));

      
      if(i < NE_EC5){
	int n_DOF=14;
	double tmptmpgCAN=pstr_paras->gCAN;
	double x[n_DOF];
	I+=Icnst_EC5;
	pstr_paras->I=I;
	if(pstr_var_EC5->Ca[i]>Ca_H)      pstr_var_EC5->hrate_can[i]+=rhigh*hstp/1000.0;
	else if(pstr_var_EC5->Ca[i]<Ca_L) pstr_var_EC5->hrate_can[i]-=rlow*hstp/1000.0;
	if (pstr_var_EC5->hrate_can[i]>2.5) pstr_var_EC5->hrate_can[i]=2.5; //check point
	else if (pstr_var_EC5->hrate_can[i]<0.3) pstr_var_EC5->hrate_can[i]=0.3;
	pstr_paras->gCAN*=pstr_var_EC5->hrate_can[i];

	x[0]=pstr_var_EC5->v[i];
	x[1]=pstr_var_EC5->n[i];
	x[2]=pstr_var_EC5->h[i];
	x[3]=pstr_var_EC5->m[i];
	x[4]=pstr_var_EC5->ahp[i];
	x[5]=pstr_var_EC5->Cal[i];
	x[6]=pstr_var_EC5->KM[i];
	x[7]=pstr_var_EC5->Napm[i];
	x[8]=pstr_var_EC5->Naph[i];
	x[9]=pstr_var_EC5->KAa[i];
	x[10]=pstr_var_EC5->KAb[i];
	x[11]=pstr_var_EC5->Ca[i];
	x[12]=pstr_var_EC5->mCAN[i];
	x[13]=pstr_var_EC5->KC[i];
	RG4(x,n_DOF,&t,pstr_paras,difffunc_EC5_E_pyr_Fransen);
	pstr_var_EC5->v[i]=x[0];
	pstr_var_EC5->n[i]=x[1];
	pstr_var_EC5->h[i]=x[2];
	pstr_var_EC5->m[i]=x[3];
	pstr_var_EC5->ahp[i]=x[4];
	pstr_var_EC5->Cal[i]=x[5];
	pstr_var_EC5->KM[i]=x[6];
	pstr_var_EC5->Napm[i]=x[7];
	pstr_var_EC5->Naph[i]=x[8];
	pstr_var_EC5->KAa[i]=x[9];
	pstr_var_EC5->KAb[i]=x[10];
	pstr_var_EC5->Ca[i]=x[11];
	pstr_var_EC5->mCAN[i]=x[12];
	pstr_var_EC5->KC[i]=x[13];

	pstr_paras->gCAN=tmptmpgCAN;
      }
      else{
	int n_DOF=3;
	double x[n_DOF];
	pstr_paras->I=I;
	
	x[0]=pstr_var_EC5->v[i];
	x[1]=pstr_var_EC5->n[i];
	x[2]=pstr_var_EC5->h[i];
	RG4(x,n_DOF,&t,pstr_paras,difffunc_CA1_PV);
	pstr_var_EC5->v[i]=x[0];
	pstr_var_EC5->n[i]=x[1];
	pstr_var_EC5->h[i]=x[2];
      }
    }

    // 1 step ahead (neurons in CA1)
    //check
    for(int i = 0; i < N_CA1; i++){
      pstr_var_CA1->gE[i].pop_front();    pstr_var_CA1->gE[i].push_back(0.0);      pstr_var_CA1->gI[i].pop_front();    pstr_var_CA1->gI[i].push_back(0.0);			
      pstr_var_CA1->gT[i].pop_front();    pstr_var_CA1->gT[i].push_back(0.0);			

      double I=-1*(pstr_var_CA1->gE[i][0]*(pstr_var_CA1->v[i]-VE)
		   +pstr_var_CA1->gI[i][0]*(pstr_var_CA1->v[i]-VI)
		   +pstr_var_CA1->gT[i][0]*(pstr_var_CA1->v[i]-VI));
      // here, VE and VI(=VT) are shared with neurons in other areas

      //check
      //      cout << pstr_var_CA1->v[i]<<" "<<pstr_var_CA1->gE[i][0] << " "<< pstr_var_CA1->gI[i][0] << " "<< pstr_var_CA1->gT[i][0] << " "<<I<<" ";

      
      if(i<NE_CA1){
	int n_DOF=5;
	// caspk n_DOF=7;
	double x[n_DOF];
	x[0]=pstr_var_CA1->v[i];
	x[1]=pstr_var_CA1->n[i];
	x[2]=pstr_var_CA1->ahp[i];
	x[3]=pstr_var_CA1->rf[i];
	x[4]=pstr_var_CA1->rs[i];
	/* caspk 
	x[5]=pstr_var_CA1->mca[i];
	x[6]=pstr_var_CA1->hca[i];
	*/
	I += Icnst_E;
	pstr_paras->I=I;
	//check point gh0.1->0.4
	RG4(x,n_DOF,&t,pstr_paras,difffunc_CA1_E_Ih);
	pstr_var_CA1->v[i]  =x[0];
	pstr_var_CA1->n[i]  =x[1];
	pstr_var_CA1->ahp[i]=x[2];
	pstr_var_CA1->rf[i] =x[3];
	pstr_var_CA1->rs[i] =x[4];
	/* caspk
	pstr_var_CA1->mca[i]=x[5];
	pstr_var_CA1->hca[i]=x[6];
	*/
      }
      else if(i<NE_CA1+NPV){
	int n_DOF=3;
	double x[n_DOF];
	x[0]=pstr_var_CA1->v[i];
	x[1]=pstr_var_CA1->n[i];
	x[2]=pstr_var_CA1->h[i];
	I += Icnst_PV;
	pstr_paras->I=I;
	RG4(x,n_DOF,&t,pstr_paras,difffunc_CA1_PV);
	pstr_var_CA1->v[i]=x[0];
	pstr_var_CA1->n[i]=x[1];
	pstr_var_CA1->h[i]=x[2];
      }
      else{
	int n_DOF=7;
	double x[n_DOF];
	x[0]=pstr_var_CA1->v[i];
	x[1]=pstr_var_CA1->n[i];
	x[2]=pstr_var_CA1->h[i];
	x[3]=pstr_var_CA1->m[i];
	x[4]=pstr_var_CA1->a[i];
	x[5]=pstr_var_CA1->b[i];
	x[6]=pstr_var_CA1->r[i];
	I += Icnst_T;
	pstr_paras->I=I;
	RG4(x,n_DOF,&t,pstr_paras,difffunc_CA1_OLM);

	pstr_var_CA1->v[i]=x[0];
	pstr_var_CA1->n[i]=x[1];
	pstr_var_CA1->h[i]=x[2];
	pstr_var_CA1->m[i]=x[3];
	pstr_var_CA1->a[i]=x[4];
	pstr_var_CA1->b[i]=x[5];
	pstr_var_CA1->r[i]=x[6];
      }
    }
  

#endif


    // for decoder
    // 1 step ahead (neurons in Dec)
    for(int i = 0; i < N_Dec; i++){
      pstr_var_Dec->gE[i].pop_front();    pstr_var_Dec->gE[i].push_back(0.0);
      pstr_var_Dec->gE_CA1[i].pop_front();    pstr_var_Dec->gE_CA1[i].push_back(0.0);
      pstr_var_Dec->gI[i].pop_front();    pstr_var_Dec->gI[i].push_back(0.0);			

      double I=-1*(pstr_var_Dec->gE[i][0]*(pstr_var_Dec->v[i]-VE)
		   +pstr_var_Dec->gI[i][0]*(pstr_var_Dec->v[i]-VI));

      if(i < NE_Dec){
	//	int n_DOF=5;
	int n_DOF=4;
	double x[n_DOF];
	I+=Icnst_Dec;
	pstr_paras->I=I;
	
	x[0]=pstr_var_Dec->v[i];
	x[1]=pstr_var_Dec->n[i];
	x[2]=pstr_var_Dec->h[i];
	x[3]=pstr_var_Dec->m[i];
	//	x[4]=pstr_var_Dec->ahp[i];
	RG4(x,n_DOF,&t,pstr_paras,difffunc_Dec_E_pyr);
	pstr_var_Dec->v[i]=x[0];
	pstr_var_Dec->n[i]=x[1];
	pstr_var_Dec->h[i]=x[2];
	pstr_var_Dec->m[i]=x[3];
	//	pstr_var_Dec->ahp[i]=x[4];
      }
      else{
	int n_DOF=3;
	double x[n_DOF];
	pstr_paras->I=I;
	
	x[0]=pstr_var_Dec->v[i];
	x[1]=pstr_var_Dec->n[i];
	x[2]=pstr_var_Dec->h[i];
	RG4(x,n_DOF,&t,pstr_paras,difffunc_CA1_PV);
	pstr_var_Dec->v[i]=x[0];
	pstr_var_Dec->n[i]=x[1];
	pstr_var_Dec->h[i]=x[2];
      }
    }


    if((int)(t*50)%50==0){
      ofsf << t <<" "<<lfp<<" ";
#ifndef readmode
      for(int i = 0; i < NE_CA1; i+=10)
	ofsf << pstr_var_CA1->v[i] << " "<<pstr_var_CA1->gE[i][0] << " "<< pstr_var_CA1->gI[i][0] << " "<<pstr_var_CA1->gT[i][0] << " " << pstr_var_CA1->olmspts[i]<<" "<<pstr_var_CA1->olmspts2[i]<<" ";
      
      for(int i = 0; i < N_EC3; i+=10)
	ofsf << pstr_var_EC3->v[i] << " "<<pstr_var_EC3->gE[i][0] << " "<< pstr_var_EC3->gI[i][0] << " ";

      //      for(int i = 0; i < NE_EC5; i+=20)	ofsf << pstr_var_EC5->v[i] << " "<<pstr_var_EC5->ahp[i] << " "<< pstr_var_EC5->mCAN[i] << " "<<pstr_var_EC5->Ca[i]<<" "<<pstr_var_EC5->hrate_can[i]<<" "<<pstr_var_EC5->KM[i]<<" "<<pstr_var_EC5->Cal[i]<<" ";
      for(int i = 0; i < NE_EC5; i+=20)	ofsf << pstr_var_EC5->v[i] << " "<<pstr_var_EC5->Ca[i]<<" "<<pstr_var_EC5->hrate_can[i]<<" "<<pstr_var_EC5->mCAN[i]<<" ";
#endif
      for(int i = 0; i < N_Dec;  i+=5)  ofsf << pstr_var_Dec->v[i] << " "<<pstr_var_Dec->gE[i][0]<< " "<<pstr_var_Dec->gI[i][0]<<" "<<pstr_var_Dec->gE_CA1[i][0]<<" ";
      ofsf<<endl;
    }
    

    //check
    //    cout<<endl;
    
  } // end of t-loop
}



void simul(int ik,int it, double c_Ach_aft, double Ach_test){
  calc(ik,it,c_Ach_aft,Ach_test);
}

int main(int argc, char **argv){
  //  srand((unsigned int)time(NULL));
  int ik = 0,it=0;
  double c_Ach_aft=0.0,Ach_test=0.0;
  int n=1;

  if(argc > n) geta = atof(argv[n]);   n+=1;
  if(argc > n) Gpv_ach = atof(argv[n]);   n+=1;
  if(argc > n) Icnst_PV = atof(argv[n]);   n+=1;
  if(argc > n) Tec3off_strt = atof(argv[n]);   n+=1;
  if(argc > n) ik        = atoi(argv[n]);   n+=1; // seed of networks
  if(argc > n) it        = atoi(argv[n]);   n+=1; // seed of timeseries
  if(argc > n) c_Ach_aft = atof(argv[n]);   n+=1;
  if(argc > n) GEE_EC5   = atof(argv[n]);   n+=1;
  if(argc > n) gCAN_def  = atof(argv[n]);   n+=1;

  if(argc > n) GEE_CA1_EC5 = atof(argv[n]);   n+=1;
  if(argc > n) GEI_EC5 = atof(argv[n]);   n+=1;
  if(argc > n) GIE_EC5 = atof(argv[n]);   n+=1;
  if(argc > n) Ach_test = atof(argv[n]);   n+=1;
  if(argc > n) Gn_Dec_IE= atof(argv[n]);   n+=1;
  if(argc > n) Icnst_Dec= atof(argv[n]);   n+=1;
#ifdef readmode
  if(argc > n) Flg_gamma= atoi(argv[n]);   n+=1;

#endif

  if(osc_flg==1){
    if(argc > n) osc_scl = atof(argv[n]);   n+=1;
  }

  simul(ik,it,c_Ach_aft,Ach_test);
  return 0;
}


