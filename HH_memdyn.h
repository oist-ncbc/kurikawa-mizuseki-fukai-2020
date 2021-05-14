

void chn_cond(double *fr, double *tau,double V, string ntype){
  double alp,beta;

  if( ntype=="OLMa" || ntype=="OLMb" ||ntype=="OLMr"){
    if(ntype=="OLMa"){
      *fr =1/(1+exp(-(V+14)/16.6));
      *tau=5;
    }
    else if(ntype=="OLMb"){
      *fr =1/(1+exp((V+71)/7.3));
      *tau=1/( 0.000009/exp((V-26.0)/18.5) + 0.014/(0.2+exp(-(V+70)/11.0))   );
    }
    else if(ntype=="OLMr"){
      *fr=1/(1+exp((V+84)/10.2));
      *tau=1/( exp(-14.59-0.086*V) + exp(-1.87+0.0701*V)  );
    }

    return ;
  }


  if(ntype=="Em"){
    alp =0.32*(V+54)/(1-exp(-(V+54)/4.));
    beta=0.28*(V+27)/(exp((V+27)/5.)-1);
  }
  else if(ntype=="En"){
    alp =0.032*(V+52)/(1-exp(-(V+52)/5.));
    beta=0.5*exp(-(V+57)/40.);
  }
  else if(ntype=="PVm"){
    alp =0.1*(V+35)/(1-exp(-(V+35)/10.));
    beta=4*exp(-(V+60)/18.);
  }
  else if(ntype=="PVh"){
    alp =0.07*exp(-(V+58)/20.);
    beta=1/(exp(-0.1*(V+28))+1);
  }
  else if(ntype=="PVn"){
    alp =0.01*(V+34)/(1-exp(-(V+34)));
    beta=0.125*exp(-(V+44)/80.0);
  }
  else if(ntype=="OLMm"){
    alp =-0.1*(V+38)/(exp(-(V+38)/10)-1);
    beta=4*exp(-(V+65)/18);
  }
  else if(ntype=="OLMh"){
    alp =0.07*exp(-(V+63)/20);
    beta=1/(exp(-0.1*(V+33))+1);
  }
  else if(ntype=="OLMn"){
    alp =0.018*(V-25)/(1-exp(-(V-25)/25.0));
    beta=0.0036*(V-35)/(exp((V-35)/12.0)-1);
  }
  else if(ntype=="EC5p"){
    alp  =1/(0.15*(1+exp(-(V+38)/6.5)));
    beta =exp(-(V+38)/6.5)/(0.15*(1+exp(-(V+38)/6.5)));
  }
  else if(ntype=="EC5m"){
    alp =-0.1*(V+23)/(exp(-0.1*(V+23))-1);
    beta=4*exp(-(V+48)/18.);
  }
  else if(ntype=="EC5h"){
    alp =0.07*exp(-(V+37)/20.);
    beta=1/(exp(-0.1*(V+7))+1);
  }
  else if(ntype=="EC5n"){
    alp =-0.01*(V+27)/(exp(-0.1*(V+27))-1);
    beta=0.125*exp(-(V+37)/80.);
  }
  else if(ntype=="EC5_caL_q"){
    alp =-0.055*(V+27)/(exp(-(V+27)/3.8)-1);
    beta=0.94*exp(-(75+V)/17.0);
  }
  else if(ntype=="EC5_caL_r"){
    alp=0.000457*exp(-(V+13)/50.0);
    beta=0.0065/(exp(-(V+15)/28)+1);
  }
  else if(ntype=="EC5_pyr_m"){   // Middleton,PNAS,2008
    alp =-0.32*(V+54)/(exp(-0.25*(V+54))-1);
    beta= 0.28*(V+27)/(exp(0.2*(V+27))  -1);
  }
  else if(ntype=="EC5_pyr_n"){    // Middleton,PNAS,2008
    alp =-0.032*(V+52)/(exp(-0.2*(V+52))-1);
    beta= 0.5*exp(-0.025*(V+57));
  }
  else if(ntype=="EC5_pyr_h"){    // Middleton,PNAS,2008
    alp =0.128*exp(-1*(V+50)/18.0);
    beta=4.0/(exp(-0.2*(V+27))+1);
  }
  else if(ntype=="EC5_pyr_ahp"){    // Middleton,PNAS,2008
    alp =1/(exp(-0.1*(V+41))+1);
    beta=500.0/(3.3*exp(0.05*(V+41))+exp(-0.05*(V+41)));  //check point original value (400 [ms]), 400->600 2017.6.30 update)
  }

  
  (*fr)=alp/(alp+beta);
  (*tau)=1/(alp+beta);
  if(ntype=="EC5_pyr_ahp"){
    (*fr)=alp;
    (*tau)=beta;
  }
  
  if(ntype=="PVh" || ntype=="PVn") *tau=*tau*0.2;
  return;
}


void difffunc_EC3_E(double x[], double dx[], STR_PARA *paras){  // Rotstein,et al., 2006, JCN
  //x[0]: membrane potential,x[1]: channel function n, x[2]: h, x[3]: m,x[4]: p  x[4]:p, x[5]: rf, x[6]: rs
  static double gNa=52.0,gK=11.0,gL=0.5,gp=0.5,gh=1.5;
  static double VNa=55.0,VK=-90.0,VL=-65.0,Vh=-20;
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,hinf,pinf,rfinf,rsinf;
  double tm,tn,th,tp,trf,trs;
  double tmpI   =paras->I;
  
  chn_cond(&ninf,&tn,x[0],"EC5n");
  chn_cond(&hinf,&th,x[0],"EC5h");
  chn_cond(&minf,&tm,x[0],"EC5m");
  chn_cond(&pinf,&tp,x[0],"EC5p");

  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;

  dx[0] =gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gp*x[4]*(VNa-x[0])+gh*(0.65*x[5]+0.35*x[6])*(Vh-x[0])+gL*(VL-x[0])+tmpI;
  
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;
  dx[3]=(minf-x[3])/tm;
  dx[4]=(pinf-x[4])/tp;
  dx[5]=(rfinf-x[5])/trf;
  dx[6]=(rsinf-x[6])/trs;

    
  return;
}



// this is an identical function 
void difffunc_EC2_stell(double x[], double dx[], STR_PARA *paras){  // Rotstein,et al., 2006, JCN
  //x[0]: membrane potential,x[1]: channel function n, x[2]: h, x[3]: m,x[4]: p  x[4]:p, x[5]: rf, x[6]: rs
  static double gNa=52.0,gK=11.0,gL=0.5,gp=0.5,gh=0.4;
  static double VNa=55.0,VK=-90.0,VL=-65.0,Vh=-20;
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,hinf,pinf,rfinf,rsinf;
  double tm,tn,th,tp,trf,trs;
  double tmpI   =paras->I;
  
  chn_cond(&ninf,&tn,x[0],"EC5n");
  chn_cond(&hinf,&th,x[0],"EC5h");
  chn_cond(&minf,&tm,x[0],"EC5m");
  chn_cond(&pinf,&tp,x[0],"EC5p");

  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;

  dx[0] =gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gp*x[4]*(VNa-x[0])+gh*(0.65*x[5]+0.35*x[6])*(Vh-x[0])+gL*(VL-x[0])+tmpI;
  
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;
  dx[3]=(minf-x[3])/tm;
  dx[4]=(pinf-x[4])/tp;
  dx[5]=(rfinf-x[5])/trf;
  dx[6]=(rsinf-x[6])/trs;

    
  return;
}




void difffunc_EC3_E_pyr(double x[], double dx[], STR_PARA *paras){  // Rotstein,et al., 2006, JCN,  Saravanan, 2015,Hipp, Middleton, 2008,PNAS
  //x[0]: membrane potential,x[1]: channel function n, x[2]: h, x[3]: m,x[4]:ahp
  static double gNa=100.0,gK=80.0,gL=0.5,gAHP=0.3; //check point in old version, gh=0.25 (2017.6.30 update)   check point gahp=0.6 -> 0.6
  static double VNa=50.0,VK=-90.0,VL=-65.0,VAHP=-100;  //check
  //  static double C=1.0　[uF/cm^2], thus, inverse of system time scale is mS/uF = 1/ms
  double minf,ninf,hinf,ahpinf,rfinf,rsinf;
  double tm,tn,th,tahp,trf,trs;
  double tmpI    =paras->I;
  double tmpgCAN =paras->gCAN;
  chn_cond(&ninf,&tn,x[0],"EC5_pyr_n");
  chn_cond(&hinf,&th,x[0],"EC5_pyr_h");
  chn_cond(&minf,&tm,x[0],"EC5_pyr_m");
  chn_cond(&ahpinf,&tahp,x[0],"EC5_pyr_ahp");

  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;

  dx[0] =gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gAHP*x[4]*(VAHP-x[0])+gL*(VL-x[0])+tmpI;
  
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;
  dx[3]=(minf-x[3])/tm;
  dx[4]=(ahpinf-x[4])/tahp;
  dx[5]=(rfinf-x[5])/trf;
  dx[6]=(rsinf-x[6])/trs;
    
  return;
}


#if VER_PROG != 9

void difffunc_EC5_E_CAN(double x[], double dx[], STR_PARA *paras){  // Rotstein,et al., 2006, JCN (we modified sttelate cells in EC2 with decreasing AHP current (see Heys's review in 2012 ),  Saravanan, 2015,Hipp
  //x[0]: membrane potential,x[1]: channel function n, x[2]: h, x[3]: m,x[4]: p  x[4]:p, x[5]: rf, x[6]: rs
  static double gNa=52.0,gK=11.0,gL=0.5,gp=0.5,gh=1.0; // original in *6.cpp.bak: gh=1.5 check point 4/7  
  static double VNa=55.0,VK=-90.0,VL=-65.0,Vh=-20,VCAN=-20;  
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,hinf,pinf,rfinf,rsinf;
  double tm,tn,th,tp,trf,trs;
  double tmpI    =paras->I;
  double tmpgCAN =paras->gCAN;
  double tmpgate_CAN=paras->gate_CAN;
  chn_cond(&ninf,&tn,x[0],"EC5n");
  chn_cond(&hinf,&th,x[0],"EC5h");
  chn_cond(&minf,&tm,x[0],"EC5m");
  chn_cond(&pinf,&tp,x[0],"EC5p");

  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;


  dx[0] =gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gp*x[4]*(VNa-x[0])+gh*(0.65*x[5]+0.35*x[6])*(Vh-x[0])+gL*(VL-x[0])+tmpI
    +tmpgCAN*tmpgate_CAN*(VCAN-x[0]);
  
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;
  dx[3]=(minf-x[3])/tm;
  dx[4]=(pinf-x[4])/tp;
  dx[5]=(rfinf-x[5])/trf;
  dx[6]=(rsinf-x[6])/trs;

    
  return;
}


void difffunc_EC5_E_pyr_CAN(double x[], double dx[], STR_PARA *paras){  // the neuron model is based on pyramidal cell in EC2 (E cell in Middleton, 2008,PNAS) + modution on Ca dynamics. (See also Fransen 2006 and 2002. the former paper builds a model of EC5 neuron which is based on EC2 neuron (Fransen 2002).
  //x[0]: membrane potential,x[1]: channel function n, x[2]: h, x[3]: m,x[4]:ahp, x[5]: rf, x[6]: rs
  static double gNa=100.0,gK=80.0,gL=0.5,gAHP=0.2,gh=0.0;   // (modification on 6/19, gahp:0.4->0.2,gh:0.25) check point 0.1->0.2
  static double VNa=50.0,VK=-100.0,VL=-65.0,VAHP=-100,Vh=-20,VCAN=-20;
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,hinf,ahpinf,rfinf,rsinf;
  double tm,tn,th,tahp,trf,trs;
  double tmpI    =paras->I;
  double tmpgCAN =paras->gCAN;
  double tmpgate_CAN=paras->gate_CAN;
  chn_cond(&ninf,&tn,x[0],"EC5_pyr_n");
  chn_cond(&hinf,&th,x[0],"EC5_pyr_h");
  chn_cond(&minf,&tm,x[0],"EC5_pyr_m");
  chn_cond(&ahpinf,&tahp,x[0],"EC5_pyr_ahp");

  
  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;


  dx[0] =gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gAHP*x[4]*(VAHP-x[0])+gh*(0.65*x[5]+0.35*x[6])*(Vh-x[0])+gL*(VL-x[0])+tmpI
    +tmpgCAN*tmpgate_CAN*(VCAN-x[0]);
  
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;
  dx[3]=(minf-x[3])/tm;
  dx[4]=(ahpinf-x[4])/tahp;
  dx[5]=(rfinf-x[5])/trf;
  dx[6]=(rsinf-x[6])/trs;

    
  return;
}
#endif


// last update 2017.6.28
void difffunc_EC5_E_pyr_Ca_dyn(double x[], double dx[], STR_PARA *paras){  // the neuron model is based on pyramidal cell in EC2 (E cell in Middleton, 2008,PNAS) + Ca dynamics (Saravanan,2015). (See also Fransen 2006 and 2002. the former paper builds a model of EC5 neuron which is based on EC2 neuron (Fransen 2002).
  // CAN channel is based on Destexhe 1994
  // Ca dynamics is based on Destexhe 1994 and Saravanan 2015
  //x[0]: membrane potential,x[1]: channel function n, x[2]: h, x[3]: m,x[4]:ahp, x[5]: rf, x[6]: rs, x[7]:qca, x[8]:rca, x[9]:Ca concentration, x[10]:m_can
  static double gNa=100.0,gK=80.0,gL=0.5,gAHP=0.1,gh=0.0,gCaL=0.01;   // [mS/cm^2] (modification on 6/19, gahp:0.4->0.2,gh:0.25)
  static double VNa=50.0,VK=-100.0,VL=-65.0,VAHP=-100,VCAN=-20,VCaL=120;  //[mV]
  static double I2C=0.51819378374737; //-k/Fd current of Ca to concentration
  static double Cainf=0.24; //[uM]
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,hinf,ahpinf,rfinf,rsinf,alpha_qca,beta_qca,alpha_rca,beta_rca,mcaninf;
  double tm,tn,th,tahp,trf,trs,t2Ca,tmcan;
  double tmpI    =paras->I;
  double tmpgCAN =paras->gCAN;
  double ICaL=0.0;
  chn_cond(&ninf,&tn,x[0],"EC5_pyr_n");
  chn_cond(&hinf,&th,x[0],"EC5_pyr_h");
  chn_cond(&minf,&tm,x[0],"EC5_pyr_m");
  chn_cond(&ahpinf,&tahp,x[0],"EC5_pyr_ahp");

  
  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;
  alpha_qca=0.055*(-27-x[0])/(exp((-27-x[0])/3.8)-1);
  beta_qca =0.94*exp(-(75+x[0])/17.0);
  alpha_rca=0.000457*exp(-(x[0]+13)/50.0);
  beta_rca =0.0065/(exp(-(x[0]+15)/28.0)+1);
  t2Ca = 250.0; // [ms]
  mcaninf=20.0*x[9]*x[9]/(20.0*x[9]*x[9]+6);  //chenge the inflection point (2000 (0.02[ms-1]*10^6) -> 100)  check point 1 8->6
  tmcan   =1.0/(20.0*x[9]*x[9]*0.002+0.003); // change from original values  0.002 -> 0.004

  ICaL=gCaL*x[7]*x[7]*x[8]*(VCaL-x[0]);
  dx[0] =gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gAHP*x[4]*(VAHP-x[0])+gL*(VL-x[0])+tmpI
    +ICaL+tmpgCAN*x[10]*x[10]*(VCAN-x[0]);
  
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;
  dx[3]=(minf-x[3])/tm;
  dx[4]=(ahpinf-x[4])/tahp;
  dx[5]=(rfinf-x[5])/trf;
  dx[6]=(rsinf-x[6])/trs;
  dx[7]=alpha_qca*(1-x[7])-beta_qca*x[7];
  dx[8]=alpha_rca*(1-x[8])-beta_rca*x[8];
  dx[9]=I2C*ICaL+(Cainf-x[9])/t2Ca;
  dx[10]=(mcaninf-x[10])/tmcan;


  return;
}


void difffunc_EC5_E_pyr_Fransen(double x[], double dx[], STR_PARA *paras){  
  // the neuron model is based on pyramidal cell in EC2 (Fransen 2006)
  // original unit system: [s], [S/m^2], [F/m^2],[V] -> [ms],[mS/cm^2],[uF/cm^2],[mV]
  //x[0]: membrane potential,x[1]: channel function n, x[2]: h, x[3]: m,x[4]:ahp, x[5]: s for CaL
  //x[6]: x for M , x[7]: m for Na p, x[8]: h for Na p, x[9]: a for K A, x[10]:b for K A
  //x[11]: concentration of Ca
  //x[12]: c for CAN
  //x[13]: c for KC


  static double gNa=100.0,gK=80.0,gL=0.5,gAHP=0.05,gKC=196,gM=3.5,gCaL=0.15,gNap=0.2,gKA=0.5;   // [mS/cm^2] 
  // original values of conductances [S/m^2] for capacitance C=0.01 [F/m^2] (unit time=[s]), here mS/cm^2, uF/cm^2, thus we just multiply 0.1 to values of conductance and set capacitance 1 (unit time=[ms]).
  static double VNa=50.0,VK=-100.0,VL=-65.0,VCAN=-20,VCa=140;  //[mV]
  static double I2C=0.51819378374737; //-k/Fd current of Ca to concentration
  double tmpI    =paras->I;
  double tmpgCAN =paras->gCAN;
  double ICaL=0.0;
  double t2Ca=250;
  double s_inf, tau_s;
  double alpha_s, beta_s;
  double tmp_s;
  double dtmp;

  dx[0] =gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0]);   // Na channel
  dx[0]+=gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0]);    // delayed K channel
  dx[0]+=gAHP*x[4]*x[4]*(VK-x[0]);          // K ahp channel
  ICaL=gCaL*x[5]*x[5]*(VCa-x[0]);         // CaL channel
  dx[0]+=ICaL;
  dx[0]+=gM*x[6]*(VK-x[0]);                 // K M channel
  dx[0]+=gNap*x[7]*x[8]*(VNa-x[0]);         // Na persistent channel
  dx[0]+=gKA*x[9]*x[10]*(VK-x[0]);          // K A channel
  dx[0]+=tmpgCAN*x[12]*x[12]*(VCAN-x[0]);   // CAN channel
  dx[0]+=gKC*x[13]*(VK-x[0]);                // K C channel
   
  dx[0]+=tmpI;

  chn_cond(&s_inf,&tau_s,x[0],"EC5_pyr_n");
  dx[1]=(s_inf-x[1])/tau_s;
  chn_cond(&s_inf,&tau_s,x[0],"EC5_pyr_h");
  dx[2]=(s_inf-x[2])/tau_s;
  chn_cond(&s_inf,&tau_s,x[0],"EC5_pyr_m");
  dx[3]=(s_inf-x[3])/tau_s;

  // ahp
  if(x[11]*1000>15){
    dtmp=(3.0+0.8*(1000*x[11]-15.0));
    alpha_s= dtmp>15.0 ? 15.0 : dtmp; //check point x[11] -> x[11]*30
  }
  else            alpha_s= 1000*x[11]*0.2;
  beta_s=1.0;
  tmp_s=x[4];
  dx[4]=0.001*(alpha_s*(1-tmp_s)-beta_s*tmp_s);
  
  // CaL
  alpha_s=1.6/(1+exp(-0.072*(x[0]-65.0)));
  beta_s =0.02*(x[0]-51.1)/(exp((x[0]-51.1)*0.2)-1.0);
  tmp_s=x[5];
  dx[5]=alpha_s*(1-tmp_s)-beta_s*tmp_s;


  //KC
  if(x[0]<50.0){
    alpha_s=exp(0.053782*x[0]-0.66835)/18.975;
    beta_s =2.0*exp((6.5-x[0])/27.0)-alpha_s;
  }
  else{
    alpha_s =2*exp((6.5-x[0])/27.0);
    beta_s  =0.0;
  }
  tmp_s=x[13];
  dx[13]=alpha_s*(1-tmp_s)-beta_s*tmp_s;

  //KM
  s_inf=1/(1+exp(-(x[0]+35.0)/5.0));
  tau_s=1000.0/(3.3*exp((x[0]+35)/40.0)+exp(-(x[0]+35)/20.0));
  tmp_s=x[6];
  dx[6]=(s_inf-tmp_s)/tau_s;

  //Nap
  s_inf=1/(1+exp(-(x[0]+48.7)/4.4));
  dtmp=x[0]+38.0;
  tau_s=1/( (0.091*dtmp/(1-exp(-dtmp/5.0))) + (-0.062*dtmp/(1-exp(dtmp/5.0))) );
  tmp_s=x[7];
  dx[7]=(s_inf-tmp_s)/tau_s;
  
  s_inf=1/(1+exp(-(x[0]+48.8)/9.98));
  tau_s=1/( (-0.00000288*(x[0]-49.1)/(1-exp(-(x[0]-49.1)/4.63))) + (0.00000694*(x[0]+44.7)/(1-exp((x[0]+44.7)/2.63))) );
  tmp_s=x[8];
  dx[8]=(s_inf-tmp_s)/tau_s;

  //KA from Traub 1991
  alpha_s=0.02*(x[0]-13.1)/(1-exp((13.1-x[0])/10.0));
  beta_s =0.175*(x[0]-40.1)/(exp((x[0]-40.1)/10.0)-1);
  tmp_s=x[9];
  dx[9]=alpha_s*(1-tmp_s)-beta_s*tmp_s;

  alpha_s=0.0016*exp(-(x[0]+13)/18.0);
  beta_s =0.05/(1+exp(10.1-x[0])/5.0);
  tmp_s=x[10];
  dx[10]=alpha_s*(1-tmp_s)-beta_s*tmp_s;


  // Ca dyn
  dx[11]=I2C*ICaL-x[11]/t2Ca;

  // CAN  check point 
  s_inf=48.0*100*x[11]*x[11]/(48.0*100*x[11]*x[11]+0.03);
  tau_s=1.0/(48.0*100*x[11]*x[11]+0.03); 
  tmp_s=x[12];
  dx[12]=(s_inf-tmp_s)/tau_s;
  
}


// last update 2017.7.25
void difffunc_EC5_E_pyr_Ca_dyn_M(double x[], double dx[], STR_PARA *paras){  // the neuron model is based on pyramidal cell in EC2 (E cell in Middleton, 2008,PNAS) + Ca dynamics (Saravanan,2015). (See also Fransen 2006 and 2002. the former paper builds a model of EC5 neuron which is based on EC2 neuron (Fransen 2002).
  // Ca-dependent K (AHP) is introduced based on Destexhe 1994
  // CAN channel is based on Destexhe 1994
  // Ca dynamics is based on Destexhe 1994 and Saravanan 2015
  // M current is based on  Saravanan, 2015 (this is model for CA1, but ECV model in Fransen 2006 used CA3 model in Traub 1992.)
  //x[0]: membrane potential,x[1]: channel function n, x[2]: h, x[3]: m,x[4]:ahp, x[5]: rf, x[6]: rs, x[7]:qca, x[8]:rca, x[9]:Ca concentration, x[10]:m_can
  static double gNa=100.0,gK=80.0,gL=0.5,gAHP=0.4,gM=0.004,gCaL=0.01;   // [mS/cm^2] (modification on 6/19, gahp:0.4->0.2)
  static double VNa=50.0,VK=-100.0,VL=-65.0,VAHP=-100,VCAN=-20,VCaL=120;  //[mV]

  static double I2C=0.51819378374737; //-k/Fd current of Ca to concentration
  static double Cainf=0.24; //[uM]
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,hinf,ahpinf,rfinf,rsinf,alpha_qca,beta_qca,alpha_rca,beta_rca,mcaninf,Minf;
  double tm,tn,th,tahp,trf,trs,t2Ca,tmcan,tM;
  double tmpI    =paras->I;
  double tmpgCAN =paras->gCAN;
  double ICaL=0.0;
  chn_cond(&ninf,&tn,x[0],"EC5_pyr_n");
  chn_cond(&hinf,&th,x[0],"EC5_pyr_h");
  chn_cond(&minf,&tm,x[0],"EC5_pyr_m");


  
  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;
  alpha_qca=0.055*(-27-x[0])/(exp((-27-x[0])/3.8)-1);
  beta_qca =0.94*exp(-(75+x[0])/17.0);
  alpha_rca=0.000457*exp(-(x[0]+13)/50.0);
  beta_rca =0.0065/(exp(-(x[0]+15)/28.0)+1);

  t2Ca = 250.0; // [ms]
  mcaninf=20.0*x[9]*x[9]/(20.0*x[9]*x[9]+6);  //chenge the inflection point (2000 (0.02[ms-1]*10^6) -> 100)  check point 1 8->6
  tmcan   =1.0/(20.0*x[9]*x[9]*0.002+0.003); // change from original values  0.002 -> 0.004

  Minf=1.0/(1.0+exp(-(x[0]+35)/10.0));
  tM  =1000/(3.3*exp((x[0]+35)/20.0)+exp(-(x[0]+35)/20.0));

  ahpinf =48.0*x[9]*x[9]/(48.0*x[9]*x[9]+90);
  tahp   =1.0/(48.0*x[9]*x[9]*0.002+0.015); // change from original values  0.002 -> 0.004

  ICaL=gCaL*x[7]*x[7]*x[8]*(VCaL-x[0]);
  dx[0] =gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gAHP*x[4]*x[4]*(VAHP-x[0])+gL*(VL-x[0])+tmpI
    +ICaL+tmpgCAN*x[10]*x[10]*(VCAN-x[0])+gM*x[11]*(VK-x[0]);
  
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;
  dx[3]=(minf-x[3])/tm;
  dx[4]=(ahpinf-x[4])/tahp;
  dx[5]=(rfinf-x[5])/trf;
  dx[6]=(rsinf-x[6])/trs;
  dx[7]=alpha_qca*(1-x[7])-beta_qca*x[7];
  dx[8]=alpha_rca*(1-x[8])-beta_rca*x[8];
  dx[9]=I2C*ICaL+(Cainf-x[9])/t2Ca;
  dx[10]=(mcaninf-x[10])/tmcan;
  dx[11]=(Minf-x[11])/tM;

  return;
}




void difffunc_EC5_E_pyr_Ca_dyn1(double x[], double dx[], STR_PARA *paras){  // the neuron model is based on pyramidal cell in EC2 (E cell in Middleton, 2008,PNAS) + Ca dynamics (Saravanan,2015). (See also Fransen 2006 and 2002. the former paper builds a model of EC5 neuron which is based on EC2 neuron (Fransen 2002).
  // Cal-dependent Kahp is introduced instead of normal AHP (Destexhe 1994)
  // CAN channel is based on Destexhe 1994
  // Ca dynamics is based on Destexhe 1994 and Saravanan 2015
  //x[0]: membrane potential,x[1]: channel function n, x[2]: h, x[3]: m,x[4]:ahp, x[5]: rf, x[6]: rs, x[7]:qca, x[8]:rca, x[9]:Ca concentration, x[10]:m_can
  static double gNa=100.0,gK=80.0,gL=0.5,gAHP=0.1,gh=0.0,gCaL=0.01;   // [mS/cm^2] (modification on 6/19, gahp:0.4->0.2,gh:0.25)
  static double VNa=50.0,VK=-100.0,VL=-65.0,VCAN=-20,VCaL=120;  //[mV]
  static double I2C=0.51819378374737; //-k/Fd current of Ca to concentration
  static double Cainf=0.24; //[uM]
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,hinf,ahpinf,rfinf,rsinf,alpha_qca,beta_qca,alpha_rca,beta_rca,mcaninf;
  double tm,tn,th,tahp,trf,trs,t2Ca,tmcan;
  double tmpI    =paras->I;
  double tmpgCAN =paras->gCAN;
  double ICaL=0.0;
  chn_cond(&ninf,&tn,x[0],"EC5_pyr_n");
  chn_cond(&hinf,&th,x[0],"EC5_pyr_h");
  chn_cond(&minf,&tm,x[0],"EC5_pyr_m");

  
  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;
  alpha_qca=0.055*(-27-x[0])/(exp((-27-x[0])/3.8)-1);
  beta_qca =0.94*exp(-(75+x[0])/17.0);
  alpha_rca=0.000457*exp(-(x[0]+13)/50.0);
  beta_rca =0.0065/(exp(-(x[0]+15)/28.0)+1);
  t2Ca = 250.0; // [ms]

  ahpinf =48.0*x[9]*x[9]/(48.0*x[9]*x[9]+90);
  tahp   =1.0/(48.0*x[9]*x[9]*0.002+0.015); // change from original values  0.002 -> 0.004
  mcaninf=20.0*x[9]*x[9]/(20.0*x[9]*x[9]+6);  //chenge the inflection point (2000 (0.02[ms-1]*10^6) -> 100)  check point 1 8->6
  tmcan   =1.0/(20.0*x[9]*x[9]*0.002+0.003); // change from original values  0.002 -> 0.004

  ICaL=gCaL*x[7]*x[7]*x[8]*(VCaL-x[0]);
  dx[0] =gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gAHP*x[4]*x[4]*(VK-x[0])+gL*(VL-x[0])+tmpI
    +ICaL+tmpgCAN*x[10]*x[10]*(VCAN-x[0]);
  
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;
  dx[3]=(minf-x[3])/tm;
  dx[4]=(ahpinf-x[4])/tahp;
  dx[5]=(rfinf-x[5])/trf;
  dx[6]=(rsinf-x[6])/trs;
  dx[7]=alpha_qca*(1-x[7])-beta_qca*x[7];
  dx[8]=alpha_rca*(1-x[8])-beta_rca*x[8];
  dx[9]=I2C*ICaL+(Cainf-x[9])/t2Ca;
  dx[10]=(mcaninf-x[10])/tmcan;


  return;
}


void difffunc_CA1_E(double x[], double dx[], STR_PARA *paras){
  //x[0]: membrane potential,x[1]: channel function n
  static double gNa=100.0,gK=80.0,gL=0.1;
  static double VNa=50.0,VK=-100.0,VL=-67.0;
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf;
  double tn;
  double hfrc;
  double dummy;
  double dtmp;
  double tmpI   =paras->I;
  
  chn_cond(&minf,&dummy,x[0],"Em");
  chn_cond(&ninf,&tn,x[0],"En");
  dtmp=(1-1.25*x[1]);
  hfrc=( dtmp> 0 ? dtmp : 0);

  dx[0]=gNa*minf*minf*minf*hfrc*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gL*(VL-x[0])+tmpI;
  dx[1]=(ninf-x[1])/tn;

  return;
}


void difffunc_CA1_E_Ih(double x[], double dx[], STR_PARA *paras){
  //x[0]: membrane potential,x[1]: channel function n, x[2]: ahp, x[3]:fast , x[4]: slow
  static double gNa=100.0,gK=80.0,gL=0.1,gAHP=0.2,gh=0.1; //check point gh0.1->0.2->0.1
  static double VNa=50.0,VK=-100.0,VL=-67.0,VAHP=-100,Vh=-20;
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,ahpinf,rfinf,rsinf;
  double tn,trf,trs,tahp;
  double hfrc;
  double dummy;
  double dtmp;
  double tmpI   =paras->I;
  
  chn_cond(&minf,&dummy,x[0],"Em");
  chn_cond(&ninf,&tn,x[0],"En");
  chn_cond(&ahpinf,&tahp,x[0],"EC5_pyr_ahp");
  dtmp=(1-1.25*x[1]);
  hfrc=( dtmp> 0 ? dtmp : 0);

  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;

 
  dx[0]=gNa*minf*minf*minf*hfrc*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gL*(VL-x[0])+gAHP*x[2]*(VAHP-x[0])+gh*(0.65*x[3]+0.35*x[4])*(Vh-x[0])+gL*(VL-x[0])+tmpI;
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(ahpinf-x[2])/tahp;
  dx[3]=(rfinf-x[3])/trf;
  dx[4]=(rsinf-x[4])/trs;

  return;
}



void difffunc_CA1_E_Ih_Caspk(double x[], double dx[], STR_PARA *paras){  // cal. spike dynamics is in Chua, 2015,front.
  //x[0]: membrane potential,x[1]: channel function n, x[2]: ahp, x[3]:fast , x[4]: slow
  static double gNa=100.0,gK=80.0,gL=0.1,gAHP=0.2,gh=0.1;
  static double VNa=50.0,VK=-100.0,VL=-67.0,VAHP=-100,Vh=-20;
  static double gca=1.0,Vca=30.0,tmca=5.0,thca=50.0; // in the paper, gca = 70 nS on a dendritic part. here, we consider point neuron(presumably soma) and set gca=1.
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,ahpinf,rfinf,rsinf,mcainf,hcainf;
  double tn,trf,trs,tahp;
  double hfrc;
  double dummy;
  double dtmp;
  double tmpI   =paras->I;
  
  chn_cond(&minf,&dummy,x[0],"Em");
  chn_cond(&ninf,&tn,x[0],"En");
  chn_cond(&ahpinf,&tahp,x[0],"EC5_pyr_ahp");
  dtmp=(1-1.25*x[1]);
  hfrc=( dtmp> 0 ? dtmp : 0);

  rfinf=1/(1+exp((x[0]+79.2)/9.78));
  trf  =0.51/(exp((x[0]-1.7)*0.1) + exp(-(x[0]+340)/52.0)) + 1;
  rsinf=1/pow((1+exp((x[0]+2.83)/15.9)),58);
  trs  =5.6/(exp((x[0]-1.7)/14) + exp(-(x[0]+260)/43.0)) + 1;

  mcainf=1/(1+exp(0.5*(x[0]+18.0)));  // original 21
  hcainf=1/(1+exp(-0.5*(x[0]+21.0))); // original 24

  
  dx[0]=gNa*minf*minf*minf*hfrc*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gL*(VL-x[0])+gAHP*x[2]*(VAHP-x[0])+gh*(0.65*x[3]+0.35*x[4])*(Vh-x[0])+gL*(VL-x[0])+gca*x[5]*x[6]*(Vca-x[0])+tmpI;
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(ahpinf-x[2])/tahp;
  dx[3]=(rfinf-x[3])/trf;
  dx[4]=(rsinf-x[4])/trs;
  dx[5]=(mcainf-x[5])/tmca;
  dx[6]=(hcainf-x[6])/thca;

  return;
}





void difffunc_CA1_OLM(double x[], double dx[], STR_PARA *paras){
  //x[0]: membrane potential,x[1]: channel function n,x[2]: channel function h, x[3]: channel function m,
  //x[4]: channel function a,x[5]: channel function b,x[6]: channel function r,

  static double gNa=40.0,gK=23.0,gL=0.05,gA=16.0,gh=6.0;  // gNa:30->40
  static double VNa=60.0,VK=-100.0,VL=-70.0,VA=-90.0,Vh=-32.0;
  static double Cinv=1/1.3;
  double minf,ninf,hinf,ainf,binf,rinf;
  double tm,tn,th,ta,tb,tr;
  double dtmp;
  double tmpI   =paras->I;
  
  chn_cond(&minf,&tm,x[0],"OLMm");
  chn_cond(&ninf,&tn,x[0],"OLMn");
  chn_cond(&hinf,&th,x[0],"OLMh");
  chn_cond(&ainf,&ta,x[0],"OLMa");
  chn_cond(&binf,&tb,x[0],"OLMb");
  chn_cond(&rinf,&tr,x[0],"OLMr");

  dx[0]=Cinv*(gNa*x[3]*x[3]*x[3]*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])
	      +gA*x[4]*x[5]*(VA-x[0])+gh*x[6]*(Vh-x[0])+gL*(VL-x[0])+tmpI);

  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;
  dx[3]=(minf-x[3])/tm;
  dx[4]=(ainf-x[4])/ta;
  dx[5]=(binf-x[5])/tb;
  dx[6]=(rinf-x[6])/tr;

  return;
}





void difffunc_CA1_PV(double x[], double dx[], STR_PARA *paras){
  //x[0]: membrane potential,x[1]: channel function n,x[2]: channel function h
  static double gNa=35.0,gK=9.0,gL=0.1;
  static double VNa=55.0,VK=-90.0,VL=-65.0;
  //  static double C=1.0;  C is unity here, so we neglect this term.
  double minf,ninf,hinf;
  double tn,th;
  double hfrc;
  double dummy;
  double dtmp;
  double tmpI   =paras->I;
  
  chn_cond(&minf,&dummy,x[0],"PVm");
  chn_cond(&ninf,&tn,x[0],"PVn");
  chn_cond(&hinf,&th,x[0],"PVh");

  dx[0]=gNa*minf*minf*minf*x[2]*(VNa-x[0])+gK*x[1]*x[1]*x[1]*x[1]*(VK-x[0])+gL*(VL-x[0])+tmpI;
  dx[1]=(ninf-x[1])/tn;
  dx[2]=(hinf-x[2])/th;

  return;
}

void RG4(double x[],int n_DOF, double *tp, STR_PARA *paras,void (*p_diff_func)(double xtmp[],double dx[],STR_PARA *paras)){
  int i0;
  double htmp,hh,h6;
  double dxm[n_DOF],dxtmp[n_DOF],xtmp[n_DOF],dxdt[n_DOF];

  htmp=hstp;
  hh=0.5*htmp;
  h6=htmp/6.0;


  (*p_diff_func)(x,dxdt,paras);                       // 1st step

  for(i0=0;i0<n_DOF;i0++)xtmp[i0]=x[i0]+hh*dxdt[i0];  
  (*p_diff_func)(xtmp,dxtmp,paras);                        // 2nd step

  for(i0=0;i0<n_DOF;i0++)xtmp[i0]=x[i0]+hh*dxtmp[i0];  
  (*p_diff_func)(xtmp,dxm,paras);                        // 3rd step

  for(i0=0;i0<n_DOF;i0++){
    xtmp[i0]=x[i0]+htmp*dxm[i0];  
    dxm[i0]+=dxtmp[i0];
  }
  (*p_diff_func)(xtmp,dxtmp,paras);                        // 4th step


  for(i0=0;i0<n_DOF;i0++) x[i0]+=h6*(dxdt[i0]+2.0*dxm[i0]+dxtmp[i0]);

  return;
}
