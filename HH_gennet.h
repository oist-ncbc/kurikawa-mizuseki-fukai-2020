#define rhebb_c 6  // 
#define rhebb_g 12
#define add_EC5 150
#define add_lrn 2
#ifdef PVinact
#define add_lrn1 2.0
#define add_lrn2 2.0
#define rhebb_gCA1_EC5 2 //check point
#define rhebb_g1 12.0 //check point 6->12
#endif


//vector< vector<double> > 
DBLMAT calc_G_CA1(string n_area){
  if(n_area!="CA1")    printf("ERROR calc_G_CA!\n");
 
  DBLMAT G;  //vector< vector<double> > 
  // * -> E   //
  for(int i = 0; i < NE_CA1; i++){   
    G.push_back(dvec);
    for(int j = 0; j < N_CA1; j++)G[i].push_back(0.0);

    
    for(int j = 0;j<NE_CA1;j++){
      if( cEE_CA1<dbl_eps) G[i][j]=0.0;
      else if(dice()<cEE_CA1 )                  G[i][j] = GEE_CA1/(NE_CA1*cEE_CA1);
    }

    for(int j = NE_CA1; j < NPV+NE_CA1; j++){
      if( cEI_CA1<dbl_eps)                      G[i][j] =0;
      else if( dice()<cEI_CA1 )                 G[i][j] = GEI_CA1/(NPV*cEI_CA1);
    }
    for(int j = NE_CA1+NPV; j < NOLM+NPV+NE_CA1; j++){
      if( cET_CA1<dbl_eps)                                G[i][j] =0;
      else if( dice()<cET_CA1 )                           G[i][j] = GET_CA1/(NOLM*cET_CA1);
    }
  }

 
  // *->PV (or I)  //
  for(int i = NE_CA1; i < NPV+NE_CA1; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_CA1; j++) G[i].push_back(0.0);
    for(int j = 0; j < NE_CA1; j++){
      if( cIE_CA1<dbl_eps)                      G[i][j] =0;
      else if( dice()<cIE_CA1 )                 G[i][j] = GIE_CA1/(NE_CA1*cIE_CA1);
    }
    for(int j = NE_CA1; j < NPV+NE_CA1; j++){
      if( cII_CA1<dbl_eps)          G[i][j] =0;
      else if( dice()<cII_CA1)      G[i][j] = GII_CA1/(NPV*cII_CA1);
    }
    for(int j = NE_CA1+NPV; j < NOLM+NPV+NE_CA1; j++){
      if( cIT_CA1<dbl_eps)                      G[i][j] =0;
      else if( dice()<cIT_CA1 )                 G[i][j] = GIT_CA1/(NOLM*cIT_CA1);
    }
  }

  // *->OLM (or T)  //
  for(int i = NE_CA1+NPV; i < NOLM+NPV+NE_CA1; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_CA1; j++) G[i].push_back(0.0);
    for(int j = 0; j < NE_CA1; j++){
      if( cTE_CA1<dbl_eps)                           G[i][j] =0;
      else if( dice()<cTE_CA1 )                      G[i][j] = GTE_CA1/(NE_CA1*cIT_CA1);
    }
    for(int j = NE_CA1; j < NPV+NE_CA1; j++){
      if( cTI_CA1<dbl_eps)                  G[i][j] =0;
      else if( dice()<cTI_CA1)            G[i][j] = GTI_CA1/(NPV*cTI_CA1);
    }
    for(int j = NE_CA1+NPV; j < NOLM+NPV+NE_CA1; j++){
      if( cTT_CA1<dbl_eps)                  G[i][j] =0;
      else if( dice()<cTT_CA1 )             G[i][j] = GTT_CA1/(NOLM*cTT_CA1);
    }
  }
  return G;
}




DBLMAT calc_G3_CA3_CA1(string n_pre_area, string n_post_area,string ver){
  char pos; // index of the cluster

  if(n_pre_area!="CA3" || n_post_area!="CA1"){
    cout<<"invalid nema!!";
    exit(1);
  }
  if( ver != "EC3_to_PV" && ver !="no_EC3toPV"){  // this case corresponding to ver is not no_EC3toPV nor EC3_to_PV
    cout<<"Error! in generation of EC3toCA1 network";
  }

  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_CA1; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_CA3; j++){
      G[i].push_back(0.0);
      if(  (i<NE_CA1_L && j<NE_CA3_L) || ((i>=NE_CA1_L && i<NE_CA1_R ) && (j>=NE_CA3_L && j<NE_CA3_R))
	   || ((i>=NE_CA1_R && i<NE_CA1_other ) && (j>=NE_CA3_R && j<NE_CA3_other))
	   || (i>=NE_CA1_other && j>=NE_CA3_other)
	   )  pos='l'; // learned connections btw. corres. clsts.
      else pos='e';  // other networks
      if(pos=='l' && dice()<cEE_CA3_CA1)G[i][j] = GE_CA3_CA1/NE_CA3;
      G[i][j]*=(dice()*2);
    }
  }
  for(int i = NE_CA1; i < NPV+NE_CA1; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_CA3; j++){
      G[i].push_back(0.0);
      if(ver == "no_EC3toPV" && dice()<cIE_CA3_CA1)G[i][j] = GI_CA3_CA1/NPV;
      G[i][j]*=(dice()*2);
    }
  }
  
  return G;
}


DBLMAT calc_G_CA3_CA1(string n_pre_area, string n_post_area,string ver){
  if(n_pre_area!="CA3" || n_post_area!="CA1"){
    cout<<"invalid nema!!";
    exit(1);
  }
  if( ver != "EC3_to_PV" && ver !="no_EC3toPV"){  // this case corresponding to ver is not no_EC3toPV nor EC3_to_PV
    cout<<"Error! in generation of EC3toCA1 network";
  }

  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_CA1; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < N_CA3; j++){
      G[i].push_back(0.0);
      if(dice()<cEE_CA3_CA1)G[i][j] = GE_CA3_CA1/NE_CA3;
      G[i][j]*=(dice()*2);
    }
  
  }
  for(int i = NE_CA1; i < NPV+NE_CA1; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_CA3; j++){
      G[i].push_back(0.0);
      if(ver == "no_EC3toPV" && dice()<cIE_CA3_CA1)G[i][j] = GI_CA3_CA1/NPV;
      G[i][j]*=(dice()*2);
    }
  }
  
  return G;
}



#ifdef NE_EC3_other
//vector< vector<double> > 
DBLMAT calc_G(string n_area){
  int Ntmp,NEtmp;
  double cEE_tmp,cEI_tmp,cIE_tmp,cII_tmp;
  double GEE_tmp,GEI_tmp,GIE_tmp,GII_tmp;
  
  if(n_area=="EC3"){
    Ntmp =N_EC3; NEtmp=NE_EC3;
    cEE_tmp=cEE_EC3;    cEI_tmp=cEI_EC3;
    cIE_tmp=cIE_EC3;    cII_tmp=cII_EC3;
    GEE_tmp=GEE_EC3;    GEI_tmp=GEI_EC3;
    GIE_tmp=GIE_EC3;    GII_tmp=GII_EC3;
  }
  else if(n_area=="EC5"){
    Ntmp =N_EC5; NEtmp=NE_EC5;
    cEE_tmp=cEE_EC5;    cEI_tmp=cEI_EC5;
    cIE_tmp=cIE_EC5;    cII_tmp=cII_EC5;
    GEE_tmp=GEE_EC5;    GEI_tmp=GEI_EC5;
    GIE_tmp=GIE_EC5;    GII_tmp=GII_EC5;
  }
  
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < Ntmp; i++)     G.push_back(dvec);
  for(int i = 0; i < NEtmp; i++){
    for(int j = 0; j < NEtmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cEE_tmp ) G[i][j] = GEE_tmp/(NEtmp*cEE_tmp);
    }
    for(int j = NEtmp; j < Ntmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cEI_tmp ) G[i][j] = GEI_tmp/((Ntmp-NEtmp)*cEI_tmp);
    }
  }
  for(int i = NEtmp; i < Ntmp; i++){
    for(int j = 0; j < NEtmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cIE_tmp ) G[i][j] = GIE_tmp/(NEtmp*cIE_tmp);
    }
    for(int j = NEtmp; j < Ntmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cII_tmp ) G[i][j] = GII_tmp/((Ntmp-NEtmp)*cII_tmp);
    }
  }
  return G;
}


// introduce strong recurrent networks within a group of neurons
DBLMAT calc_G1(string n_area){   
  int Ntmp,NEtmp,NC,NL,Nother;
  double cEE_tmp,cEI_tmp,cIE_tmp,cII_tmp,cEEin_tmp;
  double GEE_tmp,GEI_tmp,GIE_tmp,GII_tmp,GEEin_tmp;

  if(n_area=="EC3"){
    Ntmp =N_EC3; NEtmp=NE_EC3;
    NC    =NE_EC3_L;
    NL    =NE_EC3_other;
    Nother=NE_EC3;
    cEE_tmp=cEE_EC3;    cEI_tmp=cEI_EC3;
    cEEin_tmp=cEE_EC3*rhebb_c;
    cIE_tmp=cIE_EC3;    cII_tmp=cII_EC3;
    GEE_tmp=GEE_EC3;        GEI_tmp=GEI_EC3;
    GEEin_tmp=GEE_EC3*rhebb_g;
    GIE_tmp=GIE_EC3;        GII_tmp=GII_EC3;
  }
  else if(n_area=="EC5"){
    Ntmp =N_EC5; NEtmp=NE_EC5;
    NC    =NE_EC5_L;
    NL    =NE_EC5_other;
    Nother=NE_EC5;
    cEE_tmp=cEE_EC5;    cEI_tmp=cEI_EC5;
    cEEin_tmp=cEE_EC5*rhebb_c;
    cIE_tmp=cIE_EC5;    cII_tmp=cII_EC5;
    GEE_tmp=GEE_EC5;        GEI_tmp=GEI_EC5;
    GEEin_tmp=GEE_EC5*rhebb_g*add_EC5;
    GIE_tmp=GIE_EC5;        GII_tmp=GII_EC5;
  }

  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < Ntmp; i++)     G.push_back(dvec);
  for(int i = 0; i < NEtmp; i++){
    for(int j = 0; j < NEtmp; j++){
      G[i].push_back(0.0);
      if( (i<NC && j<NC) || ((i>=NC && i<NL) && (j>=NC && j<NL)) || ((i>=NL && i<Nother)&&(j>=NL && j<Nother))){
	if( i!= j && dice() < cEEin_tmp)  G[i][j] = GEEin_tmp/(NEtmp*cEE_tmp);
#ifndef PVinact
	if(i<NC)   G[i][j]=add_lrn*G[i][j];
#else   
	if(i>=NC && i<NL) G[i][j]=add_lrn1*G[i][j];
	if(i<NC)           G[i][j]=add_lrn*G[i][j];
#endif
      }
      else{
	if( i!= j && dice() < cEE_tmp)    G[i][j] = GEE_tmp/(NEtmp*cEE_tmp);
      }
    }
    for(int j = NEtmp; j < Ntmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cEI_tmp ) G[i][j] = GEI_tmp/((Ntmp-NEtmp)*cEI_tmp);
    }
  }
  for(int i = NEtmp; i < Ntmp; i++){
    for(int j = 0; j < NEtmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cIE_tmp ) G[i][j] = GIE_tmp/(NEtmp*cIE_tmp);
    }
    for(int j = NEtmp; j < Ntmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cII_tmp ) G[i][j] = GII_tmp/((Ntmp-NEtmp)*cII_tmp);
    }
  }
  return G;
}








// last update 2017.6.16:  introduce ver as an argument
DBLMAT calc_G_EC3_CA1(string n_pre_area, string n_post_area,string ver){
  
  if(n_pre_area!="EC3" || n_post_area!="CA1"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_CA1; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < N_EC3; j++){
      G[i].push_back(0.0);
      if( cEE_EC3_CA1<dbl_eps)               G[i][j] =0;
      else if( j<NE_EC3 && dice() < cEE_EC3_CA1 ) G[i][j] = GEE_EC3_CA1/(NE_EC3*cEE_EC3_CA1);
    }
  }
  if( ver == "EC3_to_PV" ){
    for(int i = NE_CA1; i < NPV+NE_CA1; i++){
      G.push_back(dvec);
      for(int j = 0; j < N_EC3; j++){
	G[i].push_back(0.0);
	if( j<NE_EC3 && dice() < cIE_EC3_CA1 ) G[i][j] = GIE_EC3_CA1/(NE_EC3*cIE_EC3_CA1);
      }
    }
  }
  else if( ver == "no_EC3toPV"){
    for(int i = NE_CA1; i < NPV+NE_CA1; i++){
      G.push_back(dvec);
      for(int j = 0; j < N_EC3; j++)
	G[i].push_back(0.0);
    }
  }
  else{
    cout<<"Error! in generation of EC3toCA1 network";
  }

  for(int i = NE_CA1+NPV; i < N_CA1; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_EC3; j++)
      G[i].push_back(0.0);
  }
  return G;
}



DBLMAT calc_G_CA1_EC5(string n_pre_area, string n_post_area){
  
  if(n_pre_area!="CA1" || n_post_area!="EC5"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC5; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_CA1; j++){
      G[i].push_back(0.0);
      if(i<NE_EC5_L                      && j<NE_CA1_L                      && dice()<cEE_CA1_EC5 ) G[i][j] = GEE_CA1_EC5/((NE_CA1_L)*cEE_CA1_EC5);
      if((i>=NE_EC5_L && i<NE_EC5_other) && (j>=NE_CA1_L && j<NE_CA1_other) && dice()<cEE_CA1_EC5 ) G[i][j] = GEE_CA1_EC5/((NE_CA1_other-NE_CA1_L)*cEE_CA1_EC5);
      if(i>=NE_EC5_other                 && j>=NE_CA1_other                 && dice()<cEE_CA1_EC5 ) G[i][j] = GEE_CA1_EC5/((NE_CA1-NE_CA1_other)*cEE_CA1_EC5);
    }
    for(int j = NE_CA1; j < N_CA1; j++)
      G[i].push_back(0.0);
  }
  for(int i = NE_EC5; i < N_EC5; i++){
    G.push_back(dvec);

#ifndef cIE_CA1_EC5
    for(int j = 0; j < N_CA1; j++)
      G[i].push_back(0.0);
#else
    for(int j = 0; j < NE_CA1; j++){
      G[i].push_back(0.0);
      if(dice()<cIE_CA1_EC5) G[i][j] = GIE_CA1_EC5/((NE_CA1)*cIE_CA1_EC5);
    }
    for(int j = NE_CA1; j < N_CA1; j++)
      G[i].push_back(0.0);
#endif	
    
  }

  return G;
}

DBLMAT calc_G_EC5_EC3(string n_pre_area, string n_post_area){
  
  if(n_pre_area!="EC5" || n_post_area!="EC3"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC3; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_EC5; j++){
      G[i].push_back(0.0);
      if(i<NE_EC3_L                      && j<NE_EC5_L                      && dice()<cEE_EC5_EC3 ) G[i][j] = GEE_EC5_EC3/((NE_EC5_L)*cEE_EC5_EC3);
      if((i>=NE_EC3_L && i<NE_EC3_other) && (j>=NE_EC5_L && j<NE_EC5_other) && dice()<cEE_EC5_EC3 ) G[i][j] = GEE_EC5_EC3/((NE_EC5_other-NE_EC5_L)*cEE_EC5_EC3);
      if(i>=NE_EC3_other                 && j>=NE_EC5_other                 && dice()<cEE_EC5_EC3 ) G[i][j] = GEE_EC5_EC3/((NE_EC5_other-NE_EC5_L)*cEE_EC5_EC3);
    }
    for(int j = NE_EC5;j<N_EC5; j++)
      G[i].push_back(0.0);
  }
  for(int i = NE_EC3; i < N_EC3; i++){
    G.push_back(dvec);
#ifndef cIE_EC5_EC3
    for(int j = 0; j < N_EC5; j++)   G[i].push_back(0.0);
#else
    for(int j = 0; j < NE_EC5; j++){
      G[i].push_back(0.0);
      if(dice()<cIE_EC5_EC3) G[i][j] = GIE_EC5_EC3/((NE_EC5)*cIE_EC5_EC3);
    }
    for(int j = NE_EC5; j < N_EC5; j++)      G[i].push_back(0.0);
#endif
  }

  return G;
}




// subnetworks between two areas which are corresponding to the identical place are strongly connected.
// 2017.6.22  fix some incorrect descriptions
DBLMAT calc_G1_EC3_CA1(string n_pre_area, string n_post_area, string ver){
  
  if(n_pre_area!="EC3" || n_post_area!="CA1"){
    cout<<"invalid nema!!";
    exit(1);
  }

  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_CA1; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_EC3; j++){
      G[i].push_back(0.0);

      if(i<NE_CA1_L){
	if(j<NE_EC3_L){
	  if( dice()<rhebb_c*cEE_EC3_CA1 ) 	G[i][j] = GEE_EC3_CA1*rhebb_g*add_lrn/((NE_EC3_L)*cEE_EC3_CA1);
	}
	else{
	  if( dice()<cEE_EC3_CA1 )            G[i][j] = GEE_EC3_CA1/((NE_EC3-NE_EC3_L)*cEE_EC3_CA1);
	}
      }
      else if(i>=NE_CA1_L && i<NE_CA1_other){
	if(j>=NE_EC3_L && j<NE_EC3_other){
	  if( dice()<rhebb_c*cEE_EC3_CA1 ) 	G[i][j] = GEE_EC3_CA1*rhebb_g/((NE_EC3_other-NE_EC3_L)*cEE_EC3_CA1);
	}
	else{
	  if( dice()<cEE_EC3_CA1 )            G[i][j] = GEE_EC3_CA1/((NE_EC3-NE_EC3_L)*cEE_EC3_CA1);
	}
      }
      else{
	if(j>=NE_EC3_other){
	  if( dice()<rhebb_c*cEE_EC3_CA1 ) 	G[i][j] = GEE_EC3_CA1*rhebb_g/((NE_EC3-NE_EC3_other)*cEE_EC3_CA1);
	}
	else{
	  if( dice()<cEE_EC3_CA1 )            G[i][j] = GEE_EC3_CA1/((NE_EC3-NE_EC3_L)*cEE_EC3_CA1);
	}
      }
    }
    for(int j = NE_EC3; j < N_EC3; j++)      G[i].push_back(0.0);
  }
  if( ver == "EC3_to_PV" ){
    for(int i = NE_CA1; i < NPV+NE_CA1; i++){
      G.push_back(dvec);
      for(int j = 0; j < N_EC3; j++){
	G[i].push_back(0.0);
	if( j<NE_EC3 && dice() < cIE_EC3_CA1 ) G[i][j] = GIE_EC3_CA1/(NE_EC3*cIE_EC3_CA1);
	G[i][j]*=(dice()*2);
      }
    }
  }
  else if( ver == "no_EC3toPV"){
    for(int i = NE_CA1; i < NPV+NE_CA1; i++){
      G.push_back(dvec);
      for(int j = 0; j < N_EC3; j++)G[i].push_back(0.0);
    }
  }
  else{
    cout<<"Error! in generation of EC3toCA1 network";
  }

  for(int i = NE_CA1+NPV; i < N_CA1; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_EC3; j++) G[i].push_back(0.0);
  }
  return G;
}





// last update 2017.4.8
DBLMAT calc_G1_CA1_EC5(string n_pre_area, string n_post_area){
  
  if(n_pre_area!="CA1" || n_post_area!="EC5"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC5; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_CA1; j++){
      G[i].push_back(0.0);
      if(i<NE_EC5_L){
	if(j<NE_CA1_L){
	  if( dice()<rhebb_c*cEE_CA1_EC5 )      G[i][j] = GEE_CA1_EC5*rhebb_g*add_lrn/((NE_CA1_L)*cEE_CA1_EC5);
	}
	else{
	  if( dice()<cEE_CA1_EC5 )            G[i][j] = GEE_CA1_EC5/((NE_CA1-NE_CA1_L)*cEE_CA1_EC5);
	}
      }
      else if(i>=NE_EC5_L && i<NE_EC5_other){
	if(j>=NE_CA1_L && j<NE_CA1_other){
#ifndef PVinact 
	  if( dice()<rhebb_c*cEE_CA1_EC5 )      G[i][j] = GEE_CA1_EC5*rhebb_g/((NE_CA1_other-NE_CA1_L)*cEE_CA1_EC5);
#else
	  if( dice()<rhebb_c*cEE_CA1_EC5 )      G[i][j] = GEE_CA1_EC5*rhebb_g*add_lrn1/((NE_CA1_other-NE_CA1_L)*cEE_CA1_EC5);

#endif
	}
	else{
	  if( dice()<cEE_CA1_EC5 )            G[i][j] = GEE_CA1_EC5/((NE_CA1-NE_CA1_L)*cEE_CA1_EC5);
	}
      }
      else{
	if(j>=NE_CA1_other){
#ifndef PVinact
	  if( dice()<rhebb_c*cEE_CA1_EC5 )      G[i][j] = GEE_CA1_EC5*rhebb_g/((NE_CA1-NE_CA1_other)*cEE_CA1_EC5);
#else
	  if( dice()<rhebb_c*cEE_CA1_EC5 )      G[i][j] = GEE_CA1_EC5*rhebb_g/((NE_CA1-NE_CA1_other)*cEE_CA1_EC5)
						  ;
#endif
	}
	else{
	  if( dice()<cEE_CA1_EC5 )            G[i][j] = GEE_CA1_EC5/((NE_CA1-NE_CA1_L)*cEE_CA1_EC5);
	}
      }
    }
    for(int j = NE_CA1; j < N_CA1; j++)      G[i].push_back(0.0);
  }
  for(int i = NE_EC5; i < N_EC5; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_CA1; j++)           G[i].push_back(0.0);
  }

  return G;
}


DBLMAT calc_G1_EC5_EC3(string n_pre_area, string n_post_area){
  
  if(n_pre_area!="EC5" || n_post_area!="EC3"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC3; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_EC5; j++){
      G[i].push_back(0.0);
      if(i<NE_EC3_L){
	if(j<NE_EC5_L){
	  if( dice()<rhebb_c*cEE_EC5_EC3 )      G[i][j] = GEE_EC5_EC3*rhebb_g*add_lrn/((NE_EC5_L)*cEE_EC5_EC3);
	}
	else{
	  if( dice()<cEE_EC5_EC3 )            G[i][j] = GEE_EC5_EC3/((NE_EC5-NE_EC5_L)*cEE_EC5_EC3);
	}
      }
      else if(i>=NE_EC3_L && i<NE_EC3_other){
	if(j>=NE_EC5_L && j<NE_EC5_other){
	  if( dice()<rhebb_c*cEE_EC5_EC3 )      G[i][j] = GEE_EC5_EC3*rhebb_g/((NE_EC5_other-NE_EC5_L)*cEE_EC5_EC3);
	}
	else{
	  if( dice()<cEE_EC5_EC3 )            G[i][j] = GEE_EC5_EC3/((NE_EC5-NE_EC5_L)*cEE_EC5_EC3);
	}
      }
      else{
	if(j>=NE_EC5_other){
	  if( dice()<rhebb_c*cEE_EC5_EC3 )      G[i][j] = GEE_EC5_EC3*rhebb_g/((NE_EC5-NE_EC5_other)*cEE_EC5_EC3);
	}
	else{
	  if( dice()<cEE_EC5_EC3 )            G[i][j] = GEE_EC5_EC3/((NE_EC5-NE_EC5_L)*cEE_EC5_EC3);
	}
      }
    }
    for(int j = NE_EC5;j<N_EC5; j++)  G[i].push_back(0.0);
  }
  for(int i = NE_EC3; i < N_EC3; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_CA1; j++)      G[i].push_back(0.0);
  }

  return G;
}
















//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// introduce strong recurrent networks within a group of neurons
// introduce right arms
// exchange notation of NC and NL, because its meaning is different
DBLMAT calc_G2(string n_area){   
  int Ntmp,NEtmp,NC,NL,NR,Nother;
  double cEE_tmp,cEI_tmp,cIE_tmp,cII_tmp,cEEin_tmp;
  double GEE_tmp,GEI_tmp,GIE_tmp,GII_tmp,GEEin_tmp;
  
  if(n_area=="EC3"){
    Ntmp =N_EC3; NEtmp=NE_EC3;
    NL    =NE_EC3_L;
    NR    =NE_EC3_R;
    NC    =NE_EC3_other;
    Nother=NE_EC3;
    cEE_tmp=cEE_EC3;    cEI_tmp=cEI_EC3;
    cEEin_tmp=cEE_EC3*rhebb_c;
    cIE_tmp=cIE_EC3;    cII_tmp=cII_EC3;
    GEE_tmp=GEE_EC3;        GEI_tmp=GEI_EC3;
    GEEin_tmp=GEE_EC3*rhebb_g;
    GIE_tmp=GIE_EC3;        GII_tmp=GII_EC3;
  }
  else if(n_area=="EC5"){
    Ntmp =N_EC5; NEtmp=NE_EC5;
    NL    =NE_EC5_L;
    NR    =NE_EC5_R;
    NC    =NE_EC5_other;
    Nother=NE_EC5;
    cEE_tmp=cEE_EC5;    cEI_tmp=cEI_EC5;
    cEEin_tmp=cEE_EC5*rhebb_c;
    cIE_tmp=cIE_EC5;    cII_tmp=cII_EC5;
    GEE_tmp=GEE_EC5;        GEI_tmp=GEI_EC5;
    GEEin_tmp=GEE_EC5*rhebb_g*add_EC5;
    GIE_tmp=GIE_EC5;        GII_tmp=GII_EC5;
  }



  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < Ntmp; i++)     G.push_back(dvec);
  for(int i = 0; i < NEtmp; i++){
    for(int j = 0; j < NEtmp; j++){
      G[i].push_back(0.0);
      if( (i<NL && j<NL) || ((i>=NL && i<NR) && (j>=NL && j<NR)) ||((i>=NR && i<NC) && (j>=NR && j<NC)) ||  ((i>=NC && i<Nother)&&(j>=NC && j<Nother))){
	if( i!= j && dice() < cEEin_tmp)  G[i][j] = GEEin_tmp/(NEtmp*cEE_tmp);
#ifndef PVinact
	if(i<NR)   G[i][j]=add_lrn*G[i][j];
#else   
	/* check point
	   if(i>=NC && i<Nother)G[i][j]=add_lrn2*G[i][j];
	   if(i>=NR && i<NC)    G[i][j]=add_lrn1*G[i][j];
	*/
	if(i<NR)             G[i][j]=add_lrn*G[i][j];
#endif
      }
      else{
	if( i!= j && dice() < cEE_tmp)    G[i][j] = GEE_tmp/(NEtmp*cEE_tmp);
      }
    }
    for(int j = NEtmp; j < Ntmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cEI_tmp ) G[i][j] = GEI_tmp/((Ntmp-NEtmp)*cEI_tmp);
    }
  }
  for(int i = NEtmp; i < Ntmp; i++){
    for(int j = 0; j < NEtmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cIE_tmp ) G[i][j] = GIE_tmp/(NEtmp*cIE_tmp);
    }
    for(int j = NEtmp; j < Ntmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cII_tmp ) G[i][j] = GII_tmp/((Ntmp-NEtmp)*cII_tmp);
    }
  }
  return G;
}



// 2017.7.14  fix some incorrect descriptions
DBLMAT calc_G2_EC3_CA1(string n_pre_area, string n_post_area, string ver){
  
  if(n_pre_area!="EC3" || n_post_area!="CA1"){
    cout<<"invalid nema!!";
    exit(1);
  }
  char pos; // index of the cluster
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_CA1; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_EC3; j++){
      G[i].push_back(0.0);

      if(  (i<NE_CA1_L && j<NE_EC3_L) || ((i>=NE_CA1_L && i<NE_CA1_R ) && (j>=NE_EC3_L && j<NE_EC3_R)))  pos='l'; // learned connections btw. corres. clsts.
      else if( ((i>=NE_CA1_R && i<NE_CA1_other) && (j>=NE_EC3_R && j<NE_EC3_other)) ||
	       (i>=NE_CA1_other && j>=NE_EC3_other))    	pos='c'; // non-learned connection btw corres. clsts.
      else if( (j>=NE_EC3_R && j<NE_EC3) )    	pos='t'; //
      else pos='e';  // other networks

      if(pos=='l'){
	if( dice()<rhebb_c*cEE_EC3_CA1 ) G[i][j] = GEE_EC3_CA1*rhebb_g*add_lrn/((NE_EC3_L)*cEE_EC3_CA1);	
      }
      else if(pos=='c'){
	if( dice()<rhebb_c*cEE_EC3_CA1 ) G[i][j] = GEE_EC3_CA1*rhebb_g1/((NE_EC3_L)*cEE_EC3_CA1);  //check point
      }
      else if(pos=='t'){
	if( dice()<rhebb_c*cEE_EC3_CA1 ) G[i][j] = GEE_EC3_CA1*rhebb_g1/((NE_EC3_L)*cEE_EC3_CA1);  //check point
      }
      else if(pos=='e'){
	if( dice()<cEE_EC3_CA1 )         G[i][j] = GEE_EC3_CA1/((NE_EC3-NE_EC3_R)*cEE_EC3_CA1); //check point
      }
    }
    for(int j = NE_EC3; j < N_EC3; j++)      G[i].push_back(0.0);
  }
  if( ver == "no_EC3toPV"){
    for(int i = NE_CA1; i < NPV+NE_CA1; i++){
      G.push_back(dvec);
      for(int j = 0; j < N_EC3; j++)G[i].push_back(0.0);
    }
  }
  else{
    cout<<"Error! in generation of EC3toCA1 network";
  }

  for(int i = NE_CA1+NPV; i < N_CA1; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_EC3; j++) G[i].push_back(0.0);
  }
  return G;
}



// last update 2017.7.14
DBLMAT calc_G2_CA1_EC5(string n_pre_area, string n_post_area){
  char pos; // index of the cluster
  if(n_pre_area!="CA1" || n_post_area!="EC5"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC5; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_CA1; j++){
      G[i].push_back(0.0);


      if(  (i<NE_EC5_L && j<NE_CA1_L) || ((i>=NE_EC5_L && i<NE_EC5_R ) && (j>=NE_CA1_L && j<NE_CA1_R)))  pos='l'; // learned connections btw. corres. clsts.
      else if( ((i>=NE_EC5_R && i<NE_EC5_other) && (j>=NE_CA1_R && j<NE_CA1_other)))                     pos='c'; 
      else if( (j>=NE_CA1_R && j<NE_CA1))                     pos='t'; 
      else if(  (i>=NE_EC5_other && j>=NE_CA1_other))   pos='d'; // non-learned connection btw corres. clsts.
      else pos='e';  // other networks
      
      if(pos=='l'){
	if( dice()<rhebb_c*cEE_CA1_EC5 ) G[i][j] = GEE_CA1_EC5*rhebb_g*add_lrn/((NE_CA1_L)*cEE_CA1_EC5);	
      }
      else if(pos=='c'){
	if( dice()<rhebb_c*cEE_CA1_EC5 ) G[i][j] = GEE_CA1_EC5*rhebb_g1/((NE_CA1_L)*cEE_CA1_EC5); //check point
      }
      else if(pos=='d'){
	if( dice()<rhebb_c*cEE_CA1_EC5 ) G[i][j] = GEE_CA1_EC5*rhebb_g1/((NE_CA1_L)*cEE_CA1_EC5);
      }
      else if(pos=='t'){
	if( dice()<rhebb_c*cEE_CA1_EC5 ) G[i][j] = GEE_CA1_EC5*rhebb_g1/((NE_CA1_L)*cEE_CA1_EC5);
      }
      else if(pos=='e'){
	if( dice()<cEE_CA1_EC5 )         G[i][j] = GEE_CA1_EC5/((NE_CA1-NE_CA1_R)*cEE_CA1_EC5);
      }

    }
    for(int j = NE_CA1; j < N_CA1; j++)      G[i].push_back(0.0);
  }
  for(int i = NE_EC5; i < N_EC5; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_CA1; j++)           G[i].push_back(0.0);
  }

  return G;
}



DBLMAT calc_G2_EC5_EC3(string n_pre_area, string n_post_area){
  char pos; // index of the cluster
  if(n_pre_area!="EC5" || n_post_area!="EC3"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC3; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_EC5; j++){
      G[i].push_back(0.0);
  
      if(  (i<NE_EC3_L && j<NE_EC5_L) || ((i>=NE_EC3_L && i<NE_EC3_R ) && (j>=NE_EC5_L && j<NE_EC5_R)))  pos='l'; // learned connections btw. corres. clsts.
      else if( ((i>=NE_EC3_R && i<NE_EC3_other) && (j>=NE_EC5_R && j<NE_EC5_other)))                     pos='c'; 
      else if(  (i>=NE_EC3_other && j>=NE_EC5_other))   pos='d'; // non-learned connection btw corres. clsts.
      else pos='e';  // other networks

      if(pos=='l'){
	if( dice()<rhebb_c*cEE_EC5_EC3 )      G[i][j] = GEE_EC5_EC3*rhebb_g*add_lrn/(NE_EC5_L*cEE_EC5_EC3);
      }
      else if(pos=='c' or pos=='d'){
	if( dice()<rhebb_c*cEE_EC5_EC3 )      G[i][j] = GEE_EC5_EC3*rhebb_g/(NE_EC5_L*cEE_EC5_EC3);
      }
      else if(pos=='e'){
	if( dice()<cEE_EC5_EC3 )            G[i][j] = GEE_EC5_EC3/(NE_EC5_L*cEE_EC5_EC3);
      }
    }
    for(int j = NE_EC5;j<N_EC5; j++)  G[i].push_back(0.0);
  }
  for(int i = NE_EC3; i < N_EC3; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_CA1; j++)      G[i].push_back(0.0);
  }

  return G;
}
#endif



/****************************************************/

/********   networks without NE_EC?_other *********/

/****************************************************/
#ifndef NE_EC3_other
DBLMAT calc_G3(string n_area){   
  int Ntmp,NEtmp,NC,NL,NR;
  double cEE_tmp,cEI_tmp,cIE_tmp,cII_tmp,cEEin_tmp;
  double GEE_tmp,GEI_tmp,GIE_tmp,GII_tmp,GEEin_tmp;
  
  if(n_area=="EC3"){
    Ntmp =N_EC3; NEtmp=NE_EC3;
    NL    =NE_EC3_L;
    NR    =NE_EC3_R;
    cEE_tmp=cEE_EC3;    cEI_tmp=cEI_EC3;
    cEEin_tmp=cEE_EC3*rhebb_c;
    cIE_tmp=cIE_EC3;    cII_tmp=cII_EC3;
    GEE_tmp=GEE_EC3;        GEI_tmp=GEI_EC3;
    GEEin_tmp=GEE_EC3*rhebb_g;
    GIE_tmp=GIE_EC3;        GII_tmp=GII_EC3;
  }
  else if(n_area=="EC5"){
    Ntmp =N_EC5; NEtmp=NE_EC5;
    NL    =NE_EC5_L;
    NR    =NE_EC5_R;
    cEE_tmp=cEE_EC5;    cEI_tmp=cEI_EC5;
    cEEin_tmp=cEE_EC5*rhebb_c;
    cIE_tmp=cIE_EC5;    cII_tmp=cII_EC5;
    GEE_tmp=GEE_EC5;        GEI_tmp=GEI_EC5;
    GEEin_tmp=GEE_EC5*rhebb_g*add_EC5;
    GIE_tmp=GIE_EC5;        GII_tmp=GII_EC5;
  }
  /*
  else if(n_area=="EC2"){
    Ntmp =N_EC2; NEtmp=NE_EC2;
    NL    =NE_EC2_L;
    NR    =NE_EC2_R;
    cEE_tmp=cEE_EC2;    cEI_tmp=cEI_EC2;
    cEEin_tmp=cEE_EC2*rhebb_c;
    cIE_tmp=cIE_EC2;    cII_tmp=cII_EC2;
    GEE_tmp=GEE_EC2;        GEI_tmp=GEI_EC2;
    GEEin_tmp=GEE_EC2*rhebb_g;
    GIE_tmp=GIE_EC2;        GII_tmp=GII_EC2;
  }
  */


  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < Ntmp; i++)     G.push_back(dvec);
  for(int i = 0; i < NEtmp; i++){
    for(int j = 0; j < NEtmp; j++){
      G[i].push_back(0.0);
      if( (i<NL && j<NL) || ((i>=NL && i<NR) && (j>=NL && j<NR)) ){
	if( i!= j && dice() < cEEin_tmp)  G[i][j] = add_lrn*GEEin_tmp/(NEtmp*cEE_tmp);
      }
      else{
	if( i!= j && dice() < cEE_tmp)    G[i][j] = GEE_tmp/(NEtmp*cEE_tmp);
      }
    }
    for(int j = NEtmp; j < Ntmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cEI_tmp ) G[i][j] = GEI_tmp/((Ntmp-NEtmp)*cEI_tmp);
    }
  }
  for(int i = NEtmp; i < Ntmp; i++){
    for(int j = 0; j < NEtmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cIE_tmp ) G[i][j] = GIE_tmp/(NEtmp*cIE_tmp);
    }
    for(int j = NEtmp; j < Ntmp; j++){
      G[i].push_back(0.0);
      if( i != j && dice() < cII_tmp ) G[i][j] = GII_tmp/((Ntmp-NEtmp)*cII_tmp);
    }
  }
  return G;
}




DBLMAT calc_G3_EC3_CA1(string n_pre_area, string n_post_area, string ver){
  
  if(n_pre_area!="EC3" || n_post_area!="CA1"){
    cout<<"invalid nema!!";
    exit(1);
  }
  char pos; // index of the cluster
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_CA1; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_EC3; j++){
      G[i].push_back(0.0);

      if(  (i<NE_CA1_L && j<NE_EC3_L) || ((i>=NE_CA1_L && i<NE_CA1_R ) && (j>=NE_EC3_L && j<NE_EC3_R)))  pos='l'; // learned connections btw. corres. clsts.
      else pos='e';  // other networks

      if(pos=='l'){
	if( dice()<rhebb_c*cEE_EC3_CA1 ) G[i][j] = GEE_EC3_CA1*rhebb_g*add_lrn/(80*cEE_EC3_CA1);	
      }
      else if(pos=='e'){
	if( dice()<cEE_EC3_CA1 )         G[i][j] = GEE_EC3_CA1/(40*cEE_EC3_CA1); // in the original ver., NE_EC3-NE_EC3_R is used instead of 40
      }
    }
    for(int j = NE_EC3; j < N_EC3; j++)      G[i].push_back(0.0);
  }
  if( ver == "no_EC3toPV"){
    for(int i = NE_CA1; i < NPV+NE_CA1; i++){
      G.push_back(dvec);
      for(int j = 0; j < N_EC3; j++)G[i].push_back(0.0);
    }
  }
  else{
    cout<<"Error! in generation of EC3toCA1 network";
  }

  for(int i = NE_CA1+NPV; i < N_CA1; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_EC3; j++) G[i].push_back(0.0);
  }
  return G;
}



DBLMAT calc_G3_CA1_EC5(string n_pre_area, string n_post_area){
  char pos; // index of the cluster
  if(n_pre_area!="CA1" || n_post_area!="EC5"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC5; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_CA1; j++){
      G[i].push_back(0.0);
      if(  (i<NE_EC5_L && j<NE_CA1_L) || ((i>=NE_EC5_L && i<NE_EC5_R ) && (j>=NE_CA1_L && j<NE_CA1_R)))  pos='l'; // learned connections btw. corres. clsts.
      else if( i<NE_EC5_R && j>=NE_CA1_R )  pos='k';
      else pos='e';  // other networks
      
      if(pos=='l'){
	if( dice()<rhebb_c*cEE_CA1_EC5 ) G[i][j] = GEE_CA1_EC5*rhebb_g*add_lrn/((NE_CA1_L)*cEE_CA1_EC5);	
      }
      else if(pos=='k'){
	if( dice()<6*cEE_CA1_EC5 )       G[i][j] = GEE_CA1_EC5*36/((NE_CA1-NE_CA1_R)*cEE_CA1_EC5);  //check point for ver6
      }
      else if(pos=='e'){
	if( dice()<cEE_CA1_EC5 )         G[i][j] = GEE_CA1_EC5/((NE_CA1-NE_CA1_R)*cEE_CA1_EC5);
      }

    }
    for(int j = NE_CA1; j < N_CA1; j++)      G[i].push_back(0.0);
  }
  for(int i=NE_EC5;i<N_EC5;i++){
    G.push_back(dvec);
    for(int j = 0; j < NE_CA1; j++){
      G[i].push_back(0.0);
      if(dice()<cIE_CA1_EC5) G[i][j] = GIE_CA1_EC5/((NE_CA1)*cIE_CA1_EC5);
    }
    for(int j = NE_CA1; j < N_CA1; j++)
      G[i].push_back(0.0);
  }

  return G;
}


DBLMAT calc_G3_EC5_EC3(string n_pre_area, string n_post_area){
  char pos; // index of the cluster
  if(n_pre_area!="EC5" || n_post_area!="EC3"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC3; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_EC5; j++){
      G[i].push_back(0.0);
  
      if(  (i<NE_EC3_L && j<NE_EC5_L) || ((i>=NE_EC3_L && i<NE_EC3_R ) && (j>=NE_EC5_L && j<NE_EC5_R)))  pos='l'; // learned connections btw. corres. clsts.
      else pos='e';  // other networks

      if(pos=='l'){
	if( dice()<rhebb_c*cEE_EC5_EC3 )      G[i][j] = GEE_EC5_EC3*rhebb_g*add_lrn/(NE_EC5_L*cEE_EC5_EC3);
      }
      else if(pos=='e'){
	if( dice()<cEE_EC5_EC3 )            G[i][j] = GEE_EC5_EC3/(NE_EC5_L*cEE_EC5_EC3);
      }
    }
    for(int j = NE_EC5;j<N_EC5; j++)  G[i].push_back(0.0);
  }
  for(int i = NE_EC3; i < N_EC3; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_EC5; j++)      G[i].push_back(0.0);
  }
  
  return G;
}



DBLMAT calc_Dec(string n_pre_area,string n_post_area){   // older version? 2019/8/13 checks
  int Ntmp,NEtmp,NC,NL,NR;
  double cEE_tmp,cEI_tmp,cIE_tmp,cII_tmp,cEEin_tmp;
  double GEE_tmp,GEI_tmp,GIE_tmp,GII_tmp,GEEin_tmp;
  
  if(n_pre_area=="Dec" && n_post_area=="Dec"){

    DBLMAT G;
    //vector< vector<double> > 
    for(int i = 0; i < N_dec; i++)     G.push_back(dvec);
    for(int i = 0; i < NE_Dec; i++){
      for(int j = 0; j < NE_Dec; j++){
	G[i].push_back(0.0);
	if( i!= j && dice() < cEE_Dec)  G[i][j] = GEE_Dec/(NE_Dec*cEE_Dec);
      }
      for(int j = NE_Dec; j < N_Dec; j++){
	G[i].push_back(0.0);
	if( (i<NE_Dec_L  && j<NI_Dec_L)   && dice() < cEI_Dec ) G[i][j] = GEI_Dec/((NI_Dec_L-NE_Dec)*cEI_Dec);
	if( (i>=NE_Dec_L && j>=NI_Dec_L)  && dice() < cEI_Dec ) G[i][j] = GEI_Dec/((NI_Dec_L-NE_Dec)*cEI_Dec);
      }
    }
    for(int i = NE_Dec; i < N_dec; i++){
      for(int j = 0; j < NE_Dec; j++){
	G[i].push_back(0.0);
	if( i != j && dice() < cIE_Dec ) G[i][j] = GIE_Dec/(NE_Dec*cIE_Dec);
      }
      for(int j = NE_Dec; j < N_dec; j++){
	G[i].push_back(0.0);
	if( i != j && dice() < cII_Dec ) G[i][j] = GII_Dec/((N_dec-NE_Dec)*cII_Dec);
      }
    }
  }
  else{n_post_area=="Dec" && n_pre_area!="Dec"){
    if(n_pre_area=="EC3"){
      Npre =N_EC3; NEpre=NE_EC3;
      NLpre    =NE_EC3_L;
      NRpre    =NE_EC3_R;
      cEE_tmp=cEE_EC3_Dec;  
      cIE_tmp=cIE_EC3_Dec;  
      GEE_tmp=GEE_EC3_Dec;  
      GIE_tmp=GIE_EC3_Dec;  
    }
    else if(n_pre_area=="CA1"){
      Npre =N_CA1; NEpre=NE_CA1;
      NLpre    =NE_CA1_L;
      NRpre    =NE_CA1_R;

      cEE_tmp=cEE_CA1_Dec;  
      cIE_tmp=cIE_CA1_Dec;  
      GEE_tmp=GEE_CA1_Dec;  
      GIE_tmp=GIE_CA1_Dec;  
    }

    Npost =N_Dec; NEpost=NE_Dec;
    NLpost    =NE_Dec_L;
    NRpost    =NE_Dec_R;
    NLIpost   =NI_Dec_L; 
    NRIpost   =NI_Dec_R; 


    DBLMAT G;
    //vector< vector<double> > 
    for(int i = 0; i < Npost; i++){  //i: post
      G.push_back(dvec);
      for(int j = 0; j < Npre; j++){
	G[i].push_back(0.0);
	
	if(  (i<NLpre && j<NLpost) || ((i>=NLpre && i<NRpre ) && (j>=NLpost && j<NRpost)))  pos='l'; // learned connections btw. corres. clsts.
	else( (i<NLpre && (j>=NEpost && j<NLIpost)) || ((i>=NLpre && i<NRpre ) && (j>=NILpost && j<NRIpost)))  pos='f'; // learned connections btw. corres. clsts.

	if(pos=='l'){
	  if( dice()<cEE_tmp )   G[i][j] = GEE_tmp*rhebb_g/(NLpre*cEE_tmp);
	}
	if(pos=='f'){
	  if( dice()<cIE_tmp )   G[i][j] = GIE_tmp/(NILpre*cIE_tmp);
	}
      }
    }
  }  

  return G;
}


/*
DBLMAT calc_G3_EC5_EC2(string n_pre_area, string n_post_area){
  char pos; // index of the cluster
  if(n_pre_area!="EC5" || n_post_area!="EC2"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC2; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_EC5; j++){
      G[i].push_back(0.0);
  
      if(  (i<NE_EC2_L && j<NE_EC5_L) || ((i>=NE_EC2_L && i<NE_EC2_R ) && (j>=NE_EC5_L && j<NE_EC5_R)))  pos='l'; // learned connections btw. corres. clsts.
      else pos='e';  // other networks

      if(pos=='l'){
	if( dice()<rhebb_c*cEE_EC5_EC2 )      G[i][j] = GEE_EC5_EC2*rhebb_g*add_lrn/(NE_EC5_L*cEE_EC5_EC2);
      }
      else if(pos=='e'){
	if( dice()<cEE_EC5_EC2 )            G[i][j] = GEE_EC5_EC2/(NE_EC5_L*cEE_EC5_EC2);
      }
    }
    for(int j = NE_EC5;j<N_EC5; j++)  G[i].push_back(0.0);
  }
  for(int i = NE_EC2; i < N_EC2; i++){
    G.push_back(dvec);
    for(int j = 0; j < N_EC5; j++)      G[i].push_back(0.0);
  }
  
  return G;
}


DBLMAT calc_G3_EC2_EC3(string n_pre_area, string n_post_area){
  char pos; // index of the cluster
  if(n_pre_area!="EC2" || n_post_area!="EC3"){
    cout<<"invalid nema!!";
    exit(1);
  }
  
  DBLMAT G;
  //vector< vector<double> > 
  for(int i = 0; i < NE_EC3; i++){  //i: post
    G.push_back(dvec);
    for(int j = 0; j < NE_EC2; j++){
      G[i].push_back(0.0);
  
      if(  (i<NE_EC3_L && j<NE_EC2_L) || ((i>=NE_EC3_L && i<NE_EC3_R ) && (j>=NE_EC2_L && j<NE_EC2_R)))  pos='l'; // learned connections btw. corres. clsts.
      else pos='e';  // other networks

      if(pos=='l'){
	if( dice()<rhebb_c*cEE_EC2_EC3 )      G[i][j] = GEE_EC2_EC3*rhebb_g*add_lrn/(NE_EC2_L*cEE_EC2_EC3);
      }
      else if(pos=='e'){
	if( dice()<cEE_EC2_EC3 )            G[i][j] = GEE_EC2_EC3/(NE_EC2_L*cEE_EC2_EC3);
      }
    }
    for(int j = NE_EC2;j<N_EC2; j++)  G[i].push_back(0.0);
  }
  for(int i = NE_EC3; i < N_EC3; i++){
    G.push_back(dvec);
    for(int j = 0; j < NE_EC2; j++){
      G[i].push_back(0.0);
      if(dice()<cIE_EC2_EC3)G[i][j]=GIE_EC2_EC3/(NE_EC2*cIE_EC2_EC3);
    }
    for(int j=NE_EC2;j>N_EC2;j++)  G[i].push_back(0.0);
  }
  
  return G;
}
*/
#endif






/***************************************************/





/************   calculation of delay    *******************/
/************************************************************/

vector< vector<int> > calc_d(int Ntmp, int NEtmp){
  vector< vector<int> > d;
  for(int i = 0; i < NEtmp; i++){
    d.push_back(ivec);
    for(int j = 0; j < NEtmp; j++)    d[i].push_back( (int)floor( ((dEmin_wti + (dEmax_wti-dEmin_wti)*dice())/hstp) ) );
    for(int j = NEtmp; j < Ntmp; j++) d[i].push_back( (int)floor( ((dImin_wti + (dImax_wti-dImin_wti)*dice())/hstp) ) );
  }
  for(int i = NEtmp; i < Ntmp; i++){
    d.push_back(ivec);
    for(int j = 0; j < NEtmp; j++)    d[i].push_back( (int)floor( ((dEmin_wti + (dEmax_wti-dEmin_wti)*dice())/hstp) ) );
    for(int j = NEtmp; j < Ntmp; j++) d[i].push_back( (int)floor( ((dImin_wti + (dImax_wti-dImin_wti)*dice())/hstp) ) );
  }
  return d;
}


vector< vector<int> > calc_d_btw(int Ntmp_pre, int NEtmp_pre,int Ntmp_post,int NEtmp_post){
  vector< vector<int> > d;
  for(int i = 0; i < NEtmp_post; i++){
    d.push_back(ivec);
    for(int j = 0; j < NEtmp_pre; j++)        d[i].push_back( (int)floor( ((dEmin_btw + (dEmax_btw-dEmin_btw)*dice())/hstp) ) );
    for(int j = NEtmp_pre; j < Ntmp_pre; j++) d[i].push_back( (int)floor( ((dImin_btw + (dImax_btw-dImin_btw)*dice())/hstp) ) );
  }
  for(int i = NEtmp_post; i < Ntmp_post; i++){
    d.push_back(ivec);
    for(int j = 0; j < NEtmp_pre; j++)        d[i].push_back( (int)floor( ((dEmin_btw + (dEmax_btw-dEmin_btw)*dice())/hstp) ) );
    for(int j = NEtmp_pre; j < Ntmp_pre; j++) d[i].push_back( (int)floor( ((dImin_btw + (dImax_btw-dImin_btw)*dice())/hstp) ) );
  }
  return d;
}


vector< vector<int> > calc_d1_btw(int Ntmp_pre, int NEtmp_pre,int Ntmp_post,int NEtmp_post, char l_type){ //l_type=S: intra area, l_type=L: inter layer
  double add_lag=2.0;
  if(l_type=='S') add_lag=1.0;
  else if(l_type=='L') add_lag=2.0;
  
  
  vector< vector<int> > d;
  for(int i = 0; i < NEtmp_post; i++){
    d.push_back(ivec);
    for(int j = 0; j < NEtmp_pre; j++)        d[i].push_back( (int)floor( ((dEmin_btw + add_lag*(dEmax_btw-dEmin_btw)*dice())/hstp) ) );
    for(int j = NEtmp_pre; j < Ntmp_pre; j++) d[i].push_back( (int)floor( ((dImin_btw + add_lag*(dImax_btw-dImin_btw)*dice())/hstp) ) );
  }
  for(int i = NEtmp_post; i < Ntmp_post; i++){
    d.push_back(ivec);
    for(int j = 0; j < NEtmp_pre; j++)        d[i].push_back( (int)floor( ((dEmin_btw + add_lag*(dEmax_btw-dEmin_btw)*dice())/hstp) ) );
    for(int j = NEtmp_pre; j < Ntmp_pre; j++) d[i].push_back( (int)floor( ((dImin_btw + add_lag*(dImax_btw-dImin_btw)*dice())/hstp) ) );
  }
  return d;
}
