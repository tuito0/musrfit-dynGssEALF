/***************************************************************************

  dynGssEALF2.cpp

  Author: Takashi U. Ito
  e-mail: tuito@post.j-parc.jp

***************************************************************************/

/***************************************************************************
 *   Copyright (C) 2024 by Takashi U. Ito                             *
 *   tuito@post.j-parc.jp                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "dynGssEALF2.h"



Double_t dynGssEALF2::glf_narrowlim(Double_t t, Double_t delta, Double_t Q, Double_t nu1, Double_t nu2, Double_t LF) const {
  // nu/Delta >= NARROW_LIM; narrowing limit

  const Double_t wL = GMU*LF;
  Double_t g_sta, g_dyn;
  Double_t nu = nu1 + nu2;

  if((Q>0.0)&&(Q<1.0)){
    // for the Q component
    if((nu/(delta*sqrt(Q)) >= NARROW_LIM) || (wL/(delta*sqrt(Q)) >= HLF_LIM)){
      g_dyn = exp( -2.0*Q*delta*delta/(wL*wL+nu*nu)/(wL*wL+nu*nu)*( (wL*wL+nu*nu)*nu*t + (wL*wL-nu*nu)*(1.0-exp(-nu*t-0.5*Q*delta*delta*t*t)*cos(wL*t)) - 2.0*nu*wL*exp(-nu*t-0.5*Q*delta*delta*t*t)*sin(wL*t) ) ); //Modified Abragam function [A. Keren, PRB 50, 10039 (1994)] for the Q component
    }
    else{
      g_dyn = dynGssEALF2::glf(t,delta*sqrt(Q),1.0,nu,0.0,LF);
    }
    // for the 1-Q component
    if((nu2/(delta*sqrt(1.0-Q)) >= NARROW_LIM)||(wL/(delta*sqrt(1.0-Q)) >= HLF_LIM)){
      g_sta = exp( -2.0*(1.0-Q)*delta*delta/(wL*wL+nu2*nu2)/(wL*wL+nu2*nu2)*( (wL*wL+nu2*nu2)*nu2*t + (wL*wL-nu2*nu2)*(1.0-exp(-nu2*t-0.5*(1.0-Q)*delta*delta*t*t)*cos(wL*t)) - 2.0*nu2*wL*exp(-nu2*t-0.5*(1.0-Q)*delta*delta*t*t)*sin(wL*t) ) ); //Modified Abragam function for the 1-Q component
    }
    else{
      g_sta = dynGssEALF2::glf(t,delta*sqrt(1.0-Q),1.0,nu2,0.0,LF);
    }
    return g_dyn*g_sta;
  }

  else if(Q==0.0){
    g_dyn = 1.0;
    if(((nu2/(delta*sqrt(1.0-Q)) >= NARROW_LIM))||(wL/(delta*sqrt(1.0-Q)) >= HLF_LIM)){
      g_sta = exp( -2.0*(1.0-Q)*delta*delta/(wL*wL+nu2*nu2)/(wL*wL+nu2*nu2)*( (wL*wL+nu2*nu2)*nu2*t + (wL*wL-nu2*nu2)*(1.0-exp(-nu2*t-0.5*(1.0-Q)*delta*delta*t*t)*cos(wL*t)) - 2.0*nu2*wL*exp(-nu2*t-0.5*(1.0-Q)*delta*delta*t*t)*sin(wL*t) ) ); //Modified Abragam function for the 1-Q component
    }
    else{
      g_sta = dynGssEALF2::glf(t,delta*sqrt(1.0-Q),1.0,nu2,0.0,LF);	
    }
    return g_dyn*g_sta;
  }

  else { // Q=1.0
    g_sta = 1.0;
    if((nu/(delta*sqrt(Q)) >= NARROW_LIM)||(wL/(delta*sqrt(Q)) >= HLF_LIM)){
      g_dyn = exp( -2.0*Q*delta*delta/(wL*wL+nu*nu)/(wL*wL+nu*nu)*( (wL*wL+nu*nu)*nu*t + (wL*wL-nu*nu)*(1.0-exp(-nu*t-0.5*Q*delta*delta*t*t)*cos(wL*t)) - 2.0*nu*wL*exp(-nu*t-0.5*Q*delta*delta*t*t)*sin(wL*t) ) ); //Modified Abragam function for the Q component
    }
    else{
      g_dyn = dynGssEALF2::glf(t,delta*sqrt(Q),1.0,nu,0.0,LF);
    }
    return g_dyn*g_sta;
  }
      
}

Double_t dynGssEALF2::glf_highLFlim(Double_t t, Double_t delta, Double_t Q, Double_t nu1, Double_t nu2, Double_t LF) const {
  // GMU*LF/Delta>=HLF_LIM; high LF limit (but not in narrowing limit)
  
  const Double_t wL = GMU*LF;
  Double_t g_sta, g_dyn;
  Double_t nu = nu1 + nu2;

  //Modified Abragam function [A. Keren, PRB 50, 10039 (1994)]
  g_dyn = exp( -2.0*Q*delta*delta/(wL*wL+nu*nu)/(wL*wL+nu*nu)*( (wL*wL+nu*nu)*nu*t + (wL*wL-nu*nu)*(1.0-exp(-nu*t-0.5*Q*delta*delta*t*t)*cos(wL*t)) - 2.0*nu*wL*exp(-nu*t-0.5*Q*delta*delta*t*t)*sin(wL*t) ) );

  //Modified Abragam function 
  g_sta = exp( -2.0*(1.0-Q)*delta*delta/(wL*wL+nu2*nu2)/(wL*wL+nu2*nu2)*( (wL*wL+nu2*nu2)*nu2*t + (wL*wL-nu2*nu2)*(1.0-exp(-nu2*t-0.5*(1.0-Q)*delta*delta*t*t)*cos(wL*t)) - 2.0*nu2*wL*exp(-nu2*t-0.5*(1.0-Q)*delta*delta*t*t)*sin(wL*t) ) );

  return g_sta*g_dyn;

}

Double_t dynGssEALF2::glf(Double_t t, Double_t delta, Double_t Q, Double_t nu1, Double_t nu2, Double_t LF) const
{ 

  Double_t a,b,c,d,e,suba,subb,subc,subd,sube;
  a=find_n(t,delta); b=find_k(Q); c=find_l(delta,nu1); d=find_h(LF,delta); e=find_m(delta,nu2);
  suba=a-floor(a); subb=b-floor(b); subc=c-floor(c); subd=d-floor(d); sube=e-floor(e);

  if(suba==0.0) suba+=0.0001; //approximation
  if(subb==0.0) subb+=0.0001; //approximation
  if(subc==0.0) subc+=0.0001; //approximation
  if(subd==0.0) subd+=0.0001; //approximation
  if(sube==0.0) sube+=0.0001; //approximation  
  if(suba==1.0) suba-=0.0001; //approximation
  if(subb==1.0) subb-=0.0001; //approximation
  if(subc==1.0) subc-=0.0001; //approximation
  if(subd==1.0) subd-=0.0001; //approximation
  if(sube==1.0) sube-=0.0001; //approximation  

  Double_t p00000;
  p00000=getTable((Int_t)floor(a),(Int_t)floor(b),(Int_t)floor(c),(Int_t)floor(d),(Int_t)floor(e));

  Double_t p10000,p01000,p00100,p00010,p00001;
  p10000=getTable((Int_t)floor(a)+1,(Int_t)floor(b),(Int_t)floor(c),(Int_t)floor(d),(Int_t)floor(e));
  p01000=getTable((Int_t)floor(a),(Int_t)floor(b)+1,(Int_t)floor(c),(Int_t)floor(d),(Int_t)floor(e));
  p00100=getTable((Int_t)floor(a), (Int_t)floor(b), (Int_t)floor(c)+1, (Int_t)floor(d),(Int_t)floor(e));
  p00010=getTable((Int_t)floor(a),(Int_t)floor(b),(Int_t)floor(c),(Int_t)floor(d)+1,(Int_t)floor(e));
  p00001=getTable((Int_t)floor(a),(Int_t)floor(b),(Int_t)floor(c),(Int_t)floor(d),(Int_t)floor(e)+1);    

  Double_t p11000,p10100,p10010,p10001,p01100,p01010,p01001,p00110,p00101,p00011;
  p11000=getTable((Int_t)floor(a)+1,(Int_t)floor(b)+1,(Int_t)floor(c),(Int_t)floor(d),(Int_t)floor(e));
  p10100=getTable((Int_t)floor(a)+1,(Int_t)floor(b),(Int_t)floor(c)+1,(Int_t)floor(d),(Int_t)floor(e));    
  p10010=getTable((Int_t)floor(a)+1,(Int_t)floor(b),(Int_t)floor(c),(Int_t)floor(d)+1,(Int_t)floor(e));
  p10001=getTable((Int_t)floor(a)+1,(Int_t)floor(b),(Int_t)floor(c),(Int_t)floor(d),(Int_t)floor(e)+1);
  p01100=getTable((Int_t)floor(a),(Int_t)floor(b)+1,(Int_t)floor(c)+1,(Int_t)floor(d),(Int_t)floor(e));
  p01010=getTable((Int_t)floor(a),(Int_t)floor(b)+1,(Int_t)floor(c),(Int_t)floor(d)+1,(Int_t)floor(e));
  p01001=getTable((Int_t)floor(a),(Int_t)floor(b)+1,(Int_t)floor(c),(Int_t)floor(d),(Int_t)floor(e)+1);
  p00110=getTable((Int_t)floor(a),(Int_t)floor(b),(Int_t)floor(c)+1,(Int_t)floor(d)+1,(Int_t)floor(e));
  p00101=getTable((Int_t)floor(a),(Int_t)floor(b),(Int_t)floor(c)+1,(Int_t)floor(d),(Int_t)floor(e)+1);
  p00011=getTable((Int_t)floor(a),(Int_t)floor(b),(Int_t)floor(c),(Int_t)floor(d)+1,(Int_t)floor(e)+1);

  Double_t p11100,p11010,p11001,p10110,p10101,p10011,p01110,p01101,p01011,p00111;
  p11100=getTable((Int_t)floor(a)+1,(Int_t)floor(b)+1,(Int_t)floor(c)+1,(Int_t)floor(d),(Int_t)floor(e));    
  p11010=getTable((Int_t)floor(a)+1,(Int_t)floor(b)+1,(Int_t)floor(c),(Int_t)floor(d)+1,(Int_t)floor(e));
  p11001=getTable((Int_t)floor(a)+1,(Int_t)floor(b)+1,(Int_t)floor(c),(Int_t)floor(d),(Int_t)floor(e)+1);
  p10110=getTable((Int_t)floor(a)+1,(Int_t)floor(b),(Int_t)floor(c)+1,(Int_t)floor(d)+1,(Int_t)floor(e));
  p10101=getTable((Int_t)floor(a)+1,(Int_t)floor(b),(Int_t)floor(c)+1,(Int_t)floor(d),(Int_t)floor(e)+1);
  p10011=getTable((Int_t)floor(a)+1,(Int_t)floor(b),(Int_t)floor(c),(Int_t)floor(d)+1,(Int_t)floor(e)+1);
  p01110=getTable((Int_t)floor(a),(Int_t)floor(b)+1,(Int_t)floor(c)+1,(Int_t)floor(d)+1,(Int_t)floor(e));
  p01101=getTable((Int_t)floor(a),(Int_t)floor(b)+1,(Int_t)floor(c)+1,(Int_t)floor(d),(Int_t)floor(e)+1);
  p01011=getTable((Int_t)floor(a),(Int_t)floor(b)+1,(Int_t)floor(c),(Int_t)floor(d)+1,(Int_t)floor(e)+1);
  p00111=getTable((Int_t)floor(a),(Int_t)floor(b),(Int_t)floor(c)+1,(Int_t)floor(d)+1,(Int_t)floor(e)+1);      

  Double_t p11110,p11101,p11011,p10111,p01111;
  p11110=getTable((Int_t)floor(a)+1,(Int_t)floor(b)+1,(Int_t)floor(c)+1,(Int_t)floor(d)+1,(Int_t)floor(e));     p11101=getTable((Int_t)floor(a)+1,(Int_t)floor(b)+1,(Int_t)floor(c)+1,(Int_t)floor(d),(Int_t)floor(e)+1);     p11011=getTable((Int_t)floor(a)+1,(Int_t)floor(b)+1,(Int_t)floor(c),(Int_t)floor(d)+1,(Int_t)floor(e)+1);     p10111=getTable((Int_t)floor(a)+1,(Int_t)floor(b),(Int_t)floor(c)+1,(Int_t)floor(d)+1,(Int_t)floor(e)+1);     p01111=getTable((Int_t)floor(a),(Int_t)floor(b)+1,(Int_t)floor(c)+1,(Int_t)floor(d)+1,(Int_t)floor(e)+1);        

  Double_t p11111;
  p11111=getTable((Int_t)floor(a)+1,(Int_t)floor(b)+1,(Int_t)floor(c)+1,(Int_t)floor(d)+1,(Int_t)floor(e)+1);    

  Double_t w00000;
  w00000=1.0/(suba*subb*subc*subd*sube);

  Double_t w10000,w01000,w00100,w00010,w00001;
  w10000=1.0/((1.0-suba)*subb*subc*subd*sube);
  w01000=1.0/(suba*(1.0-subb)*subc*subd*sube);
  w00100=1.0/(suba*subb*(1.0-subc)*subd*sube);
  w00010=1.0/(suba*subb*subc*(1.0-subd)*sube);
  w00001=1.0/(suba*subb*subc*subd*(1.0-sube));  
  
  Double_t w11000,w10100,w10010,w10001,w01100,w01010,w01001,w00110,w00101,w00011;
  w11000=1.0/((1.0-suba)*(1.0-subb)*subc*subd*sube);
  w10100=1.0/((1.0-suba)*subb*(1.0-subc)*subd*sube);
  w10010=1.0/((1.0-suba)*subb*subc*(1.0-subd)*sube);
  w10001=1.0/((1.0-suba)*subb*subc*subd*(1.0-sube));  
  w01100=1.0/(suba*(1.0-subb)*(1.0-subc)*subd*sube);
  w01010=1.0/(suba*(1.0-subb)*subc*(1.0-subd)*sube);
  w01001=1.0/(suba*(1.0-subb)*subc*subd*(1.0-sube));
  w00110=1.0/(suba*subb*(1.0-subc)*(1.0-subd)*sube);
  w00101=1.0/(suba*subb*(1.0-subc)*subd*(1.0-sube));
  w00011=1.0/(suba*subb*subc*(1.0-subd)*(1.0-sube));    
  
  Double_t w11100,w11010,w11001,w10110,w10101,w10011,w01110,w01101,w01011,w00111;
  w11100=1.0/((1.0-suba)*(1.0-subb)*(1.0-subc)*subd*sube);
  w11010=1.0/((1.0-suba)*(1.0-subb)*subc*(1.0-subd)*sube);
  w11001=1.0/((1.0-suba)*(1.0-subb)*subc*subd*(1.0-sube));  
  w10110=1.0/((1.0-suba)*subb*(1.0-subc)*(1.0-subd)*sube);
  w10101=1.0/((1.0-suba)*subb*(1.0-subc)*subd*(1.0-sube));
  w10011=1.0/((1.0-suba)*subb*subc*(1.0-subd)*(1.0-sube));    
  w01110=1.0/(suba*(1.0-subb)*(1.0-subc)*(1.0-subd)*sube);
  w01101=1.0/(suba*(1.0-subb)*(1.0-subc)*subd*(1.0-sube));
  w01011=1.0/(suba*(1.0-subb)*subc*(1.0-subd)*(1.0-sube));
  w00111=1.0/(suba*subb*(1.0-subc)*(1.0-subd)*(1.0-sube));      

  Double_t w11110,w11101,w11011,w10111,w01111;  
  w11110=1.0/((1.0-suba)*(1.0-subb)*(1.0-subc)*(1.0-subd)*sube);
  w11101=1.0/((1.0-suba)*(1.0-subb)*(1.0-subc)*subd*(1.0-sube));
  w11011=1.0/((1.0-suba)*(1.0-subb)*subc*(1.0-subd)*(1.0-sube));
  w10111=1.0/((1.0-suba)*subb*(1.0-subc)*(1.0-subd)*(1.0-sube));
  w01111=1.0/(suba*(1.0-subb)*(1.0-subc)*(1.0-subd)*(1.0-sube));        
  
  Double_t w11111;
  w11111=1.0/((1.0-suba)*(1.0-subb)*(1.0-subc)*(1.0-subd)*(1.0-sube));        

  Double_t result;
  result = (p00000*w00000+p10000*w10000+p01000*w01000+p00100*w00100+p00010*w00010+p00001*w00001+p11000*w11000+p10100*w10100+p10010*w10010+p10001*w10001+p01100*w01100+p01010*w01010+p01001*w01001+p00110*w00110+p00101*w00101+p00011*w00011+p11100*w11100+p11010*w11010+p11001*w11001+p10110*w10110+p10101*w10101+p10011*w10011+p01110*w01110+p01101*w01101+p01011*w01011+p00111*w00111+p11110*w11110+p11101*w11101+p11011*w11011+p10111*w10111+p01111*w01111+p11111*w11111)/(w00000+w10000+w01000+w00100+w00010+w00001+w11000+w10100+w10010+w10001+w01100+w01010+w01001+w00110+w00101+w00011+w11100+w11010+w11001+w10110+w10101+w10011+w01110+w01101+w01011+w00111+w11110+w11101+w11011+w10111+w01111+w11111);
  
  return result;

}

Double_t dynGssEALF2::find_n(Double_t t, Double_t delta) const {
  Int_t i=0;
  while(i<N_DT){
    if((Double_t)i*DT*DELTA0 <= t*delta) i++;
    else break;
  }
  i=i-1;
  if(i>=(N_DT-1)){
    return (Double_t)i-0.0001; //approximation
  }
  else{
    return (Double_t)i+(t*delta-(Double_t)i*DT*DELTA0)/(DT*DELTA0);
  }
}

Double_t dynGssEALF2::find_k(Double_t Q) const {
  Int_t i=0;
  while(i<N_Q)
    {
      if(meas_Q[i] <= Q) i++;
      else break;
    }
  i=i-1;
  if(i>=(N_Q-1)){
    return (Double_t)i-0.0001; //approximation
  }
  else{
    return (Double_t)i+(Q-meas_Q[i])/(meas_Q[i+1]-meas_Q[i]);
  }
}

Double_t dynGssEALF2::find_l(Double_t delta, Double_t nu1) const{
  Int_t i=0;
  while(i<(N_NU-1))
    {
      if(meas_NUonD[i] <= nu1/delta) i++;
      else break;
    }
  i=i-1;
  if(i>=N_NU){
    return (Double_t)i-0.0001; //approximation
  }
  else{
    return (Double_t)i+(nu1/delta-meas_NUonD[i])/(meas_NUonD[i+1]-meas_NUonD[i]);
  }
}

Double_t dynGssEALF2::find_h(Double_t LF, Double_t delta) const{
  Int_t i=0;
  while(i<(N_LF-1))
    {
      if(meas_LF[i] <= GMU*LF/delta) i++;
      else break;
    }
  i=i-1;
  if(i>=N_LF){
    return (Double_t)i-0.0001; //approximation
  }
  else{
    return (Double_t)i+(GMU*LF/delta-meas_LF[i])/(meas_LF[i+1]-meas_LF[i]);
  }
}

Double_t dynGssEALF2::find_m(Double_t delta, Double_t nu2) const{
  Int_t i=0;
  while(i<(N_NU2-1))
    {
      if(meas_NU2onD[i] <= nu2/delta) i++;
      else break;
    }
  i=i-1;
  if(i>=N_NU2){
    return (Double_t)i-0.0001; //approximation
  }
  else{
    return (Double_t)i+(nu2/delta-meas_NU2onD[i])/(meas_NU2onD[i+1]-meas_NU2onD[i]);
  }
}

//#################################################//

ClassImp(dynGssEALF2)  // for the ROOT dictionary

Double_t dynGssEALF2::operator()(Double_t x, const vector<Double_t> &par) const {
  assert(par.size()==5); // make sure the number of parameters handed to the function is correct

  Double_t delta=par[2];
  if(par[2] < 1.0e-5) delta=1.0e-5;

  //return 0 when parameters are out of range
  if(par[0]<0) return 0.0; //LF
  if((par[1]<0)||(par[1]>1)) return 0.0; //Q
  if((x<0)||(x*delta>DELTA0*DT*(N_DT-1))) return 0.0; //Delta
  if((par[3]<0)||(par[4]<0)) return 0.0; //nu, nu2


  if(((par[3]+par[4])>=NARROW_LIM*(delta*sqrt(par[1])))||(par[4]>=NARROW_LIM*(delta*sqrt(1.0-par[1])))){ 
    //narrowing limit
    return dynGssEALF2::glf_narrowlim(x,delta,par[1],par[3],par[4],par[0]);
  }
  else if(GMU*par[0]/delta>=HLF_LIM) {
    //high LF limit (but not in narrowing limit)
    return dynGssEALF2::glf_highLFlim(x,delta,par[1],par[3],par[4],par[0]);    
  }
  else {
    return dynGssEALF2::glf(x,delta,par[1],par[3],par[4],par[0]);
  }
      

}
