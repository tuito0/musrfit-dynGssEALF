/***************************************************************************

  dynGssEALF.cpp

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

#include "dynGssEALF.h"



Double_t dynGssEALF::glf_narrowlim(Double_t t, Double_t delta, Double_t Q, Double_t nu, Double_t LF) const {
  // nu/Delta >= NARROW_LIM; narrowing limit

  const Double_t wL = GMU*LF;
  Double_t g_sta, g_dyn;

  //Modified Abragam function [A. Keren, PRB 50, 10039 (1994)]
  g_dyn = exp( -2.0*Q*delta*delta/(wL*wL+nu*nu)/(wL*wL+nu*nu)*( (wL*wL+nu*nu)*nu*t + (wL*wL-nu*nu)*(1.0-exp(-nu*t-0.5*Q*delta*delta*t*t)*cos(wL*t)) - 2.0*nu*wL*exp(-nu*t-0.5*Q*delta*delta*t*t)*sin(wL*t) ) );

  if(wL==0){
    g_sta = 1.0/3.0+2.0/3.0*(1.0-(1.0-Q)*delta*delta*t*t)*exp(-0.5*(1.0-Q)*delta*delta*t*t); //Static KT function
  }   
  else if(wL/(delta*sqrt(1.0-Q))>=5.0){

    g_sta = exp( -2.0*(1.0-Q)*delta*delta/(wL*wL)*(1.0-exp(-0.5*(1.0-Q)*delta*delta*t*t)*cos(wL*t)) );   //Modified Abragam function with nu=0
  }
  else{
    g_sta = dynGssEALF::glf(t,delta*sqrt(1.0-Q),1.0,0.0,LF);
  }

  return g_sta*g_dyn;

}

Double_t dynGssEALF::glf_highLFlim(Double_t t, Double_t delta, Double_t Q, Double_t nu, Double_t LF) const {
  // GMU*LF/Delta>=HLF_LIM; high LF limit (but not in narrowing limit)
  
  const Double_t wL = GMU*LF;
  Double_t g_sta, g_dyn;

  //Modified Abragam function [A. Keren, PRB 50, 10039 (1994)]
  g_dyn = exp( -2.0*Q*delta*delta/(wL*wL+nu*nu)/(wL*wL+nu*nu)*( (wL*wL+nu*nu)*nu*t + (wL*wL-nu*nu)*(1.0-exp(-nu*t-0.5*Q*delta*delta*t*t)*cos(wL*t)) - 2.0*nu*wL*exp(-nu*t-0.5*Q*delta*delta*t*t)*sin(wL*t) ) );

  g_sta = exp( -2.0*(1.0-Q)*delta*delta/(wL*wL)*(1.0-exp(-0.5*(1.0-Q)*delta*delta*t*t)*cos(wL*t)) );   //Modified Abragam function with nu=0

  return g_sta*g_dyn;

}

Double_t dynGssEALF::glf(Double_t t, Double_t delta, Double_t Q, Double_t nu, Double_t LF) const
{ 

  Double_t a,b,c,d,suba,subb,subc,subd;
  a=find_n(t,delta); b=find_k(Q); c=find_l(delta,nu); d=find_h(LF,delta);
  suba=a-floor(a); subb=b-floor(b); subc=c-floor(c); subd=d-floor(d);

  if(suba==0.0) suba+=0.0001; //approximation
  if(subb==0.0) subb+=0.0001; //approximation
  if(subc==0.0) subc+=0.0001; //approximation
  if(subd==0.0) subd+=0.0001; //approximation  
  if(suba==1.0) suba-=0.0001; //approximation
  if(subb==1.0) subb-=0.0001; //approximation
  if(subc==1.0) subc-=0.0001; //approximation
  if(subd==1.0) subd-=0.0001; //approximation

  Double_t p0000,p1000,p0100,p0010,p0001;
  Double_t p1100,p1010,p1001,p0110,p0101,p0011;
  Double_t p1110,p1101,p1011,p0111,p1111;
  p0000=(Double_t)table[(Int_t)floor(a)][(Int_t)floor(b)][(Int_t)floor(c)][(Int_t)floor(d)];
  p1000=(Double_t)table[(Int_t)floor(a)+1][(Int_t)floor(b)][(Int_t)floor(c)][(Int_t)floor(d)];
  p0100=(Double_t)table[(Int_t)floor(a)][(Int_t)floor(b)+1][(Int_t)floor(c)][(Int_t)floor(d)];
  p0010=(Double_t)table[(Int_t)floor(a)][(Int_t)floor(b)][(Int_t)floor(c)+1][(Int_t)floor(d)];
  p0001=(Double_t)table[(Int_t)floor(a)][(Int_t)floor(b)][(Int_t)floor(c)][(Int_t)floor(d)+1];  
  p1100=(Double_t)table[(Int_t)floor(a)+1][(Int_t)floor(b)+1][(Int_t)floor(c)][(Int_t)floor(d)];
  p1010=(Double_t)table[(Int_t)floor(a)+1][(Int_t)floor(b)][(Int_t)floor(c)+1][(Int_t)floor(d)];    
  p1001=(Double_t)table[(Int_t)floor(a)+1][(Int_t)floor(b)][(Int_t)floor(c)][(Int_t)floor(d)+1];
  p0110=(Double_t)table[(Int_t)floor(a)][(Int_t)floor(b)+1][(Int_t)floor(c)+1][(Int_t)floor(d)];
  p0101=(Double_t)table[(Int_t)floor(a)][(Int_t)floor(b)+1][(Int_t)floor(c)][(Int_t)floor(d)+1];
  p0011=(Double_t)table[(Int_t)floor(a)][(Int_t)floor(b)][(Int_t)floor(c)+1][(Int_t)floor(d)+1];
  p1110=(Double_t)table[(Int_t)floor(a)+1][(Int_t)floor(b)+1][(Int_t)floor(c)+1][(Int_t)floor(d)];    
  p1101=(Double_t)table[(Int_t)floor(a)+1][(Int_t)floor(b)+1][(Int_t)floor(c)][(Int_t)floor(d)+1];
  p1011=(Double_t)table[(Int_t)floor(a)+1][(Int_t)floor(b)][(Int_t)floor(c)+1][(Int_t)floor(d)+1];
  p0111=(Double_t)table[(Int_t)floor(a)][(Int_t)floor(b)+1][(Int_t)floor(c)+1][(Int_t)floor(d)+1];
  p1111=(Double_t)table[(Int_t)floor(a)+1][(Int_t)floor(b)+1][(Int_t)floor(c)+1][(Int_t)floor(d)+1];    
  
  Double_t w0000,w1000,w0100,w0010,w0001;
  Double_t w1100,w1010,w1001,w0110,w0101,w0011;
  Double_t w1110,w1101,w1011,w0111,w1111;
  Double_t result;
  w0000=1.0/(suba*subb*subc*subd);
  w1000=1.0/((1.0-suba)*subb*subc*subd);
  w0100=1.0/(suba*(1.0-subb)*subc*subd);
  w0010=1.0/(suba*subb*(1.0-subc)*subd);
  w0001=1.0/(suba*subb*subc*(1.0-subd));
  w1100=1.0/((1.0-suba)*(1.0-subb)*subc*subd);
  w1010=1.0/((1.0-suba)*subb*(1.0-subc)*subd);
  w1001=1.0/((1.0-suba)*subb*subc*(1.0-subd));
  w0110=1.0/(suba*(1.0-subb)*(1.0-subc)*subd);
  w0101=1.0/(suba*(1.0-subb)*subc*(1.0-subd));
  w0011=1.0/(suba*subb*(1.0-subc)*(1.0-subd));
  w1110=1.0/((1.0-suba)*(1.0-subb)*(1.0-subc)*subd);
  w1101=1.0/((1.0-suba)*(1.0-subb)*subc*(1.0-subd));
  w1011=1.0/((1.0-suba)*subb*(1.0-subc)*(1.0-subd));
  w0111=1.0/(suba*(1.0-subb)*(1.0-subc)*(1.0-subd));
  w1111=1.0/((1.0-suba)*(1.0-subb)*(1.0-subc)*(1.0-subd));

  result= (p0000*w0000+p1000*w1000+p0100*w0100+p0010*w0010+p0001*w0001+p1100*w1100+p1010*w1010+p1001*w1001+p0110*w0110+p0101*w0101+p0011*w0011+p1110*w1110+p1101*w1101+p1011*w1011+p0111*w0111+p1111*w1111)/(w0000+w1000+w0100+w0010+w0001+w1100+w1010+w1001+w0110+w0101+w0011+w1110+w1101+w1011+w0111+w1111);
  
  return result;

}

Double_t dynGssEALF::find_n(Double_t t, Double_t delta) const {
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

Double_t dynGssEALF::find_k(Double_t Q) const {
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

Double_t dynGssEALF::find_l(Double_t delta, Double_t nu) const{
  Int_t i=0;
  while(i<(N_NU-1))
    {
      if(meas_NUonD[i] <= nu/delta) i++;
      else break;
    }
  i=i-1;
  if(i>=N_NU){
    return (Double_t)i-0.0001; //approximation
  }
  else{
    return (Double_t)i+(nu/delta-meas_NUonD[i])/(meas_NUonD[i+1]-meas_NUonD[i]);
  }
}

Double_t dynGssEALF::find_h(Double_t LF, Double_t delta) const{
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


//#################################################//

ClassImp(dynGssEALF)  // for the ROOT dictionary

Double_t dynGssEALF::operator()(Double_t x, const vector<Double_t> &par) const {
  assert(par.size()==4); // make sure the number of parameters handed to the function is correct

  Double_t delta=par[2];
  if(par[2] < 1.0e-5) delta=1.0e-5;

  //return 0 when parameters are out of range
  if(par[0]<0) return 0.0;
  if((par[1]<0)||(par[1]>1)) return 0.0;
  if((x<0)||(x*delta>DELTA0*DT*(N_DT-1))) return 0.0;
  if(par[3]<0) return 0.0;


  if((par[3]/delta)>=NARROW_LIM){
    //narrowing limit
    return dynGssEALF::glf_narrowlim(x,delta,par[1],par[3],par[0]);
  }
  else if(GMU*par[0]/delta>=HLF_LIM) {
    //high LF limit (but not in narrowing limit)
    return dynGssEALF::glf_highLFlim(x,delta,par[1],par[3],par[0]);    
  }
  else {
    return dynGssEALF::glf(x,delta,par[1],par[3],par[0]);
  }

}
