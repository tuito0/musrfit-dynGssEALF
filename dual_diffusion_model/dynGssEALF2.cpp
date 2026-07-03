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

Bool_t dynGssEALF2::cache_valid = false;
Float_t***** dynGssEALF2::table = nullptr;

void dynGssEALF2::allocateTable(){
  if(table != nullptr) {
    return;
  }
  table = new Float_t****[N_DT];
  for (Int_t i=0; i < N_DT; i++) {
    table[i] = new Float_t***[N_Q];
    for (Int_t j=0; j < N_Q; j++) {
      table[i][j] = new Float_t**[N_NU];
      for (Int_t k=0; k < N_NU; k++) {
	table[i][j][k] = new Float_t*[N_LF];
	for (Int_t l=0; l < N_LF; l++) {
	  table[i][j][k][l] = new Float_t[N_NU2];
	}
      }
    }
  }
}

void dynGssEALF2::loadTable(){
  ifstream file(TABLE_PATH, ios::binary);
  if (!file.is_open()) {
    cerr << "Table for dynGssEALF2 not detected." << endl;
    throw runtime_error("file open error");
  }
  else{
    cout << "Table for dynGssEALF2 detected." << endl;
    for (Int_t i=0; i < N_DT; i++) {
      for (Int_t j=0; j < N_Q; j++) {
	for (Int_t k=0; k < N_NU; k++) {
	  for (Int_t l=0; l < N_LF; l++) {
	    file.read(reinterpret_cast<char*>(table[i][j][k][l]), N_NU2*sizeof(Float_t));
	  }
	}
      }
    }
    file.close();
  }
}


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


Double_t dynGssEALF2::glf(Double_t t, Double_t delta, Double_t Q,
                          Double_t nu1, Double_t nu2, Double_t LF) const
{
  Double_t a = find_n(t, delta);
  Double_t b = find_k(Q);
  Double_t c = find_l(delta, nu1);
  Double_t d = find_h(LF, delta);
  Double_t e = find_m(delta, nu2);

  Int_t ia = static_cast<Int_t>(std::floor(a));
  Int_t ib = static_cast<Int_t>(std::floor(b));
  Int_t ic = static_cast<Int_t>(std::floor(c));
  Int_t id = static_cast<Int_t>(std::floor(d));
  Int_t ie = static_cast<Int_t>(std::floor(e));

  Double_t fa = a - ia;
  Double_t fb = b - ib;
  Double_t fc = c - ic;
  Double_t fd = d - id;
  Double_t fe = e - ie;

  const Double_t wa[2] = {1.0 - fa, fa};
  const Double_t wb[2] = {1.0 - fb, fb};
  const Double_t wc[2] = {1.0 - fc, fc};
  const Double_t wd[2] = {1.0 - fd, fd};
  const Double_t we[2] = {1.0 - fe, fe};

  Double_t result = 0.0;

  for (Int_t da = 0; da <= 1; ++da)
    for (Int_t db = 0; db <= 1; ++db)
      for (Int_t dc = 0; dc <= 1; ++dc)
        for (Int_t dd = 0; dd <= 1; ++dd)
          for (Int_t de = 0; de <= 1; ++de) {
            const Double_t w =
                wa[da] * wb[db] * wc[dc] * wd[dd] * we[de];

            result += w * getTable(
                ia + da, ib + db, ic + dc, id + dd, ie + de);
          }

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
