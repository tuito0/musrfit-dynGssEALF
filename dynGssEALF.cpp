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

Bool_t dynGssEALF::cache_valid = false;
Float_t dynGssEALF::table[N_DT][N_Q][N_NU][N_LF] = {};

Double_t dynGssEALF::glf_narrowlim(Double_t t, Double_t delta, Double_t Q, Double_t nu, Double_t LF) const {
  // nu/Delta >= NARROW_LIM; narrowing limit

  const Double_t wL = GMU*LF;
  Double_t g_sta, g_dyn;

  //Modified Abragam function [A. Keren, PRB 50, 10039 (1994)]
  g_dyn = exp( -2.0*Q*delta*delta/(wL*wL+nu*nu)/(wL*wL+nu*nu)*( (wL*wL+nu*nu)*nu*t + (wL*wL-nu*nu)*(1.0-exp(-nu*t-0.5*Q*delta*delta*t*t)*cos(wL*t)) - 2.0*nu*wL*exp(-nu*t-0.5*Q*delta*delta*t*t)*sin(wL*t) ) );

  if(wL==0){
    g_sta = 1.0/3.0+2.0/3.0*(1.0-(1.0-Q)*delta*delta*t*t)*exp(-0.5*(1.0-Q)*delta*delta*t*t); //Static KT function
  }   
  else if(wL>=5.0*(delta*sqrt(1.0-Q))){

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


Double_t dynGssEALF::glf(Double_t t,
                         Double_t delta,
                         Double_t Q,
                         Double_t nu,
                         Double_t LF) const
{
  const Double_t a = find_n(t, delta);
  const Double_t b = find_k(Q);
  const Double_t c = find_l(delta, nu);
  const Double_t d = find_h(LF, delta);

  const Int_t ia = static_cast<Int_t>(std::floor(a));
  const Int_t ib = static_cast<Int_t>(std::floor(b));
  const Int_t ic = static_cast<Int_t>(std::floor(c));
  const Int_t id = static_cast<Int_t>(std::floor(d));

  const Double_t fa = a - static_cast<Double_t>(ia);
  const Double_t fb = b - static_cast<Double_t>(ib);
  const Double_t fc = c - static_cast<Double_t>(ic);
  const Double_t fd = d - static_cast<Double_t>(id);

  const Double_t wa[2] = {1.0 - fa, fa};
  const Double_t wb[2] = {1.0 - fb, fb};
  const Double_t wc[2] = {1.0 - fc, fc};
  const Double_t wd[2] = {1.0 - fd, fd};

  Double_t result = 0.0;

  for (Int_t da = 0; da <= 1; ++da) {
    for (Int_t db = 0; db <= 1; ++db) {
      for (Int_t dc = 0; dc <= 1; ++dc) {
        for (Int_t dd = 0; dd <= 1; ++dd) {

          const Double_t w =
              wa[da] * wb[db] * wc[dc] * wd[dd];

          result += w * static_cast<Double_t>(
              table[ia + da][ib + db][ic + dc][id + dd]
          );
        }
      }
    }
  }

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
