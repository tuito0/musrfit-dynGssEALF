/***************************************************************************

  dynGssEALF2.h

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

#include "PUserFcnBase.h"
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define NARROW_LIM 10.0
#define HLF_LIM 10.0
// Standard delta used in the MC simulation
#define DELTA0 0.6022
// Time interval in us used in the MC table
#define DT 0.1
// Number of data points in the MC table
#define N_DT 101
#define N_Q  51
#define N_NU 14
#define N_LF 18
#define N_NU2 14
// Muon gyromagnetic ratio in 10^6 rad/s
#define GMU 851.58649

#ifndef TABLE_PATH
#define TABLE_PATH "dynGssEALF2_tbl_v0.1.0.bin"
#endif

using namespace std;

class dynGssEALF2 : public PUserFcnBase {

 public:

  // default constructor and destructor
  dynGssEALF2(){
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


  ~dynGssEALF2(){
    for (Int_t i=0; i < N_DT; i++) {
      for (Int_t j=0; j < N_Q; j++) {
	for (Int_t k=0; k < N_NU; k++) {
	  for (Int_t l=0; l < N_LF; l++) {
	    delete[] table[i][j][k][l];
	    }
	  delete[] table[i][j][k];
	  }
	delete[] table[i][j];
	}
      delete[] table[i];
      }
    delete[] table;
  }

    Bool_t NeedGlobalPart() const { return false; }
    void SetGlobalPart(vector<void *> &globalPart, UInt_t idx) { }
    Bool_t GlobalPartIsValid() const { return true; }

  // function operator
  Double_t operator()(Double_t, const vector<Double_t>&) const;

  // definition of the class for the ROOT dictionary
  ClassDef(dynGssEALF2,1)

  private:
  Float_t***** table;
  const Double_t meas_Q[N_Q]={0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58,0.60,0.62,0.64,0.66,0.68,0.70,0.72,0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1};
  const Double_t meas_NUonD[N_NU]={0,0.1,0.2,0.333,0.5,0.667,1,1.429,2,2.5,3.333,5,6.667,10.00}; //nu/Delta
  const Double_t meas_LF[N_LF]={0,0.333,0.666,1,1.333,1.666,2,2.333,2.666,3,3.5,4,5,6,7,8,9,10};  //GMU*LF/Delta
  const Double_t meas_NU2onD[N_NU2]={0,0.1,0.2,0.333,0.5,0.667,1,1.429,2,2.5,3.333,5,6.667,10.00}; //nu2/Delta
  Double_t glf_narrowlim(Double_t t, Double_t delta, Double_t Q, Double_t nu1, Double_t nu2, Double_t LF) const;
  Double_t glf_highLFlim(Double_t t, Double_t delta, Double_t Q, Double_t nu1, Double_t nu2, Double_t LF) const;  
  Double_t glf(Double_t t, Double_t delta, Double_t Q, Double_t nu1, Double_t nu2, Double_t LF) const;
  Double_t find_n(Double_t t, Double_t delta) const;
  Double_t find_k(Double_t Q) const;
  Double_t find_l(Double_t delta, Double_t nu) const;
  Double_t find_h(Double_t LF, Double_t delta) const;
  Double_t find_m(Double_t delta, Double_t nu2) const;  
  Double_t getTable(Int_t id1, Int_t id2, Int_t id3, Int_t id4, Int_t id5) const {
    return (Double_t) dynGssEALF2::table[id1][id2][id3][id4][id5];
  }
  

};
