/***************************************************************************

  dynGssEALF2.h

  Author: Takashi U. Ito
  e-mail: tuito@post.j-parc.jp

  mmap-based table access version

***************************************************************************/

#ifndef DYNGSSEALF2_H
#define DYNGSSEALF2_H

#include "PUserFcnBase.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <vector>
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
  dynGssEALF2();
  ~dynGssEALF2();

  Bool_t NeedGlobalPart() const { return false; }
  void SetGlobalPart(vector<void *> &globalPart, UInt_t idx) { }
  Bool_t GlobalPartIsValid() const { return true; }

  // function operator
  Double_t operator()(Double_t, const vector<Double_t>&) const;

  // definition of the class for the ROOT dictionary
  ClassDef(dynGssEALF2,1)

 private:

  static constexpr std::size_t tableSize =
      static_cast<std::size_t>(N_DT)  *
      static_cast<std::size_t>(N_Q)   *
      static_cast<std::size_t>(N_NU)  *
      static_cast<std::size_t>(N_LF)  *
      static_cast<std::size_t>(N_NU2);

  static constexpr std::size_t tableBytes = tableSize * sizeof(Float_t);

  static std::once_flag table_once;
  static const Float_t* table;
  static std::size_t mapped_bytes;

  static void mapTableOnce();
  static void unmapTable();

  static inline std::size_t index(Int_t id1, Int_t id2, Int_t id3,
                                  Int_t id4, Int_t id5)
  {
    return ((((static_cast<std::size_t>(id1) * N_Q +
               static_cast<std::size_t>(id2)) * N_NU +
               static_cast<std::size_t>(id3)) * N_LF +
               static_cast<std::size_t>(id4)) * N_NU2 +
               static_cast<std::size_t>(id5));
  }

  static inline Double_t getTable(Int_t id1, Int_t id2, Int_t id3,
                                  Int_t id4, Int_t id5)
  {
    return static_cast<Double_t>(
        dynGssEALF2::table[index(id1, id2, id3, id4, id5)]
    );
  }

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
};

#endif
