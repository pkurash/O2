// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "FV0Simulation/Digitizer.h"
#include "FV0Base/Geometry.h"
#include "FV0Base/Constants.h"

#include <TRandom.h>
#include <algorithm>

ClassImp(o2::fv0::Digitizer);

using namespace o2::math_utils;
using namespace o2::fv0;

void Digitizer::clear()
{
  mEventId = -1;
  mSrcId = -1;
  // for (auto& analogSignal : mPmtChargeVsTime) {
  // std::fill_n(std::begin(analogSignal), analogSignal.size(), 0);
  // }
}

//_______________________________________________________________________
void Digitizer::init()
{
  LOG(INFO) << "V0Digitizer::init -> start = ";

  mNBins = FV0DigParam::Instance().waveformNbins;      //Will be computed using detector set-up from CDB
  mBinSize = FV0DigParam::Instance().waveformBinWidth; //Will be set-up from CDB

  for (Int_t i = 0; i < Constants::nFv0Channels; i++) {
//    mPmtChargeVsTime[i].resize(mNBins);
  }

  auto const roundVc = [&](int i) -> int {
    return (i / Vc::float_v::Size) * Vc::float_v::Size;
  };
  // set up PMT response tables
  Float_t offset = -0.5f * mBinSize; // offset \in [-0.5..0.5] * mBinSize
  Int_t const nBins = roundVc(std::lround(4.0f * FV0DigParam::Instance().pmtTransitTime / mBinSize));
  for (auto& table : mPmtResponseTables) {
    table.resize(nBins);
    Float_t t = -2.0f * FV0DigParam::Instance().pmtTransitTime + offset; // t \in offset + [-2 2] * FV0DigParam::Instance().mPmtTransitTime
    for (Int_t j = 0; j < nBins; ++j) {
      table[j] = Digitizer::PmtResponse(t);
      t += mBinSize;
    }
    offset += mBinSize / Float_t(DP::NUM_PMT_RESPONSE_TABLES - 1);
  }

  TF1 scintDelayFn("fScintDelay", "gaus",
                   -6.0f * FV0DigParam::Instance().intrinsicTimeRes,
                   +6.0f * FV0DigParam::Instance().intrinsicTimeRes);
  scintDelayFn.SetParameters(1, 0, FV0DigParam::Instance().intrinsicTimeRes);
  mRndScintDelay.initialize(scintDelayFn);

  // Initialize function describing the PMT time response
  TF1 pmtResponseFn("mPmtResponse",
                    &Digitizer::PmtResponse,
                    -1.0f * FV0DigParam::Instance().pmtTransitTime,
                    +2.0f * FV0DigParam::Instance().pmtTransitTime, 0);
  pmtResponseFn.SetNpx(100);
  mPmtTimeIntegral = pmtResponseFn.Integral(-1.0f * FV0DigParam::Instance().pmtTransitTime,
                                            +2.0f * FV0DigParam::Instance().pmtTransitTime);

  // Initialize function describing PMT response to the single photoelectron
  TF1 singlePhESpectrumFn("mSinglePhESpectrum",
                          &Digitizer::SinglePhESpectrum,
                          FV0DigParam::Instance().photoelMin,
                          FV0DigParam::Instance().photoelMax, 0);
  Float_t const meansPhE = singlePhESpectrumFn.Mean(FV0DigParam::Instance().photoelMin, FV0DigParam::Instance().photoelMax);
  mRndGainVar.initialize([&]() -> float {
    return singlePhESpectrumFn.GetRandom(FV0DigParam::Instance().photoelMin, FV0DigParam::Instance().photoelMax) / meansPhE;
  });

  TF1 signalShapeFn("signalShape", "crystalball", 0, 200);
  signalShapeFn.SetParameters(FV0DigParam::Instance().shapeConst,
                              FV0DigParam::Instance().shapeMean,
                              FV0DigParam::Instance().shapeSigma,
                              FV0DigParam::Instance().shapeAlpha,
                              FV0DigParam::Instance().shapeN);
  mRndSignalShape.initialize([&]() -> float {
    return signalShapeFn.GetRandom(0, mBinSize * Float_t(mNBins));
  });

  LOG(INFO) << "V0Digitizer::init -> finished";
}

void Digitizer::process(const std::vector<o2::fv0::Hit>& hits)
{
  LOG(INFO) << "[FV0] Digitizer::process(): begin with " << hits.size() << " hits";

  std::vector<int> hitIdx(hits.size());
  std::iota(std::begin(hitIdx), std::end(hitIdx), 0);
  std::sort(std::begin(hitIdx), std::end(hitIdx), [&hits](int a, int b) { return hits[a].GetTrackID() < hits[b].GetTrackID(); });

  auto const roundVc = [&](int i) -> int {
    return (i / Vc::float_v::Size) * Vc::float_v::Size;
  };
  Int_t parentIdPrev = -10;
  Int_t nEmp = 0;
  // use ordered hits
  for (auto ids : hitIdx) {
    const auto& hit = hits[ids];
    Int_t detId = hit.GetDetectorID();
    Double_t hitEdep = hit.GetHitValue() * 1e3; //convert to MeV

    // TODO: check how big is inaccuracy if more than 1 'below-threshold' particles hit the same detector cell
    if (hitEdep < FV0DigParam::Instance().singleMipThreshold) {
      continue;
    }

    float distanceFromXc = 0;
    if (Geometry::instance()->isRing5(detId)) {
      distanceFromXc = getDistFromCellCenter(detId, hit.GetX(), hit.GetY());
    }

    int iChannelPerCell = 0;
    while (iChannelPerCell < 2) { // loop over 2 channels, into which signal from each cell in ring 5 is split
      if (Geometry::instance()->isRing5(detId)) {
        // The first channel number is located counter-clockwise from the cell center
        //   and remains identical to the detector number, the second one is clockwise and incremented by 8
        if (iChannelPerCell == 1) {
          detId += 8;
        }

        // Split signal magnitude to fractions depending on the distance of the hit from the cell center
        hitEdep = (hit.GetHitValue() * 1e3) * getSignalFraction(distanceFromXc, iChannelPerCell == 0);
        LOG(INFO) << "  detId: " << detId << "-" << iChannelPerCell << " hitEdep: " << hitEdep << " distanceFromXc: " << distanceFromXc;
        ++iChannelPerCell;
      } else {
        iChannelPerCell = 2; // not a ring 5 cell -> don't repeat the loop
      }

      Double_t const nPhotons = hitEdep * DP::N_PHOTONS_PER_MEV;
      Int_t const nPhE = SimulateLightYield(detId, nPhotons);
      Float_t const t = hit.GetTime() * 1e9 + FV0DigParam::Instance().pmtTransitTime;
      Float_t const charge = TMath::Qe() * FV0DigParam::Instance().pmtGain * mBinSize / mPmtTimeIntegral;

      //auto& analogSignal = mPmtChargeVsTime[detId];
      
      float timeHit = t;
      timeHit += mIntRecord.getTimeNS();
      LOG(INFO) << "event: " << mEventId <<  " initial time: " 
                << t << " hit number: " << ids << " hit time: " 
                << timeHit << " IR time: " << mIntRecord.getTimeNS();
      o2::InteractionTimeRecord irHit(timeHit);
      LOG(INFO) << "irHit: " << irHit.getTimeNS();
      std::array<o2::InteractionTimeRecord, NBC2Cache> cachedIR;

      int nCachedIR = 0;
      for (int i = BCCacheMin; i < BCCacheMax + 1; i++) {
        double tNS = timeHit + o2::constants::lhc::LHCBunchSpacingNS * i;
        cachedIR[nCachedIR].setFromNS(tNS);
        if (tNS < 0 && cachedIR[nCachedIR] > irHit) {
          continue; // don't go to negative BC/orbit (it will wrap)
        }
        setBCCache(cachedIR[nCachedIR++]); // ensure existence of cached container
      }
      LOG(INFO) << "mCache.size() = " << mCache.size();

      timeHit = timeHit - cachedIR[0].bc2ns();

      for (int ir = 0; ir < mCache.size(); ir ++) {
         auto bcCache = getBCCache(mCache[ir]);
         // LOG(INFO) << "ir = " << ir << " timeNS: " << (*bcCache).getTimeNS() << ", mIR.timeNS: " << mIntRecord.getTimeNS()
         // << ", difference with mIR: " << (*bcCache).getTimeNS() - mIntRecord.getTimeNS();
         for (int ich = 0; ich < Constants::nFv0Channels; ich ++) {
           (*bcCache).mPmtChargeVsTime[ich].resize(mNBins);
         }
      }

      for (int ir = 0; ir < mCache.size(); ir ++) {
         auto bcCache = getBCCache(mCache[ir]);
         // LOG(INFO) << "ir = " << ir << " timeNS: " << (*bcCache).getTimeNS() << ", mIR.timeNS: " << mIntRecord.getTimeNS()
         // << ", difference with mIR: " << (*bcCache).getTimeNS() - mIntRecord.getTimeNS();
         for (int ich = 0; ich < Constants::nFv0Channels; ich ++) {
           (*bcCache).mPmtChargeVsTime[ich].resize(mNBins);
         }
      }

      std::vector<bool> added;
      added.resize(mCache.size());
      std::fill(added.begin(), added.end(), false);

      int  parentId = hit.GetTrackID();  


      for (Int_t iPhE = 0; iPhE < nPhE; ++iPhE) {
        Float_t const tPhE = t + mRndSignalShape.getNextValue();
        Int_t const firstBin = roundVc(
          TMath::Max((Int_t)0, (Int_t)((tPhE - FV0DigParam::Instance().pmtTransitTime) / mBinSize)));
        // Int_t const lastBin = TMath::Min((Int_t)mNBins - 1,
        // (Int_t)((tPhE + 2. * FV0DigParam::Instance().pmtTransitTime) / mBinSize));
        Int_t const lastBin  = TMath::Min((Int_t)(NBC2Cache /*mCache.size()*/ * mNBins - 1),
                                      (Int_t)((tPhE + 2. * FV0DigParam::Instance().pmtTransitTime) / mBinSize)); 
        Float_t const tempT = mBinSize * (0.5f + firstBin) - tPhE;
        //Float_t* p = analogSignal.data() + firstBin;
        long iStart = std::lround((tempT + 2.0f * FV0DigParam::Instance().pmtTransitTime) / mBinSize);
        float const offset = tempT + 2.0f * FV0DigParam::Instance().pmtTransitTime - Float_t(iStart) * mBinSize;
        long const iOffset = std::lround(offset / mBinSize * Float_t(DP::NUM_PMT_RESPONSE_TABLES - 1));
        if (iStart < 0) { // this should not happen
          LOG(ERROR) << "V0Digitizer: table lookup failure";
        }
        iStart = roundVc(std::max(long(0), iStart));

        Vc::float_v workVc;
        Vc::float_v pmtVc;
        Float_t const* q = mPmtResponseTables[DP::NUM_PMT_RESPONSE_TABLES / 2 + iOffset].data() + iStart;
        Float_t const* qEnd = &mPmtResponseTables[DP::NUM_PMT_RESPONSE_TABLES / 2 + iOffset].back();

        for (int ir = firstBin / mNBins; ir <= lastBin / mNBins; ir++) {
          int localFirst = (ir == firstBin / mNBins) ? firstBin : 0;
          int localLast = (ir < lastBin / mNBins) ? mNBins : (lastBin - ir * mNBins);
    
          auto bcCache = getBCCache(mCache[ir]);
          (*bcCache).mPmtChargeVsTime[detId].resize(mNBins);
    
          auto& analogSignal = (*bcCache).mPmtChargeVsTime[detId];
          float* p = analogSignal.data() + localFirst;
    
          for (Int_t i = localFirst, iEnd = roundVc(localLast); q < qEnd && i < iEnd; i += Vc::float_v::Size) {
            pmtVc.load(q);
            q += Vc::float_v::Size;
            Vc::prefetchForOneRead(q);
            workVc.load(p);
            workVc += mRndGainVar.getNextValueVc() * charge * pmtVc;
            workVc.store(p);
            p += Vc::float_v::Size;
            Vc::prefetchForOneRead(p);
          }
         added[ir] = true;
        }
      } //photo electron loop

      // Charged particles in MCLabel
      for ( int ir = 0; ir < mCache.size(); ir ++) {
        Int_t const parentId = hit.GetTrackID();
        if (parentId != parentIdPrev) {
          auto bcCache = getBCCache(mCache[ir]);
          (*bcCache).labels.emplace_back(parentId, mEventId, mSrcId, detId);
          parentIdPrev = parentId;
        }
      }
    }
  } //hit loop
}

void Digitizer::analyseWaveformsAndStore(std::vector<fv0::BCData>& digitsBC,
                                         std::vector<fv0::ChannelData>& digitsCh,
                                         dataformats::MCTruthContainer<fv0::MCLabel>& labels)
{
  // Sum charge of all time bins to get total charge collected for a given channel

  size_t const first = digitsCh.size();
  size_t nStored = 0;
  int nbc = 0;

  LOG(INFO) << "mIntRecord: diff front: " << mIntRecord.differenceInBC(mCache.front()) 
            << " diff back: "              << mIntRecord.differenceInBC(mCache.back());

    for (auto& bc:mCache) {
   
  // LOG(INFO) << "Cache element " << nbc << "/" <<  mCache.size() << ", bc = " << bc
  // << " difference front = " << bc.differenceInBC(mCache.front())
  // << " difference back = "  << bc.differenceInBC(mCache.back());
 
      nbc ++ ;          
      
    /* if (bc.IsProcessed) { 
       continue;
     }*/

      Float_t sum_charges = 0;

      for (Int_t ipmt = 0; ipmt < Constants::nFv0Channels; ++ipmt) {
        Float_t timeCfd =   SimulateTimeCfd(bc.mPmtChargeVsTime[ipmt]);
        Float_t totalCharge = integrateCharge(bc.mPmtChargeVsTime[ipmt]);
        if (!timeCfd > 0 || !totalCharge > 0){ 
         continue;
        }
        sum_charges += totalCharge;
        //  LOG(INFO) << "nStored= " << nStored << " timeCfd = " << timeCfd << " totalCharge = " << totalCharge;
        digitsCh.emplace_back(ipmt, timeCfd, totalCharge);
        nStored ++ ;
      }

      if (sum_charges == 0) {
        continue;
      }
      
      if (bc.bc2ns() >= mIntRecord.bc2ns() + o2::constants::lhc::LHCBunchSpacingNS * NBC2Cache) {
        continue;
      }

      if (bc.bc2ns() < mIntRecord.bc2ns()) {
        mCache.pop_front();
        continue;
      }
       // Send MClabels and digitsBC to storage
      size_t const nBC = digitsBC.size();
      digitsBC.emplace_back(first, nStored, bc);
      for (auto const& lbl : bc.labels) {
        labels.addElement(nBC, lbl);
        LOG(INFO) << "lbls: EventID = " << lbl.getEventID() << " source ID = " << lbl.getSourceID();
      }

     bc.IsProcessed = true;  
     // mCache.pop_front();
    }

   // while (TMath::Abs(mCache.front().labels.front().getEventID()  - mCache.back().labels.back().getEventID())  > 0 && !mCache.empty()){
 // while (TMath::Abs(mCache.front().bc2ns() - mCache.back().bc2ns) < 250) {
    // while (mCache.back().differenceInBC(mCache.front()) > 2 * NBC2Cache + 1) {
    // mCache.pop_front();
    // }

}
//_____________________________________________________________________________
Float_t Digitizer::integrateCharge(const ChannelBCDataF& pulse)//, TH1F *h1)
{
    if (pulse.empty() ){
     return 0;
    }
    Float_t totalCharge = 0.0f;
    for (Int_t iTimeBin = 0; iTimeBin < mNBins; ++ iTimeBin) {
      totalCharge += pulse.at(iTimeBin);
      if (pulse[iTimeBin] > 0 && totalCharge > 0) {     
       /* LOG(INFO) << "iTimeBin = "        << iTimeBin
                  <<  " timeBinCharge = " << pulse[iTimeBin] * DP::INV_CHARGE_PER_ADC    
                  << " totalCharge = "    << totalCharge * DP::INV_CHARGE_PER_ADC;*/
      }
    }
    totalCharge *= DP::INV_CHARGE_PER_ADC;
    return (std::lround(totalCharge));
}

/*
void Digitizer::analyseWaveformsAndStore(std::vector<fv0::BCData>& digitsBC,
                                         std::vector<fv0::ChannelData>& digitsCh,
                                         dataformats::MCTruthContainer<fv0::MCLabel>& labels)
{
  // Sum charge of all time bins to get total charge collected for a given channel
  size_t const first = digitsCh.size();
  size_t nStored = 0;
  for (Int_t ipmt = 0; ipmt < Constants::nFv0Channels; ++ipmt) {
    Float_t totalCharge = 0.0f;
    auto const& analogSignal = mPmtChargeVsTime[ipmt];
    for (Int_t iTimeBin = 0; iTimeBin < mNBins; ++iTimeBin) {
      Float_t const timeBinCharge = mPmtChargeVsTime[ipmt][iTimeBin];
      totalCharge += timeBinCharge;
      LOG(INFO) << "ipmt = " << ipmt << " totalCharge = "   << totalCharge * DP::INV_CHARGE_PER_ADC;
      LOG(INFO) << "ipmt = " << ipmt << " timeBinCharge = " << timeBinCharge * DP::INV_CHARGE_PER_ADC; 
    }
    totalCharge *= DP::INV_CHARGE_PER_ADC;
    digitsCh.emplace_back(ipmt, SimulateTimeCfd(ipmt), std::lround(totalCharge));
    ++nStored;
  }

  // Send MClabels and digitsBC to storage
  size_t const nBC = digitsBC.size();
  digitsBC.emplace_back(first, nStored, mIntRecord);

  LOG(INFO) << "mMCLabels.size() = " << mMCLabels.size();
  int nl = 1;

  for (auto const& lbl : mMCLabels) {
    labels.addElement(nBC, lbl);
    nl ++;
  }

 while (TMath::Abs(mMCLabels.front().getEventID() - mMCLabels.back().getEventID())  > 0){ 
   mMCLabels.pop_front();
 }

  LOG(INFO) << "mMCLabels.size() after = " << mMCLabels.size();

//  mMCLabels.clear();
}
*/
//_____________________________________________________________________________
o2::fv0::Digitizer::BCCache& Digitizer::setBCCache(const o2::InteractionTimeRecord& ir)
{
   if (mCache.empty() || mCache.back() < ir) {
     mCache.emplace_back();
     auto& cb = mCache.back();
     cb = ir;
     return cb;
   }
   if (mCache.front() > ir) {
     mCache.emplace_front();
     auto& cb = mCache.front();
     cb = ir;
     return cb;
   }

   for (auto cb = mCache.begin(); cb != mCache.end(); cb++) {
     if ((*cb) == ir) {
       return *cb;
     }
     if (ir < (*cb)) {
       auto cbnew = mCache.emplace(cb); // insert new element before cb
       (*cbnew) = ir;
       return (*cbnew);
     }
   }
   return mCache.front();
}
//_____________________________________________________________________________
o2::fv0::Digitizer::BCCache* Digitizer::getBCCache(const o2::InteractionTimeRecord& ir)
{
   // get pointer on existing cache
   for (auto cb = mCache.begin(); cb != mCache.end(); cb++) {
     if ((*cb) == ir) {
       return &(*cb);
     }
   }
   return nullptr;
}
// -------------------------------------------------------------------------------
// --- Internal helper methods related to conversion of energy-deposition into ---
// --- photons -> photoelectrons -> electrical signal                          ---
// -------------------------------------------------------------------------------
Int_t Digitizer::SimulateLightYield(Int_t pmt, Int_t nPhot) const
{
  const Float_t epsilon = 0.0001f;
  const Float_t p = FV0DigParam::Instance().lightYield * FV0DigParam::Instance().photoCathodeEfficiency;
  if ((fabs(1.0f - p) < epsilon) || nPhot == 0) {
    return nPhot;
  }
  const Int_t n = Int_t(nPhot < 100
                          ? gRandom->Binomial(nPhot, p)
                          : gRandom->Gaus(p * nPhot + 0.5, TMath::Sqrt(p * (1 - p) * nPhot)));
  return n;
}

Float_t Digitizer::SimulateTimeCfd(const ChannelBCDataF& pulse) const //Int_t channel) const
{
  Float_t timeCfd = -1024.0f;

  if (pulse.empty() ){
     return timeCfd;
  }

  Int_t const binShift = TMath::Nint(FV0DigParam::Instance().timeShiftCfd / mBinSize);
  Float_t sigPrev = -pulse.at(0);//[0];
  for (Int_t iTimeBin = 1; iTimeBin < mNBins; ++iTimeBin) {
    Float_t const sigCurrent = (iTimeBin >= binShift
                                  ? 5.0f * pulse[iTimeBin - binShift] - pulse[iTimeBin]
                                  : -pulse[iTimeBin]);
    if (sigPrev < 0.0f && sigCurrent >= 0.0f) {
      timeCfd = Float_t(iTimeBin) * mBinSize;
      break;
    }
    sigPrev = sigCurrent;
  }
  //TH1F * hCfdTimes = (TH1F*)mHist->FindObject("hCfdTimes");
  if (timeCfd > 0){
   //  hCfdTimes->Fill(timeCfd);
  }

  return timeCfd;
}
/*
Float_t Digitizer::SimulateTimeCfd(Int_t channel) const
{
  Float_t timeCfd = -1024.0f;
  Int_t const binShift = TMath::Nint(FV0DigParam::Instance().timeShiftCfd / mBinSize);
  Float_t sigPrev = -mPmtChargeVsTime[channel][0];
  for (Int_t iTimeBin = 1; iTimeBin < mNBins; ++iTimeBin) {
    Float_t const sigCurrent = (iTimeBin >= binShift
                                  ? 5.0f * mPmtChargeVsTime[channel][iTimeBin - binShift] - mPmtChargeVsTime[channel][iTimeBin]
                                  : -mPmtChargeVsTime[channel][iTimeBin]);
    if (sigPrev < 0.0f && sigCurrent >= 0.0f) {
      timeCfd = Float_t(iTimeBin) * mBinSize;
      break;
    }
    sigPrev = sigCurrent;
  }
  return timeCfd;
}
*/
Double_t Digitizer::PmtResponse(Double_t* x, Double_t*)
{
  return Digitizer::PmtResponse(x[0]);
}
//_______________________________________________________________________
Double_t Digitizer::PmtResponse(Double_t x)
{
  // this function describes the PMT time response to a single photoelectron
  if (x > 2 * FV0DigParam::Instance().pmtTransitTime)
    return 0.0;
  if (x < -FV0DigParam::Instance().pmtTransitTime)
    return 0.0;
  x += FV0DigParam::Instance().pmtTransitTime;
  Double_t const x2 = x * x;
  return x2 * std::exp(-x2 * FV0DigParam::Instance().oneOverPmtTransitTime2);
}

Double_t Digitizer::SinglePhESpectrum(Double_t* x, Double_t*)
{
  // x -- number of photo-electrons emitted from the first dynode
  // this function describes the PMT amplitude response to a single photoelectron
  if (x[0] < 0.0)
    return 0.0;
  return (TMath::Poisson(x[0], FV0DigParam::Instance().pmtNbOfSecElec) +
          FV0DigParam::Instance().pmtTransparency * TMath::Poisson(x[0], 1.0));
}

// The Distance is positive for top half-sectors (when the hit position is above the cell center (has higher y))
// TODO: performance check needed
float Digitizer::getDistFromCellCenter(UInt_t cellId, double hitx, double hity)
{
  Geometry* geo = Geometry::instance();

  // Parametrize the line (ax+by+c=0) that crosses the detector center and the cell's middle point
  Point3Dsimple* pCell = &geo->getCellCenter(cellId);
  float x0, y0, z0;
  geo->getGlobalPosition(x0, y0, z0);
  double a = -(y0 - pCell->y) / (x0 - pCell->x);
  double b = 1;
  double c = -(y0 - a * x0);
  // Return the distance from hit to this line
  return (a * hitx + b * hity + c) / TMath::Sqrt(a * a + b * b);
}

float Digitizer::getSignalFraction(float distanceFromXc, bool isFirstChannel)
{
  float fraction = sigmoidPmtRing5(distanceFromXc);
  if (distanceFromXc > 0) {
    return isFirstChannel ? fraction : (1. - fraction);
  } else {
    return isFirstChannel ? (1. - fraction) : fraction;
  }
}
