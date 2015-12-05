//Modified from StMiniMcMaker.h
///////////////////////////////////////////////////////////////
#ifndef StAnaSimMaker_hh     
#define StAnaSimMaker_hh

#include "StMaker.h"
#include "TString.h"
#include <vector>
#include <utility>
#include <map>
#include "TTree.h"
#include "StThreeVectorF.hh"
#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"
//#include "StarRoot/THelixTrack.h"
#include "StTrackGeometry.h"
#include "StiMaker/StKFVertexMaker.h"
#include "StiMaker/StKFVerticesCollection.h"

//StKFVertexMaker includes
#include "TObjArray.h"
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "StEnumerations.h"
#include "TCanvas.h"
#include "TH1K.h"
class StPrimaryVertex;
class StDcaGeometry;
class KFParticle;
class StKFVerticesCollection;
class TSpectrum;

class StEvent;
class StTrack;
class TFile;
class TNtuple;
class TString;
class StMcIstHit;
class StMcPxlHit;
class StIstHit;
class StPxlHit;
class TGeoHMatrix;
class TFile;
class TTree;
class StRun;
class StRandom;
class StMcEvent;
class StMcTrack;
class StPrimaryTrack;
class StPxlDb;
class StKFVertexMaker;

const Double_t PXL_ACTIVE_X_LENGTH = 1.921;
const Double_t PXL_ACTIVE_Y_LENGTH = 1.9872;

const double IST_Width_X = 3.8016;
const double IST_Width_Z = 7.53;
const double PXL_Width_X = 0.9605*2;
const double PXL_Width_Z = 0.9936*2;

class StAnaSimMaker : public StMaker {
	public:

		StAnaSimMaker(const Char_t *outname="");     // constructor
		~StAnaSimMaker();                                 // destructor

		void Clear(Option_t *option="");    // called after every event to cleanup 
		Int_t  Init();                      // called once at the beginning of your job
		Int_t  InitRun(int);                //
		Int_t  Make();                      // invoked for every event
		Int_t  Finish();                    // called once at the end
		Int_t  FinishRun(int);              //
		//Routine to smear hit by resolution with gaussian, mean zero and width res
		double distortHit(double x, double res, double detLength);
		void Projection(StPhysicalHelixD hlx, int index, THashList* istRot);
		void GetAssHit(StEvent *ev, int mcId, int index);

		virtual const char *GetCVS() const {
			static const char cvs[]="Tag $Name:  $ $Id: StAnaSimMaker.h,v 2.4 2003/09/10 19:47:02 perev Exp $ built "__DATE__" "__TIME__ ; 
			return cvs;
		}

		typedef map<int,StMcIstHit*> mcIstHitMap;
		typedef map<int,StMcPxlHit*> mcPxlHitMap;

		//StKFVertexMaker public
		void Fit();
  		TH1F *VtxM() {return fVtxM;}
  		void SetZwindow(Double_t z = 2) {fzWindow = z;}
  		void SetDefaultTempLog(Double_t tLog = 2) {fTempLog = tLog;}
  		static Double_t AnnelingFcn(Double_t TInv=1);
  		TH1 *Vtx() {return fVtx;}
  		StKFVerticesCollection* Vertices() {return fcVertices;}
  		TObjArray &Particles() {return *fParticles;}
  		KFParticle *AddTrackAt(const StDcaGeometry *dca,Int_t kg);
  		void calculateRank(StPrimaryVertex *primV);
  		void SetCanvas(TCanvas *c1) {fc1 = c1;}
  		TCanvas *Canvas() {return fc1;}
  		TH1F *GetVtxs(Int_t pass = 0) {return fVtxs[pass];}
  		TH1K *GetVtxKs(Int_t pass = 0) {return fVtxKs[pass];}
  		TH1F *GetVtxM() {return fVtxM;}

        protected:
		StRandom* myRandom;

		StPxlDb* mPxlDb;

		Double_t resXIst1;
		Double_t resZIst1;
		Int_t mSmear; //to turn smearing on and off

		Double_t mResXPix;
		Double_t mResZPix;
		Double_t mResYPix;

	private:
		//THashList *IstSensorOnGlobal;//YF
		TString outName;           //!
		TFile* mTupFile;               //!
		TTree* trackTree;           //!
		StEvent*         mRcEvent;      //!
		StMcEvent*       mMcEvent;      //!
		rcTrackMapType*  mRcTrackMap;   //!
		mcTrackMapType*  mMcTrackMap;   //!
		const StThreeVectorF*  mRcVtxPos;
		const StThreeVectorF*  mMcVtxPos;
		mcIstHitMap istHitMap;
		mcPxlHitMap pxlHitMap;

		//StKFVertexMaker private
		TObjArray *fParticles; // KF particles
  		Int_t fNzBins;
  		Int_t fNPasses;
  		TSpectrum *fSpectrum;
  		Double_t fzWindow;
  		TH1F  *fVtxM;
  		StKFVerticesCollection **fVerticesPass;
  		static StKFVerticesCollection *fcVertices;  // current vertex collection
  		Double_t fTempLog;
  		ROOT::Math::GSLMinimizer1D *fminBrent;
  		ROOT::Math::Functor1D      *func;
  		TH1F *fVtxs[2];
  		TH1  *fVtx;
  		TH1K *fVtxKs[2];
  		Bool_t mBeamLine;
  		StPrimaryVertexOrder     mVertexOrderMethod; // will default to 0 i.e. orderByNumberOfDaughters
  		TCanvas                 *fc1;

		//  This is needed to make your maker known to root4star.
		//  It must be always the last statement in the class.
		//  Note that this is a macro, that's why the ';' is missing.
		//
		ClassDef(StAnaSimMaker,0)
};
#endif
