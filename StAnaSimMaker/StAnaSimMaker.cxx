//YF analysis maker for HFT simulation 05/08/2014
//
// - updated to match latest software  05/12/2014
// - added projected hit position      05/15/2014
// - added hit collection for residual 06/05/2014
// - added volume ID for MC hits       07/02/2014
// - added DCA vector and nGTracks     07/08/2014
// - added Global RefMult              07/11/2014
///////////////////////////////////////////////////////////

#include "StAnaSimMaker.h"
#include "StEventTypes.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"

#include "StMcEvent/StMcIstHitCollection.hh"
#include "StMcEvent/StMcIstLayerHitCollection.hh"
#include "StMcEvent/StMcIstHit.hh"
#include "StMcEvent/StMcPxlHitCollection.hh"
#include "StMcEvent/StMcPxlHit.hh"
#include "StEvent/StIstHitCollection.h"
#include "StEvent/StIstHit.h"
#include "StEvent/StPxlHitCollection.h"
#include "StEvent/StPxlHit.h"
//#include "StIstUtil/StIstDigiHit.h"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StEnumerations.h"
//#include "StEnumerations.h"
#include "StEventUtilities/StuRefMult.hh"
#include "StMcEvent/StMcTrack.hh"
#include "StBFChain.h"
#include "StMcEventTypes.hh"
#include "StMcEvent.hh"
#include "StMessMgr.h"
#include "StIstUtil/StIstConsts.h"
#include "tables/St_g2t_ist_hit_Table.h"
#include "tables/St_HitError_Table.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "StPxlDbMaker/StPxlDb.h"
#include "StIstDbMaker/StIstDbMaker.h"
#include "AnaSim.h"
#include "StarClassLibrary/StRandom.hh"
#include "TDataSet.h"
#include "TObjectSet.h"

//StKFVertexMaker includes
#include "RVersion.h"
#if ROOT_VERSION_CODE < 331013
#include "TCL.h"
#else
#include "TCernLib.h"
#endif
#include <map>
using std::map;
#include "TH1.h"
#include "TCanvas.h"
#include "StDcaGeometry.h"
#include "KFParticle.h"
#include "KFVertex.h"
#include "MTrack.h"
#include "VVertex.h"
#include "TH1K.h"
#include "StiMaker/StAnneling.h"
#include "StiMaker/StKFEvent.h"
#include "StiMaker/StKFTrack.h"
#include "StiMaker/StKFVertex.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StiMaker/StVertexP.h"
#include "StiMaker/StVertexT.h"
#include "TDirectory.h"
#include "Stypes.h"
#include "SystemOfUnits.h"
#include "StiMaker/StKFVertexMaker.h"
#include "StDetectorDbMaker/St_vertexSeedC.h"
#include "Sti/StiHit.h"
#include "Sti/StiKalmanTrack.h"
#include "Sti/StiKalmanTrackNode.h"
#include "StiMaker/StiStEventFiller.h"
#include "TRMatrix.h"
#include "TRSymMatrix.h"
#include "Sti/StiToolkit.h"
#include "TArrayI.h"
#include "StiMaker/StiDefaultToolkit.h"

StKFVerticesCollection *StAnaSimMaker::fcVertices = 0;

const StGlobalTrack*  partnerTrack(mcTrackMapType* map, StMcTrack* mT) {
	mcTrackMapIter i = map->find(mT);
	const StGlobalTrack* rT = 0; 
	if (i != map->end()) rT = (*i).second->partnerTrack();
	return rT; 
} 

const StMcTrack*  partnerMcTrack(rcTrackMapType* map, StGlobalTrack* rT) {
	rcTrackMapIter i = map->find(rT);
	const StMcTrack* mT = 0;
	if (i != map->end()) mT = (*i).second->partnerMcTrack();
	return mT;
} 

ClassImp(StAnaSimMaker)

StAnaSimMaker::StAnaSimMaker(const Char_t *outname) : StMaker("name","title"),mTupFile(0),mRcEvent(0),mMcEvent(0),trackTree(0){
	outName=outname;
}



StAnaSimMaker::~StAnaSimMaker() {}

Int_t StAnaSimMaker::Init()
{

  //for hits smearing
  int seed=time(NULL);
  myRandom=new StRandom();
  myRandom->setSeed(seed);

  mSmear=0; //smearing on/off, MC hits no need to smear, just rotate for global position.

	outName.ReplaceAll("root","rectree.root");
	mTupFile =  new TFile(outName, "RECREATE");
	//mTupFile->SetCompressionLevel(9);
	trackTree=new TTree("AnaT","trackTree");
	trackTree->SetAutoSave(10000000);

	cout << "Initialize the recotree ... " << endl;

	trackTree->Branch("mIEvt",&AnaT.mIEvt,"mIEvt/I");
	trackTree->Branch("mIRun",&AnaT.mIRun,"mIRun/I");
	trackTree->Branch("mMcMult",&AnaT.mMcMult,"mMcMult/I");
	trackTree->Branch("mRcMult",&AnaT.mRcMult,"mRcMult/I");
	trackTree->Branch("mMcVtxX",&AnaT.mMcVtxX,"mMcVtxX/F");
	trackTree->Branch("mMcVtxY",&AnaT.mMcVtxY,"mMcVtxY/F");
	trackTree->Branch("mMcVtxZ",&AnaT.mMcVtxZ,"mMcVtxZ/F");
	trackTree->Branch("mNMcPTracks",&AnaT.mNMcPTracks,"mNMcPTracks/I");
	trackTree->Branch("mRcVtxX",&AnaT.mRcVtxX,"mRcVtxX[2]/F");
	trackTree->Branch("mRcVtxY",&AnaT.mRcVtxY,"mRcVtxY[2]/F");
	trackTree->Branch("mRcVtxZ",&AnaT.mRcVtxZ,"mRcVtxZ[2]/F");
	trackTree->Branch("mNPTracks",&AnaT.mNPTracks,"mNPTracks[2]/I");
        trackTree->Branch("mNGTracks",&AnaT.mNGTracks,"mNGTracks/I");
        trackTree->Branch("mNGRefMult",&AnaT.mNGRefMult,"mNGRefMult/I");
	trackTree->Branch("MagField",&AnaT.MagField,"MagField/F");
	//mc tracks
	trackTree->Branch("mNMcTrk",&AnaT.mNMcTrk,"mNMcTrk/I");
	trackTree->Branch("mMcId",AnaT.mMcId,"mMcId[mNMcTrk]/I");
	trackTree->Branch("mGeantId",AnaT.mGeantId,"mGeantId[mNMcTrk]/I");
	trackTree->Branch("mParentMcId",AnaT.mParentMcId,"mParentMcId[mNMcTrk]/I");
	trackTree->Branch("mParentGeantId",AnaT.mParentGeantId,"mParentGeantId[mNMcTrk]/I");
	trackTree->Branch("mMcPt",AnaT.mMcPt,"mMcPt[mNMcTrk]/F");
	trackTree->Branch("mMcPz",AnaT.mMcPz,"mMcPz[mNMcTrk]/F");
	trackTree->Branch("mMcEta",AnaT.mMcEta,"mMcEta[mNMcTrk]/F");
	trackTree->Branch("mMcPhi",AnaT.mMcPhi,"mMcPhi[mNMcTrk]/F");
	trackTree->Branch("mMcMass",AnaT.mMcMass,"mMcMass[mNMcTrk]/F");
	trackTree->Branch("mMcStartX",AnaT.mMcStartX,"mMcStartX[mNMcTrk]/F");   
	trackTree->Branch("mMcStartY",AnaT.mMcStartY,"mMcStartY[mNMcTrk]/F");
	trackTree->Branch("mMcStartZ",AnaT.mMcStartZ,"mMcStartZ[mNMcTrk]/F");
	trackTree->Branch("mMcNhits",AnaT.mMcNhits,"mMcNhits[mNMcTrk]/I");
	trackTree->Branch("mMcNhitsSsd",AnaT.mMcNhitsSsd,"mMcNhitsSsd[mNMcTrk]/I");
	trackTree->Branch("mMcNhitsIst",AnaT.mMcNhitsIst,"mMcNhitsIst[mNMcTrk]/I");
	trackTree->Branch("mMcNhitsPxl2",AnaT.mMcNhitsPxl2,"mMcNhitsPxl2[mNMcTrk]/I");
	trackTree->Branch("mMcNhitsPxl1",AnaT.mMcNhitsPxl1,"mMcNhitsPxl1[mNMcTrk]/I");
        trackTree->Branch("mMcRndHitvId",AnaT.mMcRndHitvId,"mMcRndHitvId[mNMcTrk][20]/I");
	trackTree->Branch("mMcRndHitX",AnaT.mMcRndHitX,"mMcRndHitX[mNMcTrk][20]/F");
	trackTree->Branch("mMcRndHitY",AnaT.mMcRndHitY,"mMcRndHitY[mNMcTrk][20]/F");
	trackTree->Branch("mMcRndHitZ",AnaT.mMcRndHitZ,"mMcRndHitZ[mNMcTrk][20]/F");
	trackTree->Branch("mMcRndHitLX",AnaT.mMcRndHitLX,"mMcRndHitLX[mNMcTrk][20]/F");
	trackTree->Branch("mMcRndHitLY",AnaT.mMcRndHitLY,"mMcRndHitLY[mNMcTrk][20]/F");
	trackTree->Branch("mMcRndHitLZ",AnaT.mMcRndHitLZ,"mMcRndHitLZ[mNMcTrk][20]/F");
        trackTree->Branch("mMcRndHitId",AnaT.mMcRndHitId,"mMcRndHitId[mNMcTrk][20]/I");
	trackTree->Branch("mMcAssHitX",AnaT.mMcAssHitX,"mMcAssHitX[mNMcTrk][20]/F");
	trackTree->Branch("mMcAssHitY",AnaT.mMcAssHitY,"mMcAssHitY[mNMcTrk][20]/F");
	trackTree->Branch("mMcAssHitZ",AnaT.mMcAssHitZ,"mMcAssHitZ[mNMcTrk][20]/F");
	trackTree->Branch("mMcAssHitLX",AnaT.mMcAssHitLX,"mMcAssHitLX[mNMcTrk][20]/F");
	trackTree->Branch("mMcAssHitLY",AnaT.mMcAssHitLY,"mMcAssHitLY[mNMcTrk][20]/F");
	trackTree->Branch("mMcAssHitLZ",AnaT.mMcAssHitLZ,"mMcAssHitLZ[mNMcTrk][20]/F");
	trackTree->Branch("mMcAssHitId",AnaT.mMcAssHitId,"mMcAssHitId[mNMcTrk][20]/I");
	//rc tracks
	trackTree->Branch("mNRcTrk",&AnaT.mNRcTrk,"mNRcTrk/I");
	trackTree->Branch("mRcId",AnaT.mRcId,"mRcId[mNRcTrk]/I");
	trackTree->Branch("mRcIdTruth",AnaT.mRcIdTruth,"mRcIdTruth[mNRcTrk]/I");
	trackTree->Branch("mRcAssoId",AnaT.mRcAssoId,"mRcAssoId[mNRcTrk]/I");
	trackTree->Branch("mRcPt",AnaT.mRcPt,"mRcPt[mNRcTrk]/F");
	trackTree->Branch("mRcPz",AnaT.mRcPz,"mRcPz[mNRcTrk]/F");
	trackTree->Branch("mRcEta",AnaT.mRcEta,"mRcEta[mNRcTrk]/F");
	trackTree->Branch("mRcPhi",AnaT.mRcPhi,"mRcPhi[mNRcTrk]/F");
	trackTree->Branch("mRcNhits",AnaT.mRcNhits,"mRcNhits[mNRcTrk]/I");
	trackTree->Branch("mRcNhitsPoss",AnaT.mRcNhitsPoss,"mRcNhitsPoss[mNRcTrk]/I");	
	trackTree->Branch("mRcNhitsPts",AnaT.mRcNhitsPts,"mRcNhitsPts[mNRcTrk]/I");
	trackTree->Branch("mRcNhitsSsd",AnaT.mRcNhitsSsd,"mRcNhitsSsd[mNRcTrk]/I");
	trackTree->Branch("mRcNhitsIst",AnaT.mRcNhitsIst,"mRcNhitsIst[mNRcTrk]/I");
	trackTree->Branch("mRcNhitsPxl2",AnaT.mRcNhitsPxl2,"mRcNhitsPxl2[mNRcTrk]/I");
	trackTree->Branch("mRcNhitsPxl1",AnaT.mRcNhitsPxl1,"mRcNhitsPxl1[mNRcTrk]/I");
	trackTree->Branch("mRcRndHitX",AnaT.mRcRndHitX,"mRcRndHitX[mNRcTrk][20]/F");
	trackTree->Branch("mRcRndHitY",AnaT.mRcRndHitY,"mRcRndHitY[mNRcTrk][20]/F");
	trackTree->Branch("mRcRndHitZ",AnaT.mRcRndHitZ,"mRcRndHitZ[mNRcTrk][20]/F");
	trackTree->Branch("mRcRndHitLX",AnaT.mRcRndHitLX,"mRcRndHitLX[mNRcTrk][20]/F");
	trackTree->Branch("mRcRndHitLY",AnaT.mRcRndHitLY,"mRcRndHitLY[mNRcTrk][20]/F");
	trackTree->Branch("mRcRndHitLZ",AnaT.mRcRndHitLZ,"mRcRndHitLZ[mNRcTrk][20]/F");
        trackTree->Branch("mRcRndHitLId",AnaT.mRcRndHitLId,"mRcRndHitLId[mNRcTrk][20]/I");
	trackTree->Branch("mRcRndHitPX",AnaT.mRcRndHitPX,"mRcRndHitPX[mNRcTrk][20]/F");
	trackTree->Branch("mRcRndHitPY",AnaT.mRcRndHitPY,"mRcRndHitPY[mNRcTrk][20]/F");
	trackTree->Branch("mRcRndHitPZ",AnaT.mRcRndHitPZ,"mRcRndHitPZ[mNRcTrk][20]/F");
        trackTree->Branch("mRcRndHitPId",AnaT.mRcRndHitPId,"mRcRndHitPId[mNRcTrk][20]/I");
	trackTree->Branch("mRcRndHitIdTruth",AnaT.mRcRndHitIdTruth,"mRcRndHitIdTruth[mNRcTrk][20]/I");
	trackTree->Branch("mRcTrackChi2",AnaT.mRcTrackChi2,"mRcTrackChi2[mNRcTrk]/F");
        trackTree->Branch("mDca2pXY",AnaT.mDca2pXY,"mDca2pXY[mNRcTrk]/F");
        trackTree->Branch("mDcaX",AnaT.mDcaX,"mDcaX[mNRcTrk]/F");
        trackTree->Branch("mDcaY",AnaT.mDcaY,"mDcaY[mNRcTrk]/F");
        trackTree->Branch("mDcaZ",AnaT.mDcaZ,"mDcaZ[mNRcTrk]/F");
	trackTree->Branch("mHelixX",AnaT.mHelixX,"mHelixX[mNRcTrk]/F");
	trackTree->Branch("mHelixY",AnaT.mHelixY,"mHelixY[mNRcTrk]/F");
	trackTree->Branch("mHelixZ",AnaT.mHelixZ,"mHelixZ[mNRcTrk]/F");
	trackTree->Branch("TPCrightTrack",AnaT.TPCrightTrack,"TPCrightTrack[mNRcTrk]/I");
	trackTree->Branch("ISTrightTrack",AnaT.ISTrightTrack,"ISTrightTrack[mNRcTrk]/I");
	trackTree->Branch("PXLrightTrack",AnaT.PXLrightTrack,"PXLrightTrack[mNRcTrk]/I");
	trackTree->Branch("sharedTpcHits",AnaT.sharedTpcHits,"sharedTpcHits[mNRcTrk]/I");
	trackTree->Branch("percentSharedTpcHits",AnaT.percentSharedTpcHits,"percentSharedTpcHits[mNRcTrk]/F");
	trackTree->Branch("ISTsharedTpcHits",AnaT.ISTsharedTpcHits,"ISTsharedTpcHits[mNRcTrk]/I");
	trackTree->Branch("ISTpercentSharedTpcHits",AnaT.ISTpercentSharedTpcHits,"ISTpercentSharedTpcHits[mNRcTrk]/F");
	trackTree->Branch("PXLsharedTpcHits",AnaT.PXLsharedTpcHits,"PXLsharedTpcHits[mNRcTrk]/I");
	trackTree->Branch("PXLpercentSharedTpcHits",AnaT.PXLpercentSharedTpcHits,"PXLpercentSharedTpcHits[mNRcTrk]/F"); 
        //hits
        trackTree->Branch("mNHits",&AnaT.mNHits,"mNHits/I");
        trackTree->Branch("mHitIdTruth",&AnaT.mHitIdTruth,"mHitIdTruth[mNHits]/I");   
        trackTree->Branch("mHitId",&AnaT.mHitId,"mHitId[mNHits]/I");
        trackTree->Branch("mHitX",AnaT.mHitX,"mHitX[mNHits]/F");
        trackTree->Branch("mHitY",AnaT.mHitY,"mHitY[mNHits]/F");
        trackTree->Branch("mHitZ",AnaT.mHitZ,"mHitZ[mNHits]/F");
        trackTree->Branch("mHitLX",AnaT.mHitLX,"mHitLX[mNHits]/F");
        trackTree->Branch("mHitLY",AnaT.mHitLY,"mHitLY[mNHits]/F");
        trackTree->Branch("mHitLZ",AnaT.mHitLZ,"mHitLZ[mNHits]/F");  
        //rc pair
        trackTree->Branch("mInvMass",&AnaT.mInvMass,"mInvMass[200][200]/F");
	//vertex refit
	trackTree->Branch("mRcVtxRefitX",&AnaT.mRcVtxRefitX,"mRcVtxRefitX/F");
	trackTree->Branch("mRcVtxRefitY",&AnaT.mRcVtxRefitY,"mRcVtxRefitY/F");
	trackTree->Branch("mRcVtxRefitZ",&AnaT.mRcVtxRefitZ,"mRcVtxRefitZ/F");
	trackTree->Branch("mNPTracksRefit",&AnaT.mNPTracksRefit,"mNPTracksRefit/I");

	return StMaker::Init();
}

Int_t StAnaSimMaker::InitRun(int runnumber)
{
  
  //PXL Db
  TDataSet *hitErrSet = GetDataBase("Calibrations/tracker/PixelHitError");
   if (!hitErrSet)
   {
      LOG_ERROR << "StPxlSimMaker - E - could not Get Calibrations/tracker." << endm;
      return kStErr;   
   }

  TObjectSet *pxlDbDataSet = (TObjectSet*)GetDataSet("pxl_db"); //pxlDb? need update
  if (!pxlDbDataSet)
    {
      mPxlDb = 0;
      LOG_ERROR << "StPxlSimMaker - E - pxlDb  is not available" << endm;
      return kStErr;
    } else 
    {
      mPxlDb = (StPxlDb *)pxlDbDataSet->GetObject();
    }
  
   St_HitError *pxlTableSet = (St_HitError*)hitErrSet->Find("PixelHitError");

   if (!pxlTableSet)
   {
      LOG_ERROR << "StPxlFastSim - E - PixelHitError is not available" << endm;
      return kStErr;
   }

   HitError_st* pxlHitError = pxlTableSet->GetTable();

   if (!pxlHitError)
   {
      LOG_ERROR << "StPxlFastSim - E - pxl hit table is not available in PixelHitError" << endm;
      return kStErr;
   }

   // please note that what is called local Y in the PXL sensor design
   // is actually called Z in STAR coordinates convention
   if (pxlHitError->coeff[0] <= 0 || pxlHitError->coeff[3] <= 0)
   {
      LOG_ERROR << "StPxlFastSim - E - negative or corrupted PXL hits errors in DB" << endm;
      return kStErr;
   }

   mResXPix = sqrt(pxlHitError->coeff[0]); // local x
   mResZPix = sqrt(pxlHitError->coeff[3]); // local Y
   mResYPix = 0;//sqrt(pxlHitError->coeff[2]); // needs to be updated in the DB later

  TDataSet *set = GetDataBase("Calibrations/tracker");
  St_HitError *istTableSet = (St_HitError *)set->Find("ist1HitError");
  HitError_st* istHitError = istTableSet->GetTable();
  resXIst1 = sqrt(istHitError->coeff[0]);
  resZIst1 = sqrt(istHitError->coeff[3]);

	return kStOK;
}

void StAnaSimMaker::Clear(Option_t *opt)
{
	StMaker::Clear();
}

Int_t StAnaSimMaker::Finish()
{

	mTupFile->Write();
	mTupFile->Close();
	return kStOK;
}

Int_t StAnaSimMaker::FinishRun(int runnumber)
{
	//IstSensorOnGlobal = 0;
	return kStOK;
}

Int_t StAnaSimMaker::Make()
{

  //IST Db
  THashList *IstSensorOnGlobal = new THashList(144,0);
  StIstDbMaker *istDbMaker = (StIstDbMaker *)GetMaker("istDb");  
  IstSensorOnGlobal = istDbMaker->GetRotations();

	mRcEvent=(StEvent *) GetInputDS("StEvent");
	if (!mRcEvent){
		gMessMgr->Warning() << "StAnaSimMaker::Make : No StEvent" << endl;
		return kStOK;        // if no event, we're done
	}
	mMcEvent = (StMcEvent*) GetDataSet("StMcEvent");
	if(!mMcEvent) return kStErr;
	if(!mMcEvent->primaryVertex()){
		cout << "\tno primary vertex from stmcevent" << endl;
		return kFALSE;
	}
	mMcVtxPos = &mMcEvent->primaryVertex()->position();

	bool validRcVertex = true;
	int nvtx = mRcEvent->numberOfPrimaryVertices();
	//cout << "Rc vertex # = " << nvtx << endl;
	if(nvtx <= 0) validRcVertex = false;
	int ntrk = 0;
	StPrimaryVertex *mRcVtx = 0;
	for(int i=0;i<nvtx;i++) {
	  StPrimaryVertex *aVtx = mRcEvent->primaryVertex(i);
	  //aVtx->Print();
	  if(!aVtx) continue;
	  int nTracks = aVtx->numberOfDaughters();
	  
	  if(nTracks>ntrk) {
	    mRcVtx = aVtx;
	    ntrk = nTracks;
	  }        
	}
	
	if(mRcVtx) {
	  mRcVtxPos = &mRcVtx->position();
	} else {
	  mRcVtxPos = 0;
	}
	



	StAssociationMaker* assoc = 0;
	assoc = (StAssociationMaker*) GetMaker("StAssociationMaker");
	//
	// multimaps
	//
	if (assoc) {
		//mRcHitMap   = assoc->rcTpcHitMap();
		mRcTrackMap = assoc->rcTrackMap();
		mMcTrackMap = assoc->mcTrackMap();
	}
	//rcSsdHitMapType* ssdHitMap   = assoc->rcSsdHitMap();

	//cout<<"begin StMiniMcMaker modification"<<endl;
	StMcIstHitCollection* istHitCol=mMcEvent->istHitCollection();
	//assert(istHitCol);
	for(int i=0;i<(int)istHitCol->numberOfLayers();i++){
	  if(istHitCol->layer(i)){
	    StMcIstLayerHitCollection* istLHCol=istHitCol->layer(i);
	    assert(istLHCol);
	    //cout<<"layer hit collection exists"<<endl;
	    //cout<<"# of hits: "<<istLHCol->numberOfHits()<<endl;
	    //StSPtrVecMcIstHit hv=istLHCol->hits();
	    ////cout<<"got hit vector"<<endl;
	    int nh=istLHCol->hits().size();
	    //cout<<"got number of hits"<<endl;
	    for(int j=0;j<nh;j++){
	      //cout<<"in hits loop"<<endl;
	      StMcIstHit* mihit=(StMcIstHit*)istHitCol->layer(i)->hits()[j];
	      ////cout<<mihit->position()<<endl;
	      istHitMap[mihit->key()]=mihit;
	    }
	  }
	}
	StMcPxlHitCollection* pixHitCol=mMcEvent->pxlHitCollection();
	//if (!pixHitCol) { cout << "No Pixel Hit Collection" << endl; return kFALSE;}
	for(unsigned int i=0; i<pixHitCol->numberOfSectors(); i++)
	{
	  StMcPxlSectorHitCollection* sectorHitCollection = pixHitCol->sector(i);
	  for(unsigned int j=0; j<sectorHitCollection->numberOfLadders(); j++)
	    {
	      StMcPxlLadderHitCollection* ladderHitCollection = sectorHitCollection->ladder(j);
	      for(unsigned int k=0; k<ladderHitCollection->numberOfSensors(); k++)
		{
		  StMcPxlSensorHitCollection* sensorHitCollection = ladderHitCollection->sensor(k);
		  for(unsigned int l=0; l<sensorHitCollection->hits().size(); l++)
		    {
		      StMcPxlHit* mphit = sensorHitCollection->hits()[l];
		      if (mphit) {pxlHitMap[mphit->key()]=mphit;}
		    }
		}
	    }
	}	
	cout<<"Fill TrackEvent object"<<endl;

	//event level
	AnaT.mIEvt        = (Int_t) mRcEvent->id();
	AnaT.mIRun        = (Int_t) mRcEvent->runId();

        StThreeVectorF PrimVtx(-999.,-999.,-999.);

	if(validRcVertex) {
          PrimVtx.set(mRcEvent->primaryVertex(0)->position().x(),mRcEvent->primaryVertex(0)->position().y(),mRcEvent->primaryVertex(0)->position().z());
	  AnaT.mRcVtxX[0]   = PrimVtx.x();
	  AnaT.mRcVtxY[0]   = PrimVtx.y();
	  AnaT.mRcVtxZ[0]   = PrimVtx.z();
	  AnaT.mRcVtxX[1]   = mRcVtxPos->x();
	  AnaT.mRcVtxY[1]   = mRcVtxPos->y();
	  AnaT.mRcVtxZ[1]   = mRcVtxPos->z();
	  AnaT.mNPTracks[0]    = mRcEvent->primaryVertex(0)->numberOfDaughters();
	  AnaT.mNPTracks[1]    = mRcVtx->numberOfDaughters();
	} else {
	  AnaT.mRcVtxX[0]   = -999.;
	  AnaT.mRcVtxY[0]   = -999.;
	  AnaT.mRcVtxZ[0]   = -999.;
	  AnaT.mRcVtxX[1]   = -999.;
	  AnaT.mRcVtxY[1]   = -999.;
	  AnaT.mRcVtxZ[1]   = -999.;
	  AnaT.mNPTracks[0]    = 0;
	  AnaT.mNPTracks[1]    = 0;
	}

	AnaT.mMcVtxX   = mMcVtxPos->x();
	AnaT.mMcVtxY   = mMcVtxPos->y();
	AnaT.mMcVtxZ   = mMcVtxPos->z();
	AnaT.mNMcPTracks  = mMcEvent->primaryVertex()->numberOfDaughters();
	AnaT.MagField     = mRcEvent->runInfo()->magneticField();

	cout << "MC primary vertex    = " << mMcVtxPos->x()<< "    "<< mMcVtxPos->y() << "    " << mMcVtxPos->z() <<endl;
	cout << "Reco primary vertex0 = " << PrimVtx <<endl;
	cout << "Reco primary vertex1 = " << mRcVtxPos->x()<< "    "<< mRcVtxPos->y() << "    " << mRcVtxPos->z() <<endl;
	cout << "Event variables filled" << endl;

	//track level
	StSPtrVecMcTrack mctracks=mMcEvent->tracks();
	//cout<<"got track vector"<<endl;
	cout<<"number of MC tracks: "<<mctracks.size()<<endl;

        //calculate reference global tracks in TPC
        StSPtrVecTrackNode& trackNode = mRcEvent->trackNodes();
        int ntotalTracks = trackNode.size();
        Int_t nglobaltracks = 0;
        Int_t ngrefmult     = 0;
        for (int i=0; i < ntotalTracks; i++) {
          StTrackNode *node = trackNode[i];
          if (!node) continue;
          StTrack *glTrack = node->track(global); 
          if (! glTrack) continue;

          double eta = glTrack->geometry()->momentum().pseudoRapidity();
          if(fabs(eta)>1.) continue;

          int nhits = glTrack->fitTraits().numberOfFitPoints();
          int nmax  = glTrack->numberOfPossiblePoints();
          double ratio = 1.*nhits/nmax;

          double gpt = glTrack->geometry()->momentum().perp();
          if(gpt>0.1 && fabs(eta)<0.5 && nhits>9) ngrefmult ++;

          if(nhits <= 15) continue;
          if(ratio < 0.52) continue;

          double gdca = glTrack->geometry()->helix().geometricSignedDistance(PrimVtx);
          if(fabs(gdca)>3.) continue;

          nglobaltracks ++; 
        }

        AnaT.mNGTracks    = nglobaltracks;
        AnaT.mNGRefMult   = ngrefmult;
	AnaT.mNMcTrk      = 0;
	AnaT.mNRcTrk      = 0;

	Int_t irc = 0; //reco tracks id
	Int_t imc = 0;

	Int_t nMcRefMult = 0;
        Int_t nRcRefMult = 0;
	for(int xyz=0;xyz<(int) mctracks.size();xyz++){
	  //cout<<"looping through tracks"<<endl;
	  //AnaT.clear();
	  //cout<<"cleared track object"<<endl;
	  
	  if(imc>=50000 || irc>=50000) { 
	    cout << " out of array size limit, exit" << endl;
	    continue;
	  }
	  
	  StMcTrack* tr = dynamic_cast<StMcTrack *>(mctracks[xyz]);
	  if(!tr) continue;
	  
	  if(tr->key()==0 && tr->geantId()==0) continue;  // not geant tracks
	  
	  const StThreeVectorF& mcMom       = tr->momentum();

       //   if(mcMom.perp() < 0.15) continue;
       //   if(fabs(mcMom.pseudoRapidity()) > 1.1) continue;

	  if(tr->startVertex()==mMcEvent->primaryVertex() && fabs(mcMom.pseudoRapidity())<0.5 
	     &&(tr->geantId()==2||tr->geantId()==3||tr->geantId()==8||tr->geantId()==9||tr->geantId()==11||tr->geantId()==12||tr->geantId()==14||tr->geantId()==15)) nMcRefMult++;
	  
	  AnaT.mMcId[imc]            = tr->key();	  
	  AnaT.mGeantId[imc]         = tr->geantId(); 
	  
	  if(!tr->parent()) {
	    AnaT.mParentGeantId[imc]    = -999; 
	    AnaT.mParentMcId[imc]       = -999;
	  } else { 
	    AnaT.mParentGeantId[imc]    = tr->parent()->geantId();
	    AnaT.mParentMcId[imc]       = tr->parent()->key();
	  }
	  
	  AnaT.mMcPt[imc]            = mcMom.perp();
	  AnaT.mMcPz[imc]            = mcMom.z();
	  AnaT.mMcEta[imc]           = mcMom.pseudoRapidity();
	  AnaT.mMcPhi[imc]           = mcMom.phi();
	  AnaT.mMcMass[imc]          = tr->fourMomentum().m();
	  
	  if(!tr->startVertex()) {
	    AnaT.mMcStartX[imc]      = -999.; 
	    AnaT.mMcStartY[imc]      = -999.; 
	    AnaT.mMcStartZ[imc]      = -999.; 
	  } else { 
	    AnaT.mMcStartX[imc]      = tr->startVertex()->position().x();
	    AnaT.mMcStartY[imc]      = tr->startVertex()->position().y();
	    AnaT.mMcStartZ[imc]      = tr->startVertex()->position().z();
	  } 

	  //TPC mc hits
	  StPtrVecMcTpcHit& myTpcHit = tr->tpcHits();
	  int numTpcHits             = myTpcHit.size();
	  AnaT.mMcNhits[imc]         = numTpcHits;

	  //HFT mc hits
	  StPtrVecMcSsdHit&   mySsdHit      = tr->ssdHits();
	  StPtrVecMcIstHit&   myIstHit      = tr->istHits();
	  StPtrVecMcPxlHit&   myPxlHit      = tr->pxlHits();
	  Int_t numSsdHits                  = mySsdHit.size();
	  Int_t numIstHits                  = myIstHit.size();
	  Int_t numPxlHits                  = myPxlHit.size();
	  if(numSsdHits > 5) numSsdHits     = 5;
	  if(numIstHits > 5) numIstHits     = 5;
	  if(numPxlHits > 5) numPxlHits     = 5;
	  AnaT.mMcNhitsSsd[imc]      = numSsdHits;
	  AnaT.mMcNhitsIst[imc]      = numIstHits;
	  //AnaT.mMcNhitsPxl[imc]      = numPxlHits;
	  
          //cout << "number of hits on layers: " << numSsdHits << "  " << numIstHits << "  " << numPxlHits << endl;

	  for(int irndhit=0; irndhit<20; irndhit++) {
            AnaT.mMcRndHitvId[imc][irndhit] = -999;
	    AnaT.mMcRndHitX[imc][irndhit]  = -999.;
	    AnaT.mMcRndHitY[imc][irndhit]  = -999.;
	    AnaT.mMcRndHitZ[imc][irndhit]  = -999.;
            AnaT.mMcRndHitLX[imc][irndhit] = -999.;
            AnaT.mMcRndHitLY[imc][irndhit] = -999.;
            AnaT.mMcRndHitLZ[imc][irndhit] = -999.;
	    AnaT.mMcAssHitX[imc][irndhit]  = -999.;
	    AnaT.mMcAssHitY[imc][irndhit]  = -999.;
	    AnaT.mMcAssHitZ[imc][irndhit]  = -999.;
            AnaT.mMcAssHitLX[imc][irndhit] = -999.;
            AnaT.mMcAssHitLY[imc][irndhit] = -999.;
            AnaT.mMcAssHitLZ[imc][irndhit] = -999.;
	    AnaT.mMcRndHitId[imc][irndhit] = -999;
	  }
	  
	  //SSD Hits not ready yet 
	  if(numSsdHits > 0) {
	    for(int issdhit=0; issdhit<numSsdHits; issdhit++) {
	      StMcSsdHit* ssdhit = mySsdHit[issdhit];
	      AnaT.mMcRndHitX[imc][issdhit]  = ssdhit->position().x();
	      AnaT.mMcRndHitY[imc][issdhit]  = ssdhit->position().y();
	      AnaT.mMcRndHitZ[imc][issdhit]  = ssdhit->position().z();
	    }
	  }
	  
	  ///StMcIstHit contains local position, need local->global transform of Ist hits
	  //new simulator for new 1-layer design
	  float smearedX = 0., smearedY = 0., smearedZ = 0.;

	  if(numIstHits > 0) {
	    if(istHitCol->layer(0)){
	      for(int iisthit=0; iisthit<numIstHits; iisthit++) {
		StMcIstHit* isthit = myIstHit[iisthit];
                int vid    = isthit->volumeId();
                int ladder = isthit->ladder()-1;
                int sensor = isthit->wafer();
                int istid  = 1000 + ladder * 6 + sensor;
		double local[3] = {isthit->position().x(),isthit->position().y(),isthit->position().z()};
		double global[3] = {isthit->position().x(),isthit->position().y(),isthit->position().z()};
		
                if(mSmear) { // smearing on
		  LOG_DEBUG << "Smearing start... " << endm;
		  smearedX=distortHit(local[0], resXIst1, kIstSensorActiveSizeRPhi/2.0);
		  smearedZ=distortHit(local[2], resZIst1, kIstSensorActiveSizeZ/2.0);
		  
                  local[0] = smearedX;
		  local[2] = smearedZ;
		  LOG_DEBUG <<Form("Smearing done...")<< endm;
                }
                else { //smearing off
		  LOG_DEBUG << "No smearing, but discreting ... " << endm;
		  //discrete hit local position (2D structure of IST sensor pads)
		  float rPhiPos   = kIstSensorActiveSizeRPhi/2.0 - local[0];
		  float zPos      = local[2] + kIstSensorActiveSizeZ/2.0;
		  short meanColumn  = (short)floor( zPos/kIstPadPitchColumn ) + 1;
		  short meanRow     = (short)floor( rPhiPos/kIstPadPitchRow ) + 1;
		  rPhiPos = (meanRow-1) * kIstPadPitchRow + 0.5 * kIstPadPitchRow; //unit: cm
		  zPos    = (meanColumn-1) * kIstPadPitchColumn + 0.5 * kIstPadPitchColumn; //unit: cm
		  local[0] = kIstSensorActiveSizeRPhi/2.0 - rPhiPos;
		  local[2] = zPos - kIstSensorActiveSizeZ/2.0;
                }
				
                TGeoHMatrix *combI=(TGeoHMatrix *)IstSensorOnGlobal->FindObject(Form("R%04i",istid));

		combI->LocalToMaster(local, global);

                AnaT.mMcRndHitvId[imc][iisthit+5] = vid;
		AnaT.mMcRndHitX[imc][iisthit+5]  = global[0];
		AnaT.mMcRndHitY[imc][iisthit+5]  = global[1];
		AnaT.mMcRndHitZ[imc][iisthit+5]  = global[2];
                AnaT.mMcRndHitLX[imc][iisthit+5] = local[0];
                AnaT.mMcRndHitLY[imc][iisthit+5] = local[1];
                AnaT.mMcRndHitLZ[imc][iisthit+5] = local[2];
		AnaT.mMcRndHitId[imc][iisthit+5] = istid;
	      }
	    }
	  }
	  
          smearedX = 0; smearedY = 0; smearedZ = 0;
          int npxl1mchits = 0, npxl2mchits = 0;
	  if(numPxlHits > 0) {
	    for(int ipxlhit=0; ipxlhit<numPxlHits; ipxlhit++) {
	      StMcPxlHit* pxlhit = myPxlHit[ipxlhit];

               Int_t vid    = pxlhit->volumeId();
               Int_t sector = pxlhit->sector()-1;
               Int_t ladder = pxlhit->ladder()-1;
               Int_t sensor = pxlhit->sensor()-1;

               Int_t pxlid = sector * 40 + ladder * 10 + sensor + 1;

               //cout << "pxl geo: " << sector << "  " << ladder << "  " << sensor << endl;

               Double_t localPixHitPos[3] = {pxlhit->position().x(), pxlhit->position().y(), pxlhit->position().z()};

               // please note that what is called local Y in the PXL sensor design
               // is actually called Z in STAR coordinates convention and vice-versa
	       if(mSmear) {
		 smearedX = distortHit(localPixHitPos[0], mResXPix, PXL_ACTIVE_X_LENGTH / 2.0);
		 smearedZ = distortHit(localPixHitPos[2], mResZPix, PXL_ACTIVE_Y_LENGTH / 2.0);
		 if (mResYPix) smearedY = distortHit(localPixHitPos[1], mResYPix, 0.0020); // Not properly constrained yet
		 else smearedY = localPixHitPos[1];
		 localPixHitPos[0] = smearedX;
		 localPixHitPos[2] = smearedZ;
		 localPixHitPos[1] = smearedY;
	       }

               Double_t GlobaPixHitPos[3] = {0, 0, 0};
	       if(mPxlDb) {
		 TGeoHMatrix *combP = (TGeoHMatrix *)mPxlDb->geoHMatrixSensorOnGlobal(sector+1,ladder+1,sensor+1);
		 combP->LocalToMaster(localPixHitPos,GlobaPixHitPos);
	       }

	      float hitsx = GlobaPixHitPos[0];
              float hitsy = GlobaPixHitPos[1];
              float hitsz = GlobaPixHitPos[2];

	      //              if(sqrt(hitsx*hitsx+hitsy*hitsy) > 5.) {
	      if(ladder > 0) {
                AnaT.mMcRndHitvId[imc][npxl2mchits+10]  = vid;
  	        AnaT.mMcRndHitX[imc][npxl2mchits+10]  = hitsx;
	        AnaT.mMcRndHitY[imc][npxl2mchits+10]  = hitsy;
	        AnaT.mMcRndHitZ[imc][npxl2mchits+10]  = hitsz;	      
	        AnaT.mMcRndHitLX[imc][npxl2mchits+10] = localPixHitPos[0];//local position
	        AnaT.mMcRndHitLY[imc][npxl2mchits+10] = localPixHitPos[1];
	        AnaT.mMcRndHitLZ[imc][npxl2mchits+10] = localPixHitPos[2];
	        AnaT.mMcRndHitId[imc][npxl2mchits+10] = pxlid;
                npxl2mchits ++;
              } else {
              //  AnaT.mMcRndHitvId[imc][npxl2mchits+15]  = vid;
                AnaT.mMcRndHitvId[imc][npxl1mchits+15]  = vid;
                AnaT.mMcRndHitX[imc][npxl1mchits+15]  = hitsx; 
                AnaT.mMcRndHitY[imc][npxl1mchits+15]  = hitsy; 
                AnaT.mMcRndHitZ[imc][npxl1mchits+15]  = hitsz;            
                AnaT.mMcRndHitLX[imc][npxl1mchits+15] = localPixHitPos[0];//local position
                AnaT.mMcRndHitLY[imc][npxl1mchits+15] = localPixHitPos[1];
                AnaT.mMcRndHitLZ[imc][npxl1mchits+15] = localPixHitPos[2];
                AnaT.mMcRndHitId[imc][npxl1mchits+15] = pxlid;
                npxl1mchits ++;
              }
	    }
	  }

          AnaT.mMcNhitsPxl1[imc] = npxl1mchits;
          AnaT.mMcNhitsPxl2[imc] = npxl2mchits;

	  //fill associated reconstructed hits
	  GetAssHit(mRcEvent, tr->key(), imc);

	  std::pair<mcTrackMapType::iterator,mcTrackMapType::iterator> itpair3=mMcTrackMap->equal_range(tr);
	  //StMcTrack* tr = dynamic_cast<StMcTrack *>(mctracks[xyz]);		
	  //if(itpair3.first != itpair3.second) cout << "pair found" << endl;
	  float mpsht=0.;
	  int hnshist,hnshpxl;
	  float hpshist,hpshpxl;
	  int hnsht;
	  const StGlobalTrack* tMatched=0;
	  const StGlobalTrack* gtIstHit=0;
	  const StGlobalTrack* gtPxlHit=0;
	  const StGlobalTrack* tTemp=0;
	  

	  for(mcTrackMapType::iterator mit3=itpair3.first;mit3!=itpair3.second;++mit3){
	    StTrackPairInfo* tpinfo=(*mit3).second;
	    tTemp=tpinfo->partnerTrack();
	    
            /*
	    StPtrVecHit tPxlHits=tTemp->detectorInfo()->hits(kPxlId);
	    if(tPxlHits.size() >= 2 && (tPxlHits.size()<=myPxlHit.size())) {
	      if((fabs(tPxlHits[0]->position().x() - myPxlHit[0]->position().x())<0.00172) &&
		 (fabs(tPxlHits[0]->position().y() - myPxlHit[0]->position().y())<0.00172) &&
		 (fabs(tPxlHits[1]->position().x() - myPxlHit[1]->position().x())<0.00172) &&
		 (fabs(tPxlHits[1]->position().y() - myPxlHit[1]->position().y())<0.00172)) {
		gtPxlHit=tTemp;
		hnshpxl=tpinfo->commonTpcHits();
		hpshpxl=tpinfo->percentOfPairedTpcHits();
	      }
	    }
	    
	    
	    StPtrVecHit tIstHits=tTemp->detectorInfo()->hits(kIstId);
	    for(int mm=0;mm<tIstHits.size();mm++){
	      if((tIstHits.size()<=myIstHit.size()) && (tIstHits[mm]==myIstHit[mm])){ 
	        gtIstHit=tTemp;
	        hnshist=tpinfo->commonTpcHits();
	        hpshist=tpinfo->percentOfPairedTpcHits();
	      }
	    }
            */
	    //      cout << " number of ist hits on this track = " << tIstHits.size() << endl;
	    //      cout << "track common tpc hits: " << tpinfo->commonTpcHits() << endl;
	    
	    if(tpinfo->percentOfPairedTpcHits()>mpsht){
	      tMatched=tTemp;
	      mpsht=tpinfo->percentOfPairedTpcHits();
	      hnsht=tpinfo->commonTpcHits();
	    }
	    
	  }
	  
	  const StGlobalTrack* tMatched1 = 0; 
	  tMatched1 = partnerTrack(mMcTrackMap,tr);
	  
	  if(tMatched){

            /*	    
	    if(gtPxlHit){
	      if(gtPxlHit->key()==tMatched->key() && gtPxlHit->impactParameter()==tMatched->impactParameter() && gtPxlHit->length()==tMatched->length() && gtPxlHit->detectorInfo()->numberOfPoints()==tMatched->detectorInfo()->numberOfPoints() && gtPxlHit->detectorInfo()->numberOfPoints(kTpcId)==tMatched->detectorInfo()->numberOfPoints(kTpcId)) AnaT.PXLrightTrack[irc]=1;
	      else AnaT.PXLrightTrack[irc]=0;
	      AnaT.PXLsharedTpcHits[irc]=hnshpxl;
	      AnaT.PXLpercentSharedTpcHits[irc]=hpshpxl;
	    }
	    else AnaT.PXLrightTrack[irc]=0;

	    if(gtIstHit){
	      if(gtIstHit->key()==tMatched->key() && gtIstHit->impactParameter()==tMatched->impactParameter() && gtIstHit->length()==tMatched->length() && gtIstHit->detectorInfo()->numberOfPoints()==tMatched->detectorInfo()->numberOfPoints() && gtIstHit->detectorInfo()->numberOfPoints(kTpcId)==tMatched->detectorInfo()->numberOfPoints(kTpcId)) AnaT.ISTrightTrack[irc]=1;
	      else AnaT.ISTrightTrack[irc]=0;
	      AnaT.ISTsharedTpcHits[irc]=hnshist;
	      AnaT.ISTpercentSharedTpcHits[irc]=hpshist;
	    }
	    else AnaT.ISTrightTrack[irc]=0;
            */
	    
	    //      cout<<"associated track found"<<endl;

            int nTpcHits = tMatched->detectorInfo()->hits(kTpcId).size();

	    AnaT.mRcId[irc]          = tMatched->key();
	    AnaT.mRcIdTruth[irc]     = tr->key();
	    AnaT.mRcAssoId[irc]      = imc;
	    AnaT.mRcPt[irc]          = tMatched->geometry()->momentum().perp();
	    AnaT.mRcPz[irc]          = tMatched->geometry()->momentum().z();
	    AnaT.mRcEta[irc]         = tMatched->geometry()->momentum().pseudoRapidity();
	    AnaT.mRcPhi[irc]         = tMatched->geometry()->momentum().phi();
	    AnaT.mRcNhits[irc]       = nTpcHits;
	    AnaT.mRcNhitsPoss[irc]   = tMatched->numberOfPossiblePoints(kTpcId);
	    AnaT.mRcNhitsPts[irc]    = tMatched->fitTraits().numberOfFitPoints(kTpcId);
            
	    StPtrVecHit PartnerSsdHits      = tMatched->detectorInfo()->hits(kSsdId);
	    StPtrVecHit PartnerIstHits      = tMatched->detectorInfo()->hits(kIstId);
	    StPtrVecHit PartnerPxlHits      = tMatched->detectorInfo()->hits(kPxlId);
	    
	    int nPartnerSsdHits             = (int) PartnerSsdHits.size();
	    int nPartnerIstHits             = (int) PartnerIstHits.size();
	    int nPartnerPxlHits             = (int) PartnerPxlHits.size();
	    if(nPartnerSsdHits > 5) nPartnerSsdHits = 5;
	    if(nPartnerIstHits > 5) nPartnerIstHits = 5;
	    if(nPartnerPxlHits > 5) nPartnerPxlHits = 5;
	    
	    AnaT.mRcNhitsSsd[irc]    = nPartnerSsdHits;
	    AnaT.mRcNhitsIst[irc]    = nPartnerIstHits;
	    //	    AnaT.mRcNhitsPxl[irc]    = nPartnerPxlHits;
	    
	    for(int irndhit=0; irndhit<20; irndhit++) {
	      AnaT.mRcRndHitX[irc][irndhit]  = -999.;
	      AnaT.mRcRndHitY[irc][irndhit]  = -999.;
	      AnaT.mRcRndHitZ[irc][irndhit]  = -999.;
	      AnaT.mRcRndHitLX[irc][irndhit]  = -999.;
	      AnaT.mRcRndHitLY[irc][irndhit]  = -999.;
	      AnaT.mRcRndHitLZ[irc][irndhit]  = -999.;
              AnaT.mRcRndHitLId[irc][irndhit] = -1;
	      AnaT.mRcRndHitPX[irc][irndhit]  = -999.;
	      AnaT.mRcRndHitPY[irc][irndhit]  = -999.;
	      AnaT.mRcRndHitPZ[irc][irndhit]  = -999.;
              AnaT.mRcRndHitPId[irc][irndhit] = -1;
	      AnaT.mRcRndHitIdTruth[irc][irndhit] = -999;
	    }
	    if(nPartnerSsdHits > 0) {
	      for(int issdhit=0; issdhit<nPartnerSsdHits; issdhit++) {
		AnaT.mRcRndHitX[irc][issdhit]  = PartnerSsdHits[issdhit]->position().x();
		AnaT.mRcRndHitY[irc][issdhit]  = PartnerSsdHits[issdhit]->position().y();
		AnaT.mRcRndHitZ[irc][issdhit]  = PartnerSsdHits[issdhit]->position().z();
		AnaT.mRcRndHitIdTruth[irc][issdhit] = PartnerSsdHits[issdhit]->idTruth();
	      }
	    }
	    
	    if(nPartnerIstHits > 0) {
	      for(int iisthit=0; iisthit<nPartnerIstHits; iisthit++) {
		StIstHit *hit = (StIstHit *)PartnerIstHits[iisthit];
		if(!hit) continue;
                int ladder = hit->getLadder()-1;
                int sensor = hit->getSensor();
                Int_t istid = 1000 + ladder * 6 + sensor;
		AnaT.mRcRndHitX[irc][iisthit+5]  = PartnerIstHits[iisthit]->position().x();
		AnaT.mRcRndHitY[irc][iisthit+5]  = PartnerIstHits[iisthit]->position().y();
		AnaT.mRcRndHitZ[irc][iisthit+5]  = PartnerIstHits[iisthit]->position().z();
		AnaT.mRcRndHitLX[irc][iisthit+5] = hit->localPosition(0);
		AnaT.mRcRndHitLY[irc][iisthit+5] = hit->localPosition(1);
		AnaT.mRcRndHitLZ[irc][iisthit+5] = hit->localPosition(2);
                AnaT.mRcRndHitLId[irc][iisthit+5] = istid;
		AnaT.mRcRndHitIdTruth[irc][iisthit+5] = PartnerIstHits[iisthit]->idTruth();
	      }
	    }

            int npxl1hit = 0, npxl2hit = 0;
	    if(nPartnerPxlHits > 0) {
	      for(int ipxlhit=0; ipxlhit<nPartnerPxlHits; ipxlhit++) {
		StPxlHit *hit = (StPxlHit *)PartnerPxlHits[ipxlhit];
		if(!hit) continue;
                int sector = hit->sector()-1;
		int ladder = hit->ladder()-1;
                int sensor = hit->sensor()-1;
                Int_t pxlid = sector * 40 + ladder * 10 + sensor + 1;
		float pxlhitx = PartnerPxlHits[ipxlhit]->position().x();
		float pxlhity = PartnerPxlHits[ipxlhit]->position().y();
		float pxlhitz = PartnerPxlHits[ipxlhit]->position().z();
		//if(sqrt(pxlhitx*pxlhitx+pxlhity*pxlhity) > 5.) {
		if(ladder > 0) {
		  AnaT.mRcRndHitX[irc][npxl2hit+10]  = pxlhitx;
		  AnaT.mRcRndHitY[irc][npxl2hit+10]  = pxlhity;
		  AnaT.mRcRndHitZ[irc][npxl2hit+10]  = pxlhitz;
		  AnaT.mRcRndHitLX[irc][npxl2hit+10] = hit->localPosition(0);
		  AnaT.mRcRndHitLY[irc][npxl2hit+10] = hit->localPosition(1);
		  AnaT.mRcRndHitLZ[irc][npxl2hit+10] = hit->localPosition(2);
                  AnaT.mRcRndHitLId[irc][npxl2hit+10] = pxlid;
		  AnaT.mRcRndHitIdTruth[irc][npxl2hit+10] = PartnerPxlHits[ipxlhit]->idTruth();
		  npxl2hit ++;
		} else {
		  AnaT.mRcRndHitX[irc][npxl1hit+15]  = pxlhitx;
		  AnaT.mRcRndHitY[irc][npxl1hit+15]  = pxlhity;
		  AnaT.mRcRndHitZ[irc][npxl1hit+15]  = pxlhitz;
		  AnaT.mRcRndHitLX[irc][npxl1hit+15] = hit->localPosition(0);
		  AnaT.mRcRndHitLY[irc][npxl1hit+15] = hit->localPosition(1);
		  AnaT.mRcRndHitLZ[irc][npxl1hit+15] = hit->localPosition(2);
                  AnaT.mRcRndHitLId[irc][npxl1hit+15] = pxlid;
		  AnaT.mRcRndHitIdTruth[irc][npxl1hit+15] = PartnerPxlHits[ipxlhit]->idTruth();
		  npxl1hit ++;
		}
	      }	      
	    }
	    AnaT.mRcNhitsPxl2[irc] = npxl2hit;
	    AnaT.mRcNhitsPxl1[irc] = npxl1hit;

            //cout << "Ist hits: " << nPartnerIstHits << "   Pxl2 hits: " << npxl2hit << "   Pxl1 hits: " << npxl1hit << endl;
	    
	    AnaT.mRcTrackChi2[irc] = tMatched->fitTraits().chi2();
	    
	    StPhysicalHelixD PartnerHelix = tMatched->geometry()->helix();
	    StPhysicalHelixD dcaGhelix = tMatched->dcaGeometry()->helix();

	    Projection(dcaGhelix, irc, IstSensorOnGlobal);

            double thePath = dcaGhelix.pathLength(PrimVtx);
            StThreeVectorF DCAPos = dcaGhelix.at(thePath);

            AnaT.mDcaX[irc]    = DCAPos.x();
            AnaT.mDcaY[irc]    = DCAPos.y();
            AnaT.mDcaZ[irc]    = DCAPos.z();
            AnaT.mDca2pXY[irc] = dcaGhelix.geometricSignedDistance(PrimVtx.x(),PrimVtx.y());
            //AnaT.mDca2pZ[irc]  = dcaGhelix.geometricSignedDistance(PrimVtx);
	    
	    AnaT.mHelixX[irc] = PartnerHelix.origin().x();
	    AnaT.mHelixY[irc] = PartnerHelix.origin().y();
	    AnaT.mHelixZ[irc] = PartnerHelix.origin().z();    
	    
	    AnaT.TPCrightTrack[irc]        = 1;
	    AnaT.percentSharedTpcHits[irc] = mpsht;
	    AnaT.sharedTpcHits[irc]        = hnsht;
	    //cout<<"PartnerHelix.h() = "<<PartnerHelix.h()<<endl;
	    if(fabs(tMatched->geometry()->momentum().pseudoRapidity()) < 0.5 && 
               fabs(tMatched->impactParameter()) < 3.0                       &&
               !(PartnerHelix.h() == 0)                                      &&
               nTpcHits > 9 ) nRcRefMult++; 
	    irc++;
	  } else AnaT.TPCrightTrack[irc] = 0;
	  
	  imc++;
	}
        AnaT.mRcMult  = nRcRefMult; 
	AnaT.mMcMult  = /*(Int_t) mMcEvent->numberOfPrimaryTracks()*10000+*/nMcRefMult;
	AnaT.mNRcTrk = irc;
	AnaT.mNMcTrk = imc;

        cout << "McTracks: " << imc << "   RcTracks: " << irc << endl;
	cout << "RcTotalTracks: " << ntotalTracks << endl;

	//filling hit branch
        int nhits = 0;
        StIstHitCollection* istRcHitCol=mRcEvent->istHitCollection();
        for(int i_ladder=0; i_ladder<kIstNumLadders; i_ladder++)
          {
            StIstLadderHitCollection* ladderRcHitCollection = istRcHitCol->ladder(i_ladder);
            for(int i_sensor=0; i_sensor<kIstNumSensorsPerLadder; i_sensor++)
              {
		StIstSensorHitCollection* sensorRcHitCollection = ladderRcHitCollection->sensor(i_sensor);
		Int_t istid = 1000 + i_ladder * 6 + i_sensor + 1;
		for(unsigned int l=0; l<sensorRcHitCollection->hits().size(); l++)
		  {
		    StIstHit* hit=(StIstHit*)sensorRcHitCollection->hits()[l];
		    //StIstDigiHit* hit = new StIstDigiHit(*mhit);
		    const StThreeVectorF &P = hit->position();
		    AnaT.mHitX[nhits]  = P.x();
		    AnaT.mHitY[nhits]  = P.y();
		    AnaT.mHitZ[nhits]  = P.z();
		    AnaT.mHitLX[nhits] = hit->localPosition(0);
		    AnaT.mHitLY[nhits] = hit->localPosition(1);
		    AnaT.mHitLZ[nhits] = hit->localPosition(2);
		    AnaT.mHitId[nhits] = istid;
		    AnaT.mHitIdTruth[nhits] = hit->idTruth();
		    nhits ++;
		  }
	      }
	  }

	StPxlHitCollection* pixRcHitCol=mRcEvent->pxlHitCollection();
	for(unsigned int i_sector=0; i_sector<pixRcHitCol->numberOfSectors(); i_sector++)
	  {
	    StPxlSectorHitCollection* sectorRcHitCollection = pixRcHitCol->sector(i_sector);
	    for(unsigned int i_ladder=0; i_ladder<sectorRcHitCollection->numberOfLadders(); i_ladder++)
	      {
		StPxlLadderHitCollection* ladderRcHitCollection = sectorRcHitCollection->ladder(i_ladder);
		for(unsigned int i_sensor=0; i_sensor<ladderRcHitCollection->numberOfSensors(); i_sensor++)
		  {
		    StPxlSensorHitCollection* sensorRcHitCollection = ladderRcHitCollection->sensor(i_sensor);
		    Int_t pxlid = i_sector * 40 + i_ladder * 10 + i_sensor + 1;
		    for(unsigned int l=0; l<sensorRcHitCollection->hits().size(); l++)
		      {
			StPxlHit* mhit = sensorRcHitCollection->hits()[l];
			const StThreeVectorF &P = mhit->position();
			float hitsx = P.x();
			float hitsy = P.y();
			float hitsz = P.z();
			AnaT.mHitX[nhits]  = hitsx;
			AnaT.mHitY[nhits]  = hitsy;
			AnaT.mHitZ[nhits]  = hitsz;	      
			AnaT.mHitLX[nhits] = mhit->localPosition(0);
			AnaT.mHitLY[nhits] = mhit->localPosition(1);
			AnaT.mHitLZ[nhits] = mhit->localPosition(2);
			AnaT.mHitId[nhits] = pxlid;
			AnaT.mHitIdTruth[nhits] = mhit->idTruth();
			nhits ++;
		      }
		  }
	      }
	  }

	AnaT.mNHits = nhits;	

        cout << "Number of hits stored: " << nhits << endl;
	
	//pair branch
        for (int i=0; i < 200; i++) {
	  for (int j=0; j < 200; j++) {
	    AnaT.mInvMass[i][j] = -999.0;
	  }
	}
	vector<int> indxDdaughterCand;

        for (int i=0; i < ntotalTracks; i++) {
          StTrackNode *inode = trackNode[i];
          if (!inode) continue;
          StTrack *iglTrack = inode->track(global); 
          if (! iglTrack) continue;
	  StPhysicalHelixD ihelix = iglTrack->geometry()->helix();

          for (int j=i+1; j < ntotalTracks; j++) {
            StTrackNode *jnode = trackNode[j];
            if (!jnode) continue;
            StTrack *jglTrack = jnode->track(global); 
            if (! jglTrack) continue;
	    StPhysicalHelixD jhelix = jglTrack->geometry()->helix();

	    pair<double,double> pathLength = ihelix.pathLengths(jhelix);
	    StThreeVectorD idcaMom = ihelix.momentumAt(pathLength.first,mRcEvent->runInfo()->magneticField()*pow(10,-14));
	    StThreeVectorD jdcaMom = jhelix.momentumAt(pathLength.second,mRcEvent->runInfo()->magneticField()*pow(10,-14));

	    int iGeantId=0; 
	    int jGeantId=0;
	    for(int xyz=0;xyz<(int) mctracks.size();xyz++){
	       StMcTrack* tr = dynamic_cast<StMcTrack *>(mctracks[xyz]);
	       if(!tr) continue;
	       if(tr->key()==0 && tr->geantId()==0) continue;  // not geant tracks

	       std::pair<mcTrackMapType::iterator,mcTrackMapType::iterator> itpair3=mMcTrackMap->equal_range(tr);
	       float mpsht=0.;
	       const StGlobalTrack* tMatched=0;
	       const StGlobalTrack* tTemp=0;
	       for(mcTrackMapType::iterator mit3=itpair3.first;mit3!=itpair3.second;++mit3){
	       	  StTrackPairInfo* tpinfo=(*mit3).second;
	    	  tTemp=tpinfo->partnerTrack();
	    	  if(tpinfo->percentOfPairedTpcHits()>mpsht){
		     tMatched=tTemp;
	      	     mpsht=tpinfo->percentOfPairedTpcHits();
	    	  }
	       }

	       if (tMatched && tMatched->key() == iglTrack->key()) iGeantId=tr->geantId();
	       if (tMatched && tMatched->key() == jglTrack->key()) jGeantId=tr->geantId();
	    }

	    if (iGeantId==12&&jGeantId==8){
	       StLorentzVectorD LVK1, LVpi1;
	       LVK1.setPx(idcaMom.x());
	       LVK1.setPy(idcaMom.y());
	       LVK1.setPz(idcaMom.z());
	       LVK1.setE(idcaMom.massHypothesis(0.493677));
	       LVpi1.setPx(jdcaMom.x());
	       LVpi1.setPy(jdcaMom.y());
	       LVpi1.setPz(jdcaMom.z());
	       LVpi1.setE(jdcaMom.massHypothesis(0.13957018));
	       StLorentzVectorD rcD1 = LVK1+LVpi1;
	       AnaT.mInvMass[i][j] = rcD1.m();
	       if (rcD1.m()>1.8&&rcD1.m()<1.92){ 
	       	  cout << "i: " << i << "   j: " << j << "   invMass: " << rcD1.m() << endl;
	       	  indxDdaughterCand.push_back(i);
	       	  indxDdaughterCand.push_back(j);
	       }
	    }

	    if (iGeantId==8&&jGeantId==12){
	       StLorentzVectorD LVK2, LVpi2;
	       LVK2.setPx(jdcaMom.x());
	       LVK2.setPy(jdcaMom.y());
	       LVK2.setPz(jdcaMom.z());
	       LVK2.setE(jdcaMom.massHypothesis(0.493677));
	       LVpi2.setPx(idcaMom.x());
	       LVpi2.setPy(idcaMom.y());
	       LVpi2.setPz(idcaMom.z());
	       LVpi2.setE(idcaMom.massHypothesis(0.13957018));
	       StLorentzVectorD rcD2 = LVK2+LVpi2;
	       AnaT.mInvMass[j][i] = rcD2.m();
	       if (rcD2.m()>1.8&&rcD2.m()<1.92){
	       	  cout << "j: " << j << "   i: " << i << "   invMass: " << rcD2.m() << endl;
	       	  indxDdaughterCand.push_back(i);
	       	  indxDdaughterCand.push_back(j);
	       }
	    }
	  }
        }
	cout << "indxDdaughterCand.size(): " << indxDdaughterCand.size() << endl;

	//StKFVertexMaker variables
	fNzBins = 2500;
	fNPasses = 2;
	fSpectrum = 0;
	fzWindow = 2;
	fVtxM = 0;
	fVerticesPass = 0;
	fTempLog = 2;
	fminBrent = 0;
	func = 0;
	mBeamLine = kFALSE;
	fc1 = 0;

	//StKFVertexMaker constr.
  	Int_t npeaks = 100;
  	Double_t zmin = -250;
  	Double_t zmax = 250;
  	//  StKFVertex::_debug = 1;
  	for (Int_t pass = 0; pass < fNPasses; pass++) {
    	    fVtxs[pass] = new TH1F(Form("Vtx%1i",pass),Form("z-dca distribution for pass = %1i",pass),fNzBins,zmin,zmax);
    	    fVtxs[pass]->SetDirectory(0);
    	    if (pass)  fVtxs[pass]->SetLineColor(5);
    	    fVtxs[pass]->SetDefaultSumw2();
    	    fVtxs[pass]->SetStats(0);
    	    fVtxKs[pass] = new TH1K(Form("VtxK%1i",pass),Form("z-dca distribution for pass = %1i",pass),fNzBins,zmin,zmax);
    	    fVtxKs[pass]->SetDirectory(0);
    	    fVtxKs[pass]->SetStats(0);
    	    fVtxKs[pass]->SetLineColor(2);
  	}
  	fVtxM = new TH1F("VtxM","MuDst reconstructed multiplicities versus Z",fNzBins,zmin,zmax);
  	fVtxM->SetDirectory(0);
  	fSpectrum = new TSpectrum(2*npeaks);
  	func = new ROOT::Math::Functor1D(&StAnaSimMaker::AnnelingFcn);
  	fminBrent = new ROOT::Math::GSLMinimizer1D();
  	fVerticesPass = new StKFVerticesCollection *[fNPasses+1];
  	memset (fVerticesPass, 0, (fNPasses+1)*sizeof(StKFVerticesCollection *));
  	fParticles = new TObjArray();
  	fParticles->SetOwner(kTRUE);
  	mVertexOrderMethod = orderByRanking; // change ordering by ranking

	//StKFVertexMaker Make part
 	StEvent* pEvent = dynamic_cast<StEvent*> (GetInputDS("StEvent"));
  	if (! pEvent) {
    	   LOG_WARN << "StKFVertexMaker::fit: no StEvent " << endm;
    	   return kStOK;        // if no event, we're done
  	}
  	Double_t bField = 0;
  	if (pEvent->runInfo()) bField = pEvent->runInfo()->magneticField();
  	KFParticle::SetField(bField);
  	if (mBeamLine) {
    	   St_vertexSeedC* vSeed = St_vertexSeedC::instance();
    	   Double_t x0   = vSeed->x0()  ; Double_t err_x0   = vSeed->err_x0();  
    	   Double_t y0   = vSeed->y0()  ; Double_t err_y0   = vSeed->err_y0();
    	   Double_t dxdz = vSeed->dxdz(); Double_t err_dxdz = vSeed->err_dxdz();
    	   Double_t dydz = vSeed->dydz(); Double_t err_dydz = vSeed->err_dydz();
    	   Double_t weight = vSeed->weight();
    	   if (err_x0 < 0.010) err_x0 = 0.010;
    	   if (err_y0 < 0.010) err_y0 = 0.010;
    	   static Bool_t firstTime = kTRUE;
    	   if (firstTime) {
      	      firstTime = kFALSE;
      	      LOG_INFO << "BeamLine Constraint: weight =  " << weight << endm;
      	      LOG_INFO << "x(z) = (" << x0 << " +- " << err_x0 << ") + (" << dxdz << " +- " << err_dxdz << ") * z" << endm;
      	      LOG_INFO << "y(z) = (" << y0 << " +- " << err_y0 << ") + (" << dydz << " +- " << err_dydz << ") * z" << endm;
    	   }
    	   static Double_t pZ = 1000;
    	   static MTrack track;
    	   Double_t xyzP[6] = {     x0,      y0, 0.,
			      	    pZ*dxdz, pZ*dydz, pZ};
    	   Double_t CovXyzP[21] = {
      	   	    err_x0*err_x0,
      	   	    0            ,err_y0*err_y0,
      		    0            ,0              , 0,
      		    0            ,0              , 0, (err_dxdz*pZ)*(err_dxdz*pZ),
      		    0            ,0              , 0,                             0, (err_dydz*pZ)*(err_dydz*pZ)
    	   };
    	   track.SetParameters(xyzP);
    	   track.SetCovarianceMatrix(CovXyzP);
    	   track.SetNDF(1);
    	   track.SetID(0);
    	   track.SetCharge(1);
    	   KFParticle *beam = new KFParticle(track, 321);
    	   fParticles->AddAt(beam, 0);
  	}
  	StSPtrVecTrackNode& ptrackNode = pEvent->trackNodes();
  	UInt_t nTracks = ptrackNode.size();
  	StTrackNode *node=0;
  	Int_t NGoodGlobals = 0;
  	map<Int_t,StTrackNode*> TrackNodeMap;
  	for (UInt_t i=0; i < nTracks; i++) {
	   int flagDdaughterCand = 0;
	   for (int j=0; j < indxDdaughterCand.size(); j++) {
	      if (indxDdaughterCand[j] == i) flagDdaughterCand = 1;
	   }
	   if (flagDdaughterCand == 1) continue;
    	   node = ptrackNode[i]; 
    	   if (!node) continue;
    	   StGlobalTrack  *gTrack = static_cast<StGlobalTrack *>(node->track(global));
    	   if (! gTrack) continue;
    	   const StDcaGeometry* dca = gTrack->dcaGeometry();
    	   if (! dca) continue;
    	   if (gTrack->flag()     <   0) continue;     // Bad fit
    	   if (gTrack->flag()     > 700) continue;     // FTPC
    	   if (gTrack->flag()%100 == 11) continue;     // Short track pointing to EEMC
    	   if ((gTrack->isWestTpcOnly() || gTrack->isEastTpcOnly()) && gTrack->isPostXTrack()) continue; // wrong TPC side track
    	   Int_t kg = gTrack->key();
    	   TrackNodeMap[kg] = node;
    	   KFParticle *particle = AddTrackAt(dca,kg);
    	   if (Debug() > 1) {
      	      cout << Form("particle: %4i/%4i ",NGoodGlobals,kg) << *particle << endl;
      	      cout << "Add to map[" << kg << "] node = " << TrackNodeMap[kg] << endl;
    	   }
    	   NGoodGlobals++;
	   //cout << "Good Globals: " << i << endl;
     	}
  	if (NGoodGlobals < 2) return 0;
  	Fit();
  	if (! Vertices()) return 0;
  	//
  	//  In case there are no tracks left we better quit
  	//
  	StSPtrVecTrackDetectorInfo& detInfoVec = pEvent->trackDetectorInfo();
  	Int_t Nvtx = Vertices()->NoVertices();
	cout << "NGoodGlobals: " <<  NGoodGlobals << "   # Vertices: " << Nvtx << endl;
  	for (Int_t l = 0; l < Nvtx; l++) {
    	    const StKFVertex *V = Vertices()->Vertex(l);
    	    if (! V) continue;
    	    //if (Debug() > 2) 
	    V->PrintW();
    	    // Store vertex
    	    StPrimaryVertex *primV  = new StPrimaryVertex;
    	    StThreeVectorF XVertex(&V->Vertex().X());
    	    primV->setPosition(XVertex);
    	    primV->setChiSquared(V->Vertex().Chi2()/V->Vertex().GetNDF());  
    	    primV->setProbChiSquared(TMath::Prob(V->Vertex().GetChi2(),V->Vertex().GetNDF()));
    	    Float_t cov[6];
    	    TCL::ucopy(&V->Vertex().Covariance(0),cov,6);
    	    primV->setCovariantMatrix(cov); 
    	    primV->setVertexFinderId(KFVertexFinder);
    	    primV->setFlag(1); // was not set earlier by this vertex finder ?? Jan
    	    primV->setRanking(333);
    	    primV->setNumTracksUsedInFinder(V->NoTracks());
    	    Bool_t beam = kFALSE;
    	    Double_t Pars[6];
    	    TCL::ucopy(&V->Vertex().X(), Pars, 6);
    	    Double_t Cov[21];
    	    TCL::ucopy(&V->Vertex().Covariance(0), Cov, 21);
	    if (!StiToolkit::instance()) new StiDefaultToolkit;
    	    StiHit *Vertex = StiToolkit::instance()->getHitFactory()->getInstance();
    	    Vertex->setGlobal(0, 0, V->Vertex().X(), V->Vertex().Y(), V->Vertex().Z(), 0);
    	    Vertex->setError(cov);
    	    Int_t NoTracks = V->NoTracks();
	    cout << "Vertex #: " <<  l << "   # tracks: " << NoTracks << endl;
    	    TArrayI indexT(NoTracks); Int_t *indexes = indexT.GetArray();
    	    TArrayI IdT(NoTracks);    Int_t *Ids     = IdT.GetArray();
    	    for (Int_t itk = 0; itk < NoTracks; itk++) {
      	    	Ids[itk] = 999999;
      		const StKFTrack*   track = V->Track(itk);
      		if (! track) continue;
      		const KFParticle   &P = track->Particle();
      		Int_t kg = P.GetID()%100000;
      		Ids[itk] = kg;
    	    }
    	    TMath::Sort(NoTracks,Ids,indexes,0);
    	    for (Int_t i = 0; i < NoTracks; i++) {
      	    	Int_t itk = indexes[i];
      		const StKFTrack*   track = V->Track(itk);
      		if (! track) continue;
      		const KFParticle   &P = track->Particle();
      		Int_t kg = P.GetID()%100000;
      		if (kg == 0) {
		   assert(!beam);
		   beam = kTRUE;
		   continue;
      		}
      		if (Debug() > 2) {
		   const KFParticle   *PO = track->OrigParticle();
		   const KFParticle *PS[2] = {PO, &P};
		   for (Int_t m = 0; m < 2; m++) {
	  	      if (! m) cout << "Original";
	  	      else     cout << "Fitted  ";
	  	      static const Char_t *names[6] = {"x","y","z","px","py","pz"};
	  	      for (Int_t j = 0; j < 6; j++) {
	    	    	 cout << Form(" %2s: %8.3f +/- %8.3f",names[j], PS[m]->GetParameter(j), PS[m]->GetCovariance(j,j) > 0 ? TMath::Sqrt(PS[m]->GetCovariance(j,j)) : -13);
	   	      }
	  	      cout << endl;
		   }
      	    	}
      	    	node = TrackNodeMap[kg];
    	    }
    	    if (beam ) primV->setBeamConstrained();
    	    //..... add vertex to the list
    	    //if (primV->numberOfDaughters() < 1) {
      	    //   delete primV;
    	    //} else {
      	       primV->setTrackNumbers();
      	       calculateRank(primV);
      	       pEvent->addPrimaryVertex(primV,orderByRanking);
    	    //}
  	}

	bool validRcVtxRefit = true;
	int nvtxRefit = pEvent->numberOfPrimaryVertices();
	//cout << "Refit vertex # = " << nvtxRefit << endl;
	if(nvtxRefit <= 0) validRcVtxRefit = false;
        StThreeVectorF PrimVtxRefit(-999.,-999.,-999.);

	if(validRcVtxRefit) {
          PrimVtxRefit.set(pEvent->primaryVertex(0)->position().x(),pEvent->primaryVertex(0)->position().y(),pEvent->primaryVertex(0)->position().z());
	  AnaT.mRcVtxRefitX   = PrimVtxRefit.x();
	  AnaT.mRcVtxRefitY   = PrimVtxRefit.y();
	  AnaT.mRcVtxRefitZ   = PrimVtxRefit.z();
	  AnaT.mNPTracksRefit = pEvent->primaryVertex(0)->numTracksUsedInFinder();
	} else {
	  AnaT.mRcVtxRefitX   = -999.;
	  AnaT.mRcVtxRefitY   = -999.;
	  AnaT.mRcVtxRefitZ   = -999.;
	  AnaT.mNPTracksRefit = 0;
	}

	cout << "Refit rc primary vertex0 = " << PrimVtxRefit <<endl;
	cout << "Refit rc primary vertex0 #tracks = " << pEvent->primaryVertex(0)->numTracksUsedInFinder() <<endl;

	//StKFVertexMaker Clear
  	for (Int_t pass = 0; pass < fNPasses; pass++) {
    	    fVtxKs[pass]->Reset();
    	    fVtxKs[pass]->SetMaximum();
    	    fVtxs[pass]->Reset();
    	    fVtxs[pass]->SetMaximum();
  	}
  	fVtx = fVtxs[0]; // << switch between types    Vtx = fVtxKs[0];
  	fVtxM->Reset();
  	fcVertices = 0;
  	fParticles->Clear("C");
        //End of StKFVertexMaker part


	LOG_DEBUG<<"filled track object"<<endm;
	trackTree->Fill();

        mMcTrackMap->clear();
        mRcTrackMap->clear();

        //SafeDelete(mMcTrackMap);   
        //SafeDelete(mRcTrackMap);

	istHitMap.clear();
	pxlHitMap.clear();

	return kStOK;
}

double StAnaSimMaker::distortHit(double x, double res, double detLength){
  double test;
  test = x + myRandom->gauss(0,res);

  while( fabs(test) > detLength){
      test = x + myRandom->gauss(0,res);
  }

  return test;
}


void StAnaSimMaker::Projection(StPhysicalHelixD dcaG_helix, int index, THashList* istRot)
{

  int nisthit = 0;
  for(int i_ladder = 0; i_ladder<24; i_ladder++) {
    for(int i_sensor = 0; i_sensor<6; i_sensor++) {
      if(nisthit > 4) continue;
      Int_t id = 1000 + i_ladder * 6 + i_sensor + 1;
      TGeoHMatrix *comb = (TGeoHMatrix *)istRot->FindObject(Form("R%04i",id));
      Double_t *rot = comb->GetRotationMatrix();
      Double_t *tra = comb->GetTranslation();
      const StThreeVectorD normal(rot[1], rot[4], rot[7]);
      const StThreeVectorD middle(tra);
      Double_t sh = dcaG_helix.pathLength(middle, normal);
      if(sh<=0 || sh > 1e3) continue;  // dcaG geometry, projection pathLength should be positive
      StThreeVectorD xyzG = dcaG_helix.at(sh);
      Double_t xyzGPred[3] = {xyzG.x(), xyzG.y(), xyzG.z()};
      Double_t uvPred[3];
      comb->MasterToLocal(xyzGPred,uvPred);
      if (TMath::Abs(uvPred[0]) > IST_Width_X/2. ) continue;
      if (TMath::Abs(uvPred[2]) > IST_Width_Z/2. ) continue;
      
      AnaT.mRcRndHitPX[index][nisthit+5]  = uvPred[0];
      AnaT.mRcRndHitPY[index][nisthit+5]  = uvPred[1];
      AnaT.mRcRndHitPZ[index][nisthit+5]  = uvPred[2];
      AnaT.mRcRndHitPId[index][nisthit+5] = id;
      nisthit ++;
    }
  }

  int npxl1hit = 0, npxl2hit = 0;

  for(int i_sector = 0; i_sector<10; i_sector++) {
    for(int i_ladder = 0; i_ladder<4; i_ladder++) {
      for(int i_sensor = 0; i_sensor<10; i_sensor++) {
	Int_t id = i_sector * 40 + i_ladder * 10 + i_sensor + 1;
	TGeoHMatrix *comb = (TGeoHMatrix *)mPxlDb->geoHMatrixSensorOnGlobal(i_sector+1, i_ladder+1, i_sensor+1);
	Double_t *rot = comb->GetRotationMatrix();
	Double_t *tra = comb->GetTranslation();
	const StThreeVectorD normal(rot[1], rot[4], rot[7]);
	const StThreeVectorD middle(tra);
	Double_t sh = dcaG_helix.pathLength(middle, normal);
	if(sh<=0 || sh > 1e3) continue;  // dcaG geometry, projection pathLength should be positive
	StThreeVectorD xyzG = dcaG_helix.at(sh);
	Double_t xyzGPred[3] = {xyzG.x(), xyzG.y(), xyzG.z()};
	Double_t uvPred[3];
	comb->MasterToLocal(xyzGPred,uvPred);
	if (TMath::Abs(uvPred[0]) > PXL_Width_X/2. ) continue;
	if (TMath::Abs(uvPred[2]) > PXL_Width_Z/2. ) continue;
	if(i_ladder==0 && npxl1hit<5) {
	  AnaT.mRcRndHitPX[index][npxl1hit+15]  = uvPred[0];
	  AnaT.mRcRndHitPY[index][npxl1hit+15]  = uvPred[1];
	  AnaT.mRcRndHitPZ[index][npxl1hit+15]  = uvPred[2];
          AnaT.mRcRndHitPId[index][npxl1hit+15] = id;
	  npxl1hit ++;
	}
        if(i_ladder>0 && npxl2hit<5) {
	  AnaT.mRcRndHitPX[index][npxl2hit+10]  = uvPred[0];
	  AnaT.mRcRndHitPY[index][npxl2hit+10]  = uvPred[1];
	  AnaT.mRcRndHitPZ[index][npxl2hit+10]  = uvPred[2];
          AnaT.mRcRndHitPId[index][npxl2hit+10] = id; 
	  npxl2hit ++;
	}
      }
    }
  }
}


void StAnaSimMaker::GetAssHit(StEvent *ev, int mcId, int index)
{

  int nisthits = 0;

  StIstHitCollection* isthitcol=ev->istHitCollection();
  if(!isthitcol) return;
  for(int j=0; j<kIstNumLadders; j++)
    {
      StIstLadderHitCollection* istladdercol = isthitcol->ladder(j);
      for(int k=0; k<kIstNumSensorsPerLadder; k++)
	{
	  StIstSensorHitCollection* istsensorcol = istladdercol->sensor(k);
          Int_t istid = 1000 + j * 6 + k + 1;	
	  for(unsigned int l=0; l<istsensorcol->hits().size(); l++)
	    {
	      StIstHit* hit=(StIstHit*)istsensorcol->hits()[l];
              if(!hit) continue;
	      //StIstDigiHit* hit = new StIstDigiHit(*misthit);           
	      //if(misthit->idTruth() == mcId && hit && nisthits<5) {
             if(hit->idTruth() == mcId && nisthits<5) {
		StThreeVectorF P = hit->position();
		AnaT.mMcAssHitX[index][nisthits+5]  = P.x();
		AnaT.mMcAssHitY[index][nisthits+5]  = P.y();
		AnaT.mMcAssHitZ[index][nisthits+5]  = P.z();	      
		AnaT.mMcAssHitLX[index][nisthits+5] = hit->localPosition(0);
		AnaT.mMcAssHitLY[index][nisthits+5] = hit->localPosition(1);
		AnaT.mMcAssHitLZ[index][nisthits+5] = hit->localPosition(2);
                AnaT.mMcAssHitId[index][nisthits+5] = istid;
		nisthits ++;
                //cout << "Index = " << index << "   hits = " << nisthits << endl;
	      }
              //delete hit;
	    } //loop hits on sensor
	}
    }
  

  int npxl2hits = 0, npxl1hits = 0;
  
  StPxlHitCollection* pixHitCol=ev->pxlHitCollection();
  //if (!pixHitCol) { cout << "No Pixel Hit Collection" << endl; return kFALSE;}
  for(unsigned int i=0; i<pixHitCol->numberOfSectors(); i++)
    {
      StPxlSectorHitCollection* sectorHitCollection = pixHitCol->sector(i);
      for(unsigned int j=0; j<sectorHitCollection->numberOfLadders(); j++)
	{
	  StPxlLadderHitCollection* pxlladdercol = sectorHitCollection->ladder(j);
	  for(unsigned int k=0; k<pxlladdercol->numberOfSensors(); k++)
	    {
	      StPxlSensorHitCollection* pxlsensorcol = pxlladdercol->sensor(k);
              Int_t pxlid = i * 40 + j * 10 + k + 1;
	      for(unsigned int l=0; l<pxlsensorcol->hits().size(); l++)
		{
		  StPxlHit* mpxlhit = pxlsensorcol->hits()[l];

		  if(mpxlhit && mpxlhit->idTruth() == mcId) {
		    const StThreeVectorF &P = mpxlhit->position();
		    float hitsx = P.x();
		    float hitsy = P.y();
		    float hitsz = P.z();
		    //if(sqrt(hitsx*hitsx+hitsy*hitsy) > 5. && npxl2hits<5) {
		    if(j>0 && npxl2hits<5) {
		      AnaT.mMcAssHitX[index][npxl2hits+10]  = hitsx;
		      AnaT.mMcAssHitY[index][npxl2hits+10]  = hitsy;
		      AnaT.mMcAssHitZ[index][npxl2hits+10]  = hitsz;	      
		      AnaT.mMcAssHitLX[index][npxl2hits+10] = mpxlhit->localPosition(0);
		      AnaT.mMcAssHitLY[index][npxl2hits+10] = mpxlhit->localPosition(1);
		      AnaT.mMcAssHitLZ[index][npxl2hits+10] = mpxlhit->localPosition(2);
                      AnaT.mMcAssHitId[index][npxl2hits+10] = pxlid;
		      npxl2hits ++;
		    }
		    //if(sqrt(hitsx*hitsx+hitsy*hitsy) < 5. && npxl1hits<5) {
		    if(j==0 && npxl1hits<5) {
		      AnaT.mMcAssHitX[index][npxl1hits+15]  = hitsx;
		      AnaT.mMcAssHitY[index][npxl1hits+15]  = hitsy;
		      AnaT.mMcAssHitZ[index][npxl1hits+15]  = hitsz;	      
		      AnaT.mMcAssHitLX[index][npxl1hits+15] = mpxlhit->localPosition(0);
		      AnaT.mMcAssHitLY[index][npxl1hits+15] = mpxlhit->localPosition(1);
		      AnaT.mMcAssHitLZ[index][npxl1hits+15] = mpxlhit->localPosition(2);
                      AnaT.mMcAssHitId[index][npxl1hits+15] = pxlid;
		      npxl1hits ++;
		    }
		  }
		}
	    }
	}
    }	
 
}

//StKFVertexMaker functions_______________________________________________________
void StAnaSimMaker::calculateRank(StPrimaryVertex *primV) {    
  // Calculation of veretx ranks to select 'best' (i.e. triggered)  vertex
  // Simpilfied version (w/o weighting)
  Float_t rank = primV->probChiSquared();
  static Float_t Wveto = 1;
  static Float_t Wmatch = 4;
  if (primV->isBeamConstrained()) rank += Wmatch;
  rank -= Wveto*primV->numPostXTracks();
  rank += Wmatch*primV->numTracksWithPromptHit();
  rank += Wmatch*primV->numTracksCrossingCentralMembrane();
  rank += Wmatch*primV->numMatchesWithCTB()
    -     Wveto*primV->numNotMatchesWithCTB();
  rank += Wmatch*primV->numMatchesWithBTOF() 
    -     Wveto*primV->numNotMatchesWithBTOF();
  rank += Wmatch*(primV->numMatchesWithBEMC() + primV->numMatchesWithEEMC());
    -     Wveto*(primV->numNotMatchesWithBEMC() + primV->numNotMatchesWithEEMC());
  if (primV->numTracksTpcWestOnly() > 0 && primV->numTracksTpcEastOnly() > 0) 
    rank += Wmatch*TMath::Min(primV->numTracksTpcWestOnly(),primV->numTracksTpcEastOnly());
  rank += 100.0 + primV->numTracksUsedInFinder();
  primV->setRanking(rank); 
  //if (Debug()) primV->Print();
  cout << "rank = " << rank << endl;
}
//________________________________________________________________________________
KFParticle *StAnaSimMaker::AddTrackAt(const StDcaGeometry *dca, Int_t kg) {
  fParticles->AddAtAndExpand (0,kg);
  if (! dca) return 0;
  Double_t xyzp[6], CovXyzp[21];
  dca->GetXYZ(xyzp,CovXyzp);
  static MTrack track;
  track.SetParameters(xyzp);
  track.SetCovarianceMatrix(CovXyzp);
  track.SetNDF(1);
  //    track.SetChi2(GlobalTracks_mChiSqXY[k]);
  track.SetID(kg);
  Int_t q   = 1;
  Int_t pdg = 211;
  if (dca->charge() < 0) {
    q = -1;
    pdg = -211;
  } 
  track.SetCharge(q);
  KFParticle *particle = new KFParticle(track, pdg);
  particle->SetID(kg);
  fParticles->AddAt(particle,kg);
  return particle;
}
//________________________________________________________________________________
void StAnaSimMaker::Fit() {
  if (Debug() != 2)  StKFVertex::SetDebug(Debug());
  fcVertices = 0;
  for (Int_t i = 0; i < fNPasses+1; i++) {
    SafeDelete(fVerticesPass[i]);
  }
  Int_t NGoodGlobals = Particles().GetLast();
  
  Double_t TempLog = fTempLog; // default Temperature Log
  for (Int_t pass = 0; pass < fNPasses; pass++) {
    Int_t nAccepted = 0;
    Double_t dZ = fVtxs[pass]->GetBinWidth(1);
    for (Int_t k = 0; k < NGoodGlobals; k++) {
      KFParticle *particle = (KFParticle *) Particles()[k];
      if (! particle) continue;
      Double_t pT;
      Double_t dpT;
      particle->GetPt(pT,dpT);
      Double_t offset = 0.5*particle->GetPz()/pT;
      Double_t SigmaZ = TMath::Sqrt(particle->Covariance(2,2) + offset*offset);
      SigmaZ += dZ;
      Double_t Z = particle->GetZ();
      fVtxKs[pass]->Fill(Z);
      Int_t bin1 = fVtxs[pass]->FindBin(Z - 5*SigmaZ);
      if (bin1 < 1) bin1 = 1;
      Int_t bin2 = fVtxs[pass]->FindBin(Z + 5*SigmaZ);
      if (bin2 > fNzBins) bin2 = fNzBins;
      Double_t z = fVtxs[pass]->GetBinCenter(bin1);
      for (Int_t bin = bin1; bin <= bin2; bin++, z += dZ) {
	fVtxs[pass]->Fill(z,(TMath::Erfc((z - Z - fzWindow)/SigmaZ) - TMath::Erfc((z - Z + fzWindow)/SigmaZ))/2.);
      }
      nAccepted++;
    }
    Double_t F = fVtxKs[pass]->GetEntries();
    if (F < 1) continue;
    fVtxKs[pass]->SetNormFactor(F/dZ);
    fVtx = fVtxs[0]; // << switch between types    Vtx = fVtxKs[0];
    TString opt("new");
    if (! Canvas()) opt = "goff";
    Int_t nfound = fSpectrum->Search(fVtx,3,opt,TMath::Min(0.1,5./NGoodGlobals));
    if (! nfound) continue;
    if (Canvas()) {
      Canvas()->cd();
      fVtxs[0]->Draw(); fVtxKs[0]->Draw("same");
      fVtxM->Draw("same");
      if (pass)    fVtx->Draw("same");
      Canvas()->Update();
    }
    if (StKFVertex::Debug() > 1) {
      LOG_INFO << "Found " << nfound 
	   << " candidate peaks to fit with " << NGoodGlobals
	   << " good globals from with " <<  nAccepted  << " accepted" << endm;
    }
    Double_t *zOfPeaks = new Double_t[nfound];
    Int_t npeaks = 0;
#if ROOT_VERSION_CODE > 336641 /* ROOT_VERSION(5,35,1) */
    Double_t *xpeaks = fSpectrum->GetPositionX();
#else
    Float_t *xpeaks = fSpectrum->GetPositionX();
#endif
    for (Int_t p = 0; p < nfound; p++) {
#if ROOT_VERSION_CODE > 336641 /* ROOT_VERSION(5,35,1) */
      Double_t xp = xpeaks[p];
#else
      Float_t xp = xpeaks[p];
#endif
      Int_t bin = fVtx->GetXaxis()->FindBin(xp);
      Double_t yp = fVtx->GetBinContent(bin);
      Double_t ep = fVtx->GetBinError(bin);
      if (yp-1.25*ep < 0) continue;
      zOfPeaks[npeaks] = xp;
      npeaks++;
    }
    if (StKFVertex::Debug() > 1) {
      LOG_INFO << "Found " << npeaks << " useful peaks to fit" << endm;
    }
    if (! npeaks) {delete [] zOfPeaks; break; }
    if (fVerticesPass[pass]) {delete fVerticesPass[pass]; fVerticesPass[pass] = 0;}
    fVerticesPass[pass] = new StKFVerticesCollection(npeaks, zOfPeaks);
    delete [] zOfPeaks;
    fcVertices = fVerticesPass[pass];
    fcVertices->DoTrack2VertexAssociation(Particles());
    if (! fcVertices->NoVertices())                         continue;
    if (AnnelingFcn(TMath::Exp(-TempLog)) <= 0) continue;
    if (! fcVertices->NoVertices())                         continue;
    fcVertices->UniqueTracks2VertexAssociation(); // Make track associated with only vertex
    //       fcVertices->PrintV(NoMuMcVertex,NoMuMcTrack,StMuMcVertex_time,
    // 		       StMuMcVertex_xyzV_mX1,StMuMcVertex_xyzV_mX2,StMuMcVertex_xyzV_mX3,
    // 		       StMuMcVertex_NoDaughters,StMuMcVertex_IdParTrk,StMuMcTrack_gePid);
  }
  if (! fVerticesPass[0]) return;
  if (fNPasses > 1 && Canvas()) {
    Canvas()->cd();
    fVtxs[1]->Draw("same"); 
    Canvas()->Update();
  }
  Int_t N1 = fVerticesPass[0]->NoVertices();
  if (! N1) return;
  if (fVerticesPass[1]) {
    *fVerticesPass[0] += *fVerticesPass[1];
  }
  fcVertices = fVerticesPass[0];
  fcVertices->MergeDuplicatedVertices();
  if (! fcVertices->NoVertices()) return;
  // Double_t Temperature = TMath::Exp(TempLog);
  TempLog = 5;
  Double_t Temperature = TMath::Exp(TempLog);
#if 1  
  // secondary vertices
  Int_t pass = fNPasses;
  if (fVerticesPass[pass]) {delete fVerticesPass[pass]; fVerticesPass[pass] = 0;}
  fVerticesPass[pass] = new StKFVerticesCollection();
  fcVertices = fVerticesPass[pass];
  StAnneling::SetTemperature(Temperature);
  for (Int_t k = 1; k < NGoodGlobals; k++) {
    KFParticle *particleK = (KFParticle *) Particles()[k];
    if (! particleK) continue;
    if (particleK->GetID() > 100000) continue;
    StKFVertex *vtx = 0;
    for (Int_t l = k+1; l < NGoodGlobals; l++) {
      KFParticle *particleL = (KFParticle *) Particles()[l];
      if (! particleL) continue;
      if (particleL->GetID() > 100000) continue;
      Double_t dist = particleK->GetDistanceFromParticle(*particleL);
      if (dist > 5.0) continue;
      if (! vtx) {
	vtx = new StKFVertex(fcVertices->NoVertices() + 1);
	vtx->AddTrack(new StKFTrack(k,particleK));
      }
      vtx->AddTrack(new StKFTrack(l,particleL));
    }
    if (! vtx) continue;
    vtx->Fit();
    Int_t N = vtx->NoTracks();
    if (! N) {delete vtx; vtx = 0; continue;}
    Double_t X = vtx->Vertex().X();
    Double_t Y = vtx->Vertex().Y();
    Double_t R = TMath::Sqrt(X*X + Y*Y);
    if (R > 200 ) {delete vtx; vtx = 0; continue;}
    Double_t prob = TMath::Prob(vtx->Vertex().GetChi2(),vtx->Vertex().GetNDF());
    if (N > 2 || prob > 1.e-3) {// Allow V2 to share tracks
      for (Int_t i = 0; i < N; i++) {
	KFParticle *particle = vtx->Track(i)->OrigParticle();;
	Int_t ID = particle->GetID()%100000 + 100000*vtx->ID();;
	particle->SetID(ID);
      }
    }
    fcVertices->AddVertex(vtx);
  }
  if (StKFVertex::Debug() > 1) {
    LOG_INFO << "Candidate for secondary vertices: " << fcVertices->NoVertices() << endm;
  }
  if ( fcVertices->NoVertices() ) {
    //       fcVertices->PrintV(NoMuMcVertex,NoMuMcTrack,StMuMcVertex_time,
    // 		       StMuMcVertex_xyzV_mX1,StMuMcVertex_xyzV_mX2,StMuMcVertex_xyzV_mX3,
    // 		       StMuMcVertex_NoDaughters,StMuMcVertex_IdParTrk,StMuMcTrack_gePid);
    *fVerticesPass[0] += *fVerticesPass[fNPasses];
  }
  // end of loop for secondary vertices
#endif
  fcVertices = fVerticesPass[0];
  fcVertices->Compress();
  if (! fcVertices->NoVertices()) return;
  fcVertices->MergeDuplicatedVertices();
  fminBrent->SetFunction(*func,TMath::Exp(-0.5*(TempLog)),TMath::Exp(-TempLog),1);
  if (! fminBrent->Minimize(10,0.1,0.1)) {
    LOG_WARN << "Temperature fit has failed" << endm;
    Temperature = 1;
  } else {
    Temperature = 1./fminBrent->XMinimum();
  }
  StAnneling::SetTemperature(Temperature);
  fcVertices->UniqueTracks2VertexAssociation(); // Make track associated with only vertex
  fcVertices->Fit(29,Canvas(),fVtx);
  if (Canvas()) Canvas()->Update();
}
//________________________________________________________________________________
Double_t StAnaSimMaker::AnnelingFcn(Double_t TInv) {
  if (! fcVertices) return 0;
  Double_t Temperature = 1./TInv;
  StAnneling::SetTemperature(Temperature);
  Double_t Chi2 =  fcVertices->Fit();
  if (StKFVertex::Debug()) 
    LOG_INFO << "StKFVertexMaker::AnnelingFcn\tTemperature = " << Temperature << " Chi2 = " << Chi2 << endm;
  return Chi2;
}
//________________________________________________________________________________


