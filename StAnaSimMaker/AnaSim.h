
struct TupTrackEntry
{
  //Event level
  Int_t mIEvt;
  Int_t mIRun;
  Int_t mMcMult;
  Int_t mRcMult;
  Float_t mMcVtxX; 
  Float_t mMcVtxY; 
  Float_t mMcVtxZ;
  Int_t   mNMcPTracks; 
  Float_t mRcVtxX[2]; //0: highest ranking; 1: largest No. of daughters
  Float_t mRcVtxY[2];
  Float_t mRcVtxZ[2];
  Int_t   mNPTracks[2];
  Int_t   mNGTracks;
  Int_t   mNGRefMult;
  Float_t  MagField;

  //mc track quantities
  Int_t mNMcTrk;
  Int_t mMcId[50000];
  Int_t mGeantId[50000];
  Int_t mParentMcId[50000];
  Int_t mParentGeantId[50000];
  float mMcPt[50000];
  float mMcPz[50000];
  float mMcEta[50000];
  float mMcPhi[50000];
  float mMcMass[50000];
  float mMcStartX[50000];
  float mMcStartY[50000];   
  float mMcStartZ[50000];
  Int_t mMcNhits[50000];
  Int_t mMcNhitsSsd[50000];
  Int_t mMcNhitsIst[50000];
  Int_t mMcNhitsPxl2[50000]; //outer
  Int_t mMcNhitsPxl1[50000]; //inner
  Int_t mMcRndHitvId[50000][20];
  float mMcRndHitX[50000][20]; //Global hits
  float mMcRndHitY[50000][20];
  float mMcRndHitZ[50000][20];
  float mMcRndHitLX[50000][20]; //Local hits
  float mMcRndHitLY[50000][20];
  float mMcRndHitLZ[50000][20];
  Int_t mMcRndHitId[50000][20];
  float mMcAssHitX[50000][20]; //associated reconstructed hits
  float mMcAssHitY[50000][20];
  float mMcAssHitZ[50000][20];
  float mMcAssHitLX[50000][20]; //reconstructed Local hits
  float mMcAssHitLY[50000][20];
  float mMcAssHitLZ[50000][20]; 
  Int_t mMcAssHitId[50000][20];

  //properties of reco global track
  Int_t mNRcTrk;
  Int_t mRcId[20000];
  Int_t mRcIdTruth[20000]; //associated Mc Track ID
  Int_t mRcAssoId[20000];  //associated Mc Track index
  float mRcPt[20000];
  float mRcPz[20000];
  float mRcEta[20000];
  float mRcPhi[20000];
  Int_t mRcNhits[20000];
  Int_t mRcNhitsPoss[20000];
  Int_t mRcNhitsPts[20000];
  Int_t mRcNhitsSsd[20000];
  Int_t mRcNhitsIst[20000];
  Int_t mRcNhitsPxl2[20000];
  Int_t mRcNhitsPxl1[20000];
  float mRcRndHitX[20000][20];
  float mRcRndHitY[20000][20];
  float mRcRndHitZ[20000][20];
  float mRcRndHitLX[20000][20];
  float mRcRndHitLY[20000][20];
  float mRcRndHitLZ[20000][20];
  Int_t mRcRndHitLId[20000][20];
  float mRcRndHitPX[20000][20]; //project position, SSD->IST->PXL2->PXL1
  float mRcRndHitPY[20000][20];
  float mRcRndHitPZ[20000][20];
  Int_t mRcRndHitPId[20000][20];
  Int_t mRcRndHitIdTruth[20000][20];
  float mRcTrackChi2[20000];
  float mDca2pXY[20000];
  float mDcaX[20000];
  float mDcaY[20000]; 
  float mDcaZ[20000]; 
  float mHelixX[20000];
  float mHelixY[20000];
  float mHelixZ[20000];
  //is there a reco track associated to this mc track?
  int TPCrightTrack[20000];//0 if no, 1 if yes
  int ISTrightTrack[20000];
  int PXLrightTrack[20000];

  //how well does this track match to the mc track?
  int sharedTpcHits[20000];
  float percentSharedTpcHits[20000];
  int ISTsharedTpcHits[20000];
  float ISTpercentSharedTpcHits[20000];
  int PXLsharedTpcHits[20000];
  float PXLpercentSharedTpcHits[20000];
  
  //Rc Hits (IST, PXL)
  Int_t mNHits;
  Int_t mHitIdTruth[50000];
  Int_t mHitId[50000];
  float mHitX[50000];
  float mHitY[50000];
  float mHitZ[50000];
  float mHitLX[50000];
  float mHitLY[50000];
  float mHitLZ[50000];
  
  //pairs quantities
  Float_t  mInvMass[200][200];

  //vertex refit
  Float_t mRcVtxRefitX;
  Float_t mRcVtxRefitY;
  Float_t mRcVtxRefitZ;
  Int_t   mNPTracksRefit;


};
  TupTrackEntry AnaT;

