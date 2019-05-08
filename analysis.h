//head file
Double_t energy;
ULong64_t timestamp;
TTree *opt;

//output variable
Int_t xid;
Int_t yid;
Double_t de;
Float_t me;
ULong64_t dt;

void BranchOpt(){
opt->Branch("timestamp",&timestamp,"timestamp/l");
opt->Branch("dt",&dt,"dt/l");
opt->Branch("xid",&xid,"xid/I");
opt->Branch("yid",&yid,"yid/I");
opt->Branch("de",&de,"de/D");
opt->Branch("me",&me,"me/F");
}

void ResetOpt(){
xid=0;
yid=0;
de=0;
me=0;
timestamp=0;
}
