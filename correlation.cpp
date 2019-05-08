#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "correlation.h"
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string>

using namespace std;


void correlation(){

TCanvas *c1=new TCanvas("c1","c1");
TFile *fin=new TFile("sort00335.root");
TTree *tree=(TTree*)fin->Get("tree");
//tree->Print();
tree->SetBranchAddress("timestamp",&timestamp);
tree->SetBranchAddress("xid",&xid);
tree->SetBranchAddress("dt",&dt);
tree->SetBranchAddress("yid",&yid);
tree->SetBranchAddress("de",&de);
tree->SetBranchAddress("me",&me);

//tree->Scan("xstrip:ystrip:timestamp:me","","colsize=30",100,10000);
/*
tree->Draw("de>>hde(3000,0,30000)");//all
tree->Draw("de>>hdem(3000,0,30000)","me>200");//coincide with mwpc
tree->Draw("de>>hdenm(3000,0,30000)","me<200");//mwpc veto
TH1F *hde, *hdem, *hdenm;
hde=(TH1F*)gROOT->FindObject("hde");
hdem=(TH1F*)gROOT->FindObject("hdem");
hdenm=(TH1F*)gROOT->FindObject("hdenm");
hde->Draw();
hdem->SetLineColor(kGreen);
hdenm->SetLineColor(kRed);
hdem->Draw("same");
hdenm->Draw("same");
c1->Draw();
cout<<hdem->Integral(600,750)<<endl;
hdenm->GetXaxis()->SetRangeUser(5000,8000);
hdenm->Draw();
c1->Draw();

tree->Draw("yid:xid>>(48,0,48,128,0,128)","","colz");
c1->SetLogy(0);
c1->SetLogz(0);
c1->Draw();*/

//mark the decay events and implanted ions by mwpc

dssd ds;
multimap<ULong64_t, dssd> multimapimp, multimapdec;//implantaion, decay

Long64_t nentries = tree->GetEntriesFast();
cout<<nentries<<endl;
for (Long64_t jentry=0; jentry<nentries;jentry++) {
    tree->GetEntry(jentry);
	timestamp=dt;
        ds.energy = de;
        ds.xstrip = xid;
        ds.ystrip = yid;
    //cout<<timestamp<<"  "<<de<<"  "<<xid<<"  "<<yid<<"  "<<me<<endl;

    if(me>200) multimapimp.insert(pair<ULong64_t,dssd>(timestamp,ds));
        else multimapdec.insert(pair<ULong64_t,dssd>(timestamp,ds));
}

cout<<"The number of implantation/decay : "<<multimapimp.size()<<"  "<<multimapdec.size()<<endl;
Int_t cnt=0;
ULong64_t twindow = 100000000000;//100s
ULong64_t twl,twh;
TH2F *hedt=new TH2F("hedt","decayenergy vs decaytime",200,0,1e11,300,5000,8000);
multimap<ULong64_t,dssd>::iterator its,jts;

for(its=multimapimp.begin(); its!=multimapimp.end();its++) {
	twl=its->first;
	twh=its->first+twindow;
     for(jts=multimapdec.lower_bound(twl);jts->first<twh;jts++){
	if(jts==multimapdec.end()) break;
        if(abs(its->second.xstrip-jts->second.xstrip)>1) continue;
        if(abs(its->second.ystrip-jts->second.ystrip)>1) continue;
        ULong64_t decaytime = jts->first-its->first;
        hedt->Fill(decaytime,jts->second.energy);
    }

  cnt++;
  if(cnt%10000==0) cout<<"Processing ..."<<cnt<<endl;
}
hedt->Draw("colz");
c1->SetLogz();
c1->Draw();

}
