#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "analysis.h"
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string>

using namespace std;

void analysis(){
  // Open files using class TFile
  TFile *fds=new TFile("dssd.root");
  TTree *tree=(TTree*)fds->Get("tree");
  TCanvas *c1=new TCanvas("c1","c1");
  Double_t energy;
  ULong64_t timestamp;
  Int_t side,strip;
  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("timestamp",&timestamp);
  tree->SetBranchAddress("side",&side);
  tree->SetBranchAddress("strip",&strip);
  //tree->Print();
  //tree->Scan("side:strip:timestamp:energy","","colsize=15",10,1);
  //get timestamp table for front(side=0) and back(side=1) side
  struct dssd
  {
    Double_t energy;
    Int_t strip;
    Bool_t flag;
  };
  //sort the events of frant and back sides by timestamp
  multimap<ULong64_t, dssd> multimapf1, multimapb;//multimap for front and back side
  dssd ds;
  Long64_t nentries = tree->GetEntriesFast();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    tree->GetEntry(jentry);
    ds.energy = energy;
    ds.strip = strip;
    ds.flag = true;//flag used for search of correlated events
    if(side==0) multimapf1.insert(pair<ULong64_t,dssd>(timestamp,ds));
    if(side==1) multimapb.insert(pair<ULong64_t,dssd>(timestamp,ds));
    if(nentries%100000==0) cout<<jentry<<endl;
  }
  cout<<"Total number of frontside/backside events = "<<multimapf1.size()<<" / "<<multimapb.size()<<endl;    
  fds->Close();

  //Checking the events
  TString sout;
  auto ifts=multimapf1.begin();//multimap<ULong64_t, dssd>::iterator ifts=;
  auto ibts=multimapb.begin();
  TH2F *hxyraw=new TH2F("hxyraw","front-back correlation for unmatched data",3000,0,30000,3000,0,30000);
  for(int i=0; i<TMath::Min(multimapf1.size(),multimapb.size()); i++) {
    if(i<20) {
        sout.Form("%2d  %15llu  %15llu %3d    %3d   %5.1f    %5.1f %5d %5d"
                  , i, ifts->first, ibts->first, ifts->second.strip, ibts->second.strip, ifts->second.energy, ibts->second.energy,ifts->second.flag,ifts->second.flag);
       //cout<<sout<<endl;
    }
    hxyraw->Fill(ifts->second.energy,ibts->second.energy);
    ifts++;
    ibts++;
  }
  //hxyraw->Draw("colz");
  //c1->Draw();
  TH1F *hferaw1=(TH1F*) hxyraw->ProjectionY("hferaw");
  //hferaw1->Draw();
  //c1->Draw();

  //searching for time windows between front and back sides
  ULong64_t twindow=100000;//ns
  ULong64_t toffset=0;//ns
  TH1I *hdt=new TH1I("hdt","fts-bts distribution(ns)",20000,-100000,100000);
  for(ifts=multimapf1.begin(); ifts!=multimapf1.end();ifts++) {
     ibts=multimapb.lower_bound( (ifts->first-toffset)-twindow);
    while(ibts!=multimapb.end()) {
        if(ibts->first >= (ifts->first-toffset) + twindow) break;
        int dt=(ifts->first-toffset)-ibts->first;
        hdt->Fill(dt);
        ibts++;
    }   
  }
  //hdt->Draw();
  //c1->SetLogy();
  //c1->Draw();
  hdt->GetXaxis()->SetRangeUser(-1000,1000);
  //hdt->Draw();
  //c1->Draw();

  //corrected the timestamp
  toffset=-550;
  multimap<ULong64_t,dssd> multimapf;
  for(auto ifts=multimapf1.begin(); ifts!=multimapf1.end();ifts++) 
    multimapf.insert(pair<ULong64_t,dssd>(ifts->first-toffset, ifts->second));

  //energy correlation of front and back sides
  TH2F *hxyc=new TH2F("hxyc","front-back correlation within coincidence window",3000,0,30000,3000,0,30000);
  twindow=100;//ns
  Int_t fnum;
  Int_t bnum;
  for(ifts=multimapf.begin(); ifts!=multimapf.end();ifts++) {
     ibts=multimapb.lower_bound( (ifts->first)-twindow);
     auto ifts1=multimapf.lower_bound( (ifts->first)-twindow);
     fnum=0;
     while(ifts1!=multimapf.end()){
	if(ifts1->first >= (ifts->first) + twindow) break;
	fnum++;
	ifts1++;
     } 
     if(fnum>=2) continue;
     //cout<<fnum<<endl;
     while(ibts!=multimapb.end()) {

        if(ibts->first >= (ifts->first) + twindow) break;
        int dt=ifts->first-ibts->first;
        hxyc->Fill(ifts->second.energy,ibts->second.energy);
        ibts++;
    }   
  }
  //hxyc->Draw("colz");
  //c1->SetLogy(0);
  //c1->SetLogz();
  //c1->Draw();
  
  //energy correlation of neighboring strips in the back side
  Int_t xhit,yhit;
  Int_t xstrip[100],ystrip[100];
  ULong64_t xtimestamp[100],ytimestamp[100];
  Float_t xenergy[100],yenergy[100];
  TH2F *hx2=new TH2F("hx2"," x0-x1 energy correlation",500,0,10000,500,0,10000);
  TH1I *hdtx=new TH1I("hdtx"," tx0-tx1 ",20,-150,50);

  twindow=100;//ns
  ibts=multimapb.begin();
  int nx[5];
  while(ibts!=multimapb.end() ) {        
        xhit=0;
        ULong64_t t0=ibts->first;
        while (ibts!=multimapb.end() ) {
            if(ibts->first > t0+2*twindow) break;
            xtimestamp[xhit]= ibts->first;
            xstrip[xhit]= ibts->second.strip;
            xenergy[xhit]=ibts->second.energy;
            xhit++;
            ibts++;
        }
    
    if(xhit<5) nx[xhit]++;
    if(xhit==2) {
        hx2->Fill(xenergy[0],xenergy[1]);
        Int_t dt=xtimestamp[0]-xtimestamp[1];
        hdtx->Fill(dt);
    }
  }
  hx2->Draw("colz");
  c1->SetLogy(0);
  c1->Draw();
  //hdtx->Draw();
  //c1->SetLogy();
  //c1->Draw();
  //add back the energies of neighboring strips
  multimap<ULong64_t,dssd> multimapfn,multimapbn;
  ULong64_t t0;
  int mby[5],mfx[5];
  twindow=100;//ns
  ibts=multimapb.begin();
  while(ibts!=multimapb.end() ) {        
        yhit=0;
        t0=ibts->first;
        while (ibts!=multimapb.end() ) {
            if(ibts->first > t0+2*twindow) break;
            ytimestamp[yhit]= ibts->first;
            ystrip[yhit]= ibts->second.strip;
            yenergy[yhit]=ibts->second.energy;
            yhit++;
            ibts++;
        }
    if(yhit==1) {
        ds.strip=ystrip[0];
        ds.energy=yenergy[0];
        multimapbn.insert(pair<ULong64_t,dssd>(ytimestamp[0],ds));
        mby[1]++;
    }
    if(yhit==2 && abs(ystrip[0]-ystrip[1])==1) {
        ds.energy=yenergy[0]+yenergy[1];
        ds.strip=ystrip[1];
        t0=ytimestamp[1];
        if(yenergy[0]>yenergy[1]) {
            ds.strip=ystrip[0];
            t0=ytimestamp[0];
        }
        multimapbn.insert(pair<ULong64_t,dssd>(t0,ds));
        mby[2]++;
    }
    if(yhit==2 && abs(ystrip[0]-ystrip[1])>1) {    
        ds.strip=ystrip[0];
        ds.energy=yenergy[0];
        //mapbn.insert(pair<ULong64_t,dssd>(ytimestamp[0],ds));    
        ds.strip=ystrip[1];
        ds.energy=yenergy[1];
        //mapbn.insert(pair<ULong64_t,dssd>(ytimestamp[1],ds));  
        mby[3]++;
        }
  }
  cout<<mby[1]<<"  "<<mby[2]<<"  "<<mby[2]<<endl;

  twindow=100;//ns
  ifts=multimapf.begin();
  while(ifts!=multimapf.end() ) {        
        xhit=0;
        t0=ifts->first;
        while (ifts!=multimapf.end() ) {
            if(ifts->first > t0+2*twindow) break;
            xtimestamp[xhit]= ifts->first;
            xstrip[xhit]= ifts->second.strip;
            xenergy[xhit]=ifts->second.energy;
            xhit++;
            ifts++;
        }
    if(xhit==1) {
        ds.strip=xstrip[0];
        ds.energy=xenergy[0];
        multimapfn.insert(pair<ULong64_t,dssd>(xtimestamp[0],ds));
        mfx[1]++;
    }
    if(xhit==2 && abs(xstrip[0]-xstrip[1])==1) {
        ds.energy=xenergy[0]+xenergy[1];
        ds.strip=xstrip[1];
        t0=xtimestamp[1];
        if(xenergy[0]>xenergy[1]) {
            ds.strip=xstrip[0];
            t0=xtimestamp[0];
        }
        multimapfn.insert(pair<ULong64_t,dssd>(t0,ds));
        mfx[2]++;
    }
    if(xhit==2 && abs(xstrip[0]-xstrip[1])>1) {    
        ds.strip=xstrip[0];
        ds.energy=xenergy[0];
        //mapfn.insert(pair<ULong64_t,dssd>(xtimestamp[0],ds));    
        ds.strip=xstrip[1];
        ds.energy=xenergy[1];
        //mapfn.insert(pair<ULong64_t,dssd>(xtimestamp[1],ds));  
        mfx[3]++;
        }
  }
  cout<<mfx[1]<<"  "<<mfx[2]<<"  "<<mfx[3]<<endl;
//front-back correlation for corrected data
  TH2F *hxy=new TH2F("hxy","front-back correlation for corrected data",3000,0,30000,3000,0,30000);
  TH2F *hxyfbc=new TH2F("hxyfbc","front-back correlation for corrected data with f-b correlation",3000,0,30000,3000,0,30000);
  struct dssdxy
  {
    int xstrip;
    int ystrip;
    Float_t energy;
  };
  dssdxy xy;
  multimap<ULong64_t, dssdxy> multimapdssd;
  twindow=100;//ns
  for(ifts=multimapfn.begin(); ifts!=multimapfn.end();ifts++) {
    ibts=multimapbn.lower_bound( (ifts->first)-twindow);
    while(ibts!=multimapbn.end()) {
        if(ibts->first >= (ifts->first) + twindow) break;
        int dt=ifts->first-ibts->first;
        hxy->Fill(ifts->second.energy,ibts->second.energy);
        if(abs(ifts->second.energy-ibts->second.energy)<500) {
            hxyfbc->Fill(ifts->second.energy,ibts->second.energy);
            xy.xstrip=ifts->second.strip;
            xy.ystrip=ibts->second.strip;
            xy.energy=ifts->second.energy;
            multimapdssd.insert(pair<ULong64_t, dssdxy>(ifts->first, xy));
        }
        ibts++;
    }   
  }
/*
//front-back correlation for corrected data
  hxy->Draw("colz");
  c1->SetLogy(0);
  c1->SetLogz();
  c1->Draw(); 
//front-back correlation for corrected data
  hxy->GetXaxis()->SetRangeUser(6000,8000);
  hxy->GetYaxis()->SetRangeUser(6000,8000);
  hxy->Draw("colz");
  c1->Draw();
//front-back correlation for corrected data with f-b correlation
  TH1F *hfe=(TH1F*) hxy->ProjectionX("hfe");
  hfe->Draw();
  c1->Draw();
//front-back correlation for corrected data with f-b correlation
  hxyfbc->GetXaxis()->SetRangeUser(6000,8000);
  hxyfbc->GetYaxis()->SetRangeUser(6000,8000);
  hxyfbc->Draw("colz");
  c1->Draw();
  TH1F *hfefbc=(TH1F*) hxyfbc->ProjectionX("hfefbc");
  hfefbc->Draw();
  c1->Draw();
  int jj=0;
  for(auto its=multimapdssd.begin(); its!=multimapdssd.end(); its++) {
    if(jj<20) cout<<its->first<<" "<<its->second.xstrip<<" "<<its->second.ystrip<<" "<<its->second.energy<<endl;
    jj++;
  }*/

  //mwpc events//
  TFile *fmw=new TFile("mwpc.root");
  TTree *tmw=(TTree*)fmw->Get("tree");
  tmw->SetBranchAddress("energy",&energy);
  tmw->SetBranchAddress("timestamp",&timestamp);
  //tmw->Print();
  //tmw->Scan("timestamp:energy","","colsize=15",10,1);

//mutimap for mwpc events
  multimap<ULong64_t, Float_t> multimapmw1;//map for mwpc
  nentries = tmw->GetEntriesFast();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    tmw->GetEntry(jentry);
    multimapmw1.insert(pair<ULong64_t,Float_t>(timestamp,energy));
    if(nentries%100000==0) cout<<jentry<<endl;
  }
  cout<<"Total number of mwpc events = "<<multimapmw1.size()<<endl;    
  fmw->Close();

  //energy spectrum of mwpc detector
  auto imw=multimapmw1.begin();//
  TH1F *hmwe=new TH1F("hmwe","energy of mwpc",3500,0,35000);
  for(int i=0; i<multimapmw1.size(); i++) {
    if(i<10) {
        sout.Form("%2d  %15llu   %5.1f" , i, imw->first, imw->second);
       cout<<sout<<endl;
    }
    hmwe->Fill(imw->second);
    imw++;
  }
  //hmwe->Draw();
  //c1->SetLogy();
  //c1->Draw();//trheshold:1000
  twindow=5000;//ns
  TH1I *hmdt=new TH1I("hmdt","mts-fts distribution(ns)",1000,-5000,5000);
  for(auto idts=multimapdssd.begin(); idts!=multimapdssd.end();idts++) {
    auto imts=multimapmw1.lower_bound( (idts->first)-twindow);
    while(imts!=multimapmw1.end()) {
        if(imts->first >= (idts->first) + twindow) break;
        int dt=imts->first-idts->first;
        hmdt->Fill(dt);
        imts++;
    }   
  }
//time distribution of mwpc-dssd
  //hmdt->Draw();
  //c1->SetLogy();
  //c1->Draw();
  //hmdt->GetXaxis()->SetRangeUser(400,1000);
  //hmdt->Draw();
  //c1->Draw();
//subtracting offset from timestamp of mwpc
  toffset=650;
  multimap<ULong64_t, Float_t> multimapmw;
  for(auto imts=multimapmw1.begin(); imts!=multimapmw1.end();imts++) 
    multimapmw.insert(pair<ULong64_t,Float_t>(imts->first - toffset, imts->second));
//dssd events coincident with mwpc
  twindow=100;//ns
  TH1F *hmwec=new TH1F("hmwec","energy of mwpc coincided with DSSD ",3500,0,35000);
  TH1F *hdsec=new TH1F("hdsec","energy of DSSD coincided with MWPC ",3500,0,35000);
  hmdt->Reset();

  opt=new TTree("tree","tree");
  TFile *opf=new TFile("sort00335.root","RECREATE");
  BranchOpt();
  for(auto idts=multimapdssd.begin(); idts!=multimapdssd.end();idts++) {
    ResetOpt();
    auto imts=multimapmw.lower_bound( (idts->first)-twindow);
    while(imts!=multimapmw.end()) {
        if(imts->first >= (idts->first) + twindow) break;
        int dt=imts->first-idts->first;
        hmwec->Fill(imts->second);
        hdsec->Fill(idts->second.energy);
	me=imts->second;
	//cout<<me<<"  "<<de<<endl;
        imts++;
    }
    dt=idts->first;
    xid=idts->second.xstrip;
    yid=idts->second.ystrip;
    de=idts->second.energy;

    opt->Fill();
    cout<<dt<<"  "<<de<<"  "<<xid<<endl;
  }
  opf->cd();
  opt->Write();
  opf->Close();
  
//energy spectrum of mwpc coincidence with dssd
 /* hmwe->Draw();
  hmwec->SetLineColor(kRed);
  hmwec->Draw("same");
  c1->SetLogy();
  c1->Draw();*/
//energy spectrum of dssd coincidence with mwpc
  /*hferaw1->Draw();//raw
  hdsec->SetLineColor(kGreen);
  hdsec->Draw("same");//coincide with MWPC
  c1->SetLogy(0);
  c1->Draw();*/
}
