// **************************************************
//
// run with: root -l -b -q computeTriggerEffect.C++
//
// *************************************************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"



void readFiles(){

  
  const int nDatasets = 5;
  TString samplesName[nDatasets] = {"ggH125", "VBFH125", "WminusH125", "WplusH125", "ZH125"}; //, "ttH125", "bbH125"};

  // --- define input file characteristics
  TString rootInDir = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/Run2UL_22/Hadded/MC/MC2016preVFP/";
  TString rootDirName  = "ZZTree";
  TString rootTreeName = "candTree";

  TFile* inputFile[nDatasets];
  TTree* inputTree[nDatasets];

  // --- define input file branches
  Float_t ZZMass;
  vector<Short_t> *LepLepId = 0;
  vector<Float_t> *LepEta = 0;

  // --- define histos
  TH1F* hEvents_eta2p5[nDatasets];
  TH1F* hEvents_eta2p4_single[nDatasets];
  TH1F* hEvents_eta2p3_single[nDatasets];
  TH1F* hEvents_eta2p4_double[nDatasets];
  TH1F* hEvents_eta2p3_double[nDatasets];
  for(int d=0; d<nDatasets; d++){
    hEvents_eta2p5[d]        = new TH1F("hEvents_eta2p5_"+samplesName[d],"",1,0.,1.);
    hEvents_eta2p4_single[d] = new TH1F("hEvents_eta2p4_single_"+samplesName[d],"",1,0.,1.);
    hEvents_eta2p3_single[d] = new TH1F("hEvents_eta2p3_single_"+samplesName[d],"",1,0.,1.);
    hEvents_eta2p4_double[d] = new TH1F("hEvents_eta2p4_double_"+samplesName[d],"",1,0.,1.);
    hEvents_eta2p3_double[d] = new TH1F("hEvents_eta2p3_double_"+samplesName[d],"",1,0.,1.);
  }

  //--- loop over all datasets
  for(int d=0; d<nDatasets; d++){

    // --- open input file
    TString inputFileName;
    inputFileName = rootInDir + samplesName[d] + "/ZZ4lAnalysis.root";
    cout<<"Opening input file "<<inputFileName<<" ..."<<endl;
    inputFile[d] = TFile::Open(inputFileName);
    inputTree[d] = (TTree*)inputFile[d]->Get(rootDirName + "/" + rootTreeName);

    // --- set branch addresses
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("LepLepId", &LepLepId);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);

    // --- loop over input tree entries
    Long64_t entries = inputTree[d]->GetEntries();    
    cout<<"Processing file: "<< samplesName[d] << " (" << entries <<" entries) ..."<< endl;    
    for(Long64_t z=0; z<entries; z++){
      
      inputTree[d]->GetEntry(z);

      // -2.5 < eta < 2.5
      hEvents_eta2p5[d]->Fill(0.5, 1.);     

      // -2.4 < eta < 2.4
      if((TMath::Abs(LepEta->at(0))<2.4) && TMath::Abs(LepLepId->at(0))==11){
          hEvents_eta2p4_single[d]->Fill(0.5, 1.);     
      }

      if((TMath::Abs(LepEta->at(0))<2.4 && TMath::Abs(LepLepId->at(0))==11) && (TMath::Abs(LepEta->at(1))<2.4 && TMath::Abs(LepLepId->at(1))==11)){
          hEvents_eta2p4_double[d]->Fill(0.5, 1.);     
      }
      
      // -2.3 < eta < 2.3
      if((TMath::Abs(LepEta->at(0))<2.3) && TMath::Abs(LepLepId->at(0))==11){
          hEvents_eta2p3_single[d]->Fill(0.5, 1.);     
      }

      if((TMath::Abs(LepEta->at(0))<2.3 && TMath::Abs(LepLepId->at(0))==11) && (TMath::Abs(LepEta->at(1))<2.3 && TMath::Abs(LepLepId->at(1))==11)){
          hEvents_eta2p3_double[d]->Fill(0.5, 1.);     
      }


    } // end loop over entries


  } // end loop over datasets

  
  // --- save histos
  TFile* fout = new TFile("histos.root", "recreate");
  fout->cd();
  for(int d=0; d<nDatasets; d++){
    hEvents_eta2p5       [d]->Write();
    hEvents_eta2p4_single[d]->Write();
    hEvents_eta2p3_single[d]->Write();
    hEvents_eta2p4_double[d]->Write();
    hEvents_eta2p3_double[d]->Write();
  }
  fout->Close();
  cout<<"Histos saved"<<endl;


}




void computeTriggerEffect(){

  readFiles();

}
