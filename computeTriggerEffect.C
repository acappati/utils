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


#define REDOHISTOS 1



const int nDatasets = 5;
TString samplesName[nDatasets] = {"ggH125", "VBFH125", "WminusH125", "WplusH125", "ZH125"}; //, "ttH125", "bbH125"};


void makeHistos(){

  // --- define input file characteristics
  TString rootInDir = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/Run2UL_22/Hadded/MC/MC2016preVFP/";
  TString rootDirName  = "ZZTree";
  TString rootTreeName = "candTree";

  TFile* inputFile[nDatasets];
  TTree* inputTree[nDatasets];
  TH1F* hCounters[nDatasets];
  Float_t gen_sumWeights[nDatasets];
  Float_t partialSampleWeight[nDatasets];

  // --- define input file branches
  Float_t lumi = 36.3;
  Float_t overallEventWeight;
  Float_t xsec;
  Float_t ZZMass;
  vector<Short_t> *LepLepId = 0;
  vector<Float_t> *LepEta = 0;

  // --- define histos
  TH1F* hEvents_eta2p5[nDatasets];
  TH1F* hEvents_eta2p4_single[nDatasets];
  TH1F* hEvents_eta2p3_single[nDatasets];
  TH1F* hEvents_eta2p4_double[nDatasets];
  TH1F* hEvents_eta2p3_double[nDatasets];
  TH1F* hYields_eta2p5[nDatasets];
  TH1F* hYields_eta2p4_single[nDatasets];
  TH1F* hYields_eta2p3_single[nDatasets];
  TH1F* hYields_eta2p4_double[nDatasets];
  TH1F* hYields_eta2p3_double[nDatasets];
  for(int d=0; d<nDatasets; d++){
    // events
    hEvents_eta2p5[d]        = new TH1F("hEvents_eta2p5_"+samplesName[d],"",1,0.,1.);
    hEvents_eta2p4_single[d] = new TH1F("hEvents_eta2p4_single_"+samplesName[d],"",1,0.,1.);
    hEvents_eta2p3_single[d] = new TH1F("hEvents_eta2p3_single_"+samplesName[d],"",1,0.,1.);
    hEvents_eta2p4_double[d] = new TH1F("hEvents_eta2p4_double_"+samplesName[d],"",1,0.,1.);
    hEvents_eta2p3_double[d] = new TH1F("hEvents_eta2p3_double_"+samplesName[d],"",1,0.,1.);
    // yields
    hYields_eta2p5[d]        = new TH1F("hYields_eta2p5_"+samplesName[d],"",1,0.,1.);
    hYields_eta2p4_single[d] = new TH1F("hYields_eta2p4_single_"+samplesName[d],"",1,0.,1.);
    hYields_eta2p3_single[d] = new TH1F("hYields_eta2p3_single_"+samplesName[d],"",1,0.,1.);
    hYields_eta2p4_double[d] = new TH1F("hYields_eta2p4_double_"+samplesName[d],"",1,0.,1.);
    hYields_eta2p3_double[d] = new TH1F("hYields_eta2p3_double_"+samplesName[d],"",1,0.,1.);
  }

  //--- loop over all datasets
  for(int d=0; d<nDatasets; d++){

    // --- open input file
    TString inputFileName;
    inputFileName = rootInDir + samplesName[d] + "/ZZ4lAnalysis.root";
    cout<<"Opening input file "<<inputFileName<<" ..."<<endl;
    inputFile[d] = TFile::Open(inputFileName);

    hCounters[d] = (TH1F*)inputFile[d]->Get("ZZTree/Counters");
    gen_sumWeights[d] = (Long64_t)hCounters[d]->GetBinContent(40);
    partialSampleWeight[d] = lumi * 1000 / gen_sumWeights[d];
    inputTree[d] = (TTree*)inputFile[d]->Get(rootDirName + "/" + rootTreeName);

    // --- set branch addresses
    inputTree[d]->SetBranchAddress("overallEventWeight", &overallEventWeight);
    inputTree[d]->SetBranchAddress("xsec", &xsec);
    inputTree[d]->SetBranchAddress("ZZMass", &ZZMass);
    inputTree[d]->SetBranchAddress("LepLepId", &LepLepId);
    inputTree[d]->SetBranchAddress("LepEta", &LepEta);

    // --- loop over input tree entries
    Long64_t entries = inputTree[d]->GetEntries();    
    cout<<"Processing file: "<< samplesName[d] << " (" << entries <<" entries) ..."<< endl;    
    for(Long64_t z=0; z<entries; z++){
      
      inputTree[d]->GetEntry(z);

      Double_t eventWeight = 1.; 
      eventWeight = partialSampleWeight[d] *xsec *overallEventWeight;

      // -2.5 < eta < 2.5
      hEvents_eta2p5[d]->Fill(0.5, 1.);     
      hYields_eta2p5[d]->Fill(0.5, eventWeight);     

      // -2.4 < eta < 2.4
      if((TMath::Abs(LepEta->at(0))<2.4 && TMath::Abs(LepLepId->at(0))==11) || TMath::Abs(LepLepId->at(0))==13){
          hEvents_eta2p4_single[d]->Fill(0.5, 1.);     
          hYields_eta2p4_single[d]->Fill(0.5, eventWeight);     
      }

      if(((TMath::Abs(LepEta->at(0))<2.4 && TMath::Abs(LepLepId->at(0))==11) || TMath::Abs(LepLepId->at(0))==13) && ((TMath::Abs(LepEta->at(1))<2.4 && TMath::Abs(LepLepId->at(1))==11) || TMath::Abs(LepLepId->at(1))==13)){
          hEvents_eta2p4_double[d]->Fill(0.5, 1.);     
          hYields_eta2p4_double[d]->Fill(0.5, eventWeight);     
      }
      
      // -2.3 < eta < 2.3
      if((TMath::Abs(LepEta->at(0))<2.3 && TMath::Abs(LepLepId->at(0))==11) || TMath::Abs(LepLepId->at(0))==13){
          hEvents_eta2p3_single[d]->Fill(0.5, 1.);     
          hYields_eta2p3_single[d]->Fill(0.5, eventWeight);     
      }

      if(((TMath::Abs(LepEta->at(0))<2.3 && TMath::Abs(LepLepId->at(0))==11) || TMath::Abs(LepLepId->at(0))==13) && ((TMath::Abs(LepEta->at(1))<2.3 && TMath::Abs(LepLepId->at(1))==11) || TMath::Abs(LepLepId->at(1))==13)){
          hEvents_eta2p3_double[d]->Fill(0.5, 1.);     
          hYields_eta2p3_double[d]->Fill(0.5, eventWeight);     
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
    hYields_eta2p5       [d]->Write();
    hYields_eta2p4_single[d]->Write();
    hYields_eta2p3_single[d]->Write();
    hYields_eta2p4_double[d]->Write();
    hYields_eta2p3_double[d]->Write();
  }
  fout->Close();
  cout<<"Histos saved"<<endl;


}



void printYields(){

  TFile* finput = TFile::Open("histos.root");
  
  TH1F* hEvents_eta2p5       [nDatasets];
  TH1F* hEvents_eta2p4_single[nDatasets];
  TH1F* hEvents_eta2p3_single[nDatasets];
  TH1F* hEvents_eta2p4_double[nDatasets];
  TH1F* hEvents_eta2p3_double[nDatasets];
  TH1F* hYields_eta2p5       [nDatasets];
  TH1F* hYields_eta2p4_single[nDatasets];
  TH1F* hYields_eta2p3_single[nDatasets];
  TH1F* hYields_eta2p4_double[nDatasets];
  TH1F* hYields_eta2p3_double[nDatasets];

  cout<<"Retrieving histos.."<<endl;
  for(int d=0; d<nDatasets; d++){
    hEvents_eta2p5[d]       = (TH1F*)finput->Get("hEvents_eta2p5_"+samplesName[d]);
    hEvents_eta2p4_single[d]= (TH1F*)finput->Get("hEvents_eta2p4_single_"+samplesName[d]);
    hEvents_eta2p3_single[d]= (TH1F*)finput->Get("hEvents_eta2p3_single_"+samplesName[d]);
    hEvents_eta2p4_double[d]= (TH1F*)finput->Get("hEvents_eta2p4_double_"+samplesName[d]);
    hEvents_eta2p3_double[d]= (TH1F*)finput->Get("hEvents_eta2p3_double_"+samplesName[d]);
    hYields_eta2p5[d]       = (TH1F*)finput->Get("hYields_eta2p5_"+samplesName[d]);
    hYields_eta2p4_single[d]= (TH1F*)finput->Get("hYields_eta2p4_single_"+samplesName[d]);
    hYields_eta2p3_single[d]= (TH1F*)finput->Get("hYields_eta2p3_single_"+samplesName[d]);
    hYields_eta2p4_double[d]= (TH1F*)finput->Get("hYields_eta2p4_double_"+samplesName[d]);
    hYields_eta2p3_double[d]= (TH1F*)finput->Get("hYields_eta2p3_double_"+samplesName[d]);
  }

  
  // --- Events
  Float_t evt_now_ggH = hEvents_eta2p5[0]->GetBinContent(1);
  Float_t evt_now_VBF = hEvents_eta2p5[1]->GetBinContent(1);
  Float_t evt_now_VH  = hEvents_eta2p5[2]->GetBinContent(1) + hEvents_eta2p5[3]->GetBinContent(1) + hEvents_eta2p5[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t evt_now_tot = evt_now_ggH + evt_now_VBF + evt_now_VH;

  Float_t evt_singEle2p4_ggH = hEvents_eta2p4_single[0]->GetBinContent(1);
  Float_t evt_singEle2p4_VBF = hEvents_eta2p4_single[1]->GetBinContent(1);
  Float_t evt_singEle2p4_VH  = hEvents_eta2p4_single[2]->GetBinContent(1) + hEvents_eta2p4_single[3]->GetBinContent(1) + hEvents_eta2p4_single[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t evt_singEle2p4_tot = evt_singEle2p4_ggH + evt_singEle2p4_VBF + evt_singEle2p4_VH;

  Float_t evt_doubEle2p4_ggH = hEvents_eta2p4_double[0]->GetBinContent(1);
  Float_t evt_doubEle2p4_VBF = hEvents_eta2p4_double[1]->GetBinContent(1);
  Float_t evt_doubEle2p4_VH  = hEvents_eta2p4_double[2]->GetBinContent(1) + hEvents_eta2p4_double[3]->GetBinContent(1) + hEvents_eta2p4_double[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t evt_doubEle2p4_tot = evt_doubEle2p4_ggH + evt_doubEle2p4_VBF + evt_doubEle2p4_VH;

  Float_t evt_singEle2p3_ggH = hEvents_eta2p3_single[0]->GetBinContent(1);
  Float_t evt_singEle2p3_VBF = hEvents_eta2p3_single[1]->GetBinContent(1);
  Float_t evt_singEle2p3_VH  = hEvents_eta2p3_single[2]->GetBinContent(1) + hEvents_eta2p3_single[3]->GetBinContent(1) + hEvents_eta2p3_single[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t evt_singEle2p3_tot = evt_singEle2p3_ggH + evt_singEle2p3_VBF + evt_singEle2p3_VH;

  Float_t evt_doubEle2p3_ggH = hEvents_eta2p3_double[0]->GetBinContent(1);
  Float_t evt_doubEle2p3_VBF = hEvents_eta2p3_double[1]->GetBinContent(1);
  Float_t evt_doubEle2p3_VH  = hEvents_eta2p3_double[2]->GetBinContent(1) + hEvents_eta2p3_double[3]->GetBinContent(1) + hEvents_eta2p3_double[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t evt_doubEle2p3_tot = evt_doubEle2p3_ggH + evt_doubEle2p3_VBF + evt_doubEle2p3_VH;


  // print stuff
  cout<<"-------------------------------------"<<endl;
  cout<<"    "<<"Now      "<<"  "<<"singleEle |eta|<2.4  "<<"  "<<"doubleEle |eta|<2.4  "<<"  "<<"singleEle |eta|<2.3  "<<"  "<<"doubleEle |eta|<2.3  "<<endl;
  cout<<"ggH "<<evt_now_ggH<<"  "<<evt_singEle2p4_ggH     <<"  "<<evt_doubEle2p4_ggH     <<"  "<<evt_singEle2p3_ggH     <<"  "<<evt_doubEle2p3_ggH     <<endl;
  cout<<"VBF "<<evt_now_VBF<<"  "<<evt_singEle2p4_VBF     <<"  "<<evt_doubEle2p4_VBF     <<"  "<<evt_singEle2p3_VBF     <<"  "<<evt_doubEle2p3_VBF     <<endl;
  cout<<"VH  "<<evt_now_VH <<"  "<<evt_singEle2p4_VH      <<"  "<<evt_doubEle2p4_VH      <<"  "<<evt_singEle2p3_VH      <<"  "<<evt_doubEle2p3_VH      <<endl;
  cout<<"Tot "<<evt_now_tot<<"  "<<evt_singEle2p4_tot     <<"  "<<evt_doubEle2p4_tot     <<"  "<<evt_singEle2p3_tot     <<"  "<<evt_doubEle2p3_tot     <<endl;


  // compute percentage of lost events
  Float_t ploss_singEle2p4_ggH = ((evt_now_ggH - evt_singEle2p4_ggH) / evt_now_ggH) * 100.;
  Float_t ploss_singEle2p4_VBF = ((evt_now_VBF - evt_singEle2p4_VBF) / evt_now_VBF) * 100.;
  Float_t ploss_singEle2p4_VH  = ((evt_now_VH  - evt_singEle2p4_VH)  / evt_now_VH ) * 100.;
  Float_t ploss_singEle2p4_tot = ((evt_now_tot - evt_singEle2p4_tot) / evt_now_tot) * 100.;

  Float_t ploss_doubEle2p4_ggH = ((evt_now_ggH - evt_doubEle2p4_ggH) / evt_now_ggH) * 100.;
  Float_t ploss_doubEle2p4_VBF = ((evt_now_VBF - evt_doubEle2p4_VBF) / evt_now_VBF) * 100.;
  Float_t ploss_doubEle2p4_VH  = ((evt_now_VH  - evt_doubEle2p4_VH)  / evt_now_VH ) * 100.;
  Float_t ploss_doubEle2p4_tot = ((evt_now_tot - evt_doubEle2p4_tot) / evt_now_tot) * 100.;

  Float_t ploss_singEle2p3_ggH = ((evt_now_ggH - evt_singEle2p3_ggH) / evt_now_ggH) * 100.;
  Float_t ploss_singEle2p3_VBF = ((evt_now_VBF - evt_singEle2p3_VBF) / evt_now_VBF) * 100.;
  Float_t ploss_singEle2p3_VH  = ((evt_now_VH  - evt_singEle2p3_VH)  / evt_now_VH ) * 100.;
  Float_t ploss_singEle2p3_tot = ((evt_now_tot - evt_singEle2p3_tot) / evt_now_tot) * 100.;

  Float_t ploss_doubEle2p3_ggH = ((evt_now_ggH - evt_doubEle2p3_ggH) / evt_now_ggH) * 100.;
  Float_t ploss_doubEle2p3_VBF = ((evt_now_VBF - evt_doubEle2p3_VBF) / evt_now_VBF) * 100.;
  Float_t ploss_doubEle2p3_VH  = ((evt_now_VH  - evt_doubEle2p3_VH)  / evt_now_VH ) * 100.;
  Float_t ploss_doubEle2p3_tot = ((evt_now_tot - evt_doubEle2p3_tot) / evt_now_tot) * 100.;

  // print stuff
  cout<<"    "<<"singleEle|eta|<2.4  "               <<"doubleEle|eta|<2.4  "             <<"singleEle|eta|<2.3  "             <<"doubleEle|eta|<2.3  "          <<endl;
  cout<<"ggH "<<ploss_singEle2p4_ggH <<" %          "<<ploss_doubEle2p4_ggH <<" %        "<<ploss_singEle2p3_ggH <<" %        "<<ploss_doubEle2p3_ggH    <<" %  "<<endl;
  cout<<"VBF "<<ploss_singEle2p4_VBF <<" %          "<<ploss_doubEle2p4_VBF <<" %        "<<ploss_singEle2p3_VBF <<" %        "<<ploss_doubEle2p3_VBF	 <<" %  "<<endl;
  cout<<"VH  "<<ploss_singEle2p4_VH  <<" %          "<<ploss_doubEle2p4_VH  <<" %        "<<ploss_singEle2p3_VH  <<" %        "<<ploss_doubEle2p3_VH 	 <<" %  "<<endl;
  cout<<"Tot "<<ploss_singEle2p4_tot <<" %          "<<ploss_doubEle2p4_tot <<" %        "<<ploss_singEle2p3_tot <<" %        "<<ploss_doubEle2p3_tot	 <<" %  "<<endl;




  // --- Yields
  Float_t yie_now_ggH = hYields_eta2p5[0]->GetBinContent(1);
  Float_t yie_now_VBF = hYields_eta2p5[1]->GetBinContent(1);
  Float_t yie_now_VH  = hYields_eta2p5[2]->GetBinContent(1) + hYields_eta2p5[3]->GetBinContent(1) + hYields_eta2p5[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t yie_now_tot = yie_now_ggH + yie_now_VBF + yie_now_VH;

  Float_t yie_singEle2p4_ggH = hYields_eta2p4_single[0]->GetBinContent(1);
  Float_t yie_singEle2p4_VBF = hYields_eta2p4_single[1]->GetBinContent(1);
  Float_t yie_singEle2p4_VH  = hYields_eta2p4_single[2]->GetBinContent(1) + hYields_eta2p4_single[3]->GetBinContent(1) + hYields_eta2p4_single[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t yie_singEle2p4_tot = yie_singEle2p4_ggH + yie_singEle2p4_VBF + yie_singEle2p4_VH;

  Float_t yie_doubEle2p4_ggH = hYields_eta2p4_double[0]->GetBinContent(1);
  Float_t yie_doubEle2p4_VBF = hYields_eta2p4_double[1]->GetBinContent(1);
  Float_t yie_doubEle2p4_VH  = hYields_eta2p4_double[2]->GetBinContent(1) + hYields_eta2p4_double[3]->GetBinContent(1) + hYields_eta2p4_double[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t yie_doubEle2p4_tot = yie_doubEle2p4_ggH + yie_doubEle2p4_VBF + yie_doubEle2p4_VH;

  Float_t yie_singEle2p3_ggH = hYields_eta2p3_single[0]->GetBinContent(1);
  Float_t yie_singEle2p3_VBF = hYields_eta2p3_single[1]->GetBinContent(1);
  Float_t yie_singEle2p3_VH  = hYields_eta2p3_single[2]->GetBinContent(1) + hYields_eta2p3_single[3]->GetBinContent(1) + hYields_eta2p3_single[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t yie_singEle2p3_tot = yie_singEle2p3_ggH + yie_singEle2p3_VBF + yie_singEle2p3_VH;

  Float_t yie_doubEle2p3_ggH = hYields_eta2p3_double[0]->GetBinContent(1);
  Float_t yie_doubEle2p3_VBF = hYields_eta2p3_double[1]->GetBinContent(1);
  Float_t yie_doubEle2p3_VH  = hYields_eta2p3_double[2]->GetBinContent(1) + hYields_eta2p3_double[3]->GetBinContent(1) + hYields_eta2p3_double[4]->GetBinContent(1); // WplusH + WminusH + ZH
  Float_t yie_doubEle2p3_tot = yie_doubEle2p3_ggH + yie_doubEle2p3_VBF + yie_doubEle2p3_VH;


  // print stuff
  cout<<"-------------------------------------"<<endl;
  cout<<"    "<<"Now      "<<"  "<<"singleEle |eta|<2.4  "<<"  "<<"doubleEle |eta|<2.4  "<<"  "<<"singleEle |eta|<2.3  "<<"  "<<"doubleEle |eta|<2.3  "<<endl;
  cout<<"ggH "<<yie_now_ggH<<"  "<<yie_singEle2p4_ggH     <<"  "<<yie_doubEle2p4_ggH     <<"  "<<yie_singEle2p3_ggH     <<"  "<<yie_doubEle2p3_ggH     <<endl;
  cout<<"VBF "<<yie_now_VBF<<"  "<<yie_singEle2p4_VBF     <<"  "<<yie_doubEle2p4_VBF     <<"  "<<yie_singEle2p3_VBF     <<"  "<<yie_doubEle2p3_VBF     <<endl;
  cout<<"VH  "<<yie_now_VH <<"  "<<yie_singEle2p4_VH      <<"  "<<yie_doubEle2p4_VH      <<"  "<<yie_singEle2p3_VH      <<"  "<<yie_doubEle2p3_VH      <<endl;
  cout<<"Tot "<<yie_now_tot<<"  "<<yie_singEle2p4_tot     <<"  "<<yie_doubEle2p4_tot     <<"  "<<yie_singEle2p3_tot     <<"  "<<yie_doubEle2p3_tot     <<endl;


  // compute percentage of lost events
  Float_t plossy_singEle2p4_ggH = ((yie_now_ggH - yie_singEle2p4_ggH) / yie_now_ggH) * 100.;
  Float_t plossy_singEle2p4_VBF = ((yie_now_VBF - yie_singEle2p4_VBF) / yie_now_VBF) * 100.;
  Float_t plossy_singEle2p4_VH  = ((yie_now_VH  - yie_singEle2p4_VH)  / yie_now_VH ) * 100.;
  Float_t plossy_singEle2p4_tot = ((yie_now_tot - yie_singEle2p4_tot) / yie_now_tot) * 100.;

  Float_t plossy_doubEle2p4_ggH = ((yie_now_ggH - yie_doubEle2p4_ggH) / yie_now_ggH) * 100.;
  Float_t plossy_doubEle2p4_VBF = ((yie_now_VBF - yie_doubEle2p4_VBF) / yie_now_VBF) * 100.;
  Float_t plossy_doubEle2p4_VH  = ((yie_now_VH  - yie_doubEle2p4_VH)  / yie_now_VH ) * 100.;
  Float_t plossy_doubEle2p4_tot = ((yie_now_tot - yie_doubEle2p4_tot) / yie_now_tot) * 100.;

  Float_t plossy_singEle2p3_ggH = ((yie_now_ggH - yie_singEle2p3_ggH) / yie_now_ggH) * 100.;
  Float_t plossy_singEle2p3_VBF = ((yie_now_VBF - yie_singEle2p3_VBF) / yie_now_VBF) * 100.;
  Float_t plossy_singEle2p3_VH  = ((yie_now_VH  - yie_singEle2p3_VH)  / yie_now_VH ) * 100.;
  Float_t plossy_singEle2p3_tot = ((yie_now_tot - yie_singEle2p3_tot) / yie_now_tot) * 100.;

  Float_t plossy_doubEle2p3_ggH = ((yie_now_ggH - yie_doubEle2p3_ggH) / yie_now_ggH) * 100.;
  Float_t plossy_doubEle2p3_VBF = ((yie_now_VBF - yie_doubEle2p3_VBF) / yie_now_VBF) * 100.;
  Float_t plossy_doubEle2p3_VH  = ((yie_now_VH  - yie_doubEle2p3_VH)  / yie_now_VH ) * 100.;
  Float_t plossy_doubEle2p3_tot = ((yie_now_tot - yie_doubEle2p3_tot) / yie_now_tot) * 100.;

  // print stuff
  cout<<"    "<<"singleEle|eta|<2.4  "               <<"doubleEle|eta|<2.4  "             <<"singleEle|eta|<2.3  "             <<"doubleEle|eta|<2.3  "          <<endl;
  cout<<"ggH "<<plossy_singEle2p4_ggH <<" %          "<<plossy_doubEle2p4_ggH <<" %        "<<plossy_singEle2p3_ggH <<" %        "<<plossy_doubEle2p3_ggH    <<" %  "<<endl;
  cout<<"VBF "<<plossy_singEle2p4_VBF <<" %          "<<plossy_doubEle2p4_VBF <<" %        "<<plossy_singEle2p3_VBF <<" %        "<<plossy_doubEle2p3_VBF	 <<" %  "<<endl;
  cout<<"VH  "<<plossy_singEle2p4_VH  <<" %          "<<plossy_doubEle2p4_VH  <<" %        "<<plossy_singEle2p3_VH  <<" %        "<<plossy_doubEle2p3_VH 	 <<" %  "<<endl;
  cout<<"Tot "<<plossy_singEle2p4_tot <<" %          "<<plossy_doubEle2p4_tot <<" %        "<<plossy_singEle2p3_tot <<" %        "<<plossy_doubEle2p3_tot	 <<" %  "<<endl;

}




void computeTriggerEffect(){

  if(REDOHISTOS){
    makeHistos();
  }

  printYields();
  

}
