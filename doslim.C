void doslim(int i) {
  gROOT->LoadMacro("slimmer_2017.C+");
  gROOT->ProcessLine(Form("slimmer_2017(\"/data/users/abelloni/CERN_TB_Jul17/NTUP/ana_h2_tb_run00%i_EMAP-27JUL2017_Phase2-CRF_RM3.root\")",i));
}
