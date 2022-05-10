.L SegBEGeK2ESA.cxx+
TChain ch("data");
ch.Add("root_ntuples/h.regeb.20220510T093246.root");
ch.Process("Create_Energy_Spectrum.C+","20220510T093246");
