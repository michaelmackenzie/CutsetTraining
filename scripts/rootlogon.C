//------------------------------------------------------------------------------
//  rootlogon.C: a sample ROOT logon macro allowing use of ROOT script 
//               compiler in CDF RunII environment. The name of this macro file
//               is defined by the .rootrc file
//------------------------------------------------------------------------------
{
                                // the line below tells rootcling where to look for 
				// the include files

  gInterpreter->AddIncludePath("./");
  gInterpreter->AddIncludePath("../");
  gInterpreter->AddIncludePath("./include");
  // gInterpreter->AddIncludePath(gSystem->Getenv("CLHEP_INC"));
  gInterpreter->AddIncludePath(Form("%s/include",gSystem->Getenv("ROOTSYS")));
  TString hostname = gSystem->Getenv("HOSTNAME");
  const char* exec_name = gApplication->Argv(0);

  if (exec_name) {
    if (strcmp(exec_name,"root.exe") == 0) {
					// print overflows/underflows in the stat box
      gStyle->SetOptStat(11111111);
					// print fit results in the stat box
      gStyle->SetOptFit(1110);
      TArrow::SetDefaultArrowSize(0.015);

      //-----------------------------------------------------------------------------
      // report the process ID which simplifies debugging
      //-----------------------------------------------------------------------------
      printf("process ID: %i\n",gSystem->GetPid());
      TAuthenticate::SetGlobalUser(gSystem->Getenv("USER"));
      gInterpreter->ProcessLine(".! ps | grep root");
      printf("Loading CutsetTrainer\n");
      // TString cmssw = gSystem->Getenv("CMSSW_BASE");
      TString path = gSystem->Getenv("PWD");
      path += "/../";
      gSystem->Load((path+"CutsetTrainer_cc.so").Data());
    }
  }
}


