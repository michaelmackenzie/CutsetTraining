//------------------------------------------------------------------------------
//  rootlogon.C: a sample ROOT logon macro allowing use of ROOT script
//               compiler in CDF RunII environment. The name of this macro file
//               is defined by the .rootrc file
//
// assume that the environment variable MU2E_TEST_RELEASE points to the
//  Jul 08 2014 P.Murat
//------------------------------------------------------------------------------
{
                                // the line below tells rootcling where to look for
                                // the include files

  gInterpreter->AddIncludePath("./include");
  gInterpreter->AddIncludePath(gSystem->Getenv("CLHEP_INC"));
  gInterpreter->AddIncludePath(Form("%s/include",gSystem->Getenv("ROOTSYS")));


  gInterpreter->AddIncludePath(Form("%s/tex/cdfnotes",
                                    gSystem->Getenv("HOME")));

  //  gSystem->SetMakeSharedLib("cd $BuildDir ; g++ -c -g $Opt -pipe -m32 -Wall -W -Woverloaded-virtual -fPIC -pthread $IncludePath $SourceFiles ;  g++ -g $ObjectFiles -shared -Wl,-soname,$LibName.so -m32 $LinkedLibs -o $SharedLib");
//-----------------------------------------------------------------------------
// load in ROOT physics vectors and event generator libraries
//-----------------------------------------------------------------------------
  gSystem->Load("$ROOTSYS/lib/libEG.so");
  //  gSystem->Load("$ROOTSYS/lib/libPhysics.so");
  gSystem->Load("$ROOTSYS/lib/libMinuit.so");
  gSystem->Load("$ROOTSYS/lib/libFumili.so");
  //  gSystem->Load("$ROOTSYS/lib/libTree.so");
  //  gSystem->Load("$ROOTSYS/lib/libRuby.so");
//-----------------------------------------------------------------------------
//  check batch mode
//-----------------------------------------------------------------------------
  const char* opt ;
  int batch_mode = 0;

  int nargs = gApplication->Argc();

  for (int i=1; i<nargs; i++) {
    opt  = gApplication->Argv(i);
    if (strcmp(opt,"-b") == 0) {
      batch_mode = 1;
      break;
    }
  }

  printf("   batch_mode = %i\n",batch_mode);
//-----------------------------------------------------------------------------
// always need libStntuple_loop, but the other 2 libs should be loaded in
// only if we're running bare root
//-----------------------------------------------------------------------------
  const char* exec_name = gApplication->Argv(0);

  printf(" nargs = %2i exec_name = %s\n",nargs, exec_name);

  if (exec_name) {
    if (strstr(exec_name,"root.exe") != 0) {
//-----------------------------------------------------------------------------
// assume STNTUPLE  analysis job
//-----------------------------------------------------------------------------
      if (batch_mode == 1) gSystem->Load("$ROOTSYS/lib/libGui.so");
//-----------------------------------------------------------------------------
// Mu2e Offline libraries
//-----------------------------------------------------------------------------
//     //     gSystem->Load("$MU2E_BASE_RELEASE/lib/libmu2e_Mu2eInterfaces.so");
//     //     gSystem->Load("$MU2E_BASE_RELEASE/lib/libmu2e_CalorimeterGeom.so");
//
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_base.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_obj.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_loop.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_alg.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_ana.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_val.so");

      // gSystem->Load("$MU2E_BASE_RELEASE/lib/libmurat_obj.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libsu2020_ana.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libsu2020_cutset_training.so");

                                        // print overflows/underflows in the stat box
      gStyle->SetOptStat(11111111);
                                        // print fit results in the stat box
      gStyle->SetOptFit(1110);
      TArrow::SetDefaultArrowSize(0.015);
    }
    else if (strstr(exec_name,"mu2e.NNN") != 0) {
      //      gSystem->Load("$MU2E_BASE_RELEASE/lib/libmurat_obj.so");
      gSystem->Load("$MU2E_BASE_RELEASE/lib/libStntuple_val.so");
      //      gSystem->Load("$MU2E_BASE_RELEASE/lib/libmurat_plot.so");
    }
//-----------------------------------------------------------------------------
//  databases
//-----------------------------------------------------------------------------
//   gSystem->Load("libStntuple_oracle.so");

    if (gSystem->Exec("ls $HOME/root_macros/set_style.C &> /dev/null") == 0) {
      gInterpreter->ExecuteMacro("$HOME/root_macros/set_style.C");
    }
  }
//-----------------------------------------------------------------------------
// report the process ID which simplifies debugging
//-----------------------------------------------------------------------------
  printf(" process ID: %i\n",gSystem->GetPid());
  TAuthenticate::SetGlobalUser(gSystem->Getenv("USER"));
  gInterpreter->ProcessLine(".! ps | grep root");
}
