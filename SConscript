#!/usr/bin/env python
#
# Script to build the files found in this directory.
#------------------------------------------------------------------------------
import os
Import('env')
Import('mmanalysis_helper')
#------------------------------------------------------------------------------
# 2016-10-10: STNTUPLE link: add ROOT 'EG' library after 'Physics' library
#------------------------------------------------------------------------------
rootlibs = env['ROOTLIBS']
if ( not ("EG"   in rootlibs)): rootlibs.insert(rootlibs.index("Physics")+1,"EG");

helper = mmanalysis_helper();

helper.handle_dictionaries();

skip_list = []

stntuple_libs = [ 'Stntuple_val', 'Stntuple_alg', 'Stntuple_loop',
                  'Stntuple_obj', 'Stntuple_geom', 'Stntuple_base'
                  ];

mc_libs = []

if (os.getenv("STNTUPLE_MC_GEN") != None): mc_libs = ["mc_photos", "mc_base"];

libs          = stntuple_libs + mc_libs + [
                'mu2e_Mu2eUtilities',
                # 'mu2e_ConditionsService_ConditionsService_service',
                # 'mu2e_GeometryService_GeometryService_service',
                # 'mu2e_GlobalConstantsService',
                # 'mu2e_StoppingTargetGeom',
                # 'mu2e_GeneralUtilities',
                # 'mu2e_MCDataProducts',
                # 'mu2e_RecoDataProducts',
                # 'mu2e_DataProducts',
                # 'mu2e_ConfigTools',
                # 'mu2e_Mu2eInterfaces',
                # 'mu2e_GeneralUtilities',
                # 'BTrk_BbrGeom',
                # 'gsl',
                # 'art_Persistency_Common',
                # 'art_Persistency_Provenance',
                # 'art_Framework_Services_Optional_RandomNumberGenerator_service',
                # 'art_Framework_Services_Optional_TFileService_service',
                # 'art_Framework_Services_Registry',
                # 'art_Framework_Services_Optional',
                # 'art_Framework_Principal',
                # 'art_Framework_Core',
                # 'art_Utilities',
                # 'canvas',
                # 'MF_MessageLogger',
                # 'CLHEP',
                # 'HepPDT',
                'xerces-c',
                'fhiclcpp',
                'cetlib',
                'cetlib_except',
                rootlibs,
                'boost_filesystem',
                'boost_system'
                ];

list_of_cc_files = Glob('*.cc',strings=True);
helper.build_libs(list_of_cc_files,skip_list,libs);
# print "tmpdir:"+env['TMP_LIB_DIR']
