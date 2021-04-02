#!/usr/bin/env python
#
# Script to build the files found in this directory.
#------------------------------------------------------------------------------
import os
Import('env')
Import('stntuple_helper')
#------------------------------------------------------------------------------
# 2016-10-10: STNTUPLE link: add ROOT 'EG' library after 'Physics' library
#------------------------------------------------------------------------------
def local_build():

    local_env = env.Clone()
    rootlibs  = local_env['ROOTLIBS']

    if ( not ("EG"   in rootlibs)): rootlibs.insert(rootlibs.index("Physics")+1,"EG");
    if ( not ("Eve"  in rootlibs)): rootlibs.insert(rootlibs.index("EG"     )+1,"Eve");
    if ( not ("Geom" in rootlibs)): rootlibs.insert(rootlibs.index("Eve"    )+1,"Geom");
    if ( not ("RGL"  in rootlibs)): rootlibs.insert(rootlibs.index("Geom"   )+1,"RGL");
    if ( not ("TMVA" in rootlibs)): rootlibs.append("TMVA");

    debug  = False
    helper = stntuple_helper(env,debug);

    dict_skip_list = []

    helper.handle_dictionaries(".hh",dict_skip_list);

    list_of_cc_files =  Glob('*.cc', strings=True);
    skip_list        = [ ]

    stntuple_libs    = ['Stntuple_alg', 'Stntuple_obj', 'Stntuple_loop',
                        'Stntuple_gui', 'Stntuple_geom', 'Stntuple_base'];

    libs             = stntuple_libs + [ rootlibs, ];

    helper.build_libs(list_of_cc_files, skip_list,libs);

#------------------------------------------------------------------------------
local_build()
