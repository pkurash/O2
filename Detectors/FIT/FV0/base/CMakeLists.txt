# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(FV0Base
               SOURCES src/Geometry.cxx
               PUBLIC_LINK_LIBRARIES ROOT::Geom 
                                     FairRoot::Base
                                     O2::SimulationDataFormat)

o2_target_root_dictionary(FV0Base HEADERS include/FV0Base/Geometry.h)

