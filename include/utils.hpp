#pragma once

#include <vector>
#include <TTree.h>
#include <TString.h>

std::vector<Int_t> numberOfHitsPerEvent(TTree*, const TString&, bool);
std::vector<Int_t> mapChannelToLayer(TTree*, const TString&, bool);