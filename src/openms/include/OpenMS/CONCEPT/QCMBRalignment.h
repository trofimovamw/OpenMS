#ifndef OPENMS_CONCEPT_MBRAlignment_H
#define OPENMS_CONCEPT_MBRAlignment_H
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <vector>

class OPENMS_DLLAPI QCMBRalignment
{
std::vector<std::pair<OpenMS::String,OpenMS::FeatureMap>> feat_map_;
  public:
    QCMBRalignment(std::vector<std::pair<OpenMS::String,OpenMS::FeatureMap>> files):
      feat_map_(files)
      {
      };
    ~QCMBRalignment();
    int MBRAlignment(OpenMS::MzTab&) const; //std::vector<std::pair<OpenMS::String,OpenMS::FeatureMap>> FeatMaps_
  
};
#endif
