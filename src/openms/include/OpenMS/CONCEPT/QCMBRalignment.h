#ifndef OPENMS_CONCEPT_MBRAlignment_H
#define OPENMS_CONCEPT_MBRAlignment_H
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/CONCEPT/QCMetricMap.h>

class OPENMS_DLLAPI QCMBRalignment
{
  public:
    QCMBRalignment(std::vector<std::pair<OpenMS::String,OpenMS::FeatureMap>> files):
      FeatMap_(files)
      {
      };
    ~QCMBRalignment();
    int MBRAlignment(MetricMap&,MetricMap&,MetricMap&,MetricMap&,MetricMap&) const; //std::vector<std::pair<OpenMS::String,OpenMS::FeatureMap>> FeatMaps_
  protected:
    std::vector<std::pair<OpenMS::String,OpenMS::FeatureMap>> FeatMap_;
};
#endif
