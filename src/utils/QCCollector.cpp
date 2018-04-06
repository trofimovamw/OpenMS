// $Maintainer: Dragan Haberland, Leo Wurthillini
// $Authors: Dragan Haberland, Leo Wurthillini

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/QCCollector.h>

using namespace OpenMS;
using namespace std;

class QCCollector:
	public TOPPBase
{
public: 
	QCCollector():
	TOPPBase("QCCollector","Will collect several mzTabs from several utils.",false)
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
	StringList valid_in = ListUtils::create<String>("IdXML,featureXML,csv");
	registerInputFileList_("in","<files>", StringList(), "Input files");
	setValidFormats_("in", valid_in);
}
ExitCodes main_(int, const char**)
{
StringList ins = getStringList_("in");
typedef vector<FeatureMap> fxml;
fxml fvec;
typedef vector<CsvFile> csv;
csv cvec;
typedef vector<vector<PeptideIdentification>> vpep_id;
vpep_id peps_id;
typedef vector<vector<ProteinIdentification>> vprot_id;
vprot_id prot_id;
for (Size i=0;i<ins.size();++i){
	FileTypes::Type in_type = FileHandler::getType(ins[i]);
	if (in_type == FileTypes::IDXML){
		IdXMLFile il;
		vector<PeptideIdentification> pep_ids;
		vector<ProteinIdentification> prot_ids;
		il.load(ins[i], prot_ids, pep_ids);
		peps_id.push_back(pep_ids);
		prot_id.push_back(prot_ids);
		vector<PeptideHit> all_hits;
                for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
        {
          all_hits.insert(all_hits.end(), it->getHits().begin(), it->getHits().end());
        }

	}
	else if (in_type == FileTypes::FEATUREXML){
		FeatureMap features;
		FeatureXMLFile().load(ins[i], features);
		const vector<PeptideIdentification>& pep_ids = features[0].getPeptideIdentifications();
		vector<double> window;
		window.push_back(features[0].getConvexHull().getBoundingBox().minX());
		window.push_back(features[0].getConvexHull().getBoundingBox().maxX());
		double charge = features[0].getCharge();
		cout<<charge<<endl;
		vector<PeptideHit> all_hits;
		for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
        {
          all_hits.insert(all_hits.end(), it->getHits().begin(), it->getHits().end());
        }
		fvec.push_back(features);
	}
	else if (in_type == FileTypes::CSV){
		CsvFile fl(ins[i],'	',false,-1);
		cvec.push_back(fl);
	}
}
Metriken Metrik(fvec,peps_id,prot_id,cvec);
	Metrik.runAllMetrics();
return EXECUTION_OK;
}
};
int main(int argc, const char** argv)
{
QCCollector tool;
return tool.main(argc,argv);
}	
