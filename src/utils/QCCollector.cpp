// $Maintainer: Dragan Haberland, Leo Wurthillini
// $Authors: Dragan Haberland, Leo Wurthillini

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/CsvFile.h>

using namespace OpenMS;
using namespace std;

#include <OpenMS/CONCEPT/QCCollector.h>

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
	registerInputFileList_("in_ProteinQuantifier","<files>", StringList(), "Input files",false,false);
	registerInputFileList_("in_IDMapper","<files>", StringList(), "Input files",false,false);
	registerInputFileList_("in_FalseDiscoveryRate","<files>", StringList(), "Input files",false,false);
	registerInputFileList_("in_FeatureLinkerUnlabeledQT","<files>", StringList(), "Input files",false,false);
	setValidFormats_("in_ProteinQuantifier", ListUtils::create<String>("csv"));
	setValidFormats_("in_IDMapper", ListUtils::create<String>("FeatureXML"));
	setValidFormats_("in_FalseDiscoveryRate", ListUtils::create<String>("IdXML"));
	setValidFormats_("in_FeatureLinkerUnlabeledQT", ListUtils::create<String>("consensusXML"));
	registerOutputFile_("out", "<file>", "", "Output file (mzTab)", true);
  setValidFormats_("out", ListUtils::create<String>("tsv"));
}
ExitCodes main_(int, const char**)
{
StringList ins_ProteinQuantifier = getStringList_("in_ProteinQuantifier");
StringList ins_IDMapper = getStringList_("in_IDMapper");
StringList ins_FalseDiscoveryRate = getStringList_("in_FalseDiscoveryRate");
StringList ins_FeatureLinkerUnlabeledQT = getStringList_("in_FeatureLinkerUnlabeledQT");
String out = getStringOption_("out");
vector<pair<string,FeatureMap>> fvec;
vector<pair<string,CsvFile>> cvec;
vector <pair<string,ConsensusMap>> CMapVec;
vector<pair<string,string>> ivec;
typedef vector<vector<PeptideIdentification>> vpep_id;
vpep_id peps_id;
typedef vector<vector<ProteinIdentification>> vprot_id;
vprot_id prot_id;
if (ins_ProteinQuantifier.size()!=0){
		for(StringList::const_iterator it=ins_ProteinQuantifier.begin();it!=ins_ProteinQuantifier.end();++it){
				CsvFile fl(*it,'	',false,-1);
				cvec.push_back(make_pair("ProteinQuantifier",fl));
				cout<<cvec.size()<<endl;
			}
}
		/*IdXMLFile il;
		vector<PeptideIdentification> pep_ids;
		vector<ProteinIdentification> prot_ids;
		il.load(ins[i], prot_ids, pep_ids);
		peps_id.push_back(pep_ids);
		prot_id.push_back(prot_ids);
		vector<PeptideHit> all_hits;
                for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
        {
          all_hits.insert(all_hits.end(), it->getHits().begin(), it->getHits().end());
        }*/
			//ivec.push_back(make_pair(Tools[i],ins[i]));
	//}
else if (ins_IDMapper.size()!=0){
		for(StringList::const_iterator it=ins_IDMapper.begin();it!=ins_IDMapper.end();++it){
			FeatureMap features;
			FeatureXMLFile().load(*it, features);
		/*for(Size i=0;i<features.size();i++){
		const vector<PeptideIdentification>& pep_ids = features[i].getPeptideIdentifications();
		if(pep_ids.size()>0){cout<<"Vector von PeptideIdentifications hat die Größe "<<pep_ids.size()<<endl;}
		}
		const vector<PeptideIdentification>& pepit_ids = features[0].getPeptideIdentifications();
		vector<double> window;
		window.push_back(features[0].getConvexHull().getBoundingBox().minX());
		window.push_back(features[0].getConvexHull().getBoundingBox().maxX());
		double charge = features[0].getCharge();
		cout<<charge<<endl;
		vector<PeptideHit> all_hits;
		vector<DataProcessing>& test = features.getDataProcessing();
		//cout<<"Data Processing "<<test.size()<<endl;
		for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
        {
          all_hits.insert(all_hits.end(), it->getHits().begin(), it->getHits().end());
        }*/
		fvec.push_back(make_pair("IDMapper",features));
	}
}
else if (ins_FalseDiscoveryRate.size()!=0){
		for(StringList::const_iterator it=ins_FalseDiscoveryRate.begin();it!=ins_FalseDiscoveryRate.end();++it){
			ivec.push_back(make_pair("FalseDiscoveryRate",*it));
		}
	}
else if(ins_FeatureLinkerUnlabeledQT.size()!=0){
		for(StringList::const_iterator it=ins_FeatureLinkerUnlabeledQT.begin();it!=ins_FeatureLinkerUnlabeledQT.end();++it){
			ConsensusMap CMap;
			ConsensusXMLFile().load(*it,CMap);
			CMapVec.push_back(make_pair("FeatureLinkerUnlabeledQT",CMap));
		}
	}

Metriken Metrik(fvec,ivec,cvec,CMapVec,out);
	Metrik.runAllMetrics();
return EXECUTION_OK;
}
};
int main(int argc, const char** argv)
{
QCCollector tool;
return tool.main(argc,argv);
}
