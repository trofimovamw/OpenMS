#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCMS2IdentificationRate.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <boost/regex.hpp>
#include <OpenMS/FORMAT/PercolatorOutfile.h>


using namespace OpenMS;
using namespace std;

QCMS2IdentificationRate::~QCMS2IdentificationRate(){

}

int QCMS2IdentificationRate::MS2IDRateidentifier_( MzTabFile& MzTabOutputFile, String out)
{
  MzTab mztab;
  vector<String> idXMLFiles;
  Size control;
  vector<String> irawfiles;
  boost::regex rgx("[A-Za-z0-9]+[.]idXML");
  boost::regex rgx2(".idXML");
  for(vector<pair<String,String>>::const_iterator it = ivec_.begin();it!=ivec_.end();++it)
  {
    if(it->first=="Post_FalseDiscoveryRate"){
      idXMLFiles.push_back(it->second);
      boost::smatch match;
      boost::regex_search(it->first,match,rgx);
      String rawfile = boost::regex_replace(match[0].str(),rgx2,".mzML");
      irawfiles.push_back(rawfile);
    }
  }
  for(vector<String>::const_iterator it=idXMLFiles.begin();it!=idXMLFiles.end();++it){
    IdXMLFile il;
		vector<PeptideIdentification> pep_ids;
		vector<ProteinIdentification> prot_ids;
    Size spectra;
    Size chromatograms;
    String test = "/home/leo/Schreibtisch/Studium/5.Semester/SoftwarepraktikumOpenMS/openms/OpenMS/share/OpenMS/examples/BSA/BSA1.mzML";
		il.load(*it, prot_ids, pep_ids);
    MzMLFile mzmlfile;
    PeakMap exp;
    mzmlfile.getOptions().addMSLevel(2);
    mzmlfile.loadSize(test,spectra,chromatograms);
    //mzmlfile.load(test,exp);
    //vector<Int> level = mzmlfile.getOptions().getMSLevels();

    //cout<<"hat spectra: "<<spectra<<" und chromatograms: "<<chromatograms<<endl;
  }
}
