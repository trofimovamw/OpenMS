#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCMBRalignment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

QCMBRalignment::~QCMBRalignment()
{

}

int QCMBRalignment::MBRAlignment(MetricMap& raws, MetricMap& outSeq, MetricMap& outCorrRT, MetricMap& outOrigRT, MetricMap& outSpectra) const
{
      vector<String> rawFiles;
  	  vector<String> sequences;
  	  vector<float> originalRTs;
  	  vector<float> correctedRTs;
  	  vector<String> spectra;
  	  
  	  vector<FeatureMap> maps;
  	  
  	  for(vector<pair<String,FeatureMap>>::const_iterator it = FeatMap_.begin();it!=FeatMap_.end();++it){
        if(it->first=="MapRTTransformer"){maps.push_back(it->second);}
      }
      
    
      for (vector<FeatureMap>::const_iterator m_it = maps.begin(); m_it!=maps.end();m_it++) {
    		
    		//cout << m_it->getMetaValue("spectra_data") << endl;
    		
    		String rfile = m_it->getMetaValue("spectra_data");
    		
    		for (vector<Feature>::const_iterator f_it = m_it->begin(); f_it!=m_it->end();f_it++) {
    		
 				vector<PeptideIdentification> pep_id = f_it->getPeptideIdentifications();
 				
 				if (pep_id.empty()) {
 					//Nan	
 					rawFiles.push_back(rfile);
 					cout << rfile << endl;
 					float corrRT = f_it->getRT();
 					float orRT = f_it->getMetaValue("original_RT");
 					correctedRTs.push_back (corrRT);
 					cout << corrRT << endl;
 					originalRTs.push_back (orRT);
 					cout << orRT << endl;
 					sequences.push_back ("NaN");
 					cout << "NaN" << endl;
 					spectra.push_back ("0");
 					cout << 0 << endl;
 				}
 				
 				else {
 					// pep_id and sequence + spectrum ident
 					for (vector<PeptideIdentification>::iterator p_it = pep_id.begin(); p_it!=pep_id.end();p_it++) {
    					//cout << "corrected_RT: " << p_it->getRT()  << endl;
    					//cout << "original_RT: " << p_it->getMetaValue("original_RT") << endl;
    					
    					vector<PeptideHit> hit = p_it->getHits();
    					AASequence seq = hit.begin()->getSequence();
    					float corrRT = p_it->getRT();
    					float orRT = p_it->getMetaValue("original_RT");
    					String spec = p_it->getMetaValue("spectrum_reference");  //"spectrum=3462"
    					String spectrum_ref = spec.substr(9);
    					//cout << spectrum_ref << endl;
    					//String spectrumID = p_it->getIdentifier();
    					
    					rawFiles.push_back (rfile);
    					cout << rfile << endl;
 						correctedRTs.push_back (corrRT);
 						cout << corrRT << endl;
 						originalRTs.push_back (orRT);
 						cout << orRT << endl;
 						sequences.push_back (seq.toString());
 						cout << seq << endl;
 						spectra.push_back (spectrum_ref);
 						cout << spec << endl;
    					
						//for (vector<PeptideHit>::const_iterator h_it = hit.begin(); h_it!=hit.end();h_it++) {
			    		//	cout << h_it->getSequence() << endl;
    					//}
    				}
 			    }
 				
    		
 			}
 			
 	  }
 	  //MetricMap& raws, MetricMap& outSeq, MetricMap& outCorrRT, MetricMap& outOrigRT, MetricMap& outSpectra
 	  raws.pushDataString("raw_file",rawFiles);
 	  outSeq.pushDataString("sequences",sequences);
 	  outCorrRT.pushDataFloat("corrected_RT",correctedRTs);
 	  outOrigRT.pushDataFloat("original_RT",originalRTs);
 	  outSpectra.pushDataString("spectre_id",spectra);
 	  
 	  return raws.isEmpty()&&outSeq.isEmpty()&&outCorrRT.isEmpty()&&outOrigRT.isEmpty()&&outSpectra.isEmpty()?0:1;
 	  
}
