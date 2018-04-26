#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCMBRalignment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

QCMBRalignment::~QCMBRalignment()
{

}

int QCMBRalignment::MBRAlignment(MzTab& mztab) const
{
  vector<String> rawFiles;
  vector<String> sequences;
  vector<float> originalRTs;
  vector<float> correctedRTs;
  vector<String> spectra;
  	  
  vector<FeatureMap> maps;
  	  
  for(vector<pair<String,FeatureMap>>::const_iterator it = feat_map_.begin();it!=feat_map_.end();++it)
  {
    if(it->first=="MapRTTransformer"){maps.push_back(it->second);}
    }
      
      for (vector<FeatureMap>::const_iterator m_it = maps.begin(); m_it!=maps.end();m_it++) {
    		    		
        String rfile = m_it->getMetaValue("spectra_data");	
    	for (vector<Feature>::const_iterator f_it = m_it->begin(); f_it!=m_it->end();f_it++) 
    	{
    		
 		  vector<PeptideIdentification> pep_id = f_it->getPeptideIdentifications();	
 		  if (pep_id.empty()) 
 		  {
 				
 			rawFiles.push_back(rfile);
 			float corrRT = f_it->getRT();
 			float orRT = f_it->getMetaValue("original_RT");
 			correctedRTs.push_back (corrRT);
 			originalRTs.push_back (orRT);
 			sequences.push_back ("NaN");
 			spectra.push_back ("0");
 		  }
 				
 		  else 
 		  {
 		    // pep_id and sequence + spectrum ident
 			for (vector<PeptideIdentification>::iterator p_it = pep_id.begin(); p_it!=pep_id.end();p_it++) 
 			{
    					
    		  vector<PeptideHit> hit = p_it->getHits();
    		  AASequence seq = hit.begin()->getSequence();
    		  float corrRT = p_it->getRT();
    		  float orRT = p_it->getMetaValue("original_RT");
    		  String spec = p_it->getMetaValue("spectrum_reference");  //"spectrum=3462"
    		  String spectrum_ref = spec.substr(9);
    					
    					
    		  rawFiles.push_back (rfile);
    					//cout << rfile << endl;
 			  correctedRTs.push_back (corrRT);
 						//cout << corrRT << endl;
 			  originalRTs.push_back (orRT);
 						//cout << orRT << endl;
 			  sequences.push_back (seq.toString());
 						//cout << seq << endl;
 			  spectra.push_back (spectrum_ref);
 						
    		}
 		 }  
 				
    		
 	  }  
 			
 	}
 	  //vector<String> rawFiles;
  	  //vector<String> sequences;
  	  //vector<float> originalRTs;
  	  //vector<float> correctedRTs;
  	  //vector<String> spectra;
 	  
 	//int a = rawFiles.isEmpty()&&sequences.isEmpty()&&originalRTs.isEmpty()&&correctedRTs.isEmpty()&&psectra.isEmpty()?1:0;
 	
 	 int numOfcorrRT = originalRTs.size();
        
     if(numOfcorrRT>0)
     {
        MzTabPeptideSectionRows PepROWS;
        for(int i = 0; i< numOfcorrRT;i++)
        {
            MzTabPeptideSectionRow PepROW;
            MzTabString PepSeq;
            MzTabDouble corrRT;
            MzTabDouble oriRT;
                
            if (spectra[i]=="0") 
            {
              
              PepSeq.set(sequences[i]);
              corrRT.set(correctedRTs[i]);
              //oriRT.set(originalRTs[i]);
              vector<MzTabDouble> cRTs;
              cRTs.push_back (corrRT);
              //cRTs.push_back (oriRT);
              MzTabDoubleList listC;
              listC.set(cRTs);
              PepROW.sequence = PepSeq;
              PepROW.retention_time = listC;
                
                //typedef std::pair<String, MzTabString> MzTabOptionalColumnEntry
              String ori = to_string(originalRTs[i]);
              MzTabString str = MzTabString(ori);
              MzTabOptionalColumnEntry oRT = make_pair("original_retention_time",str);
              vector<MzTabOptionalColumnEntry> v;
              v.push_back (oRT);
                
              MzTabString name = MzTabString(rawFiles[i]);
              MzTabOptionalColumnEntry sraw = make_pair("raw_source_file",name);
              v.push_back (sraw);
              PepROW.opt_ = v;
                
              PepROWS.push_back(PepROW);

           }  
                
           else if (spectra[i]!="0") 
           { 
             //vector<String> rawFiles;
  	         //vector<String> sequences;
  	         //vector<float> originalRTs;
  	         //vector<float> correctedRTs;
  	         //vector<String> spectra;
  	  
             MzTabSpectraRef ref;
             ref.setSpecRef(spectra[i]);
             String out_ref = ref.getSpecRef();
             PepSeq.set(sequences[i]);
             corrRT.set(correctedRTs[i]);
             vector<MzTabDouble> cRTs;
             cRTs.push_back (corrRT);
             MzTabDoubleList listC;
             listC.set(cRTs);
             PepROW.sequence = PepSeq;
             PepROW.retention_time = listC;
             PepROW.spectra_ref = ref;
                
                //typedef std::pair<String, MzTabString> MzTabOptionalColumnEntry
             String ori = to_string(originalRTs[i]);
             MzTabString str = MzTabString(ori);
             MzTabOptionalColumnEntry oRT = make_pair("original_retention_time",str);
             vector<MzTabOptionalColumnEntry> v;
             v.push_back (oRT);
                
             MzTabString name = MzTabString(rawFiles[i]);
             MzTabOptionalColumnEntry sraw = make_pair("raw_source_file",name);
             v.push_back (sraw);
             PepROW.opt_ = v;
                
             PepROWS.push_back(PepROW);

          }  
                

       }
    mztab.setPeptideSectionRows(PepROWS);
    }
    
            
  
    return numOfcorrRT!=0 ? 1:0;
}

