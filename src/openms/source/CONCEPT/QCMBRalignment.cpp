// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Maria Trofimova $
// $Authors: Maria Trofimova $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCMBRalignment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

//Destructor
QCMBRalignment::~QCMBRalignment()
{

}

//Constructor
QCMBRalignment::QCMBRalignment(std::vector<OpenMS::FeatureMap> files):
  maps(files)
  {
      
  };
  
//Main method to write mztab peptide section data needed for MBR alignment plot (PTXQC)
int QCMBRalignment::MBRAlignment(MzTab& mztab) const
{

  //vector<FeatureMap> maps;
  MzTabPSMSectionRows rows;
  MzTabPSMSectionRows mztabRows = mztab.getPSMSectionRows();
  vector<MzTabString> unique_ids_;
  int pepIDCount = 0;
  
  //maps = feat_map_;
   
  for (Size m = 0; m < maps.size(); m++)
  {     
  
    String rfile;
    //Keep index of current map for spectra reference 
    Size run = m; 
    
    if (maps[m].metaValueExists("spectra_data")) 
    {		
      StringList rfiles = maps[m].getMetaValue("spectra_data");	
      rfile = rfiles[0];
    } 
    
    for (vector<Feature>::const_iterator f_it = maps[m].begin(); f_it!=maps[m].end();f_it++)
    {

      vector<PeptideIdentification> pep_id = f_it->getPeptideIdentifications();
      UInt64 unique_id = f_it->getUniqueId();
 		  
      if (pep_id.empty()) 
      {
        //If we want to write empty lines peptide_data_ for features with retention times
        //	
      }
 				
      else 
      {
        //Iterate over peptide hits of a feature and write data into mztab
        //In case of 2 and more hits take the first one 
 	    for (vector<PeptideIdentification>::iterator p_it = pep_id.begin(); p_it!=pep_id.end(); p_it++) 
 	    {
 	      pepIDCount++;
 	      MzTabPSMSectionRow row;
 	      MzTabString PepSeq;
          MzTabDouble correctRT;
          MzTabDouble oriRT;
 		  
 		  //Set sequence	
          vector<PeptideHit> hits = p_it->getHits(); 
          PeptideHit hit = hits[0];
          AASequence seq = hit.getSequence();
          PepSeq.set(seq.toString());
          row.sequence = PepSeq;
          
          //Set corrected RTs 
          float corrRT = f_it->getRT(); 
          correctRT.set(corrRT);
          vector<MzTabDouble> cRTs;
          cRTs.push_back (correctRT);
          MzTabDoubleList listC;
          listC.set(cRTs);
          row.retention_time = listC;
          
    	  float orRT;
    	  if (p_it->metaValueExists("original_RT")) {orRT = p_it->getMetaValue("original_RT");}
 	      else {throw "Retention time was not written in the feature.";}
 	      
 	      //Set spectrum reference
    	  String spectrum_ref = p_it->getMetaValue("spectrum_reference");
    	  const String& const_spec = spectrum_ref;
    	  MzTabSpectraRef spec_ref;
    	  Size index = run+1;
    	  spec_ref.setMSFile(index);
    	  spec_ref.setSpecRef(spectrum_ref);
    	  spec_ref.setSpecRefFile(const_spec);
    	  row.spectra_ref = spec_ref;
    	      	  
    	  //Set optional columns: original RT and source file
    	  String ori = to_string(orRT);
          MzTabOptionalColumnEntry oRT = make_pair("opt_original_retention_time",MzTabString(ori));
          vector<MzTabOptionalColumnEntry> v;
          v.push_back (oRT);
          
          UInt64 id = f_it->getUniqueId();
          MzTabOptionalColumnEntry u_id = make_pair("opt_unique_id",MzTabString(id));
          v.push_back (u_id);
          unique_ids_.push_back (MzTabString(id));
                
          MzTabOptionalColumnEntry sraw = make_pair("opt_raw_source_file",MzTabString(rfile));
          v.push_back (sraw);
          row.opt_ = v;
                
          rows.push_back(row);
    	  		
         }   
 	  }
 	    
 	}			
    		
  }
  //Write unique ids from existing mztab data structure passed to the constructor
  vector<MzTabString> ids_;
  for(vector<MzTabPSMSectionRow>::const_iterator it = mztabRows.begin();it!=mztabRows.end();++it) 
  { 
    vector<MzTabOptionalColumnEntry> opt = it->opt_;
    for (vector<MzTabOptionalColumnEntry>::const_iterator o_it = opt.begin();o_it!=opt.end();++o_it)
    {
      if(o_it->first=="opt_unique_id")
        {
          ids_.push_back(o_it->second);
        }
    }
  }
  
  //Merge new lines and existing lines. Based on unique ids (UniqueIdInterface)
  //If PSM section was not written before: append rows
  
  if (ids_.empty()) 
  {
    mztab.setPSMSectionRows(rows); 
  }
  //Else: append rows
  //If unique ids are equal: append columns (relevant for metric)
  //If not: insert unique id in dictionary and add new line to mztab
  else
  { 
    //Assign vectors for accurate merging
    vector<MzTabString> ids;
    vector<MzTabString> unique_ids;
    if (ids_.size() < unique_ids_.size())
    {
      ids = ids_;
      unique_ids = unique_ids_;   
    } 
    else 
    {
      ids = unique_ids_;
      unique_ids = ids_;
    }
    
    for (unsigned i = 0; i < unique_ids.size(); i++)
    {
      if (ids[i].toCellString().compare(unique_ids[i].toCellString())==0)
      {
        MzTabPSMSectionRow mz_r = mztabRows[i];
        MzTabPSMSectionRow r = rows[i];
        MzTabString seq = r.sequence;
        mz_r.sequence = seq;
        mz_r.retention_time = r.retention_time;
        mz_r.spectra_ref = r.spectra_ref;
        vector<MzTabOptionalColumnEntry> v = mz_r.opt_;
        for (Size i = 0; i < r.opt_.size(); i++)
        {
          v.push_back (r.opt_[i]);
        }
        mz_r.opt_ = v;
        mztabRows[i] = mz_r;
        
      }
      else 
      {
        
        MzTabPSMSectionRows split_f (mztabRows.begin(), mztabRows.begin()+i);
        MzTabPSMSectionRows split_e (mztabRows.begin()+i, mztabRows.end());
        split_f.push_back (rows[i]);
        for (unsigned t = 0; t < split_e.size(); t++)
        {
          split_f.push_back (split_e[t]);
        }
        mztabRows = split_f;
        vector<MzTabString> id_f (ids.begin(), ids.begin()+i);
        vector<MzTabString> id_e (ids.begin()+i, ids.end());
        id_f.push_back (unique_ids[i]);
        for (unsigned t = 0; t < id_e.size(); t++)
        {
          id_f.push_back (id_e[t]);
        }
        ids = id_f;
      }
    
    }
    mztab.setPSMSectionRows(mztabRows);  
  }    
  
  return rows.size()!=0 ? 1:0;
}
