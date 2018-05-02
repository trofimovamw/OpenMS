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
#include <OpenMS/CONCEPT/QCChargeDistribution.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

//Destructor
QCChargeDistribution::~QCChargeDistribution()
{

}

//Constructor
QCChargeDistribution::QCChargeDistribution(std::vector<std::pair<OpenMS::String,OpenMS::FeatureMap>> files):
  feat_map_(files)
  {
      
  };
  
//Main method to write mztab peptide section data needed for charge distribution plot (PTXQC)
int QCChargeDistribution::ChargeDistribution(MzTab& mztab) const
{

  vector<FeatureMap> maps;
  MzTabPeptideSectionRows rows;
  MzTabPeptideSectionRows mztabRows = mztab.getPeptideSectionRows();
  int pepIDCount = 0;
  
  for(vector<pair<String,FeatureMap>>::const_iterator it = feat_map_.begin();it!=feat_map_.end();++it)
  {
    maps.push_back (it->second);
  }
   
  for (unsigned long m = 0; m < maps.size(); m++)
  {     
  
    String rfile;
    
    if (maps[m].metaValueExists("spectra_data")) 
    {		
      StringList rfiles = maps[m].getMetaValue("spectra_data");	
      rfile = rfiles[0];
    } 
    
    for (vector<Feature>::const_iterator f_it = maps[m].begin(); f_it!=maps[m].end();f_it++)
    {

      vector<PeptideIdentification> pep_id = f_it->getPeptideIdentifications();	
 		  
      if (pep_id.empty()) 
      {
        //Empty lines peptide_data_ for features with retention times	
      }
 				
      else 
      {
 	    for (vector<PeptideIdentification>::iterator p_it = pep_id.begin(); p_it!=pep_id.end(); p_it++) 
 	    {
 	      pepIDCount++;
 	      MzTabPeptideSectionRow row;
 	      MzTabString PepSeq;
 	      MzTabInteger charge;
 		  
 		  //Set sequence and charge	
          vector<PeptideHit> hits = p_it->getHits(); 
          PeptideHit hit = hits[0];
          AASequence seq = hit.getSequence();
          int ch = hit.getCharge();
          PepSeq.set(seq.toString());
          charge.set(ch);
          row.charge = charge;
          row.sequence = PepSeq;
          
    	  //Set optional columns: original RT and source file
          vector<MzTabOptionalColumnEntry> v;                
          MzTabString name = MzTabString(rfile);
          MzTabOptionalColumnEntry sraw = make_pair("raw_source_file",name);
          v.push_back (sraw);
          row.opt_ = v;
                
          rows.push_back(row);
    	  		
         }   
 	  }
 	    
 	}			
    		
  }
  
  //If peptide section was written from features before: append columns
  if (pepIDCount==mztabRows.size()) 
  {
    for (unsigned i = 0; i < mztabRows.size(); i++)
    {
      MzTabPeptideSectionRow mz_r = mztabRows[i];
      MzTabPeptideSectionRow r = rows[i];
      MzTabInteger ch = r.charge;
      vector<MzTabOptionalColumnEntry> v = r.opt_;
      mz_r.charge = ch;
      vector<MzTabOptionalColumnEntry> mz_v = mz_r.opt_;
      for(unsigned j = 0; j < v.size(); j++)
      {
        mz_v.push_back(v[j]);
      }
      mz_r.opt_ = mz_v;
      mztabRows[i] = mz_r;
    }
    mztab.setPeptideSectionRows(mztabRows);
  }
  //Else: append rows
  else
  {
    for (unsigned i = 0; i < rows.size(); i++)
    {
      mztabRows.push_back(rows[i]);
    }
    mztab.setPeptideSectionRows(mztabRows);  
  }  
  
  return rows.size()!=0 ? 1:0;
}
