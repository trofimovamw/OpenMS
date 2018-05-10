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
QCChargeDistribution::QCChargeDistribution(std::vector<OpenMS::FeatureMap> files):
  maps(files)
  {
      
  };
  
//Main method to write mztab peptide section data needed for charge distribution plot (PTXQC)
int QCChargeDistribution::ChargeDistribution(MzTab& mztab) const
{
  MzTabPSMSectionRows rows;
  MzTabPSMSectionRows mztabRows = mztab.getPSMSectionRows();
  vector<MzTabString> unique_ids_;

  int pepIDCount = 0;
  int featCount = 0;
   
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
        //Empty lines of psm_data_ for features with ion charge	
        MzTabPSMSectionRow row;
        MzTabInteger charge;
        int ch = f_it->getCharge();
        charge.set(ch);
        row.charge = charge;
        rows.push_back (row);
        UInt64 id = f_it->getUniqueId();
        unique_ids_.push_back (MzTabString(id));
        featCount++;
        
      }
 				
      else 
      {
 	    for (vector<PeptideIdentification>::iterator p_it = pep_id.begin(); p_it!=pep_id.end(); p_it++) 
 	    {
 	      pepIDCount++;
 	      MzTabPSMSectionRow row;
 	      MzTabInteger charge;
 		  
 		  //Set charge and raw file
          int ch = f_it->getCharge();
          charge.set(ch);
          row.charge = charge;
          
          vector<MzTabOptionalColumnEntry> v;  
          MzTabString name = MzTabString(rfile);
          MzTabOptionalColumnEntry sraw = make_pair("opt_raw_source_file",name);
          v.push_back (sraw);
          UInt64 id = f_it->getUniqueId();
          MzTabOptionalColumnEntry u_id = make_pair("opt_unique_id",MzTabString(id));
          v.push_back (u_id);
          unique_ids_.push_back (MzTabString(id));
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
  //If PSM section was not written from features before: append rows
  if (ids_.empty()) 
  {
    mztab.setPSMSectionRows(rows); //PSM
  }
  //Else: append columns
  //If unique ids are equal: append column (charge)
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
        MzTabInteger charge = r.charge;
        mz_r.charge = charge;
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
