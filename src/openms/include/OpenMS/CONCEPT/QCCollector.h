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
// $Maintainer: Dragan Haberland, Leo Wurthillini, Mohammad el Ali
// $Authors: Dragan Haberland, Leo Wurthillini, Mohammad el Ali 
// --------------------------------------------------------------------------
#include <vector>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>

#include <regex>
#include <iostream>
#include <string>

using namespace std;
using namespace OpenMS;

    class OPENMS_DLLAPI Metriken
    {
    const vector<FeatureMap> FeatureMaps;
    const vector<vector<PeptideIdentification>> Peptides;
    const vector<vector<ProteinIdentification>> Proteins;
    const vector<CsvFile> CsvFiles;

    public:
        Metriken(const vector<FeatureMap> FM, const vector<vector<PeptideIdentification>> PeI, const vector<vector<ProteinIdentification>> PrI, const vector<CsvFile> CsvFiles):
        FeatureMaps(FM),
        Peptides(PeI),
        Proteins(PrI),
	CsvFiles(CsvFiles)
        {
        }
        void runAllMetrics();
    protected:
        int ProteinandPeptideCount_(vector<vector<Size>>&,vector<StringList>&) const;       
   };
    int Metriken::ProteinandPeptideCount_(vector<vector<Size>>& ResultVector,vector<StringList>& Metadata) const {
	for(vector<CsvFile>::const_iterator it = CsvFiles.begin();it != CsvFiles.end(); ++it){
		vector<Size> abundance;
		CsvFile fl=*it;
		StringList currentrow;
		Size i = 0;
		vector<Size> index;
		Size maxRow = fl.rowCount();
		StringList currentmetadata;
		while((fl.getRow(i,currentrow))==false){
			currentmetadata.push_back(currentrow[0]);
			i++;
		}
		while(i<maxRow){
			fl.getRow(i,currentrow);
			if((currentrow[0]=="\"peptide\"")||(currentrow[0]=="\"protein\"")){
				regex e("(abundance_)[0-9]+");
				Metadata.push_back(currentmetadata);
				for(Size i=0;i<currentrow.size();++i){
					smatch m;
        				regex_search(currentrow[i],m,e);
					if(m.size()>=1){
						index.push_back(i);
					}
				}
			i++;
			abundance.resize(index.size());
			}
			else{
				for(Size i=0;i<index.size();i++){
					abundance[i]+=stoi(currentrow[index[i]]);
				}
				i++;
			}
		}
		abundance.size()>0?ResultVector.push_back(abundance):abundance.clear();
		}
	return ResultVector.size()>0?1:0;	
    }
    
void Metriken::runAllMetrics(){
	typedef pair<vector<StringList>,vector<vector<Size>>> ProteinandPeptideResult;
	vector<vector<Size>> ProteinandPeptideCount;
	vector<StringList> Metadata;
	int a = this->ProteinandPeptideCount_(ProteinandPeptideCount,Metadata);
	cout<<Metadata[0].size()<<endl;
	ProteinandPeptideResult PAPR = make_pair(Metadata,ProteinandPeptideCount);
    }
   

    







