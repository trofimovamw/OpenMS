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
// $Maintainer: Anton Haberland, Leo Wurth, Mohammad El-Ismail$
// $Authors: Anton Haberland, Leo Wurth, Mohammad El-Ismail $
// --------------------------------------------------------------------------
#include <vector>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>

#include <regex>


using namespace std;
using namespace OpenMS;


    /*class OPENMS_DLLAPI MetrikOutput
    {
    vector<String> Metadata;
    vector<String> Heads;
    vector<String> Data;
    public
    */

    class OPENMS_DLLAPI Metriken
    {
    //das sind die eingelesenen Daten, sie können von den Metriken gelesen aber nicht umgeschrieben werden//
    const vector<FeatureMap> FeatureMaps;			//Alle FeatureXML Datein
    const vector<vector<PeptideIdentification>> PeptidesId;	//Peptide der IDXML's
    const vector<vector<ProteinIdentification>> ProteinsId;	//Proteine der IDXML's
    const vector<CsvFile> CsvFiles;   				//Alle CSV Datein
    const vector<ConsensusMap> ConsensusMaps;    //Alle ConsensusXMLFiles
 
    public:
        Metriken(const vector<FeatureMap> FM, const vector<vector<PeptideIdentification>> PeI, const vector<vector<ProteinIdentification>> PrI,const vector<CsvFile> CF, const vector<ConsensusMap> CM):
        FeatureMaps(FM),
        PeptidesId(PeI),
        ProteinsId(PrI),
        CsvFiles(CF),
        ConsensusMaps(CM)
        {
        }
        void runAllMetrics();
    protected:
        int ProteinAndPeptideCount_(vector<vector<Size>>&,vector<StringList>&) const;
	//->Hier Metriken deklarieren<-//       
	};



    	//->Hier die Metriken definieren<-//


    int Metriken::ProteinAndPeptideCount_(vector<vector<Size>>& Results,vector<StringList>& MetaData) const{
        typedef vector<Size> VSize;
        for(vector<CsvFile>::const_iterator it = CsvFiles.begin(); it!=CsvFiles.end();it++){
            StringList MetaList;
            VSize Abundances;
            VSize AbundancesPosition;
            bool headfinder = false;
            CsvFile fl = *it;
            Size line = 0;
            StringList CurrentRow;
            Size maxRow = fl.rowCount();
            while(fl.getRow(line,CurrentRow)==false){
                MetaList.push_back(CurrentRow[0]);
                line++;
            }
            cout<<line<<" <- line"<<endl;
            while(line<maxRow){
                fl.getRow(line,CurrentRow);
                if((CurrentRow[0]=="\"peptide\"")||(CurrentRow[0]=="\"protein\"")){
                    regex e("(abundance_)[0-9]+");
                    MetaData.push_back(MetaList);
                    for(int k = 2;k<CurrentRow.size();k++){
                        smatch m;
                        regex_search(CurrentRow[k],m,e);
                        if(m.size()>0){
                            AbundancesPosition.push_back(k);
                            headfinder = true;
                        }
                    }
                    Abundances.resize(AbundancesPosition.size(),0);
                }
                else if(headfinder){
                    for(int q = 0;q<Abundances.size();q++){
                        Abundances[q]+= stoi(CurrentRow[AbundancesPosition[q]]);
                    }
                }
                line++;
            }
            cout<<Abundances[0]<<" <- Abundances[0]"<<endl;
            cout<<Abundances[1]<<" <- Abundances[1]"<<endl;
            cout<<Abundances[2]<<" <- Abundances[2]"<<endl;
            headfinder==true ? Results.push_back(Abundances): Abundances.clear();
        }        
        return Results.size()>0?1:0;
    }



    
    void Metriken::runAllMetrics(){
        ////////////////Metrik1: Protein And Peptide Count /////////////////////////////////
        vector<vector<Size>> ProteinAndPeptideCount;
        vector<StringList> ProteinAndPeptideCountMetaData;
        int a = this->ProteinAndPeptideCount_(ProteinAndPeptideCount,ProteinAndPeptideCountMetaData);
	    ////////////////Metrik2: ....................../////////////////////////////////////



    }
   

    






