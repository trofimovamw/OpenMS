#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/QCMetrics.h>
#include <OpenMS/CONCEPT/QCMetricMap.h>
#include <OpenMS/CONCEPT/QCProteinAndPeptideCount.h>
#include <OpenMS/CONCEPT/QCMBRalignment.h>
#include <vector>
#include <utility>

using namespace OpenMS;
using namespace std;

Metrics::~Metrics()
 {

 }

//Wenn ihr die Metriken schreibt lasst euch bitte ein Int ausgeben 1 wenn erfolgreich, 0  wenn nicht
//Die Daten die ihr erhaltet am besten als MetricMap, die wichtigen Funktionen dafür stehen oben
void Metrics::runAllMetrics()
{
////////////////Metrik1: Protein And Peptide Count /////////////////////////////////
QCProteinAndPeptideCount papc(CFiles_);
MetricMap PeptideCountData;
MetricMap ProteinCountData;
int a = papc.ProtAndPepCount(PeptideCountData,ProteinCountData);
////////////////Metrik2: MBR Alignment /////////////////////////////////
QCMBRalignment mbra(FeatMaps_);
MetricMap raws;
MetricMap outSeq; 
MetricMap outCorrRT; 
MetricMap outOrigRT; 
MetricMap outSpectra;
int b = mbra.MBRAlignment(raws,outSeq,outCorrRT,outOrigRT,outSpectra);
////////////////Metrik3: ....................../////////////////////////////////////
////////////////Metrik4: ....................../////////////////////////////////////
////////////////Metrik5: ....................../////////////////////////////////////

//MzTab Writer:
MzTabFile MzTabOutputFile;
MzTab mztab;
if(a == 1)
{
    int numOfPeptides = PeptideCountData.size();
    int numOfProteins = ProteinCountData.size();
    vector<String> empty;
    vector<String> allPeptides = numOfPeptides > 0 ? PeptideCountData.getStringsByHead("peptide") : empty;
    vector<String> allProteins = numOfProteins > 0 ? ProteinCountData.getStringsByHead("protein") : empty;
    MzTabPeptideSectionRows PepROWS;
    MzTabProteinSectionRows ProtROWS;
    if (numOfPeptides > 0)
    {
      for( int i = 0; i < numOfPeptides ; i++ )
      {
        MzTabPeptideSectionRow PepROW;
        MzTabString PepSeq;
        PepSeq.set(allPeptides[i]);
        PepROW.sequence = PepSeq;
        PepROWS.push_back(PepROW);
      }
    }
    if( numOfProteins > 0 ){
      for (int i = 0; i < numOfProteins ; i++)
      {
        MzTabProteinSectionRow ProtROW;
        MzTabString ProtSeq;
        ProtSeq.set(allProteins[i]);
        ProtROW.description = ProtSeq;
        ProtROWS.push_back(ProtROW);
      }
    }
    mztab.setPeptideSectionRows(PepROWS);
    mztab.setProteinSectionRows(ProtROWS);
}

		if(b == 1){
            int numOfcorrRT = outCorrRT.size();
            int numOfpeptides = outSeq.size();
            int numOforigRT = outOrigRT.size();
            int numOfRaw = raws.size();
            //vector<String> empty;
            //vector<float> emptyI;
            vector<float> allCorrRT = outCorrRT.getFloatsByHead("corrected_RT");
            vector<float> allOriRT = outOrigRT.getFloatsByHead("original_RT");
            vector<String> allPeptides = outSeq.getStringsByHead("sequences");
            vector<String> allRaws = raws.getStringsByHead("raw_file");
            vector<String> spectre_ref = outSpectra.getStringsByHead("spectre_id");
            MzTabPeptideSectionRows PepROWS;
            if(numOfcorrRT>0){
              for(int i = 0; i< numOfpeptides;i++){
                MzTabPeptideSectionRow PepROW;
                MzTabString PepSeq;
                MzTabDouble corrRT;
                MzTabDouble oriRT;
                MzTabSpectraRef ref;
                
                ref.setSpecRef(spectre_ref[i]);
                String out_ref = ref.getSpecRef();
                cout << out_ref << endl;
                PepSeq.set(allPeptides[i]);
                corrRT.set(allCorrRT[i]);
                oriRT.set(allOriRT[i]);
                vector<MzTabDouble> cRTs;
                cRTs.push_back (corrRT);
                cRTs.push_back (oriRT);
                MzTabDoubleList listC;
                listC.set(cRTs);
                PepROW.sequence = PepSeq;
                PepROW.retention_time = listC;
                PepROW.spectra_ref = ref;
                
                
                //typedef std::pair<String, MzTabString> MzTabOptionalColumnEntry
                String ori = to_string(allOriRT[i]);
                MzTabString str = MzTabString(ori);
                MzTabOptionalColumnEntry oRT = make_pair("original_retention_time",str);
                vector<MzTabOptionalColumnEntry> v;
                v.push_back (oRT);
                
                // writes only one optional column?
                MzTabString name = MzTabString(allRaws[i]);
                MzTabOptionalColumnEntry sraw = make_pair("raw_source_file",name);
                v.push_back (sraw);
                PepROW.opt_ = v;
                
                
                PepROWS.push_back(PepROW);

              }
            }
            
            mztab.setPeptideSectionRows(PepROWS);
        }
MzTabOutputFile.store(out_,mztab);
}
