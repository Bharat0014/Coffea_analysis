import numpy as np
import hist
from coffea.analysis_tools import Weights
from coffea import processor
import awkward as ak
from coffea.nanoevents import  NanoEventsFactory, NanoAODSchema
from coffea.lookup_tools import extractor
from coffea.analysis_tools import PackedSelection
import warnings
warnings.filterwarnings("ignore")
import correctionlib

class NanoProcessor(processor.ProcessorABC):
    def __init__(self):
        pass
        # Define histograms with hist class- https://hist.readthedocs.io/
        
       
    @property
    ## accumulating processed results
    def accumulator(self):
        return self._accumulator

   
    def process(self, events):
      try: 
        '''
        create the dictionary contains 
        this structure is likely used within a processing framework
        to keep track of various statistics or quantities of interest
        during the processing of data.
        
        sumw: This seems to be an accumulator for summing weights or similar quantities.
        It's initialized with a default value of 0.0, indicating it will store a running total of some quantity.
        NoSel: This appears to be a nested structure with two keys: selEvents and wtEvents. 
        This structure seems to be used for accumulating values related to events that 
        do not pass a selection criterion (NoSel stands for "No Selection").
        selEvents: This seems to accumulate the count of selected events.
        wtEvents: This seems to accumulate some weighted quantity associated with the selected events.
        '''
          
        output = {
            "sumw" : processor.defaultdict_accumulator(float),
            "NoSel":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float),
                },            
            "HLT":{
                "selEvents" : processor.defaultdict_accumulator(float),
                "wtEvents" : processor.defaultdict_accumulator(float),
                },
        }
        
        #checking for data
        isRealData = not hasattr(events, "Generator")
        dataset = events.metadata["dataset"]
        if isRealData:
            output["sumw"] = len(events)
            #print("dfd")
        else:
            output["sumw"] = ak.sum(events.Generator.weight)

        
        
        ####################
        #    Selections    #
        ####################
        
        selection = PackedSelection()
        selection.add("HLT",  events.HLT.IsoMu27)
        selectionList = {
        "NoSel":{},
        "HLT":{"HLT": True},
        }
        #apply loop in selection list
        for region, cuts in selectionList.items():
            event_level = selection.require(**cuts)
            num_selected_events = float(sum(event_level))
            output[region]["selEvents"] = num_selected_events
            if sum(event_level) == 0: # if no event selected then loop endes here, it not further in code.
                print("No event selected")
                output[region]["wtEvents"] = float(sum(event_level))
                #print(output)
                if region == 'Mll':
                    return {dataset: output}
                continue
            
            
            ####################
            # Weight & Geninfo #
            ####################
            weights = Weights(sum(event_level), storeIndividual=True)# information of weights objest about dataset size it will manazing
                   
            if  not isRealData:  # Weights for MC samples
                try:
                    weights.add("genWt", weight = events[event_level].Generator.weight/fabs(events[event_level].Generator.weight))
                    #print("Added genWt")
                except:
                    #print(f"genWeight is not there for dataset {dataset}; adding 1s as genWeights")
                    weights.add("genWt", weight = np.ones(sum(event_level), dtype = float))
                try:
                    weights.add("LHEWeightSign", weight = events[event_level].LHEWeight.originalXWGTUP/fabs(events[event_level].LHEWeight.originalXWGTUP))
                    #print("Added LHEWeightSign")
                except:
                    #print(f"LHEWeight is not there for dataset {dataset}; adding +1s as LHEWeightSign")
                    weights.add("LHEWeightSign", weight = np.ones(sum(event_level), dtype = float))

                
                try:
                    weights.add("PUWt", weight = events[event_level].puWeight, weightUp = events[event_level].puWeightUp, weightDown = events[event_level].puWeightDown)
                    #print("Added PUWt")
                except:
                    #print(f"puWeight is not there for dataset {dataset}; adding 1s as puWeights")
                    weights.add("PUWt", weight = np.ones(sum(event_level), dtype = float))

            
    
            
            

                ##################################
                # getting information of weights #
                ##################################
            
            output[region]["wtEvents"] = float(sum(weights.weight()))
            
        return {dataset: output}
      except Exception as e:
        print(f"Error processing {events.metadata['filename']}: {e}")
        raise
    def postprocess(self, accumulator):
        return accumulator

    #Here processor ends

#======================================= Data Set========================================
fileset = {
    'Data':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/Data_mu/Run2017B_mu/SingleMuon/Tree_27_Feb24_Run2017B_mu/240227_035410/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/Data_mu/Run2017C_mu/SingleMuon/Tree_27_Feb24_Run2017C_mu/240227_035504/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/Data_mu/Run2017D_mu/SingleMuon/Tree_27_Feb24_Run2017D_mu/240227_035555/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/Data_mu/Run2017E_mu/SingleMuon/Tree_27_Feb24_Run2017E_mu/240227_035645/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/Data_mu/Run2017F_mu/SingleMuon/Tree_27_Feb24_Run2017F_mu/240227_035735/0000/*.root'
    ], 
    'WJets0J':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/WJetsToLNu_0J/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/Tree_14_Mar24_MCUL2017_WJetsToLNu_0J/240314_052449/0000/*.root'
    ],
    'WJets1J':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/WJetsToLNu_1J/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/Tree_14_Mar24_MCUL2017_WJetsToLNu_1J/240314_052105/0000/*.root'
    ],

    'WJets2J':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/WJetsToLNu_2J/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/Tree_14_Mar24_MCUL2017_WJetsToLNu_2J/240314_052219/0000/*.root'
    ],    
       
    'DYJetsToLL':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/DYJetsToLL/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/Tree_14_Mar24_MCUL2017_DYJetsToLL/240314_052334/0000/*.root'
    ],
    'Tchannel':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/Tchannel/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/Tree_09_Feb24_MCUL2017_Tchannel/240209_134248/0000/*.root'
    ],
    'Tbarchannel':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/Tbarchannel/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/Tree_09_Mar24_MCUL2017_Tbarchannel/240309_145624/0000/*.root'
    ],

    'Schannel':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/Schannel/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/Tree_09_Feb24_MCUL2017_Schannel/240209_133837/0000/*.root'
    ],
    
    'WWTo2L2Nu':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/WWTo2L2Nu/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/Tree_09_Feb24_MCUL2017_WWTo2L2Nu/240209_135338/0000/*.root'
    ],
    'WWTo1L1Nu2Q':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/WWTo1L1Nu2Q/WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/Tree_09_Feb24_MCUL2017_WWTo1L1Nu2Q/240209_134352/0000/*.root'
    ],
    'WZTo3LNu':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/WZTo3LNu/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/Tree_09_Feb24_MCUL2017_WZTo3LNu/240209_135237/0000/*.root'
    ],
    'WZTo2Q2L':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/WZTo2Q2L/WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/Tree_09_Feb24_MCUL2017_WZTo2Q2L/240209_135637/0000/*.root'
    ],
    'WZTo1L1Nu2Q':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/WZTo1L1Nu2Q/WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/Tree_09_Feb24_MCUL2017_WZTo1L1Nu2Q/240209_133731/0000/*.root'
    ],
    'TTWJetsToLNu':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/TTWJetsToLNu/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/Tree_09_Feb24_MCUL2017_TTWJetsToLNu/240209_134934/0000/*.root'
    ],
    'TTWJetsToQQ':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/TTWJetsToQQ/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/Tree_09_Feb24_MCUL2017_TTWJetsToQQ/240209_135939/0000/*.root'
    ],

    'TTZToLLNuNu':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/TTZToLLNuNu/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/Tree_09_Feb24_MCUL2017_TTZToLLNuNu/240209_135738/0000/*.root'
    ],
    'TTZToLL':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/TTZToLL/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/Tree_09_Feb24_MCUL2017_TTZToLL/240209_134712/0000/*.root'
    ],
    'TTZToNuNu':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/TTZToNuNu/TTZToNuNu_TuneCP5_13TeV-amcatnlo-pythia8/Tree_09_Feb24_MCUL2017_TTZToNuNu/240209_135839/0000/*.root'
    ],
    'TTZToQQ':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/TTZToQQ/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/Tree_09_Feb24_MCUL2017_TTZToQQ/240209_134608/0000/*.root'
    ],
    
    'twTop':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/tw_top/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/Tree_09_Feb24_MCUL2017_tw_top/240209_133942/0000/*.root'
    ],

    'twAntiTop':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/tw_antitop/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/Tree_09_Feb24_MCUL2017_tw_antitop/240209_134826/0000/*.root'
    ],
    
     'ZZTo4L':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/ZZTo4L/ZZTo4L_5f_TuneCP5_13TeV-madgraphMLM-pythia8/Tree_09_Feb24_MCUL2017_ZZTo4L/240209_134045/0000/*.root'
    ],

     'ZZTo2Q2L':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/ZZTo2Q2L/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/Tree_09_Feb24_MCUL2017_ZZTo2Q2L/240209_135036/0000/*.root'
    ],
    'ZZTo2L2Nu':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/ZZTo2L2Nu/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/Tree_09_Feb24_MCUL2017_ZZTo2L2Nu/240209_133627/0000/*.root'
    ],
    
    'ttbarHadronic': [
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/TTToHadronic/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/Tree_12_Mar24_MCUL2017_TTToHadronic/240312_053322/0000/*.root'
    ],

    'ttbarSemiLeptonic': [
        '/persistent/data1/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/ttbar_SemiLeptonic/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/Tree_09_Feb24_MCUL2017_ttbar_SemiLeptonic/240209_135137/0000/*.root'
    ],

   
    'ttbarFullyLeptonic': [
        '/persistent/data1/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/ttbar_FullyLeptonic/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/Tree_09_Feb24_MCUL2017_ttbar_FullyLeptonic/240209_134459/0000/*.root'
    ],

    'LQLQToTopMu':[
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-200/LQLQToTopMu_M-200_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-200/240226_054928/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-300/LQLQToTopMu_M-300_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-300/240226_054727/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-400/LQLQToTopMu_M-400_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-400/240226_054348/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-500/LQLQToTopMu_M-500_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-500/240226_054150/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-600/LQLQToTopMu_M-600_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-600/240226_053852/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-700/LQLQToTopMu_M-700_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-700/240226_054620/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-800/LQLQToTopMu_M-800_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-800/240226_054249/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-900/LQLQToTopMu_M-900_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-900/240226_053951/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-1000/LQLQToTopMu_M-1000_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-1000/240226_055127/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-1200/LQLQToTopMu_M-1200_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-1200/240226_054829/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-1400/LQLQToTopMu_M-1400_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-1400/240226_054448/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-1700/LQLQToTopMu_M-1700_TuneCP5_13TeV_pythia8/Tree_26_Feb24_MCUL2017_LQLQToTopMu_M-1700/240226_055027/0000/*.root',
        '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/LQLQToTopMu_M-2000/LQLQToTopMu_M-2000_TuneCP5_13TeV_pythia8/Tree_05_Mar24_MCUL2017_LQLQToTopMu_M-2000/240305_192140/0000/*.root'   
    ],
    
}
#=====================================Here we add wildcard path ==========================================
import glob

'''
# Define the base directory where QCD datasets are located
qcd_base_dir = '/persistent/data2/ph20059/RUN2_UL/Tree_crab/SEVENTEEN/MC/QCD/'

# Use glob to recursively search for all root files within QCD directory and its subdirectories
qcd_files = glob.glob(qcd_base_dir + '**/*.root', recursive=True)

# Create a dictionary with QCD dataset name as key and list of root files as value
qcd_datasets = {'QCD': qcd_files}

# Update the existing fileset with the QCD datasets
fileset.update(qcd_datasets)
'''

resolved_fileset = {}
for dataset, paths in fileset.items():
    resolved_paths = []
    for path in paths:
        if '*' in path:  # Check if the path contains a wildcard
            # Resolve wildcard and extend the resolved paths list
            resolved_paths.extend(glob.glob(path))
        else:
            # No wildcard, just add the path directly
            resolved_paths.append(path)
    resolved_fileset[dataset] = resolved_paths

resolved_fileset

#================================processing the processor====================================
processor_instance = NanoProcessor()
output = processor.run_uproot_job(
    resolved_fileset,
    treename='Events',
    processor_instance=processor_instance,
    executor=processor.futures_executor,
    executor_args={'workers': 40, 'schema': processor.NanoAODSchema}
)
#================printing sellected events===============================================
for dataset, result in output.items():
    print(f"Results for {dataset}:")
    for region, data in result.items():
        if region != 'sumw':
            print(f"  {region}: {data['selEvents']} selected events")
#===============printting Wtevents======================================================
for dataset, result in output.items():
    print(f"Results for {dataset}:")
    for region, data in result.items():
        if region != 'sumw':
            print(f"  {region}: {data['wtEvents']} WT events")
#============================adding Luminosity wt====================================
import MCsample_Nevent_Xsec2
import pandas as pd

# Corrected selection points, ensure all your selections are listed here
selection_points = [
    "NoSel", "HLT"
]

data_for_df = []

target_luminosity = float(49.8*1000)

# Assuming 'output' is structured as {dataset: {region: {"selEvents": value, "wtEvents": value}}}
datasets = output.keys()  # Dynamically get dataset names from the output
for dataset in datasets:
    results = [dataset]  # Start with dataset name
    if dataset in MCsample_Nevent_Xsec2.UL2017:
        num_events = float(MCsample_Nevent_Xsec2.UL2017[dataset][0])
        cross_section = float(MCsample_Nevent_Xsec2.UL2017[dataset][1])
        luminosity_weight = (target_luminosity * cross_section) / num_events
        print(f"{dataset:<20} {num_events:<15.1f} {cross_section:<10.4f} {luminosity_weight:<10.4f}")
    else:
        # For 'Data' or any dataset not found in the dictionary, set a flag value or handle accordingly
        #print(f"Skipping luminosity weight calculation for {dataset} as it's not found in MCsample_Nevent_Xsec2.UL2017")
        luminosity_weight = 1.0  # Or set to None, if you prefer

    for point in selection_points:
        # Use the correct key "wtEvents" to access weighted events count
        wt_events = output[dataset].get(point, {}).get("wtEvents", 0)  # Default to 0 if not found
        
        # Apply luminosity weight
        wt_events *= luminosity_weight
        
        # Round the weighted events count to two decimal places
        wt_events_rounded = round(wt_events, 2)
        
        results.append(wt_events_rounded)  
    
    data_for_df.append(results)

columns = ["Dataset"]
for point in selection_points:
    columns.append(point)

# Creating DataFrame with the corrected structure and column names
dfff = pd.DataFrame(data_for_df, columns=columns)

# Printing the DataFrame with rounded weighted event counts
print(dfff.to_string(index=False))


import MCsample_Nevent_Xsec2
import pandas as pd

# Corrected selection points, ensure all your selections are listed here
selection_points = [
    "NoSel", "HLT"
]

data_for_df = []

target_luminosity = float(49.8*1000)

# Assuming 'output' is structured as {dataset: {region: {"selEvents": value, "wtEvents": value}}}
datasets = output.keys()  # Dynamically get dataset names from the output
for dataset in datasets:
    results = [dataset]  # Start with dataset name
    if dataset in MCsample_Nevent_Xsec2.UL2017:
        num_events = float(MCsample_Nevent_Xsec2.UL2017[dataset][0])
        cross_section = float(MCsample_Nevent_Xsec2.UL2017[dataset][1])
        luminosity_weight = (target_luminosity * cross_section) / num_events
        #print(f"{dataset:<20} {num_events:<15.1f} {cross_section:<10.4f} {luminosity_weight:<10.4f}")
    else:
        # For 'Data' or any dataset not found in the dictionary, set a flag value or handle accordingly
        #print(f"Skipping luminosity weight calculation for {dataset} as it's not found in MCsample_Nevent_Xsec2.UL2017")
        luminosity_weight = 1.0  # Or set to None, if you prefer

    for point in selection_points:
        # Use the correct key "wtEvents" to access weighted events count
        wt_events = output[dataset].get(point, {}).get("wtEvents", 0)  # Default to 0 if not found
        
        # Apply luminosity weight
        wt_events *= luminosity_weight
        
        # Round the weighted events count to two decimal places
        wt_events_rounded = round(wt_events, 2)
        
        results.append(wt_events_rounded)  
    
    data_for_df.append(results)

columns = ["Dataset"]
for point in selection_points:
    columns.append(point)

# Creating DataFrame with the corrected structure and column names
dff = pd.DataFrame(data_for_df, columns=columns)

mc_subset = dff[~dff['Dataset'].isin(['Data', 'LQLQToTopMu'])]

# Calculate total MC excluding the specified datasets
total_mc = mc_subset.sum(numeric_only=True)
total_mc['Dataset'] = 'Total MC'

# Append the total MC row to the DataFrame
dff = pd.concat([dff, pd.DataFrame([total_mc])], ignore_index=True)

# Find the 'Data' row for the ratio calculation
data_row = dff.loc[dff['Dataset'] == 'Data'].select_dtypes(include=[float, int])

# Find the 'Total MC' row for the ratio calculation
total_mc_row = dff.loc[dff['Dataset'] == 'Total MC'].select_dtypes(include=[float, int])

# Calculate the Data/MC ratio
data_mc_ratio = data_row.values[0] / total_mc_row.values[0]
data_mc_ratio_series = pd.Series(data_mc_ratio, index=data_row.columns)
data_mc_ratio_series['Dataset'] = 'Data/Total MC'


# Append the Data/Total MC ratio row to the DataFrame
dff = pd.concat([dff, pd.DataFrame([data_mc_ratio_series])], ignore_index=True)
dff = dff.round(1)
# Print or display the updated DataFrame
#print(dff.to_string(index=False))

# Assuming 'dff' is your DataFrame after the initial processing

# Extract the rows for 'LQLQToTopMu', 'Data', 'Total MC', and 'Data/Total MC'
lqlq_row = dff[dff['Dataset'] == 'LQLQToTopMu']
data_row = dff[dff['Dataset'] == 'Data']
total_mc_row = dff[dff['Dataset'] == 'Total MC']
data_total_mc_ratio_row = dff[dff['Dataset'] == 'Data/Total MC']

# Remove these specific rows from the DataFrame to sort the rest
dff = dff[~dff['Dataset'].isin(['LQLQToTopMu', 'Data', 'Total MC', 'Data/Total MC'])]

# Sort the remaining DataFrame by 'Mll' in descending order
#dff = dff.sort_values(by='Mll', ascending=False)

# Reassemble the DataFrame in the desired order
dff = pd.concat([lqlq_row, dff, data_row, total_mc_row, data_total_mc_ratio_row], ignore_index=True)

# Print or display the updated DataFrame
print(dff.to_string(index=False))


latex_code = dff.to_latex(index=False, 
                         longtable=False,
                         float_format="{:0.1f}".format)  # Formatting floating point numbers

# Assuming df has 3 columns. Adjust '|l|l|l|' as per your DataFrame's columns.
column_format = '|l|l|l|' 

latex_document = f"""
\\documentclass[10pt]{{article}}
\\usepackage[utf8]{{inputenc}}
\\usepackage{{geometry}}
\\geometry{{a4paper}}
\\usepackage{{pdflscape}}
\\usepackage{{booktabs}}

\\begin{{document}}
\\begin{{landscape}}
\\small
\\begin{{center}}
\\begin{{tabular}}{{{column_format}}}

{latex_code}
\\end{{tabular}}
\\end{{center}}
\\end{{landscape}}
\\end{{document}}
"""

# Save the corrected LaTeX code
with open("Latex/Cutflow_Mu_For_Nosel_HLT.tex", "w") as file:
    file.write(latex_document)




