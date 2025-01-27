#===========================
Bradon R. McDonald
brmcdonald@wisc.edu
Kisat Lab
University of Wisconsin-Madison
2025
===========================#

inDir = "./brakenCountData_DominguezMcDonald2025/SampleCountsPos"
inDirNeg = "./brakenCountData_DominguezMcDonald2025/SampleCountsNeg"
classStats = "./totalFragmentTable_DominguezMcDonald2025.tsv"
pathYaml = "./pathogenData_DominguezMcDonald2025.yaml"
pathTbl = "./cultures_cohortB_DominguezMcDonald2025.tsv"
cfdnaCon = "./cfDNA_concentration_DominguezMcDonald2025.tsv"

outDir = "./analysisOutput_DominguezMcDonald2025"
dataName = "DominguezMcDonald2025"

include("./anaylsisCode_DominguezMcDonald2025_functions.jl")
include("./anaylsisCode_DominguezMcDonald2025_plots.jl")
pathogensDict = YAML.load(open(pathYaml))
dfPath = CSV.File(pathTbl) |> DataFrame
PATH = Dict{String,String}()
for r in eachrow(dfPath)
	PATH[r.Patient] = r.Pathogen
end

#===========================
Running analysis
===========================#

run(`mkdir -p $outDir`)

D01 = parseTaxonData(inDir)
annotateDepth!(D01,classStats)
annotateUnmapped!(D01,classStats)
annotateInput!(D01,cfdnaCon)
dfBKA_mp,dfBKB_mp,BK_A_mp,BK_B_mp = backgroundGenus(D01,pathogensDict,metric="microprop")

dfMaster = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict)
dfMasterPropMicro = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="propMicro")
dfMasterDepth = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="depth")
dfMasterCount = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="count")
dfMasterTotal = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="total")

dfMasterADM_mp = buildKnownTable(D01,["B251","B266","B268","B297","B304"],pathogensDict,metric="adm_mp",BK=BK_B_mp)

CSV.write("$(outDir)/positivePatientTable_$(dataName).tsv",dfMaster,delim='\t')
CSV.write("$(outDir)/positivePatientTable_$(dataName)_depth.tsv",dfMasterDepth,delim='\t')
CSV.write("$(outDir)/positivePatientTable_$(dataName)_counts.tsv",dfMasterCount,delim='\t')
CSV.write("$(outDir)/positivePatientTable_$(dataName)_totalMicrobial.tsv",dfMasterTotal,delim='\t')
CSV.write("$(outDir)/positivePatientTable_$(dataName)_ADM.tsv",dfMasterADM_mp,delim='\t')

dfTM_total = microFractionSummary(D01)
dfTM_unmapped = microFractionSummary(D01,denom="unmapped")

MI_01_SS = microFractionComparison(D01,prepType="SS")
MI_01_EMS = microFractionComparison(D01,prepType="EMS")
MI_01_STD = microFractionComparison(D01,prepType="STD")
MI_01_SSS = microFractionComparison(D01,prepType="SS_EMS")

#Boxplots of total microbial fraction
tmData_STD = vcat([MI_01_STD["mfA"],MI_01_STD["mfBPos"],MI_01_STD["mfBneg"]]...)
tmData_EMS = vcat([MI_01_EMS["mfA"],MI_01_EMS["mfBPos"],MI_01_EMS["mfBneg"]]...)
tmData_SS = vcat([MI_01_SS["mfA"],MI_01_SS["mfBPos"],MI_01_SS["mfBneg"]]...)
tmData_SSS = vcat([MI_01_SSS["mfA"],MI_01_SSS["mfBPos"],MI_01_SSS["mfBneg"]]...)

prepLabels = ["dsDNA","dsDNA\nSize selected","ssDNA","ssDNA\nSize selected"]
p0 = microFractionBoxplot2([tmData_STD,tmData_EMS,tmData_SS,tmData_SSS],prepLabels;title="Microbial Fragment Fraction");
save("$(outDir)/microbialFraction_$(dataName).pdf",p0)

#Comparing total microbial fraction vs cfDNA yield
YVT = Dict()
YVT["STD"] = [vcat([MI_01_STD["inputA"],MI_01_STD["inputBPos"],MI_01_STD["inputBneg"]]...),vcat([MI_01_STD["mfA"],MI_01_STD["mfBPos"],MI_01_STD["mfBneg"]]...)]
YVT["EMS"] = [vcat([MI_01_EMS["inputA"],MI_01_EMS["inputBPos"],MI_01_EMS["inputBneg"]]...),vcat([MI_01_EMS["mfA"],MI_01_EMS["mfBPos"],MI_01_EMS["mfBneg"]]...)]
YVT["SS"] = [vcat([MI_01_SS["inputA"],MI_01_SS["inputBPos"],MI_01_SS["inputBneg"]]...),vcat([MI_01_SS["mfA"],MI_01_SS["mfBPos"],MI_01_SS["mfBneg"]]...)]
YVT["SS_EMS"] = [vcat([MI_01_SSS["inputA"],MI_01_SSS["inputBPos"],MI_01_SSS["inputBneg"]]...),vcat([MI_01_SSS["mfA"],MI_01_SSS["mfBPos"],MI_01_SSS["mfBneg"]]...)]

p1 = scatterPlot(YVT["STD"][1],YVT["STD"][2],
	logY=true,logX=true,w=3.5,h=3.5,ptSize=8,tickSize=14,title="cfDNA vs Total Microbial Takara");
p2 = scatterPlot(YVT["EMS"][1],YVT["EMS"][2],
	logY=true,logX=true,w=3.5,h=3.5,ptSize=8,tickSize=14,title="cfDNA vs Total Microbial Takara-Size Selected");
p3 = scatterPlot(YVT["SS"][1],YVT["SS"][2],
	logY=true,logX=true,w=3.5,h=3.5,ptSize=8,tickSize=14,title="cfDNA vs Total Microbial SRSLY");
p4 = scatterPlot(YVT["SS_EMS"][1],YVT["SS_EMS"][2],
	logY=true,logX=true,w=3.5,h=3.5,ptSize=8,tickSize=14,title="cfDNA vs Total Microbial SRSLY-SS");

save("$(outDir)/cfDNA_vs_microbialFraction_std_$(dataName).pdf",p1)
save("$(outDir)/cfDNA_vs_microbialFraction_ems_$(dataName).pdf",p2)
save("$(outDir)/cfDNA_vs_microbialFraction_ss_$(dataName).pdf",p3)
save("$(outDir)/cfDNA_vs_microbialFraction_ss-ss_$(dataName).pdf",p4)

#===========================
Patient Longitudinal plots
===========================#
pats = ["B251","B266","B268","B297","B304"]
preps = ["STD","EMS","SS","SS_EMS"]
pathogens = ["Streptococcus","Staphylococcus","Klebsiella","Haemophilus","Staphylococcus"]
colorSet = ["blue","red","purple","teal"]

CTP = Dict{String,Vector{Int}}()
CTP["B251"] = [2,9]
CTP["B266"] = [1]
CTP["B268"] = [6]
CTP["B297"] = [7]
CTP["B304"] = [2]

CTN = Dict{String,Vector{Int}}()
CTN["B251"] = [1,3,8]
CTN["B266"] = [5]
CTN["B268"] = [1,5,10]
CTN["B297"] = [6,8]
CTN["B304"] = [1]

outF = "./$(outDir)/posPatients_Lon_PathYield_sig_microProp"
plotPatientLonADM(dfMasterPropMicro,dfMasterADM_mp,["STD","EMS","SS","SS_EMS"],outF,yaxis="microprop")

#===========================
Negative culture patient analysis
===========================#

DN = parseTaxonData(inDirNeg)
annotateDepth!(DN,classStats)
annotateUnmapped!(DN,classStats)

MIN_SS = microFractionComparison(DN,prepType="SS")
MIN_EMS = microFractionComparison(DN,prepType="EMS")
MIN_STD = microFractionComparison(DN,prepType="STD")
MIN_SSS = microFractionComparison(DN,prepType="SS_EMS")

prepLabels = ["dsDNA\nSize selected","ssDNA\nSize selected"]
p0 = microFractionBoxplot2([MIN_EMS["mfBneg"],MIN_SSS["mfBneg"]],prepLabels;title="Microbial Fragment Fraction");
save("$(outDir)/microbialFraction_$(dataName)_negativePatients.pdf",p0)

#===========================
Stats
===========================#

dfM = sort(CSV.File(classStats) |> DataFrame)
dfMP = filter(:PatientType=>x->x == "PosCulture",dfM)
dsDNA = filter(:PrepType=>x->x == "dsDNA",dfMP)
dsDNAss = filter(:PrepType=>x->x == "dsDNA-SS",dfMP)
ssDNA = filter(:PrepType=>x->x == "ssDNA",dfMP)
ssDNAss = filter(:PrepType=>x->x == "ssDNA-SS",dfMP)

#mean unmapped fraction in dsDNA vs ssDNA, enrichment fold
unmapProp_dsDNA =[r.GoodUnmapped/r.GoodFragments for r in eachrow(dsDNA)]
unmapProp_ssDNA = [r.GoodUnmapped/r.GoodFragments for r in eachrow(ssDNA)]
unmapProp_dsDNAss = [r.GoodUnmapped/r.GoodFragments for r in eachrow(dsDNAss)]
unmapProp_ssDNAss = [r.GoodUnmapped/r.GoodFragments for r in eachrow(ssDNAss)]
enrich_ssDNA = [(unmapProp_ssDNA[i]-unmapProp_dsDNA[i])/unmapProp_dsDNA[i] for i in 1:length(unmapProp_dsDNA)]
enrich_ssDNAss = [(unmapProp_ssDNAss[i]-unmapProp_dsDNA[i])/unmapProp_dsDNA[i] for i in 1:length(unmapProp_dsDNA)]

#Mean micro fragments in all preps, fold enrichment from dsDNA to ssDNA-SS
microProp_dsDNA = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(dsDNA)]
microProp_ssDNA = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(ssDNA)]
microProp_dsDNAss = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(dsDNAss)]
microProp_ssDNAss = [sum(values(D01[x.SampleName].phylumData)) / D01[x.SampleName].totalFrags for x in eachrow(ssDNAss)]

enrich_microProp_ds_ss = [(microProp_ssDNA[i]-microProp_dsDNA[i])/microProp_dsDNA[i] for i in 1:length(microProp_dsDNA)]
enrich_microProp_ds_dsss = [(microProp_dsDNAss[i]-microProp_dsDNA[i])/microProp_dsDNA[i] for i in 1:length(microProp_dsDNA)]
enrich_microProp_ss_ssss = [(microProp_ssDNAss[i]-microProp_ssDNA[i])/microProp_ssDNA[i] for i in 1:length(microProp_ssDNA)]
enrich_microProp_ds_ssss = [(microProp_ssDNAss[i]-microProp_dsDNA[i])/microProp_dsDNA[i] for i in 1:length(microProp_dsDNA)]

dfMN = filter(:PatientType=>x->x == "NegCulture",dfM)
neg_dsDNAss = filter(:PrepType=>x->x == "dsDNA-SS",dfMN)
neg_ssDNAss = filter(:PrepType=>x->x == "ssDNA-SS",dfMN)

microPropN_dsDNAss = [sum(values(DN[x.SampleName].phylumData)) / DN[x.SampleName].totalFrags for x in eachrow(neg_dsDNAss)]
microPropN_ssDNAss = [sum(values(DN[x.SampleName].phylumData)) / DN[x.SampleName].totalFrags for x in eachrow(neg_ssDNAss)]

mean(microPropN_dsDNAss)
mean(microPropN_ssDNAss)



