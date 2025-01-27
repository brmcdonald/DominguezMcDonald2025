#===========================
Bradon R. McDonald
brmcdonald@wisc.edu
Murtaza Lab
University of Wisconsin-Madison
2024
==========================#

using Pkg
Pkg.activate("./env_DominguezMcDonald2025")
using CSV,DataFrames,YAML
using Statistics,HypothesisTests,StatsBase
using CodecZlib
using DataStructures
using CairoMakie,Colors


mutable struct MRMSample
	name::String
	patient::String
	tp::Int
	prepType::String
	cohort::String
	totalFrags::Real
	unmappedFrags::Real
	microbialFrags::Real
	phylumData::DefaultDict{String, Real, Int64}
	genusData::DefaultDict{String, Real, Int64}
	speciesData::DefaultDict{String, Real, Int64}
	cfdna::Real
end

#From Kraken output, exclude human
function getkrakenProp(inDir::String)
	files = filter(x->occursin(".k2.out.gz",x),readdir(inDir))
	sL = String[]
	mL = Int[]
	tL = Int[]
	for f in files
		n = replace(f,".k2.out.gz"=>"")
		df =CSV.File(joinpath(inDir,f),header=false) |> DataFrame
		dfC = filter(:Column3=>x->x!=9606,df)
		filter!(:Column1=>x->x=="C",dfC)
		totalUnmapped = size(df,1)
		totalClass = size(dfC,1)

		push!(sL,n)
		push!(mL,totalClass)
		push!(tL,totalUnmapped)
	end
	return DataFrame(SampleName=sL,GoodClassified=mL,TotalFragments=tL)
end


function parseTaxonData(inDir::String)
	D = Dict{String,MRMSample}()
	genusData = readdir("$inDir/G_taxonomy")
	filter!(x->occursin("DS_Store",x)==false,genusData)

	for f in genusData
		getName = r"TRS22\-([^-_]+)"
		getTPa = r"_T(\d+)"
		getTPb = r"-T(\d+)"

		sampleName = replace(f,".taxoncounts.G.tsv"=>"")

		patient = match(getName,sampleName)[1]


		tp1 = match(getTPa,f)
		tp2 = match(getTPb,f)
		if isnothing(tp1) == false
			tp = parse(Int,match(getTPa,f)[1])
		elseif isnothing(tp2) == false
			tp = parse(Int,match(getTPb,f)[1])
		else
			tp = 0
		end

		if occursin("-SS-EMS",f) || occursin("-SRS-BP",f)
			prepType = "SS_EMS"
		elseif occursin("-SS",f) || occursin("-SRS",f)
			prepType = "SS"
		elseif occursin("-EMS",f) || occursin("-TAK-MDBP",f)
			prepType = "EMS"
		else
			prepType = "STD"
		end

		if startswith(patient,"B")
			cohort = "B"
		elseif startswith(patient,"A")
			cohort = "A"
		elseif startswith(patient,"NCE")
			cohort = "ctrl"
		else
			cohort = "NA"
		end
		TXG = DefaultDict{String,Real}(0)
		df = CSV.File(joinpath(inDir,"G_taxonomy",f)) |> DataFrame
		for r in eachrow(df)
			TXG[r.Taxon] = r.Fragments
		end

		TXP = DefaultDict{String,Real}(0)
		pfilename = replace(f,"G.tsv"=>"P.tsv")
		df = CSV.File(joinpath(inDir,"P_taxonomy",pfilename)) |> DataFrame
		for r in eachrow(df)
			TXP[r.Taxon] = r.Fragments
		end

		TXS = DefaultDict{String,Real}(0)
		pfilename = replace(f,"G.tsv"=>"S.tsv")
		try
			df = CSV.File(joinpath(inDir,"S_taxonomy",pfilename)) |> DataFrame
			for r in eachrow(df)
				TXS[r.Taxon] = r.Fragments
			end
		catch
			x = 1
		end

		g = MRMSample(sampleName,patient,tp,prepType,cohort,0,0,0,TXP,TXG,TXS,-1)
		D[sampleName] = g
	end
	return D
end

function annotateDepth!(D::Dict{String,MRMSample},tbl::String;datafield="goodfragments")
	df = CSV.File(tbl) |> DataFrame
	T = Dict{String,Real}()
	for r in eachrow(df)
		if datafield == "totalfragments"
			T[r.SampleName] = r.TotalFragments
		elseif datafield == "goodfragments"
			T[r.SampleName] = r.GoodFragments
		end
	end
	for i in keys(D)
		if haskey(T,i) == false
			continue
		end
		D[i].totalFrags = T[i]
	end
end

function annotateUnmapped!(D::Dict{String,MRMSample},tbl::String;denom="goodunmapped")
	df = CSV.File(tbl) |> DataFrame
	T = Dict{String,Real}()
	M = Dict{String,Real}()
	for r in eachrow(df)
		if denom == "goodunmapped"
			T[r.SampleName] = r.GoodUnmapped
		else
			T[r.SampleName] = r.TotalUnmapped
		end
		M[r.SampleName] = r.GoodClassified
	end
	for i in keys(D)
		D[i].unmappedFrags = T[i]
		D[i].microbialFrags = M[i]
	end
end

function annotateInput!(D::Dict{String,MRMSample},tbl::String)
	df = CSV.File(tbl) |> DataFrame
	I = Dict{String,Real}()
	for r in eachrow(df)
		I[r.Sample] = r.cfDNA
	end

	for s in keys(D)
		for i in keys(I)
			if occursin(i,s)
				if D[s].cfdna == -1
					D[s].cfdna = I[i]
				else
					print("Sample found already: $i\n")
				end
			end
		end
	end
end

function sortDict(D::Dict{X,Y};rev::Bool=false) where {X<:Any,Y<:Real}
	return sort(collect(D), by=x->x[2],rev=rev)
end

function buildKnownTable(D,pats,PATH;metric="prop",preps=["STD","EMS","SS","SS_EMS"],BK=nothing)
	PROP = Dict{String,Vector{Real}}()
	for prep in preps
		PROP[prep] = Real[]
	end
	patL = String[]
	tpL = Int[]
	pathL = String[]

	for (idx,p) in enumerate(pats)
		for path in PATH[p]
			samples = [x for x in keys(D) if D[x].patient == p]
			allTP = unique([D[x].tp for x in samples])
			for tp in allTP
				tpS = [x for x in samples if D[x].tp == tp]
				push!(patL,p)
				push!(tpL,tp)
				push!(pathL,path)

				for prep in preps
					s = [x for x in tpS if D[x].prepType == prep]
					if length(s) > 0 && metric == "prop"
						push!(PROP[prep],D[s[1]].genusData[path] / D[s[1]].totalFrags)
					elseif length(s) > 0 && metric == "propMicro"
						push!(PROP[prep],D[s[1]].genusData[path] / sum(values(D[s[1]].phylumData)))
					elseif length(s) > 0 && metric == "count"
						push!(PROP[prep],D[s[1]].genusData[path])
					elseif length(s) > 0 && metric == "depth"
						push!(PROP[prep],D[s[1]].totalFrags)
					elseif length(s) > 0 && metric == "total"
						push!(PROP[prep],sum(values(D[s[1]].phylumData)) / D[s[1]].totalFrags)
					elseif length(s) > 0 && metric == "adm"
						prop = D[s[1]].genusData[path] / D[s[1]].totalFrags * D[s[1]].cfdna * 1000
						adm = round(getADM(prop,BK[prep][path]),digits=3)
						push!(PROP[prep],adm)
					elseif length(s) > 0 && metric == "adm_mp"
						prop = D[s[1]].genusData[path] / sum(values(D[s[1]].phylumData))
						adm = round(getADM(prop,BK[prep][path]),digits=3)
						push!(PROP[prep],adm)
					elseif length(s) > 0 && metric == "pgml"
						val = D[s[1]].genusData[path] / D[s[1]].totalFrags * D[s[1]].cfdna * 1000
						push!(PROP[prep],val)
					else
						push!(PROP[prep],0.0)
					end
				end #For each prep

			end #For each timepoint
		end #For each pathogen
	end #For each patient

	df = DataFrame(Patient=patL,Timepoint=tpL,Pathogen=pathL)
	for p in preps
		df[!,Symbol(p)] .= PROP[p]
	end
	sort!(df,[:Patient,:Timepoint])
	return df
end


##################################################################################################

#Comparing total microbial fraction between infection/no infection across preps
function microFractionComparison(D;denom="total",prepType="STD")
	DO = Dict{String,Vector{Float64}}()
	sepPats = Set(["B268","B251","B266","B297"])
	cohortBClean = Set(["B287","B289","B296","B306","B309"])

	ssA = filter(x->D[x].cohort=="A" && D[x].prepType==prepType && D[x].cfdna >0,keys(D))
	ssB = filter(x->D[x].patient in sepPats && D[x].prepType==prepType && D[x].cfdna >0,keys(D))
	ssC = filter(x->D[x].patient in cohortBClean && D[x].prepType==prepType && D[x].cfdna >0,keys(D))
	if denom == "total"
		totalA = [sum(values(D[x].phylumData)) / D[x].totalFrags for x in ssA]
		totalB = [sum(values(D[x].phylumData)) / D[x].totalFrags for x in ssB]
		totalC = [sum(values(D[x].phylumData)) / D[x].totalFrags for x in ssC]
	elseif denom == "unmapped"
		totalA = [sum(values(D[x].phylumData)) / D[x].unmappedFrags for x in ssA]
		totalB = [sum(values(D[x].phylumData)) / D[x].unmappedFrags for x in ssB]
		totalC = [sum(values(D[x].phylumData)) / D[x].unmappedFrags for x in ssC]
	end
	inputA = [D[x].cfdna for x in ssA]
	inputB = [D[x].cfdna for x in ssB]
	inputC = [D[x].cfdna for x in ssC]

	DO["mfA"] = totalA
	DO["inputA"] = inputA
	DO["mfBPos"] = totalB
	DO["inputBPos"] = inputB
	DO["mfBneg"] = totalC
	DO["inputBneg"] = inputC
	return DO
end


function microFractionSummary(D;denom="total",prepTypes=["STD","EMS","SS","SS_EMS"])
	sepPats = Set(["B268","B251","B266","B297"])
	cohortBClean = Set(["B287","B289","B296","B306","B309"])
	pL = String[]
	aL = Real[]
	bL = Real[]
	cL = Real[]
	pvalL = Real[]
	for p in prepTypes
		push!(pL,p)
		ssA = filter(x->D[x].cohort=="A" && D[x].prepType==p,keys(D))
		ssB = filter(x->D[x].patient in sepPats && D[x].prepType==p,keys(D))
		ssC = filter(x->D[x].patient in cohortBClean && D[x].prepType==p,keys(D))
		if denom == "total"
			totalA = [D[x].microbialFrags / D[x].totalFrags for x in ssA]
			totalB = [D[x].microbialFrags / D[x].totalFrags for x in ssB]
			totalC = [D[x].microbialFrags / D[x].totalFrags for x in ssC]
		elseif denom == "unmapped"
			totalA = [D[x].microbialFrags / D[x].unmappedFrags for x in ssA]
			totalB = [D[x].microbialFrags / D[x].unmappedFrags for x in ssB]
			totalC = [D[x].microbialFrags / D[x].unmappedFrags for x in ssC]
		end

		if length(totalA) > 0 
			push!(aL,median(totalA))
		else
			push!(aL,-1.0)
		end
		if length(totalB) > 0 
			push!(bL,median(totalB))
		else
			push!(bL,-1.0)
		end
		if length(totalC) > 0
			push!(cL,median(totalC))
		else
			push!(cL,-1.0)
		end
		if length(totalA) > 0 && length(totalB) > 0
			push!(pvalL,pvalue(MannWhitneyUTest(totalA,totalB)))
		else
			push!(pvalL,1.0)
		end
	end
	dfTM = DataFrame(Prep=pL,CohortA=aL,CohortBPos=bL,CohortBNeg=cL,Pval=pvalL)
	return dfTM
end


#Background levels of each pathogen in non-sepsis patients
function backgroundGenus(D,PATH;prepTypes=["STD","EMS","SS","SS_EMS"],metric="pgml")
	POS = Dict{String,Dict{String,Vector{Real}}}()
	BND_A = Dict{String,Dict{String,Vector{Real}}}()
	BND_B = Dict{String,Dict{String,Vector{Real}}}()
	bPats = [x for x in keys(D) if haskey(PATH,D[x].patient)]
	allPath = collect(values(PATH))
	allPath = unique(vcat(allPath...))

	for prep in prepTypes
		POS[prep] = Dict{String,Vector{Real}}()
		BND_A[prep] = Dict{String,Vector{Real}}()
		BND_B[prep] = Dict{String,Vector{Real}}()
		for path in allPath

			ssA = filter(x-> D[x].cohort=="A" && D[x].prepType==prep,keys(D))
			ssB = filter(x-> !(path in PATH[D[x].patient]) && D[x].prepType==prep,bPats)
			ssPos = filter(x-> path in PATH[D[x].patient] && D[x].prepType==prep,bPats)

			if metric == "pgml"
				totalA = [D[x].genusData[path] / D[x].totalFrags * D[x].cfdna * 1000 for x in ssA]
				totalB = [D[x].genusData[path] / D[x].totalFrags * D[x].cfdna * 1000 for x in ssB]
				totalPos = [D[x].genusData[path] / D[x].totalFrags * D[x].cfdna * 1000 for x in ssPos]
			elseif metric == "microprop"
				totalA = [D[x].genusData[path] / sum(values(D[x].phylumData)) for x in ssA]
				totalB = [D[x].genusData[path] / sum(values(D[x].phylumData)) for x in ssB]
				totalPos = [D[x].genusData[path] / sum(values(D[x].phylumData)) for x in ssPos]
			end
			POS[prep][path] = totalPos
			BND_A[prep][path] = totalA
			BND_B[prep][path] = totalB
		end
	end

	prepL = String[]
	pathL = String[]
	posCount = Int[]
	totalCount = Int[]
	medianProp = Float64[]
	stdProp = Float64[]
	for prep in prepTypes
		for path in allPath
			push!(prepL,prep)
			push!(pathL,path)
			pc = [x for x in BND_A[prep][path] if x > 0]
			tc = length(BND_A[prep][path])
			if length(pc) > 0
				mp = median(pc)
				stdev = std(pc)
			else
				mp = -1
				stdev = -1
			end
			push!(posCount,length(pc))
			push!(totalCount,tc)
			push!(medianProp,mp)
			push!(stdProp,stdev)
		end
	end
	dfBNDA = DataFrame(Prep=prepL,Pathogen=pathL,PosSamples=posCount,TotalSamples=totalCount,MedianProp=medianProp,stdev=stdProp)
	dfBNDA[!,:FalsePos] .= [round(r.PosSamples / r.TotalSamples,digits=3) for r in eachrow(dfBNDA)]

	prepL = String[]
	pathL = String[]
	posCount = Int[]
	totalCount = Int[]
	medianProp = Float64[]
	stdProp = Float64[]
	for prep in prepTypes
		for path in allPath
			push!(prepL,prep)
			push!(pathL,path)
			pc = [x for x in BND_B[prep][path] if x > 0]
			tc = length(BND_B[prep][path])
			if length(pc) > 0
				mp = median(pc)
				stdev = std(pc)
			else
				mp = -1
				stdev = -1
			end
			push!(posCount,length(pc))
			push!(totalCount,tc)
			push!(medianProp,mp)
			push!(stdProp,stdev)
		end
	end
	dfBNDB = DataFrame(Prep=prepL,Pathogen=pathL,PosSamples=posCount,TotalSamples=totalCount,MedianProp=medianProp,stdev=stdProp)
	dfBNDB[!,:FalsePos] .= [round(r.PosSamples / r.TotalSamples,digits=3) for r in eachrow(dfBNDB)]

	return (dfBNDA,dfBNDB,BND_A,BND_B,POS)
end


##################################################################################################
#Utils
##################################################################################################

function getPatientGenus(D,pat,genus,prep)
	allS = filter(x->D[x].patient==pat && D[x].prepType==prep,keys(D))
	vals = Float32[]
	for s in allS
		push!(vals,(D[s].genusData[genus] /  D[s].totalFrags) * D[s].cfdna * 1000)
	end
	return vals
end

function getADM(x::Real,dist)
	return (x - median(dist)) / StatsBase.mad(dist)
end

function sigMAD(x::Real,dist;cutoff=2.5)
	if abs(getADM(x,dist)) >= cutoff
		return true
	else
		return false
	end
end
