#===========================
Bradon R. McDonald
brmcdonald@wisc.edu
Murtaza Lab
University of Wisconsin-Madison
2024
==========================#

function natureAxis!(ax)
	ax.xgridvisible = false
	ax.ygridvisible = false
	ax.xticklabelsize = 6
	ax.yticklabelsize = 6
	ax.xlabelsize = 8
	ax.ylabelsize = 8
	ax.titlesize = 9
	ax.titlefont = :regular
	ax.spinewidth = 0.75
	ax.xtickwidth = 0.25
	ax.ytickwidth = 0.25
	ax.xticksize = 2.5
	ax.yticksize = 2.5
	ax.xminortickwidth = 0.15
	ax.yminortickwidth = 0.15
	ax.xminorticksize = 1.5
	ax.yminorticksize = 1.5
	hidespines!(ax, :t, :r)
	return
end

function natureMarkers!(l)
	l.strokewidth = 0.25
	l.markersize = 5
	return
end

function natureBoxplot!(l)
	l.whiskerlinewidth = 0.75
	l.medianlinewidth = 0.75
	l.outlierstrokewidth = 0.25
	l.xticklabelsize=7
	return
end

function natureColors()
	COL = Dict()
	COL["red"] = colorant"#db2929"
	COL["blue"] = colorant"#3b76d6"
	COL["purple"] = colorant"#9438e0"
	COL["teal"] = colorant"#2ba8a2"
	COL["orange"] = colorant"#ffa500"
	return COL
end

function scatterPlot(dataX,dataY;w::Real=3.65,h::Real=3.65,title="",
	xlab="",ylab="",logX=false,logY=false,ptSize=4,tickSize=9,titleSize=12)

	fig1 = Figure(size = (w*72,h*72));
	ax, l1 = scatter(fig1[1,1],dataX, dataY,markersize=ptSize);
	ax.yminorticksvisible = true
	ax.xlabel = xlab
	ax.ylabel = ylab
	ax.yticklabelsize = tickSize
	ax.xticklabelsize = tickSize
	ax.title= title
	ax.titlesize=titleSize

	natureAxis!(ax)
	natureMarkers!(l1)

	if logX
		ax.xscale = log10
	end
	if logY
		ax.yscale = log10
	end

#	if ymax == -1
#		ymax = maximum(dataY)+(maximum(dataY)*1.05)
	#end
	#ax.limits = (xmin,xmax,ymin,ymax)

	return fig1
end


function scatterPlot(dataX::Vector{Vector{N}},dataY::Vector{Vector{N}};w::Real=3.65,h::Real=3.65,title="",
	xlab="",ylab="",logX=false,logY=false,ptSize=4,tickSize=9,titleSize=12) where N <: Real
	size_inches = (w,h)
	size_pt = 72 .* size_inches

	colors = ["blue","red","teal"]
	fig1 = Figure(size = (w*72,h*72));
	ax, l1 = scatter(fig1[1,1],dataX[1], dataY[1],markersize=ptSize,color=colors[1]);

	for idx in 2:length(dataX)
		scatter!(dataX[idx], dataY[idx],markersize=ptSize,color=colors[idx]);
	end
	ax.yminorticksvisible = true
	ax.xlabel = xlab
	ax.ylabel = ylab
	ax.yticklabelsize = tickSize
	ax.xticklabelsize = tickSize
	ax.title= title
	ax.titlesize=titleSize
	hidespines!(ax, :t, :r)

	natureAxis!(ax)
	natureMarkers!(l1)

	if logX
		ax.xscale = log10
	end
	if logY
		ax.yscale = log10
	end

#	if ymax == -1
#		ymax = maximum(dataY)+(maximum(dataY)*1.05)
	#end
	#ax.limits = (xmin,xmax,ymin,ymax)

	return fig1
end


#########################
#	!Boxplot
#########################

function microFractionBoxplot2(D,xLabels;logY=false,title="")
	if length(xLabels) != length(D)
		print("Labels dont match input data\n")
		return
	end

	COL = natureColors()

	xData = Float32[]
	yData = Float32[]
	for i in 1:length(xLabels)
		for val in D[i]
			push!(xData,i)
			if val == 0 && logY == true
				push!(yData,log10(1e-16))
				#push!(yData,1e-16)
			else
				push!(yData,log10(val))
				#push!(yData,val)
			end
		end
	end

	w,h = (3.65,3.65)

	fig1 = Figure(size = (w*72,h*72));
	ax, l1 = boxplot(fig1[1,1],xData,yData,color=(COL["blue"],0.75));
	natureAxis!(ax)
	natureBoxplot!(l1)
	natureMarkers!(l1)
	ax.xticks = (collect(1:length(xLabels)),xLabels)
	#ax.yticks = [1e-4,1e-3,1e-2]
	ax.xlabel = "Library Prep Type"
	ax.ylabel = "Total Microbial Fraction"
	ax.title = title
	hidespines!(ax, :t, :r)
	if logY
		ax.yscale = log10
	end
	return fig1
end





#########################
#	!Longitudinal plots
#########################



function plotPatientLonADM(dfCountA,dfADM,preps,outF;yaxis="pgml")
	COL = Dict()
	COL["red"] = colorant"#db2929"
	COL["blue"] = colorant"#3b76d6"
	COL["purple"] = colorant"#9438e0"
	COL["teal"] = colorant"#2ba8a2"
	COL["orange"] = colorant"#b76d2c"
	colorSet = ["blue","red","purple","teal"]

	prepNames = [Symbol(x) for x in preps]
	dfSig = copy(dfADM)[!,prepNames]
	sigNames = [Symbol(x*"_ADM") for x in preps]
	rename!(dfSig,sigNames)
	dfCount = copy(dfCountA)
	dfCount = hcat(dfCount,dfSig)
	allPats = sort(unique(dfCount.Patient))

	w,h = (2.75 * length(allPats),1.5*length(preps))
	fig4 = Figure(size = (w*72,h*72));
	colIdx = 1
	rowIdx = 0

	for dfPat in groupby(dfCount,:Patient)
		rowIdx += 1
		colIdx = 0

		sort!(dfPat,:Timepoint)
		patN = dfPat.Patient[1]
		
		for (prepIdx,prepType) in enumerate(preps)
			minVal = 1
			maxVal = 0
			colIdx += 1
			ax = Axis(fig4[colIdx,rowIdx])
			natureAxis!(ax)
			ax.yticklabelrotation=45
			ax.xticks=(1:10)
			ax.title = "Patient $patN - $prepType"
			ax.xlabel= "Timepoint"
			if yaxis == "pgml"
				ax.ylabel = "Pathogen DNA (pg/mL)"
			elseif yaxis == "microprop"
				ax.ylabel = "Pathogen / Microbial DNA"
			end

			allPathData = Vector{Vector{Point2f}}()
			allColors = Vector{Vector{Any}}()
			allYvals = Vector{Real}()
			for (idx,dfPath) in enumerate(groupby(dfPat,:Pathogen))
				pathN = dfPath.Pathogen[1]
				vals = Vector{Point2f}()
				colors = Any[]
				pathCol= COL[colorSet[idx]]
				for r in eachrow(dfPath)
					push!(vals,(r.Timepoint+1,r[prepType]))
					adm = r[Symbol(prepType*"_ADM")]
					if isnan(adm) == false && adm >= 2.5 && r[prepType] > 0
						push!(colors,pathCol)
					elseif r[prepType] > 0
						push!(colors,"lightgrey")
					else
						push!(colors,"white")
					end
				end
				push!(allYvals,[r[prepType] for r in eachrow(dfPath)]...)
				push!(allPathData,vals)
				push!(allColors,colors)
			end

			if minimum(allYvals) < minVal
				minVal = minimum(allYvals)
			end
			if maximum(allYvals) > maxVal
				maxVal = maximum(allYvals)
			end
			if maxVal == 0
				maxVal = 0.01
			end

			for (idx,pathData) in enumerate(allPathData)
				plt = lines!(pathData,color=COL[colorSet[idx]],linewidth=1)
				plt = scatter!(pathData,markersize=5,color=allColors[idx],strokewidth = 0.5,strokecolor=:black)
			end
			ax.limits = (0,11,minVal-(maxVal / 30),maxVal+(maxVal/20))
		end
	end

	rowgap!(fig4.layout, 10)
	colgap!(fig4.layout, 10)
	save("$(outF).pdf",fig4, pt_per_unit = 1)
end


