using DataFrames, Dates, PlotlyJS, ExcelReaders, ColorSchemes

""" Define styles """
def_style = let
	axis = attr(showgrid = true, gridcolor="#e2e2e2", gridwidth=0.5, zeroline=false)
	layout = Layout(xaxis = axis, yaxis=axis)
	Style(layout=layout)
end

slides_def_noleg = let
	layout = Layout(plot_bgcolor="#fafafa", paper_bgcolor="#fafafa",
		width=1920*0.45, height=1080*0.45, font_size=16, font_family="Lato",
		)
	Style(def_style, layout=layout)
end

slides_def = Style(slides_def_noleg, layout=Layout(legend = attr(orientation = "h", x=0.05)))

dark_bg = let
	axis = attr(gridcolor="#1b1b1b")
	layout = Layout(plot_bgcolor="#020202", paper_bgcolor="#020202", font_color="white", xaxis=axis,yaxis=axis)
	Style(layout=layout)
end
slides_dark = Style(slides_def, dark_bg)

paper = let
	layout = Layout(width = 1920 * 0.5, height = 1080 * 0.35, font_size=16, font_family = "Linux Libertine",
		legend = attr(orientation = "h", x=0.05))
	Style(def_style, layout=layout)
end


plot_borrowing(;style::Style=slides_def) = plot_data([3,8,9,10,11,12], "Borrowing rates", style=style)
plot_deposit(;style::Style=slides_def) = plot_data([3,4,5,6,7], "Deposit rates", style=style)

function plot_data(varvec, title; style::Style=slides_def_noleg)
	dataraw = readxlsheet("../Data/IMFIFS/Interest_Rates.xls", "Interest Rates")

	varlabels = convert(Vector{String}, dataraw[8:end,1])
	varnames = convert(Vector{String}, dataraw[8:end,3])
	vals = dataraw[8:end,5:end-1]
	vals[typeof.(vals).==String].=NaN

	vals = convert(Matrix{Float64}, vals')

	dates = Date.(replace.(dataraw[7,5:end-1], "M"=>"-"), "yyyy-mm")

	df = DataFrame(vals, varnames)
	df.date = dates

	layout = Layout(height = 0.5 * 1080, width=0.5*1920, yaxis_title="%", title=title)

	totcol = length(varvec)

	data = [
	scatter(x=df.date, y=df[!,names(df)[varvec[1]]], name=varlabels[varvec[1]], line_width = 3)
	[scatter(x=df.date, y=df[!,names(df)[jj]], name=varlabels[jj], line_color=get(ColorSchemes.lajolla, (jn+1)/totcol)) for (jn,jj) in enumerate(varvec[2:end])]
	]

	plot(data, layout, style=style)
end