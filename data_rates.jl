using DataFrames, Dates, PlotlyJS, ExcelReaders, ColorSchemes

""" Define styles """
sand() = "#F5F3F1"
darkbgd() = "#272929"
lightgrid() = "#353535"
darkgrid() = "#e2e2e2"
gridcol(dark=false) = ifelse(dark, lightgrid(), darkgrid())

q_axis(dark) = attr(showgrid=true, gridcolor=gridcol(dark), gridwidth=0.5, zeroline=false)
bgcol(slides, dark) = ifelse(slides, ifelse(dark, darkbgd(), sand()), "white")
qleg() = attr(orientation="h", x=0.05, xanchor="left")

qwidth(slides) = 864
qheight(slides) = ceil(Int, qwidth(slides) * ifelse(slides, 10 / 16, 7 / 16))

function qtemplate(; dark=false, slides=!dark)
    axis = q_axis(dark)
    width = 864 #1920 * 0.45
    l = Layout(
        xaxis=axis, yaxis=axis,
        width=width,
        height=width * ifelse(slides, 10 / 16, 7 / 16),
        font=attr(
            family=ifelse(slides, "Lato", "Linux Libertine"),
            size=16, color=ifelse(dark, sand(), darkbgd())
        ),
        paper_bgcolor=bgcol(slides, dark), plot_bgcolor=bgcol(slides, dark),
        legend=qleg(),
    )
    return Template(layout=l)
end


plot_borrowing(; slides=true, dark=slides, template::Template=qtemplate(slides=slides, dark=dark)) = plot_data([3, 8, 9, 10, 11, 12], "Borrowing rates", template=template)
plot_deposit(; slides=true, dark=slides, template::Template=qtemplate(slides=slides, dark=dark)) = plot_data([3, 4, 5, 6, 7], "Deposit rates", template=template)

function plot_data(varvec, title; slides=true, dark=slides, template::Template=qtemplate(slides=slides, dark=dark))
    dataraw = readxlsheet("../Data/IMFIFS/Interest_Rates.xls", "Interest Rates")

    varlabels = convert(Vector{String}, dataraw[8:end, 1])
    varnames = convert(Vector{String}, dataraw[8:end, 3])
    vals = dataraw[8:end, 5:end-1]
    vals[typeof.(vals).==String] .= NaN

    vals = convert(Matrix{Float64}, vals')

    dates = Date.(replace.(dataraw[7, 5:end-1], "M" => "-"), "yyyy-mm")

    df = DataFrame(vals, varnames)
    df.date = dates

    layout = Layout(template=template, height=0.5 * 1080, width=0.5 * 1920, yaxis_title="%", title=title)

    totcol = length(varvec)

    data = [
        scatter(x=df.date, y=df[!, names(df)[varvec[1]]], name=varlabels[varvec[1]], line_width=3)
        [scatter(x=df.date, y=df[!, names(df)[jj]], name=varlabels[jj], line_color=get(ColorSchemes.lajolla, (jn + 1) / totcol)) for (jn, jj) in enumerate(varvec[2:end])]
    ]

    plot(data, layout)
end