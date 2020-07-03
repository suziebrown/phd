### Create generic ternary diagrams

using(Plots)

### Create a blank ternary plot (corner labels + border)

function ternaryplot(sidelength::Float64=1.0, labels::Array{String,1}=["x","y","z"])
    # utterly blank plot
    myplot = plot(ticks=nothing, border=:none, legend=:none, aspect_ratio=1)
    # draw outline of triangle
    plot!(myplot, [0,sidelength,2*sidelength,0], [0,sqrt(3)*sidelength,0,0], line=1, color=:black)
    # label the corners
    annotate!(myplot, -0.02, 0, Plots.text(labels[3], 12, :right))
    annotate!(myplot, 2.02, 0, Plots.text(labels[1], 12, :left))
    annotate!(myplot, 1, sqrt(3)+0.02, Plots.text(labels[2], 12, :bottom))
    return(myplot)
end

### Convert x-y coordinates to points on ternary plot, and v.v.

function xytoternary(x::Array{Float64,1}, y::Array{Float64,1}, sidelength::Float64=1.0)
    terny = y * sqrt(3) * sidelength
    ternx = (2 * x .+ y) * sidelength
    return(ternx, terny)
end

function ternarytoxy(ternx::Array{Float64,1}, terny::Array{Float64,1}, sidelength::Float64=1.0)
    y = terny / (sqrt(3) * sidelength)
    x = (ternx .- y) / (2*sidelength)
    return(x=x,y=y)
end

### Test adding a parametric curve to the plot

# ternaryplot()
# xvals = Vector(0:0.01:1)
# plot!(xytoternary(xvals, 1 .- xvals))
# plot!(xytoternary(xvals, (1 .- xvals).^3))
# plot!(xytoternary(xvals, abs.(sin.(10 .* xvals)) .* (1 .-xvals)))

### Find centres of grid triamgles

function gridcentres(resolution::Int64, sidelength::Float64=1.0)
    ntriangle = resolution^2
    centresx = Array{Float64, 1}(undef, ntriangle)
    centresy = Array{Float64, 1}(undef, ntriangle)
    upwards = Array{Bool, 1}(undef, ntriangle)
    for row in 1:resolution
        rowlength = 2(resolution - row) +1
        for i in 1:rowlength
            triangleindex = (row - 1)*(2*resolution - row + 1) + i
            centresx[triangleindex] = (row + i - 1) * sidelength/resolution # x plotting coord of centre
            yoffset = iseven(i) ? (sqrt(3)/3) : (2*sqrt(3)/3)
            centresy[triangleindex] = (row*sqrt(3) - yoffset) * sidelength/resolution
            upwards[triangleindex] = !iseven(i)
        end
    end
    return(x=centresx, y=centresy, upwards=upwards)
end

### Shade the grid cells according to a specified function

function ternaryheatmap(func::Function, resolution::Int64, sidelength::Float64=1.0, labels::Array{String,1}=["x","y","z"])
    myplot = ternaryplot(sidelength, labels)
    centres = gridcentres(resolution, sidelength)
    ntriangle = length(centres.x)
    xycentres = ternarytoxy(centres.x, centres.y)
    data = func(xycentres.x, xycentres.y)
    mindata = minimum(data)
    maxdata = maximum(data)
    for i in 1:ntriangle
        if centres.upwards[i]
            cornersx = [centres.x[i] - sidelength/resolution, centres.x[i] + sidelength/resolution, centres.x[i], centres.x[i] - sidelength/resolution]
            cornersy = [centres.y[i] - sqrt(3)*sidelength/(3*resolution), centres.y[i] - sqrt(3)*sidelength/(3*resolution), centres.y[i] + 2*sqrt(3)*sidelength/(3*resolution), centres.y[i] - sqrt(3)*sidelength/(3*resolution)]
        else
            cornersx = [centres.x[i] - sidelength/resolution, centres.x[i] + sidelength/resolution, centres.x[i], centres.x[i] - sidelength/resolution]
            cornersy = [centres.y[i] + sqrt(3)*sidelength/(3*resolution), centres.y[i] + sqrt(3)*sidelength/(3*resolution), centres.y[i] - 2*sqrt(3)*sidelength/(3*resolution), centres.y[i] + sqrt(3)*sidelength/(3*resolution)]
        end
        colourval = (data[i] - mindata)/(maxdata - mindata)
        cellcolour = cgrad(:default)[colourval]
        plot!(myplot, cornersx,cornersy, fill=(0, cellcolour), linecolor=cellcolour)
    end
    display(myplot)
end
