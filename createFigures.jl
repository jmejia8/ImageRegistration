using PyPlot

function plotScatter(x, y, z, w)
    subplot(2,2,1)
    plot(x[:,1], x[:,2], marker=:o, lw = 0, markersize=4, color=:b)
    title("A")
    grid("on")
    
    subplot(2,2,2)
    plot(y[:,1], y[:,2], marker=:o, lw = 0, markersize=4, color=:r)
    title("B")
    grid("on")

    subplot(2,2, 3)
    title("ECA-MP")
    plot(x[:,1], x[:,2], marker=:o, lw = 0, markersize=4, color=:b)
    plot(z[:,1], z[:,2], marker=:o, lw = 0, markersize=4, color=:r)
    grid("on")

    subplot(2,2, 4)
    title("MSAC")
    plot(x[:,1], x[:,2], marker=:o, lw = 0, markersize=4, color=:b)
    plot(w[:,1], w[:,2], marker=:o, lw = 0, markersize=4, color=:r)
    grid("on")

end