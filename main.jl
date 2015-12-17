import MyPack
p = MyPack.Plot
g = p.Gnuplot()

typealias Pos Tuple{Int, Int}
typealias PosVal Tuple{Pos, Float64}

# Init grid, boundary, stick
function init(g, b, stick)
    const c = round(Int, (size(g)[1])/2)
    const r = c-2
    # seed
    g[c, c] = 0
    b[c, c] = true
    push!(stick, (c, c))
    # DEBUG
#    g[c+1, c] = 0
#    b[c+1, c] = true
#    push!(stick, (c+1, c))
#    g[c+2, c] = 0
#    b[c+2, c] = true
#    push!(stick, (c+2, c))
#    g[c+3, c] = 0
#    b[c+3, c] = true
#    push!(stick, (c+3, c))
#    g[c+3, c+1] = 0
#    b[c+3, c+1] = true
#    push!(stick, (c+3, c+1))
#    g[c+3, c+2] = 0
#    b[c+3, c+2] = true
#    push!(stick, (c+3, c+2))
#    g[c+3, c+3] = 0
#    b[c+3, c+3] = true
#    push!(stick, (c+3, c+3))
#    g[c+2, c+3] = 0
#    b[c+2, c+3] = true
#    push!(stick, (c+2, c+3))
#    g[c+1, c+3] = 0
#    b[c+1, c+3] = true
#    push!(stick, (c+1, c+3))
#    g[c, c+1] = 0
#    b[c, c+1] = true
#    push!(stick, (c+1, c))
#    g[c, c+2] = 0
#    b[c, c+2] = true
#    push!(stick, (c+2, c))
    for j in 1:size(g,2)
        for i in 1:size(g, 1)
            r2 = (j-c)^2 + (i-c)^2
            if r2 >= r^2
                g[i, j] = 1
                b[i, j] = true
            end
        end
    end
end

# one step SOR method
function solve(g, b)
    const omega = 1.5
    for j in 2:size(g,2)-1
        for i in 2:size(g, 1)-1
            if !b[i,j]
#                g[i, j] = (g[i+1, j] + g[i-1, j] + g[i, j+1] + g[i, j-1]) / 4
                g[i, j] = g[i, j] + omega * ((g[i+1, j] + g[i-1, j] + g[i, j+1] + g[i, j-1]) / 4 - g[i, j])
            end
        end
    end
end

# 
function get_perimeter(g, b, stick, r)
    i,j = r
    nbhd_r = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
    reduce((res, r) -> get(b, r, true)? res:push!(res, (r, get(g, r, nothing))), Set{PosVal}(), nbhd_r)
end


function perimeter(g, b, stick)
    reduce((res, x) -> union(res, get_perimeter(grid, b, stick, x)), Set{PosVal}(), stick)
end

function plist(peri)
    # eliminate zero potential
    ita = 1.0
    peri = filter(x -> x[2]>1.0e-10, peri)
    C = sum(x->x[2]^ita, peri)
    map(x->(x[1], x[2]^ita/C), peri)
end

# BUG
# 1. forではなくwhileで実装していた時に、どこかおかしかった(おそらく条件の評価)
# 2. よってfor文で書きなおした。
# 3. その後、if p >= s となってて、plの先頭を選択していたことでうまくいっていなかった
# 4. さらに、s += pl[i][2]がifの後ろについていたので、plの先頭が必ず選ばれていなかった
# 5. 下のDebug code(#DEBUG)からplistとplistからのサンプルのヒストグラムが右に１つずれていたことから
#    上記のBugを発見
# 6. 枝のループの検出はおそらく必要ない。ネットのコードにもなかった
function select(pl)
    p = rand()
    println(p)
    s = 0
    for i in 1:length(pl)
        s += pl[i][2]
        if s >= p
            return pl[i]
        end
    end
end

function select_debug(pl)
    p = rand()
    println(p)
    s = 0
    for i in 1:length(pl)
        s += pl[i][2]
        if s >= p
            return (i, pl[i])
        end
    end
end


function add!(g, b, stick, p)
    println(p)
    b[p...] = true
    push!(stick, p)
    g[p...] = 0
end

function detect_loop(g, b, stick, p)
    get_perimeter(g, b, stick, p)
end

detect_loop(grid, b, stick, (200, 200))

function print_matrix(mat, f)
    row,col = size(mat)
    for i in 1:row
        for j in 1:col
            print(f(mat[i,j]))
        end
        println(" ")
    end
end


SIZE = 200
grid = zeros(Float64, (SIZE, SIZE))
b = zeros(Bool, (SIZE, SIZE))
stick = Set{Pos}()
init(grid, b, stick)
stick
@time for i in 1:1000
    solve(grid, b)
end

p.com(g, "set size square")
p.com(g, "set grid")
p.xrange(g, (SIZE/2-20, SIZE/2+20))
p.yrange(g, (SIZE/2-20, SIZE/2+20))

for j in 1:500
    particle = select(plist(perimeter(grid, b, stick)))
    add!(grid, b, stick, particle[1])
    for i in 1:50
        solve(grid, b)
    end
end
p.splot(g, map( x -> x? 0:1, b), "matrix with image")
#p.splot_heatmap(g, grid)
#print_matrix(b[int(SIZE/2):int(SIZE/2)+2, int(SIZE/2):int(SIZE/2)+1], x->x? "X":"o")
#print_matrix(grid[int(SIZE/2):int(SIZE/2)+2, int(SIZE/2):int(SIZE/2)+1], x->string(x, " "))
#sort(plist(perimeter(grid, b, stick)), lt=(x,y)->x[2]<y[2])



# DEBUG
##pl = Vector{Float64}(map(x->x[2], sort(plist(perimeter(grid, b, stick)), lt=(x,y)->x[2]<y[2])))
#peri = perimeter(grid, b, stick)
#pl = plist(peri)
#count = Vector{Int}(zeros(length(pl)))
#SAMPLE = 50000
#for i in 1:SAMPLE
#    s = select_debug(pl)[1]
#    count[s] = count[s] + 1
#end
#p.com(g, "set style fill solid")
#pl
#count
#p.plot(g, Vector{Float64}(map(x->x[2], pl)), "using 0:1:xtic(1) w boxes")
#p.plot(g, map(x->x/SAMPLE, count), "using 0:1:xtic(1) w boxes")

print_matrix(b[200:203, 200:203], x->x? "@":" ")
print_matrix(grid[200:203, 200:203], x->string(round(x, 2), " "))

p.splot_heatmap(g, grid)

p.reset(g)
p.com(g, "set palette grey")
p.splot(g, map(x->ceil(Int, x), grid), "matrix with image")
p.splot_heatmap(g, grid)

p.reset(g)
#p.xrange(g, (0, 40))
#p.yrange(g, (0, 40))
#p.splot(g, grid, "matrix w l")
