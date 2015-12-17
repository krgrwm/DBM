import MyPack
p = MyPack.Plot
g = p.Gnuplot()

typealias Pos Tuple{Int, Int}
typealias PosVal Tuple{Pos, Float64}

type DBM
    size::Int
    grid::Array{Float64, 2}
    b::Array{Bool, 2}
    stick::Set{Pos}
    function DBM(size::Int)
        grid = zeros(Float64, (SIZE, SIZE))
        b = zeros(Bool, (SIZE, SIZE))
        stick = Set{Pos}()
        new(size, grid, b, stick)
    end
end


# Init grid, boundary, stick
function init(dbm::DBM)
    const c = round(Int, dbm.size/2)
    const r = c-2
    # seed
    dbm.grid[c, c] = 0
    dbm.b[c, c] = true
    push!(dbm.stick, (c, c))
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
    for j in 1:dbm.size
        for i in 1:dbm.size
            r2 = (j-c)^2 + (i-c)^2
            if r2 >= r^2
                dbm.grid[i, j] = 1
                dbm.b[i, j] = true
            end
        end
    end
end

# one step SOR method
function solve!(g, b)
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

function calc_potential(dbm::DBM, n)
    for i in 1:n
        solve!(dbm.grid, dbm.b)
    end
end

# 
function get_perimeter(dbm::DBM, r::Pos)
    i,j = r
    nbhd_r = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
    reduce((res, r) -> get(dbm.b, r, true)? res:push!(res, (r, get(dbm.grid, r, nothing))), Set{PosVal}(), nbhd_r)
end


function perimeter(dbm::DBM)
    reduce((res, x) -> union(res, get_perimeter(dbm, x)), Set{PosVal}(), dbm.stick)
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


function add!(dbm::DBM, p)
    println(p)
    dbm.b[p...] = true
    push!(dbm.stick, p)
    dbm.grid[p...] = 0
end

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
dbm = DBM(SIZE)
init(dbm)

@time calc_potential(dbm, 1000)

p.com(g, "set size square")
p.com(g, "set grid")
p.xrange(g, (SIZE/2-20, SIZE/2+20))
p.yrange(g, (SIZE/2-20, SIZE/2+20))
#p.com(g, "set palette grey")

for j in 1:200
    particle = select(plist(perimeter(dbm)))
    add!(dbm, particle[1])
    calc_potential(dbm, 100)
end
p.splot(g, map( x -> x? 0:1, dbm.b), "matrix with image")

p.splot_heatmap(g, dbm.grid)


# Fast Calculation
function rij(ri, rj)
    delta = map(-, ri, rj)
    sqrt(reduce((acc, v)->acc+v^2, 0, delta))
end

function phi_j2i(ri, rj)
    R1 = 0.5
    1 - R1/rij(ri, rj)
end

function phi_i(ri::Pos, charges)
    sum(map(rj->phi_j2i(ri, rj), charges))
end

SIZE = 200
dbm = DBM(SIZE)
init(dbm)
dbm.stick
perimeter(dbm)
get_perimeter(dbm, (30, 40))


function get_perimeterDICT(dbm::DBM, r::Pos)
    i,j = r
    nbhd_r = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
    reduce((res, r) -> get(dbm.b, r, true)? res:push!(res, (r, get(dbm.grid, r, nothing))), Set{PosVal}(), nbhd_r)
end


function perimeterDICT(dbm::DBM)
    reduce((res, x) -> union(res, get_perimeterDICT(dbm, x)), Set{PosVal}(), dbm.stick)
end


for peri in perimeter(dbm)
    updated = (peri[1], phi_i(peri[1], dbm.stick))
end
