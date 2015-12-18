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
    perimeters::Dict{Pos, Float64}
    function DBM(size::Int)
        grid = zeros(Float64, (SIZE, SIZE))
        b = zeros(Bool, (SIZE, SIZE))
        stick = Set{Pos}()
        new(size, grid, b, stick, Dict{Pos, Float64}())
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
    for j in 1:dbm.size
        for i in 1:dbm.size
            r2 = (j-c)^2 + (i-c)^2
            if r2 >= r^2
                dbm.grid[i, j] = 1
                dbm.b[i, j] = true
            end
        end
    end
    # 初期のperimeterは系の中心にあるとするので、periがboundaryかどうかの
    # checkはしない
    # eight site!!!
    perimeter_site = [(c-1, c), (c+1, c), (c, c-1), (c, c+1)]
#    perimeter_site = [(c+i, c+j) for i in -1:1, j in -1:1]
    for peri in perimeter_site
        dbm.perimeters[peri] = 0.0
    end
    delete!(dbm.perimeters, (c, c))
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
# 7. ループ内部はポテンシャルが境界値と等しくなるので確率は０となる
#    よってループの検出は不必要
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

for j in 1:100
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

# calc potential at new perimeter site
function calc_perimeters(dbm::DBM, perimeters)
    for (peri, _) in perimeters
        new_phi = phi_i(peri, dbm.stick)
        perimeters[peri] = new_phi
    end
    return perimeters
end

function get_new_perimeters(dbm::DBM, r::Pos)
    function set_element(dict, dbm, r)
#        isboundary = get(dbm.b, r, true)
        isboundary = dbm.b[r...]
        if (!isboundary && !(r in keys(dbm.perimeters)))
            dict[r] = 0.0
        end
        return dict
    end
    i,j = r
    peri_pos = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
# eight site!!!
#    peri_pos = [(i+_i, j+_j) for _i in -1:1, _j in -1:1]
    new_peri = Dict{Pos, Float64}()
    reduce((d, pos)->set_element(d, dbm, pos), new_peri, peri_pos)
end

function update_perimeters(dbm::DBM, new_charge_pos)
    delete!(dbm.perimeters, new_charge_pos)
    for (peri, phi) in dbm.perimeters
        dbm.perimeters[peri] = phi + phi_j2i(new_charge_pos, peri)
    end
end

function add_charge_and_perimeters!(dbm::DBM, p)
#    println(p)
    dbm.b[p...] = true
    push!(dbm.stick, p)
    dbm.grid[p...] = 0
    update_perimeters(dbm, p)
    new_perimeters = calc_perimeters(dbm, get_new_perimeters(dbm, p))
#    println(dbm.perimeters)
#    println(new_perimeters)
    merge!(dbm.perimeters, new_perimeters)
end

#n = get_new_perimeters(dbm, (int(dbm.size/2+1), int(dbm.size/2+1)))
#dbm.perimeters

function plist2(peri)
    ita = 4.0
    phi_min = min(values(peri)...)
    phi_max = max(values(peri)...)
    delta = phi_max-phi_min
    map(x->(x[1], ((x[2]-phi_min)/delta)^ita), peri)
    C = sum(x->x[2], peri)
    map(x->(x[1], x[2]/C), peri)
end

function iterate(dbm::DBM, n::Int)
    for i in 1:n
        pl = plist2(dbm.perimeters)
        new_charge = select(pl)
        add_charge_and_perimeters!(dbm, new_charge[1])
    end
end

SIZE = 200
dbm = DBM(SIZE)
init(dbm)
dbm.stick
# init perimeters
dbm.perimeters = calc_perimeters(dbm, dbm.perimeters)
# add new charge


function writef(f, stick)
    for s in stick
        write(f, string(s[1], ", ", s[2], "\n"))
    end
    write(f, "\n")
end

pl = plist2(dbm.perimeters)
#new_charge = select(pl)
#add_charge_and_perimeters!(dbm, new_charge[1])
iterate(dbm, 100)
vec = Pos[s for s in dbm.stick]
#p._plot_data(g, "plot", writef, vec, "w p  ps 2.0 pt 5")
p.splot(g, map( x -> x? 0:1, dbm.b), "matrix with image")
#p.com(g, "set size square")
#p.splot_heatmap(g, dbm.grid)

pl = plist2(dbm.perimeters)
dbm.perimeters
p.reset(g)

iterate(dbm, 1)
p.plot(g, Float64[phi for (p, phi) in dbm.perimeters], "w lp")
