module SpecialMatrices

function selection_matrix(m::Int)
    cd = round(Int, m*(m+1)/2)
    rd = m*m

    D = zeros(Int, rd, cd)
    @inbounds for j = 1:m
        for i = 1:j
            r_ij = round(Int, (j*j-j)/2 + i)
            h_ij = round(Int, m*(j-1) + i)
            h_ji = round(Int, m*(i-1) + j)
            D[h_ij, r_ij] = 1
            D[h_ji, r_ij] = 1
        end
    end
    D
end

function Hmat(i, j, m, n)
    H = zeros(m,n)
    H[i,j] = 1
    H
end

function commutation_matrix(m::Int, n::Int = m)
    K = zeros(m*n, m*n)
    for i = 1:m
        for j = 1:n
            Hij = Hmat(i, j, m, n)
            K += kron(Hij, Hij')
        end
    end
    K
end

function Nm(m::Int)
    Km = commutation_matrix(m,m)
    for j = 1:m^2
        Km[j,j] = Km[j,j] + 1.0
    end
    Km./2
end

function Emat(i, j, m)
    E = zeros(m,m)
    E[i,j] = 1
    E
end

function Tmat(i, j, m)
    if i ≠ j
        Emat(i, j, m) + Emat(j, i, m)
    else
        Emat(i,i,m)
    end
end

function uvec(i, j, m)
    r = round(Int, m*(m+1)/2)
    u = zeros(r, 1)
    idx = round(Int, (j-1)*m + i - j*(j-1)/2)
    u[idx,1] = 1
    u
end

function Dm(m::Int)
    D = zeros(m^2, round(Int, m*(m+1)/2))
    for j = 1:m
        for i = j:m
            D += vec(Tmat(i, j, m))*uvec(i,j,m)'
        end
    end
    D
end

function Lmt(m::Int)
    D = Dm(m)
    for j = 1:(m-1)
        D[j*m+(1:j),:] = 0.0
    end
    D
end

function Lm(m::Int)
    Lmt(m)'
end

export Lm, Dm, Nm, selection_matrix, commutation_matrix

end # module
