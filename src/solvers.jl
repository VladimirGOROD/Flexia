function cros_step(func::Function,jac::Function, mass::Matrix{T}, time_step::T, u_cur::Vector{T})::Vector{T} where T <: Real
    ζ = (mass - (1+im)/2 * time_step * jac(u_cur))\func(u_cur);

    return u_cur + time_step * real.(ζ)
end

function cros!(sol::Matrix{T}, u0::Vector{T}, mass::Matrix{T}, func::Function, jac::Function, time_step::T) where T<: Real
    @assert size(sol, 1) == size(u0, 1)
    @assert size(sol, 1) == size(func(u0),1)
    sol[:, 1] = u0
    for i in 2:size(sol, 2)
        sol[:, i] = cros_step(func, jac, mass, time_step, sol[:, i-1])
    end
end

function cros!(sol::Matrix{T}, u0::Vector{T}, mass::Matrix{T}, sys::MBSystem2D, time_step::T) where T<:Real
    @assert size(sol, 1) == size(u0, 1)
    @assert size(sol, 1) == size(sys.rhs(u0), 1)
    sol[:, 1] = u0
    for i in 2:size(sol, 2)
        #cros_step(sys.func, sys.jac, mass, time_step, sol[:, i-1])
        u_cur = sol[:, i-1]
        sys.prestep(u_cur)
        update_targets!(sys)
        f_i = sys.rhs(u_cur)
        ζ = (mass - (1 + im) / 2 * time_step * sys.jacobian(u_cur)) \ f_i
        sol[:, i] = u_cur + time_step * real.(ζ);
    end
end

function simulate(sys::MBSystem2D, initial::Vector{Float64}, time_span)

    sol = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
    mass = sys.mass
    cros!(sol, initial, mass, sys, step(time_span))
    return sol;
end

function vector_norm(x::AbstractVector, p::Real=2)
    if p == 1
        return sum(abs, x)
    elseif p == 2
        return sqrt(sum(abs2, x))
    elseif p == Inf
        return maximum(abs, x)
    elseif p == -Inf
        return minimum(abs, x)
    elseif p > 0
        return (sum(abs(x).^p))^(1/p)
    else
        throw(ArgumentError("Параметр p должен быть положительным или Inf"))
    end
end

function matrix_norm(A::AbstractMatrix, norm_type::Symbol=:frobenius)
    if norm_type === :one
        # 1-норма матрицы (максимальная сумма по столбцам)
        return maximum(sum(abs, A, dims=1)[:])
    elseif norm_type === :inf
        # ∞-норма матрицы (максимальная сумма по строкам)
        return maximum(sum(abs, A, dims=2)[:])
    elseif norm_type === :frobenius
        # Фробениусова норма
        return sqrt(sum(abs2, A))
    elseif norm_type === :two
        # 2-норма матрицы (спектральная норма)
        return svd(A).S[1]  # Наибольшее сингулярное значение
    elseif norm_type === :max
        # Максимальная абсолютная величина элемента
        return maximum(abs, A)
    elseif norm_type === :nuclear
        # Ядерная норма (сумма сингулярных значений)
        return sum(svd(A).S)
    else
        throw(ArgumentError("Неизвестный тип нормы: $norm_type. Допустимые значения: :one, :two, :inf, :frobenius, :max, :nuclear"))
    end
end

function norm(A::AbstractVector, p::Real=2)
    return vector_norm(A, p)
end

function norm(A::AbstractMatrix, norm_type::Symbol=:frobenius)
    return matrix_norm(A, norm_type)
end

function newton_step(func::Function, jac::Function, u_cur::Vector{T}, max_iter = 1000, tol_e = 1e-2) where T<: Real

    u = copy(u_cur)
    n = length(u)
    
    history = [copy(u)]
    residuals = [norm(func(u))]

    for iter in max_iter
     Fx = func(u)
        current_residual = norm(Fx)
        push!(residuals, current_residual)
        
        # Проверка сходимости
        if current_residual < tol_e
            return u
        end
        
        try
            Ju = jac(u)
            
            # Решаем систему J*dx = -F
            dx = Ju \ (-Fx)
            
            # Обновляем решение
            u .+= dx
            
            push!(history, copy(u))
            
            # Проверка сходимости по изменению решения
            if norm(dx) < tol_e * max(1.0, norm(u))
                return x, iter, history, residuals
            end
            
        catch e
        end
    end
    @warn "Метод не сошёлся"

    return u
end

function static_solver!(sol::Matrix{T}, u0::Vector{T}, func::Function, jac::Function) where T<: Real
    @assert size(sol, 1) == size(u0, 1)
    @assert size(sol, 1) == size(func(u0),1)
    sol[:, 1] = u0
    
    for i in 2:size(sol, 2)
        sol[:, i] = newton_step(func, jac, sol[:, i-1])
    end

end

#

function find_natural_freqs(sys::MBSystem2D, position)

    bodies_in_order = sys.bodies            # берём из sys, чтобы точно знать порядок
    nb   = last_body_dof(sys)               # = 6 * число тел
    nλ   = last_lm_dof(sys) - nb
    ntot = nb + nλ

    pos_idx = Int[]
    vel_idx = Int[]
    for bd in bodies_in_order
        ix, iy, it = get_body_position_dofs(sys, bd)
        append!(pos_idx, (ix, iy, it))
        vx, vy, wt = get_body_velocity_dofs(sys, bd)
        append!(vel_idx, (vx, vy, wt))
    end
    λ_idx = (nb+1):ntot

    npos = length(pos_idx)

    func = sys.rhs

    jacoby = (x) -> ForwardDiff.jacobian(func, x)


    J = jacoby(position)

    E = Flexia.get_mass_matrix(sys)

    # Блоки якобиана силовой части
    mK = -J[vel_idx, pos_idx]   # = -K_eff / m (включая геом. члены от λ)
    mM = E[vel_idx, vel_idx]

    # Матрица связей (однозначно)
    C = J[λ_idx, pos_idx]

    F_qr = qr(Matrix(C'))

    Qfull = F_qr.Q * Matrix(I, npos, npos)    # принудительно расширяем до 24×24
    r = rank(C)
    V_n = Qfull[:, r+1:end]


    M̃ = V_n' * mM * V_n
    K̃ = V_n' * mK    * V_n

    # Симметризуем
    M̃_sym = Symmetric((M̃ + M̃') / 2)
    K̃_sym = Symmetric((K̃ + K̃') / 2)

    eig = eigen(Matrix(K̃_sym), Matrix(M̃_sym))
    ω² = real.(eig.values)
    order = sortperm(ω²)
    ω²    = ω²[order]
    return sqrt.(ω²) / (2π)
end