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