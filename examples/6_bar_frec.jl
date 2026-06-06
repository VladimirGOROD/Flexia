# =============================================================================
# МОДАЛЬНЫЙ АНАЛИЗ — версия v3, учитывающая реальную структуру Flexia.
#
# После чтения исходников Flexia выяснилось:
#
#   rhs[pos_i]   = q̇_i                                  (кинематическое тождество)
#   rhs[vel_i]   = F_ext(q, q̇) / m_i  +  (C^T λ)_i      (СИЛЫ поделены на массу,
#                                                         а множители λ — НЕТ!)
#   rhs[λ_idx]   = Φ(q)                                 (уравнения связей)
#
# Т.е. уравнения в Flexia фактически записаны как
#   q̈ = F/m + C^T λ       (а не M q̈ = F + C^T λ)
#
# Поэтому блоки якобиана имеют смысл:
#   J[vel, pos] = (1/m) ∂F/∂q  +  ∂(C^T λ)/∂q     = −K/m + геом. члены
#   J[vel, vel] = (1/m) ∂F/∂q̇                      = −D/m
#   J[vel, λ]   = C^T                               (без деления)
#   J[λ,   pos] = ∂Φ/∂q = C                         (матрица связей)
#
# Чтобы восстановить ФИЗИЧЕСКУЮ K (для задачи Mv̈ + Kv = 0, где M — настоящая
# матрица масс), нужно:
#   K = −mM · J[vel, pos]      (умножить на массу, т.к. сила была поделена)
#   D = −mM · J[vel, vel]
#   C = J[λ, pos]
#
# Проверка: J[vel, λ] должен быть C^T (без массы). Это и есть настоящий
# санити-чек.
# =============================================================================

using Pkg;
Pkg.activate("./examples");

include("dynamic_calculations_of_six_bar.jl")
# include("real_5bar.jl")

using LinearAlgebra
using Printf

# -----------------------------------------------------------------------------
# Шаг 0. Равновесие (с доводкой Ньютоном)
# -----------------------------------------------------------------------------
println("="^70)
println("ШАГ 0. Равновесие")
println("="^70)

# -----------------------------------------------------------------------------
# Шаг 1. Индексы (теперь точно известные из исходников Flexia)
# -----------------------------------------------------------------------------
# Раскладка: для каждого тела 6 dof идут подряд как [x, y, θ, ẋ, ẏ, θ̇].
# Порядок тел в state соответствует порядку add!(sys, body).

bodies_in_order = sys.bodies            # берём из sys, чтобы точно знать порядок
nb = last_body_dof(sys)               # = 6 * число тел
nλ = last_lm_dof(sys) - nb
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
nvel = length(vel_idx)
@printf "\nТел: %d,  позиций: %d,  скоростей: %d,  λ: %d,  всего: %d\n" length(bodies_in_order) npos nvel nλ ntot

# -----------------------------------------------------------------------------
# Шаг 2. Санити-чек раскладки
# -----------------------------------------------------------------------------
println("\n" * "="^70)
println("ШАГ 2. Санити-чек раскладки и блочной структуры rhs")
println("="^70)

J = jacoby(initial)
spy(J)
# Блок позиционных уравнений (должно быть q̇ = q̇)
check_pp = norm(J[pos_idx, pos_idx])
check_pv = norm(J[pos_idx, vel_idx] - I)
check_pλ = norm(J[pos_idx, λ_idx])
@printf "Позиционные уравнения q̇=q̇:\n"
@printf "  ‖J[pos, pos]‖     = %.3e  (ждём 0)\n" check_pp
@printf "  ‖J[pos, vel] − I‖ = %.3e  (ждём 0)\n" check_pv
@printf "  ‖J[pos, λ]‖       = %.3e  (ждём 0)\n" check_pλ

# Блок алгебраических уравнений (Φ зависит только от q)
check_λp = rank(J[λ_idx, pos_idx])
check_λv = norm(J[λ_idx, vel_idx])
check_λλ = norm(J[λ_idx, λ_idx])
@printf "\nАлгебраические уравнения Φ(q)=0:\n"
@printf "  rank(J[λ, pos])   = %d / %d  (матрица связей C)\n" check_λp nλ
@printf "  ‖J[λ, vel]‖       = %.3e  (ждём 0)\n" check_λv
@printf "  ‖J[λ, λ]‖         = %.3e  (ждём 0)\n" check_λλ

# -----------------------------------------------------------------------------
# Шаг 3. Физическая матрица масс (из свойств тел)
# -----------------------------------------------------------------------------
println("\n" * "="^70)
println("ШАГ 3. Физические матрицы M, K, D, C")
println("="^70)


E = Flexia.get_mass_matrix(sys)

# Блоки якобиана силовой части
mK = -J[vel_idx, pos_idx]   # = -K_eff / m (включая геом. члены от λ)
mC = -J[vel_idx, vel_idx]   # = -D / m
mM = E[vel_idx, vel_idx]
J_fλ = J[vel_idx, λ_idx]     # = C^T (без массы!)

# Матрица связей (однозначно)
C = J[λ_idx, pos_idx]

# ГЛАВНЫЙ САНИТИ-ЧЕК: согласованность C из двух блоков
# J[vel, λ] должно равняться C^T (без деления на массу)
check_Ct = J_fλ - C'
@printf "\nСанити-чек  ‖J[vel, λ] − C^T‖ = %.3e  (ждём ~0)\n" norm(check_Ct)
# spy(mK)

# Смотрим симметрию K (теперь ДОЛЖНА быть почти идеальной, если всё правильно)
K_s = (mK + mK') / 2
K_a = (mK - mK') / 2
asym_frac_K = 100 * norm(K_a) / max(norm(mK), 1e-12)
@printf "‖K‖   = %.3e,  ‖K_asym‖ = %.3e  (доля антисимметрии = %.2f %%)\n" norm(mK) norm(K_a) asym_frac_K

# -----------------------------------------------------------------------------
# Шаг 4. Нуль-пространство матрицы связей через QR
# -----------------------------------------------------------------------------
println("\n" * "="^70)
println("ШАГ 4. Базис ker(C) через QR")
println("="^70)

F_qr = qr(Matrix(C'))
# ВАЖНО: Matrix(F_qr.Q) возвращает "тонкую" Q размера 24×23 — последний
# столбец (как раз базис нуль-пространства) теряется. Запрашиваем полную:
Qfull = F_qr.Q * Matrix(I, npos, npos)    # принудительно расширяем до 24×24
r = rank(C)
V_n = Qfull[:, r+1:end]
nmin = size(V_n, 2)

@printf "rank(C) = %d / %d\n" r size(C, 1)
@printf "Минимальных координат (истинных dof): %d\n" nmin
@printf "‖C · V_n‖ = %.3e\n" norm(C * V_n)

# -----------------------------------------------------------------------------
# Шаг 5. Проецируем и решаем обобщённую задачу на СЗ
# -----------------------------------------------------------------------------
println("\n" * "="^70)
println("ШАГ 5. Редуцированная задача и собственные частоты")
println("="^70)

M̃ = V_n' * mM * V_n
K̃ = V_n' * mK * V_n

asym_Ktilde = norm((K̃ - K̃') / 2) / max(norm(K̃), 1e-12)
@printf "После проекции: доля антисимметрии K̃ = %.2f %%\n" 100 * asym_Ktilde

# Симметризуем
M̃_sym = Symmetric((M̃ + M̃') / 2)
K̃_sym = Symmetric((K̃ + K̃') / 2)

eig = eigen(Matrix(K̃_sym), Matrix(M̃_sym))
ω² = real.(eig.values)
order = sortperm(ω²)
ω² = ω²[order]
V̂ = eig.vectors[:, order]

V = V_n * 0.01 * V̂



println("\n┌───────┬──────────────────┬────────────────┬─────────────────┐")
println("│ Мода  │      ω², 1/с²    │     ω, рад/с   │      f, Гц      │")
println("├───────┼──────────────────┼────────────────┼─────────────────┤")
for i in 1:nmin
    v = ω²[i]
    if v >= 0
        ω = sqrt(v)
        f = ω / (2π)
        @printf "│ %-5d │ %16.4f │ %14.4f │ %15.4f │\n" i v ω f
    else
        ω = sqrt(-v)
        f = ω / (2π)
        @printf "│ %-5d │ %16.4f │ %13.4fi │ %14.4fi │ НЕУСТОЙЧ\n" i v ω f
    end
end
println("└───────┴──────────────────┴────────────────┴─────────────────┘")
##
# -----------------------------------------------------------------------------
# Шаг 6. Формы колебаний в исходных координатах
# -----------------------------------------------------------------------------
println("\n" * "="^70)
println("ШАГ 6. Формы колебаний")
println("="^70)

modes_q = V_n * V̂
for j in 1:nmin
    mx = maximum(abs, modes_q[:, j])
    mx > 0 && (modes_q[:, j] ./= mx)
end

nshow = min(5, nmin)
# Порядок тел в sys.bodies: bd1,bd2,bd3,bd4,bd5,bd01,bd02,bd03 (по add!)
body_labels = ["bd1", "bd2", "bd3", "bd4", "bd5", "bd01", "bd02", "bd03"]
for j in 1:nshow
    v = ω²[j]
    ω_j = v >= 0 ? sqrt(v) : NaN
    f_j = ω_j / (2π)
    @printf "\nМода %d  (ω = %.3f рад/с, f = %.3f Гц):\n" j ω_j f_j
    for (k, lab) in enumerate(body_labels[1:length(bodies_in_order)])
        base = 3 * (k - 1)
        @printf "  %-4s:  Δx = %+.3f,  Δy = %+.3f,  Δθ = %+.3f\n" lab modes_q[base+1, j] modes_q[base+2, j] modes_q[base+3, j]
    end
end

println("\n" * "="^70)
println("Результат: ω_rad = sqrt.(ω²),  формы в modes_q (размер $npos × $nmin)")
println("="^70)


state = zeros(number_of_dofs(sys), 100)
state[pos_idx, :] .= V[:, 1]
time_span = 1:100
animate(sys, state, time_span, "triv_dyn.mp4"; framerate=60, limits=(-0.1, 0.5, -0.1, 0.5))
sys.forces
##


function apply_stiffnesses!(sys, k, gr)

    for j in eachindex(gr)
        for force in sys.forces[gr[j]]
            force.stiffness = k[j]
        end
    end
    if (!assemble!(sys))
        error("Assembling failed!")
    end

    # вызвать assemble!(sys)
end

function compute_sensitivity_matrix(sys, position, groups)
    initial_freqs = Flexia.find_natural_freqs(sys, position)
    sensitivity_matrix = Matrix{Float64}(undef, length(initial_freqs), length(groups))

    for (i, group) in enumerate(groups)
        for force in sys.forces[group]
            force.stiffness += 1e-6
        end
        assemble!(sys)
        freqs_i = Flexia.find_natural_freqs(sys, position)
        sensitivity_matrix[:, i] = (freqs_i - initial_freqs) / 1e-6
        for force in sys.forces[group]
            force.stiffness -= 1e-6
        end
    end
    assemble!(sys)
    return sensitivity_matrix
end


function update_stiffnesses!(sys, position, ω_exp, initial_stiffness, stiffness_groups; α=0.3, k_min=1e-8, tol=1e-6, max_iter=50)
    L = length(stiffness_groups)
    N = length(ω_exp)
    apply_stiffnesses!(sys, initial_stiffness, stiffness_groups)
    new_stiffness = copy(initial_stiffness)
    for m in 1:max_iter
        ω_m = Flexia.find_natural_freqs(sys, position)[1:N]
        Δω = ω_m - ω_exp
        if norm(Δω) < tol
            return new_stiffness, "сошлось за $m шагов"
        end
        T = compute_sensitivity_matrix(sys, position, stiffness_groups)
        T_inv = pinv(T[:, :])
        new_stiffness -= T_inv[:, 1:N] * Δω

        #ограничеваем шаг
        # for j in 1:L
        #    if abs(Δk[j]) > α * abs(k[j])
        #     Δk[j] = sign(Δk[j]) * α * abs(k[j])
        #    end
        # end
        #что бы не было отрицательных
        for j in 1:L
            new_stiffness[j] = max(new_stiffness[j], k_min)
        end
        apply_stiffnesses!(sys, new_stiffness, stiffness_groups)
    end
    return new_stiffness, "не сошлось за $max_iter"
end

Flexia.find_natural_freqs(sys, initial)
# k = [0.2, 10]
k = [0.02, 7500.0]
apply_stiffnesses!(sys, k, gr)
Flexia.find_natural_freqs(sys, initial)
ω_exp = [6, 30]# эксп первая и вторая частоты
# ω_exp = [1.640, 2.4675]# эксп первая и вторая частоты
gr1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
gr2 = [11, 12]
# gr1 = [2, 5, 3, 4]
# gr2 = [1, 6]
gr = [gr1, gr2]
apply_stiffnesses!(sys, k, gr)

k_result, status = update_stiffnesses!(sys, initial, ω_exp, k, gr;
                                            α=0.3, max_iter=4500, tol=1e-3)
begin
    k = [0.02, 0.02]
    apply_stiffnesses!(sys, k_result, gr)
    Flexia.find_natural_freqs(sys, initial)
end
##




println("Результат: ", k_result)
println("Статус: ", status)
# println("Отн. ошибка: ", abs.(k_result .- k_true) ./ abs.(k_true))


k_true_synth = [0.5, 1.0]
apply_stiffnesses!(sys, k_true_synth, gr)
ω_synth = Flexia.find_natural_freqs(sys, initial)[1:2]
println("Синтетический эксперимент (Гц): ", ω_synth ./ (2π))

k_start = [0.1, 0.3]
k_rec, status = update_stiffnesses!(sys, initial, ω_synth, k_start, gr;
    α=0.3, max_iter=100, tol=1e-7)
println("Истинные:        ", k_true_synth)
println("Восстановленные: ", k_rec)
println("Статус:          ", status)