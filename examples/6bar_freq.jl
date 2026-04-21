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
#   K = −M_phys · J[vel, pos]      (умножить на массу, т.к. сила была поделена)
#   D = −M_phys · J[vel, vel]
#   C = J[λ, pos]
#
# Проверка: J[vel, λ] должен быть C^T (без массы). Это и есть настоящий
# санити-чек.
# =============================================================================
using Pkg
using Pkg; Pkg.activate("./examples")
Pkg.add(url="https://github.com/kutsjuice/Flexia.git", rev="diag_stat_bars")
include("real_5bar.jl")

using LinearAlgebra
using Printf

# -----------------------------------------------------------------------------
# Шаг 0. Равновесие (с доводкой Ньютоном)
# -----------------------------------------------------------------------------
println("="^70)
println("ШАГ 0. Равновесие")
println("="^70)

sol_static = Matrix{Float64}(undef, number_of_dofs(sys), length(time_span))
static_solver!(sol_static, initial, func, jacoby)
q₀ = sol_static[:, end]

@printf "После static_solver!:  ‖func(q₀)‖ = %.3e\n" norm(func(q₀))

# println("Доводка Ньютоном...")
# for iter in 1:30
#     r = func(q₀)
#     nr = norm(r)
#     nr < 1e-10 && (println("  сошлось за $(iter-1) итераций"); break)
#     Jq = jacoby(q₀)
#     δ = (Jq + 1e-12*I) \ r
#     q₀ .-= δ
#     iter % 5 == 0 && @printf "  iter %2d:  ‖func‖ = %.3e\n" iter nr
# end
# @printf "Итог:  ‖func(q₀)‖ = %.3e\n" norm(func(q₀))

# -----------------------------------------------------------------------------
# Шаг 1. Индексы (теперь точно известные из исходников Flexia)
# -----------------------------------------------------------------------------
# Раскладка: для каждого тела 6 dof идут подряд как [x, y, θ, ẋ, ẏ, θ̇].
# Порядок тел в state соответствует порядку add!(sys, body).

bodies_in_order = sys.bodies            # берём из sys, чтобы точно знать порядок
nb   = last_body_dof(sys)               # = 6 * число тел
nλ   = number_of_dofs(sys) - nb
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

J = jacoby(q₀)

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

m_phys = zeros(npos)
for (k, bd) in enumerate(bodies_in_order)
    base = 3*(k-1)
    m_phys[base+1] = bd.mass
    m_phys[base+2] = bd.mass
    m_phys[base+3] = bd.inertia
end
M_phys = Diagonal(m_phys)
@printf "Массы/инерции тел (первые 6): %s ...\n" string(round.(m_phys[1:6], digits=2))

# Блоки якобиана силовой части
J_fq = J[vel_idx, pos_idx]   # = -K_eff / m (включая геом. члены от λ)
J_fv = J[vel_idx, vel_idx]   # = -D / m
J_fλ = J[vel_idx, λ_idx]     # = C^T (без массы!)

# Матрица связей (однозначно)
C = J[λ_idx, pos_idx]

# ГЛАВНЫЙ САНИТИ-ЧЕК: согласованность C из двух блоков
# J[vel, λ] должно равняться C^T (без деления на массу)
check_Ct = J_fλ - C'
@printf "\nСанити-чек  ‖J[vel, λ] − C^T‖ = %.3e  (ждём ~0)\n" norm(check_Ct)

# Восстанавливаем физические K и D умножением на M_phys
K = -M_phys * J_fq
D = -M_phys * J_fv

# Смотрим симметрию K (теперь ДОЛЖНА быть почти идеальной, если всё правильно)
K_s = (K + K') / 2
K_a = (K - K') / 2
asym_frac_K = 100*norm(K_a)/max(norm(K), 1e-12)
@printf "‖K‖   = %.3e,  ‖K_asym‖ = %.3e  (доля антисимметрии = %.2f %%)\n" norm(K) norm(K_a) asym_frac_K

# Симметрия D аналогично
D_a = (D - D') / 2
asym_frac_D = 100*norm(D_a)/max(norm(D), 1e-12)
@printf "‖D‖   = %.3e,  ‖D_asym‖ = %.3e  (доля антисимметрии = %.2f %%)\n" norm(D) norm(D_a) asym_frac_D

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

@printf "rank(C) = %d / %d\n" r size(C,1)
@printf "Минимальных координат (истинных dof): %d\n" nmin
@printf "‖C · V_n‖ = %.3e\n" norm(C * V_n)

# -----------------------------------------------------------------------------
# Шаг 5. Проецируем и решаем обобщённую задачу на СЗ
# -----------------------------------------------------------------------------
println("\n" * "="^70)
println("ШАГ 5. Редуцированная задача и собственные частоты")
println("="^70)

M̃ = V_n' * M_phys * V_n
K̃ = V_n' * K    * V_n
D̃ = V_n' * D    * V_n

asym_Ktilde = norm((K̃ - K̃') / 2) / max(norm(K̃), 1e-12)
@printf "После проекции: доля антисимметрии K̃ = %.2f %%\n" 100*asym_Ktilde

# Симметризуем
M̃_sym = Symmetric((M̃ + M̃') / 2)
K̃_sym = Symmetric((K̃ + K̃') / 2)

eig = eigen(Matrix(K̃_sym), Matrix(M̃_sym))
ω² = real.(eig.values)
order = sortperm(ω²)
ω²    = ω²[order]
V̂    = eig.vectors[:, order]

println("\n┌───────┬──────────────────┬────────────────┬─────────────────┐")
println("│ Мода  │      ω², 1/с²    │     ω, рад/с   │      f, Гц      │")
println("├───────┼──────────────────┼────────────────┼─────────────────┤")
for i in 1:nmin
    v = ω²[i]
    if v >= 0
        ω = sqrt(v);  f = ω / (2π)
        @printf "│ %-5d │ %16.4f │ %14.4f │ %15.4f │\n" i v ω f
    else
        ω = sqrt(-v); f = ω / (2π)
        @printf "│ %-5d │ %16.4f │ %13.4fi │ %14.4fi │ НЕУСТОЙЧ\n" i v ω f
    end
end
println("└───────┴──────────────────┴────────────────┴─────────────────┘")

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
body_labels = ["bd1","bd2","bd3","bd4","bd5","bd01","bd02","bd03"]
for j in 1:nshow
    v = ω²[j]
    ω_j = v >= 0 ? sqrt(v) : NaN
    f_j = ω_j / (2π)
    @printf "\nМода %d  (ω = %.3f рад/с, f = %.3f Гц):\n" j ω_j f_j
    for (k, lab) in enumerate(body_labels[1:length(bodies_in_order)])
        base = 3*(k-1)
        @printf "  %-4s:  Δx = %+.3f,  Δy = %+.3f,  Δθ = %+.3f\n" lab modes_q[base+1,j] modes_q[base+2,j] modes_q[base+3,j]
    end
end

println("\n" * "="^70)
println("Результат: ω_rad = sqrt.(ω²),  формы в modes_q (размер $npos × $nmin)")
println("="^70)