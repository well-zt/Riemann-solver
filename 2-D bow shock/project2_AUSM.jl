#/usr/local/bin/julia
#coded based on julia v1.9
#copyright © ZHANG Ting, polyU, Hong Kong
#e-mail ting123.zhang@connect.polyu.hk 
#updated at 05/Dec/2024
using CSV
using DataFrames
using Plots
using PyCall
using LinearAlgebra
using Colors
using Printf
@pyimport plot3d as p3d
using Base.Threads: @threads, nthreads

# 确认线程数量
println("Number of threads: ", nthreads())
# 读取数据
block = p3d.read_plot3D("Cylinder.dat", binary = false)

IMAX, JMAX, KMAX = size(block[1].X)
X_coordinate = block[1].X
Y_coordinate = block[1].Y
U_bar = zeros(Float64, IMAX+3, JMAX+3, 4)


gamma = 1.4
R = 8.314
M_IG = 0.029
c_v = 717
R_G=R/M_IG
# 定义流动条件和初始条件
M_inf = 8.1
T_inf = 63.73
p_inf = 370.6
rho_inf=p_inf/(R_G*T_inf)
a_inf=sqrt(gamma*R_G*T_inf)
u_inf=M_inf*a_inf

rho_after_shock=0.116
p_after_shock=28465
M_after_shock=0.392
T_after_shock=849.75
a_after_shock=sqrt(gamma*R_G*T_after_shock)
u_after_shock=a_after_shock*M_after_shock
# p_b=101325
p_b=p_inf

T_total=T_inf*(1+((gamma-1)/2)*M_inf^2)
p_total=p_inf*(1+((gamma-1)/2)*M_inf^2)^(gamma/(gamma-1))

# 初始化
function initial_conditions()
    rho = fill(rho_after_shock, IMAX+3, JMAX+3)
    u = fill(0,IMAX+3, JMAX+3) #M_inf * sqrt(gamma * R_G * T_inf)
    v = zeros(IMAX+3, JMAX+3)
    p = fill(p_after_shock, IMAX+3, JMAX+3)
    return rho, u, v, p
end


function U_flux_decompose(U)
    rho = U[1]
    u = U[2] / U[1]
    v = U[3] / U[1]
    p = (U[4] - 0.5 * rho * (u^2 + v^2)) * (gamma - 1)
    return rho, u, v, p
end

function area(I, J)
    I=I-2
    J=J-2
    p1 = [X_coordinate[I, J][1], Y_coordinate[I, J][1],0]
    p2 = [X_coordinate[I, J+1][1], Y_coordinate[I, J+1][1],0]
    p3 = [X_coordinate[I+1, J+1][1], Y_coordinate[I+1, J+1][1],0]
    p4 = [X_coordinate[I+1, J][1], Y_coordinate[I+1, J][1],0]
    r_1 = p1 - p3
    r_2 = p2 - p4
    cross_product = cross(r_1,r_2)
    magnitude = norm(cross_product)
    return magnitude
end

function length(a, b, c, d)
    a=a-2
    b=b-2
    c=c-2
    d=d-2
    p1 = [X_coordinate[a, b][1], Y_coordinate[a, b][1]]
    p2 = [X_coordinate[c, d][1], Y_coordinate[c, d][1]]
    p = p1 - p2
    mag = norm(p)
    return mag
end

function AUSM(U_L, U_R, n)
    rho_L, u_L, v_L, p_L = U_flux_decompose(U_L)
    rho_R, u_R, v_R, p_R = U_flux_decompose(U_R)
    if p_L/rho_L*p_R/rho_R < 0
        println("rho=",rho_L, "u=",u_L, "v=",v_L, "p=",p_L)
        println("rho=",rho_R, "u=",u_R, "v=",v_R, "p=",p_R)
        println(U_L)
        println(U_R)
    end 
    a_L = sqrt(gamma *  p_L/rho_L)
    a_R = sqrt(gamma *  p_R/rho_R)
    a = (a_L + a_R) / 2
    M_L = (u_L*n[1]+v_L*n[2]) / a_L
    M_R = (u_R*n[1]+v_R*n[2]) / a_R

    M_plus = M_L > 1 ? M_L : M_L < -1 ? 0 : (M_L + 1)^2 / 4
    M_minus = M_R > 1 ? 0 : M_R < -1 ? M_R : -(M_R - 1)^2 / 4

    u = (M_plus + M_minus) * a
    if u >= 0
        PHI =U_L
        PHI[end]=p_L * (gamma / (gamma - 1)) + (rho_L * 0.5 * (u_L^2 + v_L^2))
    else
        PHI =U_R
        PHI[end]=p_R* (gamma / (gamma - 1)) + (rho_R * 0.5 * (u_R^2 + v_R^2))
    end

    F_c = u * PHI

    p_plus = M_L > 1 ? 1 : M_L < -1 ? 0 : (M_L + 1)^2 * (2 - M_L) / 4
    p_minus = M_R > 1 ? 0 : M_R < -1 ? (M_R - 1)^2 * (2 + M_R) / 4 : 1

    p_ASUM = p_plus * p_L + p_minus * p_R

    F_p = [0, p_ASUM*n[1], p_ASUM*n[2], 0]

    # M_L = v_L / a_L
    # M_R = v_R / a_R

    # M_plus = M_L > 1 ? M_L : M_L < -1 ? 0 : (M_L + 1)^2 / 4
    # M_minus = M_R > 1 ? 0 : M_R < -1 ? M_R : -(M_R - 1)^2 / 4

    # v = (M_plus + M_minus) * a
    # if v >= 0
    #     PHI =U_L
    #     PHI[end]=p_L * (gamma / (gamma - 1)) * (rho_L * 0.5 * (u_L^2 + v_L^2))
    # else
    #     PHI =U_R
    #     PHI[end]=p_R * (gamma / (gamma - 1)) * (rho_R * 0.5 * (u_R^2 + v_R^2))
    # end

    # G_c = v * PHI

    # p_plus = M_L > 1 ? 1 : M_L < -1 ? 0 : (M_L + 1)^2 * (2 - M_L) / 4
    # p_minus = M_R > 1 ? 0 : M_R < -1 ? (M_R - 1)^2 * (2 + M_R) / 4 : 1

    # p_ASUM = p_plus * p_L + p_minus * p_R

    # G_p = [0, 0, p_ASUM, 0]

    return F_c , F_p #*n[1] +G_p*n[2]
end

function normal_vector(a, b, c, d)
    a=a-2
    b=b-2
    c=c-2
    d=d-2
    p1 = [X_coordinate[a, b][1], Y_coordinate[a, b][1]]
    p2 = [X_coordinate[c, d][1], Y_coordinate[c, d][1]]
    p = p1 - p2
    magnitude = norm(p) 
    unit_vector = p / magnitude
    return unit_vector
end

function main(dt, t_end)
    # initialize
    rho, u, v, p = initial_conditions()
    U_bar = cat(rho, rho .* u, rho .* v, 0.5 .* rho .* (u.^2 .+ v.^2) .+ p ./ (gamma - 1); dims=3)
    for t in 0:dt:t_end
        println("t=$t")
        for i in 3:IMAX+1
            @threads for j in 3:JMAX+1
                U = U_bar[i, j, :]
                omega = area(i, j)
                # println(i,j)
                U_L = U_bar[i, j+1,:] .+ 0.5 .* (U_bar[i, j+1,:] .- U_bar[i, j+2,:])
                U_R = U_bar[i, j,:] .+ 0.5 .* (U_bar[i, j,:] .- U_bar[i, j-1,:])
                # println(U_bar[i, j,:])
                # println(U_bar[i, j-1,:])
                # left face
                S_left = length(i, j+1, i+1, j+1)
                n_left = normal_vector(i, j+1, i, j)
                F_c, F_p = AUSM(U_R, U_L, n_left)
                F_left = F_c .+ F_p
                # if j==JMAX
                #     println(F_left)
                # end
                
                U_L = U_bar[i, j,:] .+ 0.5 .* (U_bar[i, j,:] .- U_bar[i, j+1,:])
                U_R = U_bar[i, j-1,:] .+ 0.5 .* (U_bar[i, j-1,:] .- U_bar[i, j-2,:])
                # println("right face")
                # right face
                S_right = length(i, j, i+1, j)
                n_right = normal_vector(i, j, i, j+1)
                F_c, F_p = AUSM(U_L, U_R, n_right)
                F_right = F_c .+ F_p
                # if j==JMAX
                #     println(F_right)
                # end

                U_L = U_bar[i, j,:] .+ 0.5 .* (U_bar[i, j,:] .- U_bar[i-1, j,:])
                U_R = U_bar[i+1, j,:] .+ 0.5 .* (U_bar[i+1, j,:] .- U_bar[i+2, j,:])
                # println(U_R)
                # upper face
                S_upper = length(i+1, j, i+1, j+1)
                n_upper = normal_vector(i+1, j, i, j)
                F_c, F_p = AUSM(U_L, U_R, n_upper)
                F_upper = F_c .+ F_p
                
                U_L = U_bar[i-1, j,:] .+ 0.5 .* (U_bar[i-1, j,:] .- U_bar[i-2, j,:])
                U_R = U_bar[i, j,:] .+ 0.5 .* (U_bar[i, j,:] .- U_bar[i+1, j,:])
                # println(U_R)
                # down face
                S_down = length(i, j, i, j+1)
                n_down = normal_vector(i, j, i+1, j)
                F_c, F_p = AUSM(U_R, U_L, n_down)
                F_down = F_c .+ F_p
                
                U_temp = U .- dt / omega .* (S_left .* F_left .+ S_right .* F_right .+ S_upper .* F_upper .+ S_down .* F_down)
                U_bar[i, j,:] = U_temp
                # if j==JMAX
                #     println(F_left)
                #     println(dt / omega .*(S_left .* F_left .+ S_right .* F_right .+ S_upper .* F_upper .+ S_down .* F_down))
                # end
            end
        end

        #define some monitor
        println("rho=",U_bar[80,40,1],"u=",U_bar[80,40,2]/U_bar[80,40,1],"p=",(U_bar[80,40,4] - 0.5 * U_bar[80,40,1] * ((U_bar[80,40,2]/U_bar[80,40,1])^2 + (U_bar[80,40,3]/U_bar[80,40,1])^2)) * (gamma - 1))
        println("rho=",U_bar[80,end,1],"u=",U_bar[80,end,2]/U_bar[80,end,1],"p=",(U_bar[80,end,4] - 0.5 * U_bar[80,end,1] * ((U_bar[80,end,2]/U_bar[80,end,1])^2 + (U_bar[80,end,3]/U_bar[80,end,1])^2)) * (gamma - 1))
        println("rho=",U_bar[80,end-3,1],"u=",U_bar[80,end-3,2]/U_bar[80,end-3,1],"p=",(U_bar[80,end-3,4] - 0.5 * U_bar[80,end-3,1] * ((U_bar[80,end-3,2]/U_bar[80,end-3,1])^2 + (U_bar[80,end-3,3]/U_bar[80,end-3,1])^2)) * (gamma - 1))

        # boundary conditions
        for i in 3:IMAX+1
            # U_bar[i, 2,:] = U_bar[i, 3,:]
            # Inviscid wall
            n_right = normal_vector(i, 3, i, 4)
            n_upper = normal_vector(i+1, 3, i, 3)
            A = [1 0 0 0; 0 n_right[1] n_right[2] 0; 0 n_upper[1] n_upper[2] 0; 0 0 0 1]
            B1 = [1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1]
            B2 = [1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1]
            # println(A)
            C = [1 0 0 0; 0 n_right[1] n_upper[1] 0; 0 n_right[2] n_upper[2] 0; 0 0 0 1]
            # coe_matrix = C * B1 * A
            # println(U_bar[i, 3,:])
            # println(B2 *A*U_bar[i, 3,:])
            U_bar[i, 2,:]= C * B2 * A*U_bar[i, 3,:]
            U_bar[i, 1,:]= C * B2 * A*U_bar[i, 4,:]
            # P_U_bar = coe_matrix * U_bar[i, 3,:]
            # U_bar[i, 1,:] = 3 .* U_bar[i, 2,:] - 2 .* P_U_bar
        end
        # Supersonic inflow 
        for i in 3:IMAX+1
            rho_in = p_inf  / (T_inf* R_G)
            u_in = M_inf * sqrt(gamma * R_G * T_inf)
            v_in = 0
            p_in = p_inf
            U_bar[i, end,1]=U_bar[i, end-1,1]=rho_in
            U_bar[i, end,2]=U_bar[i, end-1,2]=rho_in*u_in
            U_bar[i, end,3]=U_bar[i, end-1,3]=rho_in*v_in
            U_bar[i, end,4]=U_bar[i, end-1,4]=0.5 * rho_in * (u_in^2 + v_in^2) + p_in / (gamma - 1)
            # println(U_bar[i, end,:])
        end
        # trans-sonic outlet boundary condition
        for j in 3:JMAX+1
            U_bar[1, j,:] = U_bar[3, j,:]
            U_bar[2, j,:] = U_bar[3, j,:]
            U_bar[end, j,:] = U_bar[IMAX+1, j,:]
            U_bar[end-1, j,:] = U_bar[IMAX+1, j,:]
        #=
            rho_out, u_out, v_out, p_out=U_flux_decompose(U_bar[IMAX+1, j,:])
            a_out=sqrt(gamma *  p_out/rho_out)
            M_out=sqrt(u_out^2+v_out^2)/a_out
            n=normal_vector(IMAX+2,j,IMAX+1,j)
            if M_out >= 1
                U_bar[end, j,:] = U_bar[IMAX+1, j,:]
                U_bar[end-1, j,:] = U_bar[IMAX+1, j,:]
            else
                p=p_b
                rho=rho_out+(p_b-p_out)/a_out^2
                u=u_out+(p_b-p_out)*n[1]/(rho_out*a_out)
                v=v_out+(p_b-p_out)*n[2]/(rho_out*a_out)
                P_U_bar=[rho,rho*u,rho*v, 0.5* rho * (u^2 + v^2) + p / (gamma - 1)]
                U_bar[end-1, j,:]=2 .*P_U_bar .- U_bar[IMAX+1, j,:]
                U_bar[end-1, j,:]=4 .*P_U_bar .- 3 .*U_bar[IMAX+1, j,:]
            end
            rho_out, u_out, v_out, p_out=U_flux_decompose(U_bar[3, j,:])
            a_out=sqrt(gamma *  p_out/rho_out)
            M_out=sqrt(u_out^2+v_out^2)/a_out
            n=normal_vector(3,j,4,j)
            if M_out >= 1
                U_bar[1, j,:] = U_bar[3, j,:]
                U_bar[2, j,:] = U_bar[3, j,:]
            else
                p=p_b
                rho=rho_out+(p_b-p_out)/a_out^2
                u=u_out+(p_b-p_out)*n[1]/(rho_out*a_out)
                v=v_out+(p_b-p_out)*n[2]/(rho_out*a_out)
                P_U_bar=[rho,rho*u,rho*v, 0.5* rho * (u^2 + v^2) + p / (gamma - 1)]
                U_bar[2, j,:]=2 .* P_U_bar .- U_bar[3, j,:]
                U_bar[1, j,:]=4 .* P_U_bar .- 3 .*U_bar[3, j,:]
            end
            =#
        end
    
    end
    
    return U_bar
end

# 调用main函数
U_bar = main(1e-6, 1e-1)

rho_final = U_bar[:, :,1]
u_final = U_bar[:, :,2] ./ rho_final
v_final = U_bar[:, :,3] ./ rho_final
p_final = (U_bar[:, :,4] .- 0.5 .* rho_final .* (u_final.^2 .+ v_final.^2)) .* (gamma - 1)
# println(rho_final)

# 定义绘制彩色方块的函数
function plot_colored_square(x_coords, y_coords, color)
    plot!(x_coords, y_coords, seriestype = :shape, fillcolor = color, linecolor = :transparent)
end
# 辅助函数：获取颜色
function get_color(value, min_value, delta)
    return cgrad(:rainbow)[(value - min_value) / delta]
end
# 绘制rho
# normalized_rho = (rho_final .- minimum(rho_final)) ./ (maximum(rho_final) - minimum(rho_final))
p_rho = plot(figsize = (800, 800),legend=false)
count_greater_than_10 = 0  # 初始化计数器
for i in 3:IMAX+1
    for j in 3:JMAX+1
        X = X_coordinate
        Y = Y_coordinate
        # 获取单元(cell)的四个顶点坐标
        I_coordinate=i-2
        J_coordinate=j-2
        cell_x = [X[I_coordinate, J_coordinate], X[I_coordinate+1, J_coordinate], X[I_coordinate+1, J_coordinate+1], X[I_coordinate, J_coordinate+1], X[I_coordinate, J_coordinate]]
        cell_y = [Y[I_coordinate, J_coordinate], Y[I_coordinate+1, J_coordinate], Y[I_coordinate+1, J_coordinate+1], Y[I_coordinate, J_coordinate+1], Y[I_coordinate, J_coordinate]]
        
        Delta = maximum(rho_final[3:IMAX+1,3:JMAX+1]) - minimum(rho_final[3:IMAX+1,3:JMAX+1]) == 0 ? 1 : maximum(rho_final[3:IMAX+1,3:JMAX+1]) - minimum(rho_final[3:IMAX+1,3:JMAX+1])
        color = get_color(rho_final[i, j], minimum(rho_final[3:IMAX+1,3:JMAX+1]), Delta)
        
        # 绘制单元
        plot_colored_square(cell_x, cell_y, color)
        global count_greater_than_10
        # 检查 rho[i, j] 是否大于10，并更新计数器
        if rho_final[i, j] > 1
            count_greater_than_10 += 1
        end
    end
end

# 输出大于10的数量
println("Number of rho values greater than 10: ", count_greater_than_10)

colorbar_ticks = range(minimum(rho_final[3:IMAX+1,3:JMAX+1]), stop=maximum(rho_final[3:IMAX+1,3:JMAX+1]), length=11)
plot!(p_rho, color=:rainbow, colorbar=:right, colorbar_ticks=colorbar_ticks, colorbar_tick_labels=[@sprintf("%.2f", x) for x in colorbar_ticks])
# Colorbar(p_rho, pltobj)
title!(p_rho, "rho, rhomin=$( @sprintf("%.3f", minimum(rho_final[3:IMAX+1,3:JMAX+1])) ), rhomax=$( @sprintf("%.3f", maximum(rho_final[3:IMAX+1,3:JMAX+1])) )")
savefig(p_rho,"rho.png")


# 绘制p
p_p=plot(figsize = (800, 800),legend=false)
for i in 3:IMAX+1
    for j in 3:JMAX+1
        X = X_coordinate
        Y = Y_coordinate
        I_coordinate=i-2
        J_coordinate=j-2
        cell_x = [X[I_coordinate, J_coordinate], X[I_coordinate+1, J_coordinate], X[I_coordinate+1, J_coordinate+1], X[I_coordinate, J_coordinate+1], X[I_coordinate, J_coordinate]]
        cell_y = [Y[I_coordinate, J_coordinate], Y[I_coordinate+1, J_coordinate], Y[I_coordinate+1, J_coordinate+1], Y[I_coordinate, J_coordinate+1], Y[I_coordinate, J_coordinate]]
        
        Delta = maximum(p_final[3:IMAX+1,3:JMAX+1]) - minimum(p_final[3:IMAX+1,3:JMAX+1]) == 0 ? 1 : maximum(p_final[3:IMAX+1,3:JMAX+1]) - minimum(p_final[3:IMAX+1,3:JMAX+1])
        color = get_color(p_final[i, j], minimum(p_final[3:IMAX+1,3:JMAX+1]), Delta)
        
        # 绘制单元
        plot_colored_square(cell_x, cell_y, color)
    end
end
title!(p_p,"p, pmin=$( @sprintf("%.3f", minimum(p_final[3:IMAX+1,3:JMAX+1])) ), pmax=$( @sprintf("%.3f", maximum(p_final[3:IMAX+1,3:JMAX+1])) )")
savefig(p_p,"p.png")


# # 绘制u
p_u=plot(figsize = (800, 800),legend=false)
for i in 3:IMAX+1
    for j in 3:JMAX+1
        X = X_coordinate
        Y = Y_coordinate
        I_coordinate=i-2
        J_coordinate=j-2
        cell_x = [X[I_coordinate, J_coordinate], X[I_coordinate+1, J_coordinate], X[I_coordinate+1, J_coordinate+1], X[I_coordinate, J_coordinate+1], X[I_coordinate, J_coordinate]]
        cell_y = [Y[I_coordinate, J_coordinate], Y[I_coordinate+1, J_coordinate], Y[I_coordinate+1, J_coordinate+1], Y[I_coordinate, J_coordinate+1], Y[I_coordinate, J_coordinate]]
        
        Delta = maximum(u_final[3:IMAX+1,3:JMAX+1]) - minimum(u_final[3:IMAX+1,3:JMAX+1]) == 0 ? 1 : maximum(u_final[3:IMAX+1,3:JMAX+1]) - minimum(u_final[3:IMAX+1,3:JMAX+1])
        color = get_color(u_final[i, j], minimum(u_final[3:IMAX+1,3:JMAX+1]), Delta)
        
        # 绘制单元
        plot_colored_square(cell_x, cell_y, color)
    end
end
title!(p_u,"u, umin=$(@sprintf("%.3f", abs(minimum(u_final[3:IMAX+1,3:JMAX+1])))), umax=$(@sprintf("%.3f", abs(maximum(u_final[3:IMAX+1,3:JMAX+1]))))")
savefig(p_u,"u.png")
#=
=#