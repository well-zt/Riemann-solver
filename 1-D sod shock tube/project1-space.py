import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.linalg import eig
gamma=1.4

# Initialize
def initial_conditions(nx, x):
    rho = np.ones(nx)
    u = np.zeros(nx)
    p = np.ones(nx)
    
    # 左侧（高密度，高压）
    rho[:nx//2] = 2.0
    p[:nx//2] = 200000
    
    # 右侧（低密度，低压）
    rho[nx//2:] = 1
    p[nx//2:] = 100000
    
    return rho, u, p

# 定义通量函数
def flux_F(rho, u, p):
    F1 = rho * u
    F2 = rho * u**2 + p
    F3 = (p / (gamma - 1)) + 0.5 * rho * u**2 + p * u
    return np.array([F1, F2, F3])

def flux_U(rho, u, p):
    U1 = rho
    U2 = rho * u
    U3 = (p / (gamma - 1)) + 0.5 * rho * u**2
    return np.array([U1, U2, U3])

# 定义Steger-Warming通量分裂方法
def steger_warming_flux_split(rho, u, p):
    c = np.sqrt(gamma * p / rho)
    lam = np.array([[u - c, 0,0],
                   [0, u, 0],
                   [0, 0, u + c]]
                   )
    lam_plus=0.5*(lam+np.abs(lam))
    lam_minus=0.5*(lam-np.abs(lam))
    # 计算特征值和特征向量

    # Jocobian Metrix
    A1=np.array([0, 1, 0])
    A2=np.array([((gamma-3)/2)*u**2, (3-gamma)*u, (gamma-1)])
    A3=np.array([(gamma-1)*u**3-gamma*u*(p/((gamma-1)*rho)+0.5*u**2),gamma*(p/((gamma-1)*rho)+0.5*u**2)-(3*(gamma-1)/2)*u**2, gamma*u])
    A= np.array([A1,A2,A3])
    # eigenvalues, eigenvectors = eig(A)

    H=c**2/(gamma-1)+0.5*u**2
    P1=np.array([1,1,1])
    P2=np.array([u-c,u,u+c])
    P3=np.array([H-u*c,0.5*u**2,H+u*c])

    P= np.array([P1,P2,P3])

    A_plus=P @ lam_plus @ np.linalg.inv(P)
    A_minus=P @ lam_minus @ np.linalg.inv(P)
    # print(lam_plus,lam_minus)

    return A,A_plus, A_minus

# 定义Van Leer通量限制器
def van_leer_limiter(x,U,direction):
    if direction==-1:
        numerator = U[:,x+1] - U[:,x]
        denominator = U[:,x] - U[:,x-1]
    if direction==1:
        numerator = U[:,x-1] - U[:,x]
        denominator = U[:,x] - U[:,x+1]
    if np.any(denominator == 0):
        r=np.array([1,1,1])
    else:
        r=numerator/denominator
    if np.any(r == 0):
        r=np.array([1,1,1])
    phi = (r + np.abs(r)) / (1 + np.abs(r))
    return r,phi

# 定义superBee通量限制器
def superBee(x,U,direction):
    if direction==-1:
        numerator = U[:,x+1] - U[:,x]
        denominator = U[:,x] - U[:,x-1]
    if direction==1:
        numerator = U[:,x-1] - U[:,x]
        denominator = U[:,x] - U[:,x+1]
    if np.any(denominator == 0):
        r=np.array([1,1,1])
    else:
        r=numerator/denominator
    if np.any(r == 0):
        r=np.array([1,1,1])
        
    phi1=max(0, min(2 * r[0], 1), min(r[0], 2))
    phi2=max(0, min(2 * r[1], 1), min(r[1], 2))
    phi3=max(0, min(2 * r[2], 1), min(r[2], 2))
    phi=np.array([phi1,phi2,phi3])
    return r,phi



# 定义求解Sod激波管问题的主函数
def solve_sod_shock_tube_upwind(nx, L, dt, t_end):
    dx = L / (nx - 1)
    # dt = t_end / (nt - 1)
    nt=int((t_end/dt)+1)
    CFL=dt/dx
    x_list = np.linspace(0, L, nx)
    
    # 初始化变量
    rho, u, p = initial_conditions(nx, x_list)
    U= flux_U(rho, u, p)
    U_prime=U
    # 时间步进循环
    for n in range(nt):
        for x in range(nx-2):#
            # print("Iteration report nx={},nt={}".format(x,n))
            if x ==0 :
                continue
            # 使用中心差分计算每个单元界面的通量
            A,A_plus, A_minus = steger_warming_flux_split(rho[x], u[x], p[x])

            r_x_positive,phi_x_positive=van_leer_limiter(x,U,1)
            r_x_negtive,phi_x_negtive=van_leer_limiter(x,U,-1)
            r_x_minus1,phi_x_minus1=van_leer_limiter(x-1,U,-1)
            r_x_plus1,phi_x_plus1=van_leer_limiter(x+1,U,1)

            #positive eigenvalue
            positive_term=-1*CFL*((np.array([1,1,1])+0.5*phi_x_positive-0.5*phi_x_minus1/r_x_minus1)*np.dot(A_plus,(U[:,x]-U[:,x-1])))
            #negtive eigenvalue            
            negtive_term=CFL*((1+0.5*phi_x_negtive-0.5*phi_x_plus1/r_x_plus1))*np.dot(A_minus,(U[:,x]-U[:,x+1]))

            # print(phi_x_minus1,r_x_minus1)

            U_prime[:,x]=positive_term+negtive_term+U[:,x]

            U[:,x]=U_prime[:,x]
            rho[x]=U[0,x]
            u[x]=U[1,x]/U[0,x]
            p[x]=(U[2,x]-0.5*rho[x]*u[x]*u[x])*(gamma-1)

            # print(rho[x], u[x], p[x])
    return x_list, rho,u,p

# 定义求解Sod激波管问题的主函数
def solve_sod_shock_tube_SW(nx, L, dt, t_end):
    dx = L / (nx - 1)
    # dt = t_end / (nt - 1)
    nt=int((t_end/dt)+1)
    CFL=dt/dx
    x_list = np.linspace(0, L, nx)
    
    # 初始化变量
    rho, u, p = initial_conditions(nx, x_list)
    U= flux_U(rho, u, p)
    U_prime=U
    # 时间步进循环
    for n in range(nt):
        for x in range(nx-2):#
            # print("Iteration report nx={},nt={}".format(x,n))
            if x ==0 :
                continue
            # 使用中心差分计算每个单元界面的通量
            A,A_plus, A_minus = steger_warming_flux_split(rho[x], u[x], p[x])


            #positive propagating
            positive_term= -1*CFL/2*np.dot(A_plus,(U[:,x+1]-U[:,x-1]))
            #negtive propagating            
            negtive_term= CFL/2*np.dot(A_minus,(U[:,x-1]-U[:,x+1]))
        
            # print(A_minus+A_plus,A)
            #positive propagating
            Linear_term= -1*CFL/2*np.dot(A,(U[:,x+1]-U[:,x-1]))
            A_square=A @ A
            square_term= CFL**2/2*np.dot(A_square,(U[:,x+1]-2*U[:,x]+U[:,x-1]))
            
            # square_term= CFL**2/2*np.dot(A_square,(phi_x_plus1*(U[:,x+1]-U[:,x])+phi_x_positive*(U[:,x]-U[:,x-1])))


            U_prime[:,x]=positive_term+negtive_term+square_term+U[:,x]

            U[:,x]=U_prime[:,x]
            rho[x]=U[0,x]
            u[x]=U[1,x]/U[0,x]
            p[x]=(U[2,x]-0.5*rho[x]*u[x]*u[x])*(gamma-1)

            # print(U[:,x+1]-U[:,x-1],A,np.dot((U[:,x+1]-U[:,x-1]),A),Linear_term)
    return x_list, rho,u,p


# 求解Sod激波管
x_list, rho_upwind,u_upwind,p_upwind = solve_sod_shock_tube_upwind(nx=200, L=2.0, dt=1e-5, t_end=0.001)
x_list, rho_SW,u_SW,p_SW = solve_sod_shock_tube_SW(nx=200, L=2.0, dt=1e-5, t_end=0.001)

##读精确解
data = np.loadtxt('exact.dat', skiprows=1)

X_POS_PLOT = data[:, 0]
DENSITY = data[:, 1]
VELOCITY_X = data[:, 2]
PRESSURE = data[:, 3]

# 画出求解结果
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(X_POS_PLOT[:-2], DENSITY[:-2],linestyle='solid', label='Density_exact', color='blue')
plt.plot(x_list[:-2], rho_upwind[:-2],linestyle='dashed', label='Density_upwind', color='orange')
plt.plot(x_list[:-2], rho_SW[:-2], linestyle='dashdot',label='Density_LW', color='green')
plt.xlabel('Position')
plt.ylabel('Density')
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(X_POS_PLOT[:-2], VELOCITY_X[:-2],linestyle='solid', label='Velocity_exact',color='blue')
plt.plot(x_list[:-2], u_upwind[:-2], label='Velocity_upwind', linestyle='dashed',color='orange')
plt.plot(x_list[:-2], u_SW[:-2], linestyle='dashdot',label='Velocity_LW', color='green')
plt.xlabel('Position')
plt.ylabel('Velocity')
plt.legend()

plt.subplot(3, 1, 3)
plt.plot(X_POS_PLOT[:-2], PRESSURE[:-2], label='Pressure_exact',linestyle='solid', color='blue')
plt.plot(x_list[:-2], p_upwind[:-2], label='Pressure_upwind',linestyle='dashed', color='orange')
plt.plot(x_list[:-2], p_SW[:-2],linestyle='dashdot', label='Pressure_LW', color='green')
plt.xlabel('Position')
plt.ylabel('Pressure')
plt.legend()

plt.tight_layout()
plt.show()
