'''
AUTHOR:         B ANEESH BHAT
VERSION:        4
REFERENCE:      E. F. Toro, “Riemann Solvers and Numerical Methods for Fluid Dynamics 3rd EDITION,” Springer Verlag, Berlin, 2009
'''

#%%
import numpy as np
import matplotlib.pyplot as plt
import sys

#%%---------------FUNCTIONS DEFINITIONS---------------#
#
def GUESS_P(P_L, U_L, RHO_L, P_R, U_R, RHO_R, G1, G2, G3, G4, G5, G6, G7, G8):

    Q_USER = 2.0

#----------COMPUTE GUESS PRESSURE FROM PVRS-------------#
    #-------PVRS = PRIMITIVE VARIABLE RIEMANN SOLVER--------#
    A_UP = 0.25*(RHO_L +  RHO_R) * (A_L + A_R)
    P_PV = 0.5 * (P_L + P_R) + 0.5 * (U_L - U_R) * A_UP
    P_PV = np.maximum(0.0, P_PV)
    P_MIN = np.minimum(P_L, P_R)
    P_MAX = np.maximum(P_L, P_R)
    Q_MAX = P_MAX / P_MIN

    if (Q_MAX <= Q_USER and (P_MIN <= P_PV and P_PV <= P_MAX)):
        #-------SELECT PVRS GUESS---------------#
        P_M = P_PV

    elif (P_PV <= P_MIN):
        #-------SELECT TWO-RAREFACTION RIEMANN SOLVER---------#
        P_Q = (P_L/ P_R)**G1
        U_M = ((P_Q * U_L / A_L) + (U_R / A_R) + G4*(P_Q - 1.0)) / ((P_Q / A_L) + (1.0 / A_R))
        P_TL = 1.0 + (G7 * (U_L - U_M) / A_L)
        P_TR = 1.0 + (G7 * (U_M - U_R) / A_R)
        P_M = 0.5 * ((P_L * P_TL**G3) + (P_R * P_TR**G3))
    
    else:
        #-------SELECT TWO-SHOCK RIEMANN SOLVER----------#
        #-------WITH PVRS AS AN ESTIMATE-----------------#
        G_EL = np.sqrt((G5 / RHO_L) / (G6 * P_L + P_PV))
        G_ER = np.sqrt((G5 / RHO_R) / (G6 * P_R + P_PV))
        P_M = ((G_EL * P_L) + (G_ER * P_R) - (U_R - U_L)) / (G_EL + G_ER)
    
    print("\n=>   P_START = ", P_M, "\n")

    return P_M

#---------------PRESSURE FUCNTION AND DERIVATIVES----------------#
def PRESS_FUNC(P, RHO_K, P_K, A_K):

    if P <= P_K:
        #-------------RAREFACTION WAVE--------------#
        P_RATIO = P / P_K
        F = G4 * A_K * (P_RATIO**G1 - 1.0)
        F_DASH = (1.0 / (RHO_K * A_K)) * P_RATIO**(-G2)

    else:
        #-------------SHOCK WAVE---------------------#
        A_K = G5 / RHO_K
        B_K = G6 * P_K
        SQRT = np.sqrt(A_K / (B_K +P))
        F = (P - P_K) * SQRT
        F_DASH = (1.0 - 0.5 * (P - P_K) / (B_K + P)) * SQRT

    return F, F_DASH

#---------------SOLVER FOR P* & U*---------------------------#
def P_STAR_SOLVER(P_L, U_L, RHO_L, P_R, U_R, RHO_R, G1, G2, G3, G4, G5, G6, G7, G8, MAX_ITER, A_L, A_R):

    P_START = GUESS_P(P_L, U_L, RHO_L, P_R, U_R, RHO_R, G1, G2, G3, G4, G5, G6, G7, G8)
    P_OLD = P_START
    U_DIFF = U_R - U_L

    TOL = 1e-6

    for i in range(1, MAX_ITER+1):

        print("------------ ITERATION = ", i, " ---------------")
        
        F_L, F_L_DASH = PRESS_FUNC(P_OLD, RHO_L, P_L, A_L)
        print("=>   F_L      =  ", F_L)
        print("=>   F_L_DASH =  ", F_L_DASH)

        F_R, F_R_DASH = PRESS_FUNC(P_OLD, RHO_R, P_R, A_R)
        print("=>   F_R      =  ", F_R)
        print("=>   F_R_DASH =  ", F_R_DASH)

        P = P_OLD - ((F_L + F_R + U_DIFF) / (F_L_DASH + F_R_DASH))

        DELTA_P = 2 * np.abs((P - P_OLD) / (P + P_OLD))
        print("=>   RESIDUAL = ", DELTA_P, "\n")
        
        if DELTA_P <= TOL:
            break
        if P < 0.0:
            P = TOL

        P_OLD = P

    P_STAR = P
    #------------------VELOCITY IN STAR REGION---------------#
    U_STAR = 0.5 * (U_L + U_R + F_R - F_L)
    print("------P* = ", P_STAR, " ------")
    print("------U* = ", U_STAR, " ------\n")

    return P_STAR, U_STAR

#--------SAMPLING USING P* AND U* TO GET COMPLETE SOLUTION--------#
def SAMPLING_SOL(P_STAR, U_STAR, S, P_L, U_L, RHO_L, P_R, U_R, RHO_R, \
                                                 G1, G2, G3, G4, G5, G6, G7, G8, A_L, A_R):
    GAMMA = 1.4
    if (S <= U_STAR):
        #---LEFT OF CONTACT DISCONTINUITY---
        if (P_STAR <= P_L):
            #---LEFT RAREFACTION---
            S_HL = U_L - A_L            # SPEED OF HEAD OF RAREFACTION
            if (S <= S_HL):
                #---LEFT DATA STATE---
                RHO = RHO_L
                U = U_L
                P = P_L
            else:
                A_ML = A_L * (P_STAR / P_L)**G1
                S_TL = U_STAR - A_ML    # SPEED OF TAIL OF RAREFACTION
                if (S > S_TL):
                    #---LEFT STAR STATE---
                    RHO = RHO_L * (P_STAR / P_L)**(1/GAMMA)
                    U = U_STAR
                    P = P_STAR
                else:
                    #---INSIDE LEFT FAN---
                    U = G5 * (A_L + G7 * U_L + S)
                    A = G5 * (A_L + G7 * (U_L - S))
                    RHO = RHO_L * (A / A_L)**G4
                    P = P_L * (A / A_L)**G3
            
        else:
            #---LEFT SHOCK---
            P_STARL = P_STAR / P_L
            S_L = U_L - A_L * np.sqrt(G2 * P_STARL + G1)
            if (S <= S_L):
                #---LEFT DATA STATE---
                RHO = RHO_L
                U = U_L
                P = P_L
            else:
                #---LEFT STAR STATE---
                RHO = RHO_L * (P_STARL + G6) / (P_STARL * G8 + 1.0)
                U = U_STAR
                P = P_STAR

    else:
        #---RIGHT OF THE CONTACT DISCONTINUITY---
        if (P_STAR > P_R):
            #---RIGHT SHOCK---
            P_STARR = P_STAR / P_R
            S_R = U_R + A_R * np.sqrt(G2 * P_STARR + G1)        # SHOCK SPEED
            if (S >= S_R):
                #---RIGHT DATA STATE---
                RHO = RHO_R
                U = U_R
                P = P_R
            else:
                #---RIGHT STAR STATE
                RHO = RHO_R * (P_STARR + G6) / (P_STARR * G6 + 1.0)
                U = U_STAR
                P = P_STAR
        
        else:
            #--- RIGHT RAREFACTION---
            S_HR = U_R + A_R
            if (S >= S_HR):
                #---RIGHT DATA STATE---
                RHO = RHO_R
                U = U_R
                P = P_R
            else:
                A_MR = A_R * (P_STAR / P_R)**G1
                S_TR = U_STAR + A_MR
                if (S <= S_TR):
                    #---RIGHT STAR STATE---
                    RHO = RHO_R * (P_STAR / P_R)**(1.0 / GAMMA)
                    U = U_STAR
                    P = P_STAR
                else:
                    #---INSIDE LEFT FAN---
                    U = G5 * (-A_R + G7 * U_R + S)
                    A = G5 * (A_R - G7 * (U_R - S))
                    RHO = RHO_R * (A / A_R)**G4
                    P = P_R * (A / A_R)**G3
                
    return RHO, U, P

#%%
#------------VARIABLES-------------------------------#
D_LEN = 2.0
DIA_POS = 1
N_CELL = 200
GAMMA = 1.4
SIM_TIME = 0.001


STATE_L = np.array([2.0, 0.0, 200000])
STATE_R = np.array([1, 0.0, 100000])


RHO_L = STATE_L[0]
U_L = STATE_L[1]
P_L = STATE_L[2]

RHO_R = STATE_R[0]
U_R = STATE_R[1]
P_R = STATE_R[2]

MAX_ITER = 25

#------------GAMMA RELATED CONSTANTS-----------------#
G1 = (GAMMA - 1.0) / (2.0 * GAMMA)
G2 = (GAMMA + 1.0) / (2.0 * GAMMA)
G3 = 2.0 * GAMMA / (GAMMA - 1.0)
G4 = 2.0 / (GAMMA - 1.0)
G5 = 2.0 / (GAMMA + 1.0)
G6 = (GAMMA - 1.0) / (GAMMA + 1.0)
G7 = (GAMMA - 1.0) / 2.0
G8 = GAMMA - 1.0

#-------------COMPUTE SONIC SPEEDS (RIGHT & LEFT)------------------------#
A_L = np.sqrt(GAMMA * P_L / RHO_L)
A_R = np.sqrt(GAMMA * P_R / RHO_R)

#------------PRESSURE POSITIVITY CONDITION CHECK------------------------#
if (G4 * (A_L + A_R) <= (U_R - U_L)):
    print("!- INITIAL DATA IS SUCH THAT VACUUM IS GENERATED -!")
    sys.exit()

#---------EXACT SOLUTION FOR PRESSURE AND VELOCITY IN STAR REGION----------#
P_STAR, U_STAR = P_STAR_SOLVER(P_L, U_L, RHO_L, P_R, U_R, RHO_R, G1, G2, G3, G4, G5, G6, G7, G8, MAX_ITER, A_L, A_R)

#%%
DENSITY = []
PRESSURE = []
VELOCITY_X = []
INT_ENERGY = []

DX = D_LEN / N_CELL

for i in range(1, N_CELL+1):

    XPOS = (i - 0.5) * DX
    S = (XPOS - DIA_POS) / SIM_TIME

    RHO, U, P = SAMPLING_SOL(P_STAR, U_STAR, S, P_L, U_L, RHO_L, P_R, U_R, RHO_R, \
                                                 G1, G2, G3, G4, G5, G6, G7, G8, A_L, A_R)
    E = P / (RHO * (GAMMA - 1.0))

    #------------TRACKING SOLUTIONS-----------------
    # print("\n=====At location x = ", XPOS,"=====")
    # print("-------DENSITY  = ", RHO, "-------")
    # print("-------VELOCITY = ", U, "-------")
    # print("-------PRESSURE = ", P, "-------")
    #-----------------------------------------------
    DENSITY.append(RHO)
    VELOCITY_X.append(U)
    PRESSURE.append(P)
    INT_ENERGY.append(E)

# %%
#----------------PLOTTING THE RESULTS------------------

X_POS_PLOT = np.zeros(N_CELL)
for i in range(1, N_CELL+1):
    X_POS_PLOT[i-1] = (i - 0.5) * DX
plt.figure(figsize=(12, 8))

data = np.column_stack((X_POS_PLOT, DENSITY, VELOCITY_X, PRESSURE))
np.savetxt('exact.dat', data, header='X_POS_PLOT DENSITY VELOCITY_X PRESSURE', comments='')