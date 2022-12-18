import numpy as  np
from base import solve_theta1, solve_theta2, solve_theta3, inv_bTa, solve_theta4_theta5_theta6, T_rotate_along_X, T_rotate_along_Y, T_rotate_along_Z
import math

def rot_xyz_fix_angle(thetaX, thetaY, thetaZ):
    
    return np.matmul(T_rotate_along_Z(thetaZ), np.matmul(T_rotate_along_Y(thetaY), T_rotate_along_X(thetaX)) )

def T_xyz_fix_angle(X, Y, Z, thetaX, thetaY, thetaZ):
    
    _T = rot_xyz_fix_angle(thetaX, thetaY, thetaZ)
    _T[0,-1], _T[1,-1], _T[2,-1] = X, Y, Z
    
    return _T.astype(np.float32)



if __name__ == "__main__":
    np.set_printoptions(precision=5, suppress=True)
    
    _6TC = np.array([
        [0,0,1,0],
        [0,-1,0,0],
        [1,0,0,206],
        [0,0,0,1]
    ]).astype(np.float32)
    
    # _0Tc0 = T_xyz_fix_angle(550,270,79.5,0,0, 35/180*math.pi)
    # _0T6_p = np.matmul(_0Tc0, inv_bTa(_6TC))
    
    #Q1
    _0Tc0 = T_xyz_fix_angle(630,364,20,0,0, 0/180*math.pi)
    _0T6_p0 = np.matmul(_0Tc0, inv_bTa(_6TC))
    print("Q1 Ans: ")
    print(_0T6_p0)
    print("\n")
    
    #Q2
    _0Tc1 = T_xyz_fix_angle(630,304,220,60/180*math.pi, 0, 0)
    _0T6_p1 = np.matmul(_0Tc1, inv_bTa(_6TC))
    print("Q2 Ans: ")
    print(_0T6_p1)
    print("\n")
    
    #Q3
    _0Tc2 = T_xyz_fix_angle(630,220,24,180/180*math.pi, 0, 0)
    _0T6_p2 = np.matmul(_0Tc2, inv_bTa(_6TC))
    print("Q3 Ans: ")
    print(_0T6_p2)
    print("\n")
    
    ################################################
    #DH table
    alpha0, a0, d1, theta1 = 0,0,0,None
    alpha1, a1, d2, theta2 = -math.pi/2,-30,0,None
    alpha2, a2, d3, theta3 = 0,340,0,None
    alpha3, a3, d4, theta4 = -math.pi/2,-40,338,None
    alpha4, a4, d5, theta5 = math.pi/2,0,0,None
    alpha5, a5, d6, theta6 = -math.pi/2,0,0,None
    
    # _0T6 = np.matmul(T_xyz_fix_angle(330,372,367,0,-60/180*math.pi, 0), inv_bTa(_6TC))
    # print(_0T6)
    # px, py, pz = _0T6[0,-1], _0T6[1,-1], _0T6[2,-1]
    # theta1_v1, theta1_v2 = solve_theta1(px=px, py=py, d2=d2, d3=d3)
    # print(f"theta1: {theta1_v1/math.pi*180} or {theta1_v2/math.pi*180}")
    # print('\n')

    # theta3_v1, theta3_v2 = solve_theta3(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d2=d2, d3=d3, d4=d4, theta1=theta1_v1)
    # theta3_v1 = theta3_v1 - 2*math.pi
    # theta3_v2 = theta3_v2 + 2*math.pi
    # print(f"theta3: {theta3_v1/math.pi*180} or {theta3_v2/math.pi*180}")
    
    #Q4-5
    print("Q4-5")
    px, py, pz = _0T6_p0[0,-1], _0T6_p0[1,-1], _0T6_p0[2,-1]
    theta1_v1, theta1_v2 = solve_theta1(px=px, py=py, d2=d2, d3=d3)
    print(f"theta1: {theta1_v1/math.pi*180} or {theta1_v2/math.pi*180}")
    print('\n')
    
    theta1 = theta1_v1

    theta3_v1, theta3_v2 = solve_theta3(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d2=d2, d3=d3, d4=d4, theta1=theta1)
    print(f"theta3: {theta3_v1/math.pi*180} or {theta3_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2_v1 = solve_theta2(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v1)
    theta2_v2 = solve_theta2(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v2)
    print(f"theta2: {theta2_v1/math.pi*180} or {theta2_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2, theta3 = theta2_v1, theta3_v1
    
    theta4, theta5, theta6 = solve_theta4_theta5_theta6(alpha0=alpha0, 
                                                        a0=a0, 
                                                        d1=d1, 
                                                        theta1=theta1, 
                                                        alpha1=alpha1, 
                                                        a1=a1, 
                                                        d2=d2, 
                                                        theta2=theta2, 
                                                        alpha2=alpha2, 
                                                        a2=a2, 
                                                        d3=d3, 
                                                        theta3=theta3,
                                                        alpha3=alpha3, 
                                                        a3=a3, 
                                                        d4=d4, 
                                                        _0T6=_0T6_p0.copy()
                                                        )

    print(f"theta4: {theta4/math.pi*180} , theta5: {theta5/math.pi*180} , theta6: {theta6/math.pi*180}")
    print("##################################################")
    print('\n')
    
    #Q6-7
    print("Q6-7")
    px, py, pz = _0T6_p1[0,-1], _0T6_p1[1,-1], _0T6_p1[2,-1]
    theta1_v1, theta1_v2 = solve_theta1(px=px, py=py, d2=d2, d3=d3)
    print(f"theta1: {theta1_v1/math.pi*180} or {theta1_v2/math.pi*180}")
    print('\n')
    
    theta1 = theta1_v1

    theta3_v1, theta3_v2 = solve_theta3(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d2=d2, d3=d3, d4=d4, theta1=theta1)
    print(f"theta3: {theta3_v1/math.pi*180} or {theta3_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2_v1 = solve_theta2(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v1)
    theta2_v2 = solve_theta2(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v2)
    print(f"theta2: {theta2_v1/math.pi*180} or {theta2_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2, theta3 = theta2_v1, theta3_v1
    
    theta4, theta5, theta6 = solve_theta4_theta5_theta6(alpha0=alpha0, 
                                                        a0=a0, 
                                                        d1=d1, 
                                                        theta1=theta1, 
                                                        alpha1=alpha1, 
                                                        a1=a1, 
                                                        d2=d2, 
                                                        theta2=theta2, 
                                                        alpha2=alpha2, 
                                                        a2=a2, 
                                                        d3=d3, 
                                                        theta3=theta3,
                                                        alpha3=alpha3, 
                                                        a3=a3, 
                                                        d4=d4, 
                                                        _0T6=_0T6_p1.copy()
                                                        )

    print(f"theta4: {theta4/math.pi*180} , theta5: {theta5/math.pi*180} , theta6: {theta6/math.pi*180}")
    print("##################################################")
    print('\n')
    
    #Q8-9
    print("Q8-9")
    px, py, pz = _0T6_p2[0,-1], _0T6_p2[1,-1], _0T6_p2[2,-1]
    theta1_v1, theta1_v2 = solve_theta1(px=px, py=py, d2=d2, d3=d3)
    print(f"theta1: {theta1_v1/math.pi*180} or {theta1_v2/math.pi*180}")
    print('\n')
    
    theta1 = theta1_v1

    theta3_v1, theta3_v2 = solve_theta3(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d2=d2, d3=d3, d4=d4, theta1=theta1)
    theta3_v1 = theta3_v1 - 2*math.pi
    print(f"theta3: {theta3_v1/math.pi*180} or {theta3_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2_v1 = solve_theta2(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v1)
    theta2_v2 = solve_theta2(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v2)
    print(f"theta2: {theta2_v1/math.pi*180} or {theta2_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2, theta3 = theta2_v1, theta3_v1
    
    theta4, theta5, theta6 = solve_theta4_theta5_theta6(alpha0=alpha0, 
                                                        a0=a0, 
                                                        d1=d1, 
                                                        theta1=theta1, 
                                                        alpha1=alpha1, 
                                                        a1=a1, 
                                                        d2=d2, 
                                                        theta2=theta2, 
                                                        alpha2=alpha2, 
                                                        a2=a2, 
                                                        d3=d3, 
                                                        theta3=theta3,
                                                        alpha3=alpha3, 
                                                        a3=a3, 
                                                        d4=d4, 
                                                        _0T6=_0T6_p2.copy()
                                                        )

    print(f"theta4: {theta4/math.pi*180} , theta5: {theta5/math.pi*180} , theta6: {theta6/math.pi*180}")
    print("##################################################")
    print('\n')
    
    #Q10 - Q14
    setup_table = np.array([
        [0,630,364,20,0,0,0],
        [3,630,304,220,60,0,0],
        [7,630,220,24,180,0,0]
    ]).astype(np.float32)
    
    vel_table = setup_table[1:,:] - setup_table[:-1,:]
    vel_table[0,0] = vel_table[0,0] - 0.5/2
    vel_table[-1,0] = vel_table[-1,0] - 0.5/2
    vel_table = vel_table[:,1:] / vel_table[:,0:1]
    vel_table  = np.concatenate([np.zeros((1,6)),vel_table],axis=0)
    vel_table  = np.concatenate([vel_table,np.zeros((1,6))],axis=0)
    print("Velocity table: ")
    print(vel_table)
    
    acc_table = (vel_table[1:,:] - vel_table[:-1,:])/0.5
    print("Acceleration table: ")
    print(acc_table)
    
    #Q15
    print("Q15: ")
    pos = setup_table[1:2,1:] + vel_table[2:3, :]*(5-3)
    print(pos)