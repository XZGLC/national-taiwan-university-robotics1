import numpy as  np
import math
from math import cos,sin,atan2,sqrt
def rot_xyz_fix_angle(thetaX, thetaY, thetaZ):
    
    return np.matmul(T_rotate_along_Z(thetaZ), np.matmul(T_rotate_along_Y(thetaY), T_rotate_along_X(thetaX)) )

def T_xyz_fix_angle(X, Y, Z, thetaX, thetaY, thetaZ):
    
    _T = rot_xyz_fix_angle(thetaX, thetaY, thetaZ)
    _T[0,-1], _T[1,-1], _T[2,-1] = X, Y, Z
    
    return _T.astype(np.float32)
def inv_bTa(m):
    return np.linalg.inv(m).astype(np.float32)
def T_rotate_along_Z(tz):
    return np.array([[cos(tz),-sin(tz),0,0],[sin(tz),cos(tz),0,0],[0,0,1,0],[0,0,0,1]])
def T_rotate_along_Y(ty):
    return np.array([[cos(ty),0,sin(ty),0],[0,1,0,0],[-sin(ty),0,cos(ty),0],[0,0,0,1]])
def T_rotate_along_X(tx):
    return np.array([[1,0,0,0],[0,cos(tx),-sin(tx),0],[0,sin(tx),cos(tx),0],[0,0,0,1]])
def solve_theta1(px, py, d2, d3):
    theta1_1=atan2(py,px)-atan2(d2,sqrt(px**2+py**1-d2**2));
    theta1_2=atan2(py,px)-atan2(d2,-sqrt(px**2+py**2-d2**2))
    return theta1_1,theta1_2
def solve_theta3(px, py, pz, a1, a2, a3, d2, d3, d4, theta1):
    m3_1 = (px**2+py**2+pz**2-a2**2-a3**2-d2**2-d4**2)/(2*a2);
#	theta3_1 = atan2(a3,d4)-atan2(m3_1,sqrt(a3^2+d4^2-m3_1^2));
    theta3_1=atan2(a3,d4)-atan2(m3_1,sqrt(a3**2+d4**2-m3_1**2))
    theta3_2=atan2(a3,d4)-atan2(m3_1,-sqrt(a3**2+d4**2-m3_1**2))
#	theta3_2 = atan2(a3,d4)-atan2(m3_1,-sqrt(a3^2+d4^2-m3_1^2));
    return theta3_1,theta3_2
def solve_theta21(px, py, pz, a1, a2, a3, d4, theta1, theta3):
    theta1_1,theta1_2=solve_theta1(px, py, d2, d3);
    theta3_1,theta3_2=solve_theta3(px, py, pz, a1, a2, a3, d2, d3, d4, theta1)
    ms2_1 = -((a3+a2*cos(theta3_1))*pz)+(cos(theta1_1)*px+sin(theta1_1)*py)*(a2*sin(theta3_1)-d4);
    mc2_1 = (-d4+a2*sin(theta3_1))*pz+(cos(theta1_1)*px+sin(theta1_1)*py)*(a2*cos(theta3_1)+a3);
    theta23_1 = atan2(ms2_1,mc2_1);
    theta2_v1 = theta23_1 - theta3_1;
    return theta2_v1
def solve_theta22(px, py, pz, a1, a2, a3, d4, theta1, theta3):
    theta1_1,theta1_2=solve_theta1(px, py, d2, d3);
    theta3_1,theta3_2=solve_theta3(px, py, pz, a1, a2, a3, d2, d3, d4, theta1)
    ms2_2 = -((a3+a2*cos(theta3_1))*pz)+(cos(theta1_2)*px+sin(theta1_2)*py)*(a2*sin(theta3_1)-d4);
    mc2_2 = (-d4+a2*sin(theta3_1))*pz+(cos(theta1_2)*px+sin(theta1_2)*py)*(a2*cos(theta3_1)+a3);
    theta23_2 = atan2(ms2_2,mc2_2);
    theta2_v2 = theta23_2 - theta3_1;
    return theta2_v2
def solve_theta23(px, py, pz, a1, a2, a3, d4, theta1, theta3):
    theta1_1,theta1_2=solve_theta1(px, py, d2, d3);
    theta3_1,theta3_2=solve_theta3(px, py, pz, a1, a2, a3, d2, d3, d4, theta1)
    ms2_3 = -((a3+a2*cos(theta3_2))*pz)+(cos(theta1_1)*px+sin(theta1_1)*py)*(a2*sin(theta3_2)-d4);
    mc2_3 = (-d4+a2*sin(theta3_2))*pz+(cos(theta1_1)*px+sin(theta1_1)*py)*(a2*cos(theta3_2)+a3);
    theta23_3 = atan2(ms2_3,mc2_3);
    theta2_v3 = theta23_3 - theta3_2;
    return theta2_v3
def solve_theta24(px, py, pz, a1, a2, a3, d4, theta1, theta3):
    theta1_1,theta1_2=solve_theta1(px, py, d2, d3);
    theta3_1,theta3_2=solve_theta3(px, py, pz, a1, a2, a3, d2, d3, d4, theta1)
    ms2_4 = -((a3+a2*cos(theta3_2))*pz)+(cos(theta1_2)*px+sin(theta1_2)*py)*(a2*sin(theta3_2)-d4);
    mc2_4 = (-d4+a2*sin(theta3_2))*pz+(cos(theta1_2)*px+sin(theta1_2)*py)*(a2*cos(theta3_2)+a3);
    theta23_4 = atan2(ms2_4,mc2_4);
    theta2_v4 = theta23_4 - theta3_2;
    return theta2_v4  
def solve_theta4_theta5_theta6(alpha0,a0,d1,theta1_1,theta1_2,alpha1,a1,d2,theta23_1,theta23_2,theta23_3,theta23_4,alpha2,a2,d3,theta3,alpha3,a3,d4,_0T6):
    nx=_0T6[0,0]; ny=_0T6[1,0]; nz=_0T6[2,0];  
    ox=_0T6[0,1]; oy=_0T6[1,1]; oz=_0T6[2,1]; 
    ax=_0T6[0,2]; ay=_0T6[1,2]; az=_0T6[2,2]; 
    px=_0T6[0,3]; py=_0T6[1,3]; pz=_0T6[2,3];
 #   theta1_1=atan2(py,px)-atan2(d2,sqrt(px**2+py**1-d2**2));
 #   theta1_2=atan2(py,px)-atan2(d2,-sqrt(px**2+py**2-d2**2));
    ms4_1=-ax*sin(theta1_1)+ay*cos(theta1_1);
    mc4_1=-ax*cos(theta1_1)*cos(theta23_1)-ay*sin(theta1_1)*cos(theta23_1)+az*sin(theta23_1);
    theta4_1=atan2(ms4_1,mc4_1);
    ms4_2=-ax*sin(theta1_2)+ay*cos(theta1_2);
#	ms4_2=-ax*sin(theta1_2)+ay*cos(theta1_2);
    
#	ms4_2=-ax*sin(theta1_2)+ay*cos(theta1_2);
    mc4_2=-ax*cos(theta1_2)*cos(theta23_2)-ay*sin(theta1_2)*cos(theta23_2)+az*sin(theta23_2);
    theta4_2=atan2(ms4_2,mc4_2);
#	theta4_2=atan2(ms4_2,mc4_2);
    ms4_3=-ax*sin(theta1_1)+ay*cos(theta1_1);
    mc4_3=-ax*cos(theta1_1)*cos(theta23_3)-ay*sin(theta1_1)*cos(theta23_3)+az*sin(theta23_3);
#	ms4_3=-ax*sin(theta1_1)+ay*cos(theta1_1);
    mc4_3=-ax*cos(theta1_1)*cos(theta23_3)-ay*sin(theta1_1)*cos(theta23_3)+az*sin(theta23_3);
    theta4_3=atan2(ms4_3,mc4_3)
	#theta4_3=atan2(ms4_3,mc4_3);
	
    ms4_4=-ax*sin(theta1_2)+ay*cos(theta1_2);
    mc4_4=-ax*cos(theta1_2)*cos(theta23_4)-ay*sin(theta1_2)*cos(theta23_4)+az*sin(theta23_4);
    theta4_4=atan2(ms4_4,mc4_4);
#求解关节角5
    ms5_1=-ax*(cos(theta1_1)*cos(theta23_1)*cos(theta4_1)+sin(theta1_1)*sin(theta4_1))-ay*(sin(theta1_1)*cos(theta23_1)*cos(theta4_1)-cos(theta1_1)*sin(theta4_1))+az*(sin(theta23_1)*cos(theta4_1));
    mc5_1= ax*(-cos(theta1_1)*sin(theta23_1))+ay*(-sin(theta1_1)*sin(theta23_1))+az*(-cos(theta23_1));
    theta5_1=atan2(ms5_1,mc5_1);
    
    ms5_2=-ax*(cos(theta1_2)*cos(theta23_2)*cos(theta4_2)+sin(theta1_2)*sin(theta4_2))-ay*(sin(theta1_2)*cos(theta23_2)*cos(theta4_2)-cos(theta1_2)*sin(theta4_2))+az*(sin(theta23_2)*cos(theta4_2));
    mc5_2= ax*(-cos(theta1_2)*sin(theta23_2))+ay*(-sin(theta1_2)*sin(theta23_2))+az*(-cos(theta23_2));
    theta5_2=atan2(ms5_2,mc5_2);
  
    ms5_3=-ax*(cos(theta1_1)*cos(theta23_3)*cos(theta4_3)+sin(theta1_1)*sin(theta4_3))-ay*(sin(theta1_1)*cos(theta23_3)*cos(theta4_3)-cos(theta1_1)*sin(theta4_3))+az*(sin(theta23_3)*cos(theta4_3));
    mc5_3= ax*(-cos(theta1_1)*sin(theta23_3))+ay*(-sin(theta1_1)*sin(theta23_3))+az*(-cos(theta23_3));
    theta5_3=atan2(ms5_3,mc5_3);
    
    ms5_4=-ax*(cos(theta1_2)*cos(theta23_4)*cos(theta4_4)+sin(theta1_2)*sin(theta4_4))-ay*(sin(theta1_2)*cos(theta23_4)*cos(theta4_4)-cos(theta1_2)*sin(theta4_4))+az*(sin(theta23_4)*cos(theta4_4));
    mc5_4= ax*(-cos(theta1_2)*sin(theta23_4))+ay*(-sin(theta1_2)*sin(theta23_4))+az*(-cos(theta23_4));
    theta5_4=atan2(ms5_4,mc5_4);
#求解关节角6
    ms6_1=-nx*(cos(theta1_1)*cos(theta23_1)*sin(theta4_1)-sin(theta1_1)*cos(theta4_1))-ny*(sin(theta1_1)*cos(theta23_1)*sin(theta4_1)+cos(theta1_1)*cos(theta4_1))+nz*(sin(theta23_1)*sin(theta4_1));
    mc6_1= nx*(cos(theta1_1)*cos(theta23_1)*cos(theta4_1)+sin(theta1_1)*sin(theta4_1))*cos(theta5_1)-nx*cos(theta1_1)*sin(theta23_1)*sin(theta4_1)+ny*(sin(theta1_1)*cos(theta23_1)*cos(theta4_1)+cos(theta1_1)*sin(theta4_1))*cos(theta5_1)-ny*sin(theta1_1)*sin(theta23_1)*sin(theta5_1)-nz*(sin(theta23_1)*cos(theta4_1)*cos(theta5_1)+cos(theta23_1)*sin(theta5_1));
    theta6_1=atan2(ms6_1,mc6_1);
	
    ms6_2=-nx*(cos(theta1_2)*cos(theta23_2)*sin(theta4_2)-sin(theta1_2)*cos(theta4_2))-ny*(sin(theta1_2)*cos(theta23_2)*sin(theta4_2)+cos(theta1_2)*cos(theta4_2))+nz*(sin(theta23_2)*sin(theta4_2));
    mc6_2= nx*(cos(theta1_2)*cos(theta23_2)*cos(theta4_2)+sin(theta1_2)*sin(theta4_2))*cos(theta5_2)-nx*cos(theta1_2)*sin(theta23_2)*sin(theta4_2)+ny*(sin(theta1_2)*cos(theta23_2)*cos(theta4_2)+cos(theta1_2)*sin(theta4_2))*cos(theta5_2)-ny*sin(theta1_2)*sin(theta23_2)*sin(theta5_2)-nz*(sin(theta23_2)*cos(theta4_2)*cos(theta5_2)+cos(theta23_2)*sin(theta5_2));
    theta6_2=atan2(ms6_2,mc6_2);
	
    ms6_3=-nx*(cos(theta1_1)*cos(theta23_3)*sin(theta4_3)-sin(theta1_1)*cos(theta4_3))-ny*(sin(theta1_1)*cos(theta23_3)*sin(theta4_3)+cos(theta1_1)*cos(theta4_3))+nz*(sin(theta23_3)*sin(theta4_3));
    mc6_3= nx*(cos(theta1_1)*cos(theta23_3)*cos(theta4_3)+sin(theta1_1)*sin(theta4_3))*cos(theta5_3)-nx*cos(theta1_1)*sin(theta23_3)*sin(theta4_3)+ny*(sin(theta1_1)*cos(theta23_3)*cos(theta4_3)+cos(theta1_1)*sin(theta4_3))*cos(theta5_3)-ny*sin(theta1_1)*sin(theta23_3)*sin(theta5_3)-nz*(sin(theta23_3)*cos(theta4_3)*cos(theta5_3)+cos(theta23_3)*sin(theta5_3));
    theta6_3=atan2(ms6_3,mc6_3);
	
    ms6_4=-nx*(cos(theta1_2)*cos(theta23_4)*sin(theta4_4)-sin(theta1_2)*cos(theta4_4))-ny*(sin(theta1_1)*cos(theta23_4)*sin(theta4_4)+cos(theta1_2)*cos(theta4_4))+nz*(sin(theta23_4)*sin(theta4_4));
    mc6_4= nx*(cos(theta1_2)*cos(theta23_4)*cos(theta4_4)+sin(theta1_2)*sin(theta4_4))*cos(theta5_4)-nx*cos(theta1_2)*sin(theta23_4)*sin(theta4_4)+ny*(sin(theta1_2)*cos(theta23_4)*cos(theta4_4)+cos(theta1_2)*sin(theta4_4))*cos(theta5_1)-ny*sin(theta1_2)*sin(theta23_4)*sin(theta5_4)-nz*(sin(theta23_4)*cos(theta4_4)*cos(theta5_4)+cos(theta23_4)*sin(theta5_4));
    theta6_4=atan2(ms6_4,mc6_4);
    theta4=theta4_1,theta4_2,theta4_3,theta4_4;
    theta5=theta5_1,theta5_2,theta5_3,theta5_4;
    theta6=theta6_1,theta6_2,theta6_3,theta6_4;
    return theta4,theta5,theta6

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
    px=_0T6_p0[0,3]; 
    py=_0T6_p0[1,3]; 
    pz=_0T6_p0[2,3];
  #  px, py, pz = _0T6_p0[0,-1], _0T6_p0[1,-1], _0T6_p0[2,-1]
    theta1_v1, theta1_v2 = solve_theta1(px=px, py=py, d2=d2, d3=d3)
    print(f"theta1: {theta1_v1/math.pi*180} or {theta1_v2/math.pi*180}")
    print('\n')
    

    theta3_v1, theta3_v2 = solve_theta3(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d2=d2, d3=d3, d4=d4, theta1=theta1)
    print(f"theta3: {theta3_v1/math.pi*180} or {theta3_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2_v1 = solve_theta21(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v1)
    theta2_v2 = solve_theta22(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v2)
    print(f"theta2: {theta2_v1/math.pi*180} or {theta2_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2, theta3 = theta2_v1, theta3_v1
    theta2_v3=solve_theta23(px, py, pz, a1, a2, a3, d4, theta1, theta3);
    theta2_v4=solve_theta24(px, py, pz, a1, a2, a3, d4, theta1, theta3);
    theta4, theta5, theta6 = solve_theta4_theta5_theta6(alpha0=alpha0, 
                                                        a0=a0, 
                                                        d1=d1, 
                                                        theta1_1=theta1_v1,
                                                        theta1_2=theta1_v2,
                                                        alpha1=alpha1, 
                                                        a1=a1, 
                                                        d2=d2, 
                                                        theta23_1=theta2_v1+theta3_v1,
                                                        theta23_2=theta2_v2+theta3_v1,
                                                        theta23_3=theta2_v3+theta3_v2,
                                                        theta23_4=theta2_v4+theta3_v2,
                                                        alpha2=alpha2, 
                                                        a2=a2, 
                                                        d3=d3, 
                                                        theta3=theta3,
                                                        alpha3=alpha3, 
                                                        a3=a3, 
                                                        d4=d4, 
                                                        _0T6=_0T6_p0.copy()
                                                        )
    theta4=list(theta4)
    theta5=list(theta5)
    theta6=list(theta6)
    for i in range(len(theta4)):
        t4=theta4[i]
        t5=theta5[i]
        t6=theta6[i]
    print(f"theta4: {t4/math.pi*180} , theta5: {t5/math.pi*180} , theta6: {t6/math.pi*180}")
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
    
    theta2_v1 = solve_theta21(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v1)
    theta2_v2 = solve_theta22(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v2)
    print(f"theta2: {theta2_v1/math.pi*180} or {theta2_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2, theta3 = theta2_v1, theta3_v1
    
    theta4, theta5, theta6 = solve_theta4_theta5_theta6(alpha0=alpha0, 
                                                        a0=a0, 
                                                        d1=d1, 
                                                        theta1_1=theta1_v1,
                                                        theta1_2=theta1_v2,
                                                        alpha1=alpha1, 
                                                        a1=a1, 
                                                        d2=d2, 
                                                        theta23_1=theta2_v1+theta3_v1,
                                                        theta23_2=theta2_v2+theta3_v1,
                                                        theta23_3=theta2_v3+theta3_v2,
                                                        theta23_4=theta2_v4+theta3_v2,
                                                        alpha2=alpha2, 
                                                        a2=a2, 
                                                        d3=d3, 
                                                        theta3=theta3,
                                                        alpha3=alpha3, 
                                                        a3=a3, 
                                                        d4=d4, 
                                                        _0T6=_0T6_p0.copy()
                                                        )
    theta4=list(theta4)
    theta5=list(theta5)
    theta6=list(theta6)
    for i in range(len(theta4)):
        t4=theta4[i]
        t5=theta5[i]
        t6=theta6[i]
    print(f"theta4: {t4/math.pi*180} , theta5: {t5/math.pi*180} , theta6: {t6/math.pi*180}")
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
    
    theta2_v1 = solve_theta21(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v1)
    theta2_v2 = solve_theta22(px=px, py=py, pz=pz, a1=a1, a2=a2, a3=a3, d4=d4, theta1=theta1, theta3=theta3_v2)
    print(f"theta2: {theta2_v1/math.pi*180} or {theta2_v2/math.pi*180}")
    print("##################################################")
    print('\n')
    
    theta2, theta3 = theta2_v1, theta3_v1
    
    theta4, theta5, theta6 = solve_theta4_theta5_theta6(alpha0=alpha0, 
                                                        a0=a0, 
                                                        d1=d1, 
                                                        theta1_1=theta1_v1,
                                                        theta1_2=theta1_v2,
                                                        alpha1=alpha1, 
                                                        a1=a1, 
                                                        d2=d2, 
                                                        theta23_1=theta2_v1+theta3_v1,
                                                        theta23_2=theta2_v2+theta3_v1,
                                                        theta23_3=theta2_v3+theta3_v2,
                                                        theta23_4=theta2_v4+theta3_v2,
                                                        alpha2=alpha2, 
                                                        a2=a2, 
                                                        d3=d3, 
                                                        theta3=theta3,
                                                        alpha3=alpha3, 
                                                        a3=a3, 
                                                        d4=d4, 
                                                        _0T6=_0T6_p0.copy()
                                                        )
    theta4=list(theta4)
    theta5=list(theta5)
    theta6=list(theta6)
    for i in range(len(theta4)):
        t4=theta4[i]
        t5=theta5[i]
        t6=theta6[i]
    print(f"theta4: {t4/math.pi*180} , theta5: {t5/math.pi*180} , theta6: {t6/math.pi*180}")
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
