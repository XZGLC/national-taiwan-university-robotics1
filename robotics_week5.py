import numpy as np
import math

def solve_theta1_theta2_theta3(l1,l2,x,y,phi):

    a = l1**2 + l2**2 - (x**2 + y**2)
    b = -2*l1*l2
    theta2 = np.arccos( a/b )

    a = l1**2 + (x**2 + y**2) - l2**2
    b = 2*l1*np.sqrt(x**2+y**2)
    gamma = np.arccos( a/b )

    if theta2 < 0:
        theta1 = np.arctan2(y,x) + gamma
    else:
        theta1 = np.arctan2(y,x) - gamma

    theta3 = phi - theta1 - theta2

    return theta1, theta2, theta3


if __name__ == "__main__":

    delta_t1 = 2
    delta_t2 = 2
    delta_t3 = 5 

    T12_mat = np.array(
        [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, delta_t1, delta_t1**2, delta_t1**3, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, delta_t2, delta_t2**2, delta_t2**3, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, delta_t3, delta_t3**2, delta_t3**3],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2*delta_t3, 3*(delta_t3**2)],
        [0, 1, 2*delta_t1, 3*(delta_t1**2), 0, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 2, 6*delta_t1, 0, 0, -2, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 2*delta_t2, 3*(delta_t2**2), 0, -1, 0, 0],
        [0, 0, 0, 0, 0, 0, 2, 6*delta_t2, 0, 0, -2, 0]
        ]
    )

    np.set_printoptions(precision=4, suppress=True)
    print(T12_mat)
    print("###########################################################")
    print("\n")

    inv_T12 = np.linalg.inv(T12_mat)

    # theta = np.array(
    #     [[-4, 0, math.pi/2],
    #     [0, 3, math.pi/4],
    #     [0, 3, math.pi/4],
    #     [3, 3, math.pi/6],
    #     [3, 3, math.pi/6],
    #     [4, 0, 0],
    #     [0, 0, 0],
    #     [0, 0, 0],
    #     [0, 0, 0],
    #     [0, 0, 0],
    #     [0, 0, 0],
    #     [0, 0, 0]
    #     ]
    # )

    theta = np.array(
        [[-4, 0, 120],
        [-5, 5, 45],
        [-5, 5, 45],
        [2, 3, 30],
        [2, 3, 30],
        [2, -3, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
        ]
    )

    a = np.matmul(inv_T12, theta)

    print("ans:")
    print(a)
    print("###########################################################")
    print("\n")

    #get joint space mat
    l1, l2 = 5, 3
    theta1, theta2, theta3 = solve_theta1_theta2_theta3(l1=l1, l2=l2, x=-4, y=0, phi=120/180*math.pi)
    t0_row = [theta1, theta2, theta3]

    theta1, theta2, theta3 = solve_theta1_theta2_theta3(l1=l1, l2=l2, x=-5, y=5, phi=45/180*math.pi)
    t1_row = [theta1, theta2, theta3]

    theta1, theta2, theta3 = solve_theta1_theta2_theta3(l1=l1, l2=l2, x=2, y=3, phi=30/180*math.pi)
    t2_row = [theta1, theta2, theta3]

    theta1, theta2, theta3 = solve_theta1_theta2_theta3(l1=l1, l2=l2, x=2, y=-3, phi=0/180*math.pi)
    t3_row = [theta1, theta2, theta3]

    theta = np.array(
        [t0_row,
        t1_row,
        t1_row,
        t2_row,
        t2_row,
        t3_row,
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
        ]
    )

    a = np.matmul(inv_T12, theta)

    print("ans:")
    print(a)
    print("###########################################################")
    print("\n")


