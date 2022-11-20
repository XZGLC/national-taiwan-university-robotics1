import numpy as np
import math

def get_time_space_info(info):

    spatial_info = info[:, 1:].copy()
    time_info = info[:,0:1].copy()

    return time_info, spatial_info

def parabolic_blend(info, accelration_time=0.5):

    time_info, spatial_info = get_time_space_info(info) # rx1 , rxc

    #time
    time_interval = time_info[1:] - time_info[:-1]
    time_interval[0] = time_interval[0] - accelration_time/2
    time_interval[-1] = time_interval[-1] - accelration_time/2

    #velocity
    velocity = (spatial_info[1:,:] - spatial_info[:-1,:]) / time_interval
    velocity = np.concatenate([
        np.zeros((1,spatial_info.shape[-1])),
        velocity,
        np.zeros((1,spatial_info.shape[-1]))
    ], axis=0)

    #acceleration
    acceleration = (velocity[1:] - velocity[:-1]) / accelration_time

    return time_interval, velocity, acceleration

if __name__ == "__main__":

    np.set_printoptions(precision=4, suppress=True)
    degree = True
    path_table = np.array(
        [ [0, -4, 0, 120/180*math.pi],
        [2, -5, 5, 45/180*math.pi],
        [4, 2, 3, 30/180*math.pi],
        [9, 5, -3, 0] ]
    )

    time_interval, velocity, acceleration = parabolic_blend(path_table)

    print("time interval:")
    print(time_interval)
    print("Velocity:")
    if degree:
        velocity[:,-1] = velocity[:,-1]/math.pi*180
    print(velocity)
    print("Acceleration:")
    if degree:
        acceleration[:,-1] = acceleration[:,-1]/math.pi*180
    print(acceleration)

    print("##################################################")
    print('\n')

    #Q8
    t = 4
    time_info, spatial_info = get_time_space_info(path_table)
    spatial_info[:,-1] = spatial_info[:,-1]/math.pi*180
    xytheta = spatial_info[1:2,:] + velocity[2:3,:]*(t-2) + 0.5*acceleration[2:3,:]*np.square(t-3.75)
    print("Q8 Ans:")
    print(xytheta)
    print("##################################################")
    print('\n')

    #checking
    path_table = np.array(
        [ [0, -4, 0, 90/180*math.pi],
        [2, 0, 3, 45/180*math.pi],
        [4, 3, 3, 30/180*math.pi],
        [7, 4, 0, 0] ]
    )

    time_interval, velocity, acceleration = parabolic_blend(path_table)

    print("time interval:")
    print(time_interval)
    print("Velocity:")
    if degree:
        velocity[:,-1] = velocity[:,-1]/math.pi*180
    print(velocity)
    print("Acceleration:")
    if degree:
        acceleration[:,-1] = acceleration[:,-1]/math.pi*180
    print(acceleration)

    print("##################################################")
    print('\n')

   



