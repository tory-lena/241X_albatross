#This file reading data

import numpy as np
import matplotlib.pyplot as plt
def RotatioMatrix(roll,pitch,yaw):
    '''
    compute rotation matrix, R, project from NED to body framework
    :param row:
    :param pitch:
    :param yaw:
    :return:
    '''
    Rz = np.array([[np.cos(yaw),-np.sin(yaw),0.0],[np.sin(yaw),np.cos(yaw),0.0],[0,0,1]])
    Ry = np.array([[np.cos(pitch),0.0,np.sin(pitch)],[0,1.0,0],[-np.sin(pitch), 0, np.cos(pitch)]])
    Rx = np.array([[1, 0.0, 0.0],[0, np.cos(roll), - np.sin(roll)],[0, np.sin(roll),np.cos(roll)]])
    return np.dot(Rz, np.dot(Ry,Rx))
    #return np.dot(Rx, np.dot(Ry,Rz))

def computeDragLift(IMU_AccX,IMU_AccY,IMU_AccZ,GPS_VelN,GPS_VelE,GPS_VelD,roll, pitch, yaw, g,mass):
    R = RotatioMatrix(roll, pitch, yaw)

    gravity = np.dot(R,np.array([0,0,1]))*g*mass

    thrust = np.array([0.0,0.0,0.0])

    F_tot = np.array([IMU_AccX,IMU_AccY,IMU_AccZ])*mass


    vel = np.array([GPS_VelN,GPS_VelE,GPS_VelD])
    v_dir = np.dot(R,vel/np.linalg.norm(vel))
    #drag*v + thrust*v + gravity *v = F_tot*v
    drag = (np.dot(F_tot,v_dir) - np.dot(gravity,v_dir)  - np.dot(thrust,v_dir) )* v_dir
    lift = F_tot - gravity - thrust - drag
    return np.linalg.norm(drag), np.linalg.norm(lift)

def computeAOA(GPS_VelN,GPS_VelE,GPS_VelD,roll,pitch,yaw):
    R = RotatioMatrix(roll,pitch,yaw)
    vel = np.array([GPS_VelN,GPS_VelE,GPS_VelD])
    #velocity direction in body coordinate
    v_dir = np.dot(R, vel/np.linalg.norm(vel))

    alpha = 180*np.arctan(v_dir[2]/v_dir[0])/np.pi

    return alpha

def get_max_DL(lift, drag, throttle, b, S, rho, GPS_VelN,GPS_VelE,GPS_VelD):
    e=0.85 #Oswald factor ~0.8-0.9
    Cd0=np.zeros(len(lift))
    for i in range(len(lift)):
        if throttle>10**(-8.): continue
        
        vel = np.norm(np.array([GPS_VelN[i],GPS_VelE[i],GPS_VelD[i]]))
        Cd=2*drag[i]/(rho*vel*S)
        Cl=2*lift[i]/(rho*vel*S)
        Cd0[i]=Cd-Cl**2./(np.pi*S*e)
        
    #get avg: might choose a different model here!
    Cd0_avg=np.sum(Cd0)/np.count_nonzero(Cd0)
    #compute max_DL
    max_DL=.5*np.sqrt(np.pi*b**2./S/Cd0_avg)
    
    return max_DL




if __name__  == "__main__":

    log = np.genfromtxt('log002.csv', delimiter=",",filling_values = 0.0, skip_header=4)

    """
    skip_header=5,
    """
    g = 9.8
    mass = 1.031


    roll,pitch,yaw = log[:,4],log[:,5],log[:,6]
    IMU_AccX,IMU_AccY,IMU_AccZ = log[:,21],log[:,22],log[:,23]
    GPS_VelN,GPS_VelE,GPS_VelD = log[:,74],log[:,75],log[:,76]
    throttle = log[:,106]
    print(roll[0],IMU_AccX[0],GPS_VelN[0])

    len = len(roll)
    len = int(len/5)
    drag,lift,alpha = np.zeros(len),np.zeros(len),np.zeros(len)
    for i in range(len):

        drag[i],lift[i] = computeDragLift(IMU_AccX[i],IMU_AccY[i],IMU_AccZ[i],GPS_VelN[i],GPS_VelE[i],GPS_VelD[i],roll[i],pitch[i],yaw[i],g,mass)
        alpha[i] = computeAOA(GPS_VelN[i],GPS_VelE[i],GPS_VelD[i],roll[i],pitch[i],yaw[i])
        #if(throttle[i] > 1e-8):
         #   drag[i],lift[i],alpha[i] = 0.0, 0.0, 0.0
            
    for i in range(len):
        v_min, L=(200.,200.) #default values
        vel = np.array([GPS_VelN[i],GPS_VelE[i],GPS_VelD[i]])
        if throttle[i]<10**(-8.) and np.linalg.norm(vel)<v_min:
            v_min=vel
            L=lift[i]
            
    rho=1.1455 #kg.m-3
    S= #wing area
    Cl_max=2*L/(rho*vel**2.0 *S)
    
    #find stall speed
    v_stall=np.sqrt(2*mass*g/(rho*S*Cl_max))
    
    max_DL= get_max_DL(lift, drag, throttle, b, S, rho, GPS_VelN,GPS_VelE,GPS_VelD)

    plt.figure(1)
    plt.plot(alpha,lift,'ro')
    plt.xlabel('alpha')
    plt.ylabel('lift')
    plt.figure(2)
    plt.plot(alpha,drag,'ro')
    plt.xlabel('alpha')
    plt.ylabel('drag')

    plt.show()

