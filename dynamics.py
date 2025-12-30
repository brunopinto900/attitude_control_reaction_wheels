"""
Dynamics functions for spacecraft attitude control.
Includes rigid-body dynamics and numerical integration.
"""
import numpy as np
from math_functions import skew


def dynamics(x, torque, I_inv, I):
    """
    Spacecraft rigid-body dynamics using quaternion representation.
    
    Parameters:
    -----------
    x : array
        State vector [q0, q1, q2, q3, wx, wy, wz]
    torque : array
        Applied torque vector
    I_inv : array
        Inverse of spacecraft inertia matrix
    I : array
        Spacecraft inertia matrix
    
    Returns:
    --------
    array : State derivative
    """
    q, w = x[:4], x[4:]
    wdot = I_inv @ (torque - np.cross(w, I @ w))
    Omega = np.zeros((4,4))
    Omega[0,1:] = -w
    Omega[1:,0] = w
    Omega[1:,1:] = -skew(w)
    qdot = 0.5 * Omega @ q
    return np.hstack([qdot, wdot])


def rk4(x, torque, dt, I_inv, I):
    """
    4th order Runge-Kutta integrator.
    
    Parameters:
    -----------
    x : array
        Current state vector
    torque : array
        Applied torque
    dt : float
        Time step
    I_inv : array
        Inverse of spacecraft inertia matrix
    I : array
        Spacecraft inertia matrix
    
    Returns:
    --------
    array : Updated state vector
    """
    k1 = dynamics(x, torque, I_inv, I)
    k2 = dynamics(x + dt/2*k1, torque, I_inv, I)
    k3 = dynamics(x + dt/2*k2, torque, I_inv, I)
    k4 = dynamics(x + dt*k3, torque, I_inv, I)
    return x + dt/6*(k1 + 2*k2 + 2*k3 + k4)
