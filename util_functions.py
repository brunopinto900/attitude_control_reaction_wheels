"""
Utility functions for spacecraft attitude control.
Includes ramp functions, alignment errors, and reaction wheel updates.
"""
import numpy as np
from math_functions import dcm_to_rotvec


def ramp_wd_per_axis(C_err, t, t_settle=3.0, w_max=np.array([0.2,0.2,0.2])):
    """
    Generate commanded angular rate per body axis (3-element array)
    with ramp-up, peak, ramp-down profile.
    wd = 0 when C_err = identity.
    
    Parameters:
    -----------
    C_err : array
        DCM error matrix
    t : float
        Current time
    t_settle : float
        Settling time
    w_max : array
        Maximum angular rates per axis
    
    Returns:
    --------
    array : Commanded angular rate
    """
    rot_vec = dcm_to_rotvec(C_err)  # 3-element error
    # normalize to [-1,1] per axis to keep relative contributions
    max_err = np.max(np.abs(rot_vec))
    if max_err < 1e-6:
        return np.zeros(3)
    norm_vec = rot_vec / max_err
    # Half-sine ramp from 0 -> max -> 0 over t_settle
    ramp = np.sin(np.pi * t / t_settle)
    wd = norm_vec * ramp * w_max
    # Clip to limits
    wd = np.clip(wd, -w_max, w_max)
    return wd


def wd_linear_ramp(C_err, t_since_command, t_settle=3.0, w_max=np.array([0.6,0.6,0.6])):
    """
    Linear ramp angular rate per axis:
    - Starts at 0
    - Peaks at t_settle/2
    - Returns to 0 at t_settle
    - Returns zeros if t_since_command < 0 or > t_settle
    
    Parameters:
    -----------
    C_err : array
        DCM error matrix
    t_since_command : float
        Time since command was issued
    t_settle : float
        Settling time
    w_max : array
        Maximum angular rates per axis
    
    Returns:
    --------
    array : Commanded angular rate
    """
    t_settle = t_settle
    if t_since_command < 0 or t_since_command > t_settle:
        return np.zeros(3)
    
    # Linear ramp factor
    if t_since_command <= t_settle/2:
        ramp = 2 * t_since_command / t_settle  # 0 → 1
    else:
        ramp = 2 * (1 - t_since_command / t_settle)  # 1 → 0
    
    rot_vec = dcm_to_rotvec(C_err)  # 3-element error
    max_err = np.max(np.abs(rot_vec))
    if max_err < 1e-6:
        return np.zeros(3)
    norm_vec = rot_vec / max_err

    wd = ramp * w_max * norm_vec
    return wd