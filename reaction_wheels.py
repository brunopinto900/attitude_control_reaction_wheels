"""
Reaction wheel functions for spacecraft attitude control.
Includes alignment error simulation and wheel dynamics updates.
"""
import numpy as np


def addAlignmentErrors(rw_axes_control):
    """
    Add alignment errors to reaction wheel axes.
    
    Parameters:
    -----------
    rw_axes_control : array
        Nominal reaction wheel axes (4x3)
    
    Returns:
    --------
    array : Reaction wheel axes with alignment errors
    """
    # Slightly off-axis axes for dynamics propagation
    theta_off = np.deg2rad(10.0)  # e.g., 10° tilt
    np.random.seed(0)
    perturb = np.random.randn(4,3)
    perturb = perturb / np.linalg.norm(perturb, axis=1)[:,None]

    rw_axes_dyn = np.zeros_like(rw_axes_control)
    for i in range(4):
        axis = perturb[i]
        cos_t = np.cos(theta_off)
        sin_t = np.sin(theta_off)
        cross = np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
        R = cos_t*np.eye(3) + sin_t*cross + (1-cos_t)*np.outer(axis,axis)
        rw_axes_dyn[i] = R @ rw_axes_control[i]
        rw_axes_dyn[i] /= np.linalg.norm(rw_axes_dyn[i])
    
    return rw_axes_dyn


def update_reaction_wheels(w_rw, tau_rw_cmd, J_rw, w_rw_max, dt, failure=False):
    """
    Update reaction wheel speeds with saturation-aware torque.
    
    Parameters:
    -----------
    w_rw : array
        Current wheel speeds
    tau_rw_cmd : array
        Commanded wheel torques
    J_rw : float
        Wheel rotor inertia
    w_rw_max : float
        Maximum wheel speed
    dt : float
        Time step
    failure : bool
        Whether to simulate RW failure
    
    Returns:
    --------
    tuple : (w_rw_new, tau_rw_actual) - updated speeds and actual applied torques
    """
    tau_rw_actual = np.zeros_like(tau_rw_cmd)
    w_rw_new = w_rw.copy()
    
    for i in range(4):
        if np.abs(w_rw[i]) < w_rw_max:
            # Wheel can accelerate
            tau_rw_actual[i] = tau_rw_cmd[i]
        else:
            # Wheel saturated, do not apply torque
            tau_rw_actual[i] = 0
        
        if failure:
            tau_rw_actual[0] = 0  # failure in RW #1
            # tau_rw_actual[1] = 0  # failure in RW #2

        # Update wheel speed
        w_rw_new[i] += (tau_rw_actual[i] / J_rw) * dt
        
        # Clip for safety
        w_rw_new[i] = np.clip(w_rw_new[i], -w_rw_max, w_rw_max)
    
    return w_rw_new, tau_rw_actual
