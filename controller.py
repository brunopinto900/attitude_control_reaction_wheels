"""
Attitude controller functions for spacecraft control.
Includes cascaded attitude control architecture.
"""
import numpy as np


def cascaded_attitude_controller(
    q_err_vec,
    omega,
    J,
    Ts=6.0,
    zeta=1/np.sqrt(2),
    omega_max=1.0,
    torque_max=None
):
    """
    Cascaded attitude controller with outer loop (quaternion error -> rate command)
    and inner loop (rate error -> torque).
    
    Parameters:
    -----------
    q_err_vec : array
        Quaternion error vector part
    omega : array
        Current angular velocity
    J : array
        Spacecraft inertia matrix
    Ts : float
        Settling time
    zeta : float
        Damping ratio
    omega_max : float
        Maximum angular rate
    torque_max : float, optional
        Maximum torque
    
    Returns:
    --------
    tuple : (tau, omega_cmd) - torque command and rate command
    """
    # --- Natural frequency from settling time ---
    omega_n = 4.0 / (zeta * Ts)

    # --- Outer loop: quaternion error -> rate command ---
    k_q = -2.0 * omega_n
    omega_cmd = k_q * q_err_vec

    # Rate saturation
    norm = np.linalg.norm(omega_cmd)
    if norm > omega_max:
        omega_cmd = omega_cmd / norm * omega_max

    # --- Inner loop: rate error -> torque ---
    omega_rate = 5.0 * omega_n
    K_omega = J * omega_rate

    tau = -K_omega @ (omega - omega_cmd)

    # Optional torque saturation
    if torque_max is not None:
        tau_norm = np.linalg.norm(tau)
        if tau_norm > torque_max:
            tau = tau / tau_norm * torque_max

    return tau, omega_cmd
