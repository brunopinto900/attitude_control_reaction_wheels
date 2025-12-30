"""
Mathematical functions for spacecraft attitude control.
Includes quaternion operations, matrix operations, and attitude controller.
"""
import numpy as np


def euler_to_quat(roll, pitch, yaw):
    """Convert Euler angles to quaternion representation."""
    cr = np.cos(roll/2)
    sr = np.sin(roll/2)
    cp = np.cos(pitch/2)
    sp = np.sin(pitch/2)
    cy = np.cos(yaw/2)
    sy = np.sin(yaw/2)

    q0 = cr*cp*cy + sr*sp*sy
    q1 = sr*cp*cy - cr*sp*sy
    q2 = cr*sp*cy + sr*cp*sy
    q3 = cr*cp*sy - sr*sp*cy
    return np.array([q0,q1,q2,q3])


def skew(w):
    """Create skew-symmetric matrix from vector."""
    return np.array([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])


def quat_mult(q1, q2):
    """Multiply two quaternions."""
    w1, v1 = q1[0], q1[1:]
    w2, v2 = q2[0], q2[1:]
    return np.hstack([w1*w2-np.dot(v1,v2),
                      w1*v2+w2*v1+np.cross(v1,v2)])


def quat_conj(q):
    """Compute quaternion conjugate."""
    return np.array([q[0],-q[1],-q[2],-q[3]])


def quat_norm(q):
    """Normalize a quaternion."""
    return q/np.linalg.norm(q)


def quat_to_dcm(q):
    """Convert quaternion to direction cosine matrix (DCM)."""
    q0, q1, q2, q3 = q
    return np.array([
        [1-2*(q2*q2+q3*q3), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3), 1-2*(q1*q1+q3*q3), 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1*q1+q2*q2)]
    ])


def dcm_to_rotvec(C_err):
    """Extract rotation vector from DCM error."""
    return 0.5 * np.array([
        C_err[2,1]-C_err[1,2],
        C_err[0,2]-C_err[2,0],
        C_err[1,0]-C_err[0,1]
    ])


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
