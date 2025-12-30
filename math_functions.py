"""
Mathematical functions for spacecraft attitude control.
Includes quaternion operations and matrix operations.
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
