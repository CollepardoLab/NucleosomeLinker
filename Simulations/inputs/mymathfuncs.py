import numpy as np

def quat_axis_angle(axis, angle):
    qw = np.cos(0.5 * angle)
    qx = axis[0] * np.sin(0.5 * angle)
    qy = axis[1] * np.sin(0.5 * angle)
    qz = axis[2] * np.sin(0.5 * angle)
    return np.array([qw, qx, qy, qz])

def quat_mul(qa, qb):
    r0 = qa[0] * qb[0] - qa[1] * qb[1] - qa[2] * qb[2] - qa[3] * qb[3]
    r1 = qa[0] * qb[1] + qa[1] * qb[0] + qa[2] * qb[3] - qa[3] * qb[2]
    r2 = qa[0] * qb[2] - qa[1] * qb[3] + qa[2] * qb[0] + qa[3] * qb[1]
    r3 = qa[0] * qb[3] + qa[1] * qb[2] - qa[2] * qb[1] + qa[3] * qb[0]
    return np.array([r0, r1, r2, r3])

def quat_norm(q):
    m = np.linalg.norm(q)
    return np.array([q[0] / m, q[1] / m, q[2] / m, q[3] / m])

def q_to_exyz(q):
    ex = (q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3],
          2.0 * (q[1] * q[2] + q[0] * q[3]),
          2.0 * (q[1] * q[3] - q[0] * q[2]))

    ey = (2.0 * (q[1] * q[2] - q[0] * q[3]),
          q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3],
          2.0 * (q[2] * q[3] + q[0] * q[1]))

    ez = (2.0 * (q[1] * q[3] + q[0] * q[2]),
          2.0 * (q[2] * q[3] - q[0] * q[1]),
          q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3])

    return np.array(ex), np.array(ey), np.array(ez)


def exyz_to_q(ex, ey, ez):
    """ taken from LAMMPS source code, converts ex,ey,ez frame vectors to quaternion orientations """

    # squares of quaternion components

    q0sq = 0.25 * (ex[0] + ey[1] + ez[2] + 1.0)
    q1sq = q0sq - 0.5 * (ey[1] + ez[2])
    q2sq = q0sq - 0.5 * (ex[0] + ez[2])
    q3sq = q0sq - 0.5 * (ex[0] + ey[1])

    q = np.array([0.0, 0.0, 0.0, 0.0])
    # some component must be greater than 1/4 since they sum to 1
    # compute other components from it

    if q0sq >= 0.25:
        q[0] = np.sqrt(q0sq)
        q[1] = (ey[2] - ez[1]) / (4.0 * q[0])
        q[2] = (ez[0] - ex[2]) / (4.0 * q[0])
        q[3] = (ex[1] - ey[0]) / (4.0 * q[0])
    elif q1sq >= 0.25:
        q[1] = np.sqrt(q1sq)
        q[0] = (ey[2] - ez[1]) / (4.0 * q[1])
        q[2] = (ey[0] + ex[1]) / (4.0 * q[1])
        q[3] = (ex[2] + ez[0]) / (4.0 * q[1])
    elif q2sq >= 0.25:
        q[2] = np.sqrt(q2sq)
        q[0] = (ez[0] - ex[2]) / (4.0 * q[2])
        q[1] = (ey[0] + ex[1]) / (4.0 * q[2])
        q[3] = (ez[1] + ey[2]) / (4.0 * q[2])
    elif q3sq >= 0.25:
        q[3] = np.sqrt(q3sq)
        q[0] = (ex[1] - ey[0]) / (4.0 * q[3])
        q[1] = (ez[0] + ex[2]) / (4.0 * q[3])
        q[2] = (ez[1] + ey[2]) / (4.0 * q[3])

    norm = np.linalg.norm(q)
    q = q / norm

    return q

def mag(arg):
    """
    magnitude of 3d vector
    """
    return np.sqrt(arg[0] * arg[0] + arg[1] * arg[1] + arg[2] * arg[2])

def unit_vec(vector):
    """
    returns vector converted to a unit vector
    """
    m = mag(vector)
    return np.array([vector[0] / m, vector[1] / m, vector[2] / m])

def rotation(xin, axis, angle):
    """
    y = rotate "x" by "angle" about "axis"  axis needs to be unit vector
    uses Rodrigue's Rotation Formula:
    y = x cos(t) + (k cross x)sin(t) + k(k.x)(1-cos(t))
    where t = angle, k = axis
    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    """
    axis = unit_vec(axis)
    kcrossx = np.cross(axis, xin)
    kdotx = np.dot(axis, xin)

    # print kcrossx
    # print kdotx
    yout = np.array([0.0, 0.0, 0.0])

    yout[0] = xin[0] * np.cos(angle) + kcrossx[0] * np.sin(angle) \
              + axis[0] * kdotx * (1.0 - np.cos(angle))
    yout[1] = xin[1] * np.cos(angle) + kcrossx[1] * np.sin(angle) \
              + axis[1] * kdotx * (1.0 - np.cos(angle))
    yout[2] = xin[2] * np.cos(angle) + kcrossx[2] * np.sin(angle) \
              + axis[2] * kdotx * (1.0 - np.cos(angle))
    return yout


def quatinv(q):
    return np.array([q[0],-q[1],-q[2],-q[3]])/np.linalg.norm(q)

def get_com(data):
    return np.mean(data, axis=0)
