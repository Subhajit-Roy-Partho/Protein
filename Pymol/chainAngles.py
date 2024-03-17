from pymol import cmd, stored
import numpy as np
import sys

cmd.load(str(sys.argv[1]),"mol")

def angle_between_vectors(v1, v2):
    # Calculate the dot product of the vectors
    dot_product = np.dot(v1, v2)

    # Calculate the norms (magnitudes) of the vectors
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)

    # Calculate the cosine of the angle
    cos_angle = dot_product / (norm_v1 * norm_v2)

    # Calculate the angle in radians and then convert to degrees
    angle_radians = np.arccos(np.clip(cos_angle, -1.0, 1.0))
    angle_degrees = np.degrees(angle_radians)

    return angle_degrees


def bestfit(chain):
    points = []
    i = 1
    for i in range(1, 100):
        # Select the residue
        cmd.select("A", "chain "+chain+" and resi " + str(i))
        stored.coords = []
        cmd.iterate_state(1, "A", "stored.coords.append((x,y,z))")
        # print(len(stored.coords))
        if len(stored.coords) == 0:
            continue
        points.append(np.mean(stored.coords, axis=0))
    datamean = np.mean(points, axis=0);
    uu, dd, vv = np.linalg.svd(points - datamean)
    return vv[0];

best1 = bestfit(str(sys.argv[2]));
best2 = bestfit(str(sys.argv[3]));

# print(best2);
print(angle_between_vectors(best1, best2));