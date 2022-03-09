from cmath import sqrt
import numpy as np
import matplotlib.pyplot as plt
import csv
from colour import Color

# %% Classes and functions


class Point(object):
    """
        Takes in x and y coordinates and holds number of neighbors inside a given distance.
        Also adds id to the point
    """

    def __init__(self, ids, x, y):
        self.x = int(x)
        self.y = int(y)
        self.p = np.array([int(x), int(y)])
        self.neigh = []  # neighbours within eps
        self.cneigh = []  # core neighbours within eps
        self.id = int(ids)
        self.type = 2
        self.cluster = -1

    def __repr__(self) -> np.ndarray:
        return '{:>2}:, [{:>2},{:>2}], t={}, cluster={}'.format(self.id, self.x, self.y, self.type, self.cluster)

    def __eq__(self, __o: object) -> bool:
        if isinstance(__o, Point):
            return self.id == __o.id

    def dist(self, other):
        if isinstance(other, Point):
            return np.linalg.norm(self.p - other.p)


def dbscan(ps: list, eps: float, mp: int):
    # Find all core points
    core = []
    border = []
    noise = []
    # Part 1: Label core, border and noise
    # Label core
    for i in range(len(ps)):
        nc = 0  # Neighbour count
        for j in range(len(ps)):
            if i == 8 and j==13:
                a=0
            if ps[i].dist(ps[j]) <= eps:
                nc += 1  # Also counts itself.
                if i != j:
                    ps[i].neigh.append(ps[j].id)  # Add neigbour id to point
                # print(ps[i], ps[j])
                if nc >= mp and ps[i].id not in core:
                    ps[i].type = 0  # Changes it to be core point
                    core.append(ps[i].id)
    # Label border
    for cid in core:
        for nid in ps[cid].neigh:
            if ps[nid].type == 2:
                ps[nid].type = 1  # Changes noise point to border
                border.append(ps[nid].id)
            elif ps[nid].type == 0:
                ps[cid].cneigh.append(ps[nid].id)  # create edge to core neighbour

    # The rest stay in noise
    for p in ps:
        if p.type == 2:
            noise.append(p.id)

    # Part 2: Find clusters
    clusters = []
    c = -1
    for cid in core:
        if ps[cid].cluster == -1:
            ct = -1
            for cn in ps[cid].cneigh:
                if ps[cid].cluster != -1:
                    ct = ps[cid].cluster
                    break
            if ct != -1:
                ps[cid].cluster = ct
            else:
                c += 1
                clusters.append(c)
                ps[cid].cluster = c
        for n in ps[cid].neigh:
            ps[n].cluster = ps[cid].cluster

    return (core, border, noise, clusters)


# %% Read data
f = open('redwood.dat', newline='')
dt = csv.reader(f)
# head = next(dt)

ids = []
x = []
y = []
points = []
i = 0
for l in dt:
    ids.append(str(i))
    # ids.append(str((l[0])))
    x.append(int(l[0]))
    y.append(int(l[1]))
    points.append(Point(i, l[1], l[2]))
    i += 1

f.close()


# %% Run
cpid, bpid, npid, clusters = dbscan(points, 4, 3)


# %% Plot points and core points
startC = Color('blue')
c = list(startC.range_to(Color('red'), len(ids)+1))

fig, ax = plt.subplots()

for i in range(0, len(ids)):
    plt.plot(x[i], y[i], marker='o', color=c[i].hex)
    ax.add_patch(plt.Circle((x[i], y[i]), 4, color=c[i].hex, fill=False))
    plt.text(x[i]+0.1, y[i]+0.1, ids[i], color=c[i].hex)
plt.scatter(x, y)
plt.xlim(0, 15)
plt.ylim(0, 13)
ax.set_aspect('equal')
plt.savefig("3_visualizePoints.png")
plt.show()


# Core points
fig, ax = plt.subplots()
startC = Color('blue')
cm = ['bo', 'go', 'ro']
c = list(startC.range_to(Color('green'), 2))
c.append(Color("black"))
t = ['Core', 'Border', 'Noise']

for p in points:
    plt.plot(int(p.x), int(p.y), color=c[int(p.type)].hex, marker='o')
    plt.text(p.x+0.1, p.y+0.1, p.id, color=c[int(p.type)].hex)
    if p.type == 0:  # Core
        ax.add_patch(plt.Circle((p.x, p.y), 4, color=c[int(p.type)].hex, fill=False))


plt.xlim(0, 15)
plt.ylim(0, 13)
ax.set_aspect('equal')
plt.savefig("3_corePoints.png")
plt.show()


# %% Plot resulting DBSCAN
fig, ax = plt.subplots()
startC = Color('blue')
cm = ['bo', 'go', 'ro']
ms = [14, 6, 6]
clustCol = list(startC.range_to(Color('red'), len(clusters)))
clustCol.append(Color('black'))
t = ['Core', 'Border', 'Noise']

for p in points:
    plt.plot(int(p.x), int(p.y),
             color=clustCol[p.cluster].hex,
             marker='o',
             markersize=ms[p.type])
    plt.text(p.x+0.3, p.y+0.3, p.id,
             color=clustCol[p.cluster].hex)
    # if p.type == 0:  # Core
    #     ax.add_patch(plt.Circle((p.x, p.y), 4, color=clustCol[p.cluster].hex, fill=False))


plt.xlim(0, 15)
plt.ylim(0, 13)
ax.set_aspect('equal')
plt.savefig("3_visualizeResult.png")
plt.show()
