import math
import numpy as np
import operator

###

def apply_matrix(tf, points):
    tf = np.array(tf)
    r = np.dot(points, np.transpose(tf[:, :3])) ##dot is OK; from numpy import dot as matrix_multiply
    np.add(r, tf[:, 3], r)
    return r

# -----------------------------------------------------------------------------
# Transpose the rotation part.
#
def transpose_matrix(tf):

  return ((tf[0][0], tf[1][0], tf[2][0], tf[0][3]),
          (tf[0][1], tf[1][1], tf[2][1], tf[1][3]),
          (tf[0][2], tf[1][2], tf[2][2], tf[2][3]))

# -----------------------------------------------------------------------------
#
def multiply_matrices(*mlist):

  if len(mlist) == 2:
    m1, m2 = mlist
    p = [[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0]]
    for r in range(3):
      for c in range(4):
        p[r][c] = m1[r][0]*m2[0][c] + m1[r][1]*m2[1][c] + m1[r][2]*m2[2][c]
        p[r][3] += m1[r][3]
    p = tuple(map(tuple, p))
  else:
    p = multiply_matrices(*mlist[1:])
    p = multiply_matrices(mlist[0], p)
  return p

# -----------------------------------------------------------------------------
#
def identity_matrix():
    return ((1.0, 0, 0, 0), (0, 1.0, 0, 0), (0, 0, 1.0, 0))

# -----------------------------------------------------------------------------

def translation_matrix(shift):
    tf = np.array(((1.0, 0, 0, shift[0]),
                   (0, 1.0, 0, shift[1]),
                   (0, 0, 1.0, shift[2])))
    return tf
# -----------------------------------------------------------------------------

def coordinate_transform_list(tflist, ctf):

  ctfinv = invert_matrix(ctf)
  return [multiply_matrices(ctfinv, tf, ctf) for tf in tflist]

# -----------------------------------------------------------------------------

def recenter_symmetries(tflist, center):

    if center == (0,0,0):
        return tflist
    ctf = translation_matrix([-x for x in center])
    return coordinate_transform_list(tflist, ctf)

# -----------------------------------------------------------------------------

class Icosahedron(object):

    def __init__(self, circumscribed_radius = 1, orientation = '2n5', center=(0, 0, 0)):
        """point = np.array([0, 1, PHI]"""
        self.circumscribed_radius = circumscribed_radius
        self.orientation = orientation
        self.center=center
        self.e= math.sqrt(2 - 2 / math.sqrt(5))  # Triangle edge length of unit icosahedron.
        self.vertices = None  # icosahedron vertices
        self.triangles = None # icosahedron faces
        self._3foldAxis = [] #
        self._2foldAxis = [] #
        self.edges = None # icosahedron edges
        # ---------------------------------------------------------------------------------------------
        # Coordinates systems.
        # '222'         2-fold symmetry along x, y, and z axes.
        # '222r'        '222' with 90 degree rotation around z.
        # '2n5'         2-fold symmetry along x and 5-fold along z.
        # '2n5r'        '2n5' with 180 degree rotation about y.
        # 'n25'         2-fold symmetry along y and 5-fold along z.
        # 'n25r'        'n25' with 180 degree rotation about x.
        # '2n3'         2-fold symmetry along x and 3-fold along z.
        # '2n3r'        '2n3' with 180 degree rotation about y.
        #
        self.coordinate_system_names = ('222', '222r', '2n5', '2n5r', 'n25', 'n25r', '2n3', '2n3r')
        self.icosahedron_edge_length = self.icosahedron_edge_length(circumscribed_radius)
        self.angle23, self.angle25, self.angle35 = self.icosahedron_angles()
        self.icosahedron_geometry()

    # ---------------------------------------------------------------------------------------------
    # 60 icosahedral symmetry matrices.
    #
    def icosahedral_symmetry_matrices(self):
        t = self.icosahedral_matrix_table()
        tflist = recenter_symmetries(t.get(self.orientation, None), self.center)
        return tflist


    # -----------------------------------------------------------------------------
    # Edge length of icosahedron with a certain radio of the circumscribed sphere:
    # According to Radio of the circumscribed sphere = 1 (vertices), the icosahedron_edge_length should be 1.0515

    def icosahedron_edge_length(self, circumscribed_radius):
        return 4 * circumscribed_radius / math.sqrt(10 + 2 * math.sqrt(5))

    # -----------------------------------------------------------------------------
    # Matrices for mapping between different icosahedron coordinate frames.
    #
    def coordinate_system_transform(self, from_cs, to_cs):

        self.cst = {}
        if self.cst:
            return self.cst[(from_cs, to_cs)]

        transform = self.cst

        s25 = self.e / 2  # Sin/Cos for angle between 2-fold and 5-fold axis
        c25 = math.sqrt(1 - s25 * s25)
        s35 = self.e / math.sqrt(3)  # Sin/Cos for angle between 3-fold and 5-fold axis
        c35 = math.sqrt(1 - s35 * s35)

        transform[('2n5', '222')] = ((1, 0, 0, 0),
                                     (0, c25, -s25, 0),
                                     (0, s25, c25, 0))
        transform[('2n5', '2n3')] = ((1, 0, 0, 0),
                                     (0, c35, s35, 0),
                                     (0, -s35, c35, 0))

        # Axes permutations.
        transform[('222', '222r')] = ((0, 1, 0, 0),  # 90 degree rotation about z
                                      (-1, 0, 0, 0),
                                      (0, 0, 1, 0))
        transform[('2n3', '2n3r')] = \
            transform[('2n5', '2n5r')] = ((-1, 0, 0, 0),  # 180 degree rotation about y
                                          (0, 1, 0, 0),
                                          (0, 0, -1, 0))
        transform[('n25', 'n25r')] = ((1, 0, 0, 0),  # 180 degree rotation about x
                                      (0, -1, 0, 0),
                                      (0, 0, -1, 0))
        transform[('n25', '2n5')] = ((0, 1, 0, 0),  # x <-> y and z -> -z
                                     (1, 0, 0, 0),
                                     (0, 0, -1, 0))

        # Extend to all pairs of transforms.
        tlist = []
        while len(transform) > len(tlist):

            tlist = transform.keys()

            # Add inverse transforms
            for f, t in tlist:
                if not (t, f) in transform:
                    transform[(t, f)] = transpose_matrix(transform[(f, t)])

            # Use transitivity
            for f1, t1 in tlist:
                for f2, t2 in tlist:
                    if f2 == t1 and f1 != t2 and not (f1, t2) in transform:
                        transform[(f1, t2)] = multiply_matrices(transform[(f2, t2)],
                                                                transform[(f1, t1)])

        i = identity_matrix()
        for s in self.coordinate_system_names:
            transform[(s, s)] = i

        return transform[(from_cs, to_cs)]

    # ---------------------------------------------------------------------------------------
    # Compute icosahedral transformation matrices for different coordinate systems.
    #
    icos_matrices = {}  # Maps orientation name to 60 matrices.

    def icosahedral_matrix_table(self):
        global icos_matrices
        if icos_matrices:
            return icos_matrices

        c = math.cos(2 * math.pi / 5)  # .309016994
        c2 = math.cos(4 * math.pi / 5)  # -.809016994

        icos_matrices['222'] = (

            ((1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0)),

            ((c2, -0.5, c, 0.0),
             (-0.5, c, c2, 0.0),
             (c, c2, -0.5, 0.0)),

            ((0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0)),

            ((-c2, -0.5, -c, 0.0),
             (-0.5, -c, c2, 0.0),
             (c, -c2, -0.5, 0.0)),

            ((0.5, c, c2, 0.0),
             (-c, c2, -0.5, 0.0),
             (c2, 0.5, -c, 0.0)),

            ((-c, c2, -0.5, 0.0),
             (c2, 0.5, -c, 0.0),
             (0.5, c, c2, 0.0)),

            ((c2, 0.5, -c, 0.0),
             (0.5, c, c2, 0.0),
             (-c, c2, -0.5, 0.0)),

            ((c2, -0.5, -c, 0.0),
             (0.5, -c, c2, 0.0),
             (c, c2, 0.5, 0.0)),

            ((-c, -c2, -0.5, 0.0),
             (c2, -0.5, -c, 0.0),
             (-0.5, c, -c2, 0.0)),

            ((0.5, -c, c2, 0.0),
             (-c, -c2, -0.5, 0.0),
             (-c2, 0.5, c, 0.0)),

            ((0.0, 0.0, -1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0)),

            ((-0.5, -c, c2, 0.0),
             (c, -c2, -0.5, 0.0),
             (-c2, -0.5, -c, 0.0)),

            ((-0.5, c, c2, 0.0),
             (c, c2, -0.5, 0.0),
             (c2, -0.5, c, 0.0)),

            ((-c, c2, -0.5, 0.0),
             (-c2, -0.5, c, 0.0),
             (-0.5, -c, -c2, 0.0)),

            ((c2, 0.5, -c, 0.0),
             (-0.5, -c, -c2, 0.0),
             (c, -c2, 0.5, 0.0)),

            ((0.5, c, c2, 0.0),
             (c, -c2, 0.5, 0.0),
             (-c2, -0.5, c, 0.0)),

            ((-0.5, c, c2, 0.0),
             (-c, -c2, 0.5, 0.0),
             (-c2, 0.5, -c, 0.0)),

            ((0.0, 0.0, -1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0)),

            ((-0.5, -c, c2, 0.0),
             (-c, c2, 0.5, 0.0),
             (c2, 0.5, c, 0.0)),

            ((0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0)),

            ((c2, 0.5, c, 0.0),
             (0.5, c, -c2, 0.0),
             (c, -c2, -0.5, 0.0)),

            ((-c2, 0.5, -c, 0.0),
             (0.5, -c, -c2, 0.0),
             (c, c2, -0.5, 0.0)),

            ((-c, -c2, -0.5, 0.0),
             (-c2, 0.5, c, 0.0),
             (0.5, -c, c2, 0.0)),

            ((0.5, -c, c2, 0.0),
             (c, c2, 0.5, 0.0),
             (c2, -0.5, -c, 0.0)),

            ((c2, -0.5, -c, 0.0),
             (-0.5, c, -c2, 0.0),
             (-c, -c2, -0.5, 0.0)),

            ((-c, c2, 0.5, 0.0),
             (c2, 0.5, c, 0.0),
             (-0.5, -c, c2, 0.0)),

            ((-c, -c2, 0.5, 0.0),
             (-c2, 0.5, -c, 0.0),
             (-0.5, c, c2, 0.0)),

            ((1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0)),

            ((c, -c2, -0.5, 0.0),
             (-c2, -0.5, -c, 0.0),
             (-0.5, -c, c2, 0.0)),

            ((c, c2, -0.5, 0.0),
             (c2, -0.5, c, 0.0),
             (-0.5, c, c2, 0.0)),

            ((-1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0)),

            ((-c2, 0.5, -c, 0.0),
             (-0.5, c, c2, 0.0),
             (-c, -c2, 0.5, 0.0)),

            ((0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0)),

            ((c2, 0.5, c, 0.0),
             (-0.5, -c, c2, 0.0),
             (-c, c2, 0.5, 0.0)),

            ((-0.5, -c, -c2, 0.0),
             (-c, c2, -0.5, 0.0),
             (-c2, -0.5, c, 0.0)),

            ((c, -c2, 0.5, 0.0),
             (c2, 0.5, -c, 0.0),
             (-0.5, -c, -c2, 0.0)),

            ((-c2, -0.5, c, 0.0),
             (0.5, c, c2, 0.0),
             (c, -c2, 0.5, 0.0)),

            ((-c2, 0.5, c, 0.0),
             (0.5, -c, c2, 0.0),
             (-c, -c2, -0.5, 0.0)),

            ((c, c2, 0.5, 0.0),
             (c2, -0.5, -c, 0.0),
             (0.5, -c, c2, 0.0)),

            ((-0.5, c, -c2, 0.0),
             (-c, -c2, -0.5, 0.0),
             (c2, -0.5, -c, 0.0)),

            ((0.0, 0.0, 1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0)),

            ((0.5, c, -c2, 0.0),
             (c, -c2, -0.5, 0.0),
             (c2, 0.5, c, 0.0)),

            ((0.5, -c, -c2, 0.0),
             (c, c2, -0.5, 0.0),
             (-c2, 0.5, -c, 0.0)),

            ((c, -c2, 0.5, 0.0),
             (-c2, -0.5, c, 0.0),
             (0.5, c, c2, 0.0)),

            ((-c2, -0.5, c, 0.0),
             (-0.5, -c, -c2, 0.0),
             (-c, c2, -0.5, 0.0)),

            ((-0.5, -c, -c2, 0.0),
             (c, -c2, 0.5, 0.0),
             (c2, 0.5, -c, 0.0)),

            ((0.5, -c, -c2, 0.0),
             (-c, -c2, 0.5, 0.0),
             (c2, -0.5, c, 0.0)),

            ((0.0, 0.0, 1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0)),

            ((0.5, c, -c2, 0.0),
             (-c, c2, 0.5, 0.0),
             (-c2, -0.5, -c, 0.0)),

            ((0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0)),

            ((-c2, -0.5, -c, 0.0),
             (0.5, c, -c2, 0.0),
             (-c, c2, 0.5, 0.0)),

            ((c2, -0.5, c, 0.0),
             (0.5, -c, -c2, 0.0),
             (-c, -c2, 0.5, 0.0)),

            ((c, c2, 0.5, 0.0),
             (-c2, 0.5, c, 0.0),
             (-0.5, c, -c2, 0.0)),

            ((-0.5, c, -c2, 0.0),
             (c, c2, 0.5, 0.0),
             (-c2, 0.5, c, 0.0)),

            ((-c2, 0.5, c, 0.0),
             (-0.5, c, -c2, 0.0),
             (c, c2, 0.5, 0.0)),

            ((c, -c2, -0.5, 0.0),
             (c2, 0.5, c, 0.0),
             (0.5, c, -c2, 0.0)),

            ((c, c2, -0.5, 0.0),
             (-c2, 0.5, -c, 0.0),
             (0.5, -c, -c2, 0.0)),

            ((-1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0)),

            ((-c, c2, 0.5, 0.0),
             (-c2, -0.5, -c, 0.0),
             (0.5, c, -c2, 0.0)),

            ((-c, -c2, 0.5, 0.0),
             (c2, -0.5, c, 0.0),
             (0.5, -c, -c2, 0.0)),

        )

        for cs in self.coordinate_system_names:
            if cs != '222':
                t = self.coordinate_system_transform(cs, '222')
                tinv = self.coordinate_system_transform('222', cs)
                icos_matrices[cs] = [multiply_matrices(tinv, m, t)
                                     for m in icos_matrices['222']]
        return icos_matrices

    # -----------------------------------------------------------------------------

    def icosahedron_angles(self):
        #Angles between 2-fold, 3-fold and 5-fold symmetry axes of an icosahedron.
        angle25 = math.asin(self.e/2)
        angle35 = math.asin(self.e/math.sqrt(3))
        angle23 = math.asin(self.e/(2*math.sqrt(3)))
        return angle23, angle25, angle35

    def icosahedron_geometry(self):
        a = 2 * self.angle25  # Angle spanned by edge from center
        # 5-fold symmetry axis along z
        c5 = math.cos(2 * math.pi / 5)
        s5 = math.sin(2 * math.pi / 5)
        tf5 = ((c5, -s5, 0, 0),
               (s5, c5, 0, 0),
               (0, 0, 1, 0))

        # 2-fold symmetry axis along x
        tf2 = ((1, 0, 0, 0),
               (0, -1, 0, 0),
               (0, 0, -1, 0))

        p = tuple(map(operator.add, self.center, (0, 0, self.circumscribed_radius)))
        p50 = tuple(map(operator.add, self.center, (0, self.circumscribed_radius * math.sin(a), self.circumscribed_radius * math.cos(a))))
        p51 = apply_matrix(tf5, p50)
        p52 = apply_matrix(tf5, p51)
        p53 = apply_matrix(tf5, p52)
        p54 = apply_matrix(tf5, p53)
        self.vertices = [p, p50, p51, p52, p53, p54]
        self.vertices.extend([apply_matrix(tf2, q) for q in self.vertices])
        if self.orientation != '2n5':
            self.tf = self.coordinate_system_transform('2n5', self.orientation)
            self.vertices = [apply_matrix(self.tf, p) for p in self.vertices]

            #
            # Vertex numbering
            #
            #  Top   1          Bottom
            #      2   5        9   10
            #        0            6
            #      3   4        8   11
            #                     7
        # 20 triangles composing icosahedron.
        #
        self.triangles = ((0,1,2), (0,2,3), (0,3,4), (0,4,5), (0,5,1),
                         (6,7,8), (6,8,9), (6,9,10), (6,10,11), (6,11,7),
                         (1,9,2), (2,9,8), (2,8,3), (3,8,7), (3,7,4),
                         (4,7,11), (4,11,5), (5,11,10), (5,10,1), (1,10,9))

        for triangle in self.triangles:
            self._3foldAxis.append( [(item1 + item2 + item3) /3. \
                                             for item1,  item2,  item3 in zip(self.vertices[triangle[0]],
                                                                              self.vertices[triangle[1]],
                                                                              self.vertices[triangle[2]])])

        self.edges = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
                      (1, 2), (2, 3), (3, 4), (4, 5), (5, 1),
                      (6, 7), (6, 8), (6, 9), (6, 10), (6, 11),
                      (1, 9), (2, 8), (2, 9), (3, 7), (3, 8),
                      (7, 8), (8, 9), (9, 10), (10, 11), (11, 7),
                      (4, 7), (4, 11), (5, 10), (5, 11), (1, 10))

        for edge in self.edges:
            self._2foldAxis.append([(item1 + item2 ) / 2. \
                                          for item1, item2 in zip(self.vertices[edge[0]],
                                                                  self.vertices[edge[1]])])

    def get_vertices(self):
        return self.vertices

    def get_triangles(self):
        return self.triangles

    def get_edges(self):
        return self.edges

    def get_3foldAxis(self):
        return self._3foldAxis

    def get_2foldAxis(self):
        return self._2foldAxis





"""
    def getVertices(self):
        return self.vertices

    def getFaces(self):
        return self.faces

    def getPairs(self):
        return self.pairs

    def icosahedron_vertices(self):
        
        generate icosahedron vertices
       
        # used to permute for all 12 vertices
        point = (self.point)

        # normalize
        point = point/distance(point)

        # construct 12 vertices by cyclically permuting [0, +/-1, +/- PHI]
        pqueue = collections.deque(point)

        Vertice.clear()

        for _ in range(3):
            pqueue.rotate()
            for i in [-1, 1]:
                for j in [-1, 1]:
                    for k in [-1, 1]:
                        Vertice.make(np.multiply(pqueue, [i, j, k]))

    def icosahedron_faces(self):
        
        generate icosahedron faces
        
        self.pairs = collections.defaultdict(set)
        for vert in Vertice.vertices.values():
            min_dist = min([vert.distance(other_vert) \
                for other_vert in Vertice.vertices.values() if vert.distance(other_vert) > 0])

            for other_vert in Vertice.vertices.values():
                if vert.distance(other_vert) == min_dist:
                    self.pairs[vert.key].add(other_vert.key)
                    self.pairs[other_vert.key].add(vert.key)

        Face.clear()

        for vert, verts in self.pairs.iteritems():
            for v1 in verts:
                if v1 == vert:
                    continue
                for v2 in verts:
                    if v2 == vert or v2 == v1:
                        continue
                    if v2 in self.pairs[v1] and v1 in self.pairs[v2] and vert in self.pairs[v1] and vert in self.pairs[v2]:
                        key = tuple(sorted([vert, v1, v2]))

                        face = Face.make([np.array(v, dtype=np.float64) for v in key])

                        Vertice.vertices.get(vert).face = face
                        Vertice.vertices.get(v1).face = face
                        Vertice.vertices.get(v2).face = face

    def icosahedron_bisect(self,level):
        
        bisect the faces "level" times
        

        l = 1
        while l < level:
            l += 1

            faces = Face.faces

            Vertice.clear()
            Face.clear()

            for fkey, face in faces.iteritems():
                a = midpoint(face.v1, face.v2)
                b = midpoint(face.v2, face.v3)
                c = midpoint(face.v3, face.v1)

                Face.make([a, face.v1, c])
                Face.make([a, face.v2, b])
                Face.make([a, b, c])
                Face.make([b, face.v3, c])

    def icosahedron_geojson(self):
        
        return geojson
                json = {'type': 'FeatureCollection', 'features': []}

        for fkey, face in Face.faces.iteritems():
            vert = face.vertex
            json['features'].append(dict(
                type='Feature',
                geometry=dict(
                    type='Polygon',
                    coordinates=[face.coordinates]
                ),
                properties=dict(
                    x=vert[0], y=vert[1], z=vert[2]
                )
            ))

        return json

    def icosahedron_dual(self):
        
        convert icosahedron into its dual
        

        # build a mapping of vertices to the connecting faces
        vertFaces = collections.defaultdict(set)
        for vert in Vertice.vertices.values():
            for fkey in vert.face:
                vertFaces[vert.key].add(fkey)

        # store the vertices and faces
        vertices = Vertice.vertices
        faces = Face.faces

        # clear the vertices and faces for the dual
        Vertice.clear()
        Face.clear()

        # build the new vertices and faces
        for vkey, fkeys in vertFaces.iteritems():
            # store the first face
            face = faces[fkeys.pop()]
            _faces = [face]
            l = len(fkeys)

            # loop through the other faces
            # grabbing the nearest face (by its corresponding vertex)
            while len(_faces) <= l:
                fkey = sorted(fkeys, key = lambda fkey: distance(face.vertex, faces[fkey].vertex))[0]
                face = faces[fkey]
                fkeys.remove(fkey)
                _faces.append(faces[fkey])

            # build a new face from the plane -> vertices of the original faces
            Face.make([face.vertex for face in _faces])

def main():
    
    run the program
    
    icosahedron = Icosahedron( np.array([0, 1, PHI], dtype=np.float64) )
    #for key, value in icosahedron.getVertices().iteritems():
    #    print type(key), type(value)
    #for key, value in icosahedron.getFaces().iteritems():
    #    print key, value
    for key, value in icosahedron.getPairs().iteritems():
        print key, value

if __name__ == '__main__':
    main()
"""