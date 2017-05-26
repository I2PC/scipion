"""
Example usage:
python icosahedron.py --json --dual --level 3 > icosahedron.json && ogr2ogr -f PostgreSQL PG:"user=root password=secret" -overwrite -nln icosahedron icosahedron.json
"""

import json
import math
import collections
import numpy as np

ORIGIN = np.array([0, 0, 0], dtype=np.float64)

# golden ratio
PHI = (1 + math.sqrt(5))/2

def distance(v1, v2=ORIGIN):
    """
    return distance between two vectors
    """
    return np.linalg.norm(v1-v2)


def midpoint(v1, v2):
    """
    return midpoint between two vectors
    """
    v = (v1 + v2)
    mp = v/np.linalg.norm(v)
    return mp


class Vertice(object):
    vertices = {}

    def __init__(self, vertice, face=None):
        self.vertice = vertice
        self.key = tuple(self.vertice)
        self._faces = set()

    @property
    def face(self):
        return self._faces

    @face.setter
    def face(self, face):
        if face:
            self._faces.add(face.key)

    @classmethod
    def make(cls, vertice, face=None):
        vert = Vertice(vertice, face)
        if vert.key in Vertice.vertices:
            vert = cls.vertices[vert.key]
            return vert
        else:
            cls.vertices[vert.key] = vert
            vert.id = len(cls.vertices)
            return vert


    def distance(self, other):
        return distance(self.vertice, other.vertice)


    def latlon(self):
        """
        return lon, lat of point
        """
        x, y, z = tuple(self.vertice.tolist())
        th = 25 * math.pi / 180
        x1 = math.cos(th)*x - math.sin(th)*y
        y1 = math.sin(th)*x + math.cos(th)*y
        return [math.atan2(y1, x1)*180/math.pi, math.asin(z)*180/math.pi]

    @property
    def coordinates(self):
        return self.latlon()

    @classmethod
    def clear(cls):
        cls.vertices = {}

    def __repr__(self):
        return 'V(%d)' %(self.id)

class Face(object):

    faces = {}

    def __init__(self, vertices, id):

        self.id = id
        self.vertices = [Vertice.make(vert, self) for vert in vertices]
        self.key = tuple(sorted([vert.key for vert in self.vertices]))

        for vertice in self.vertices:
            vertice.face = self

        # convert face to a vertex of the polyhedral dual

        # geometrically point x lies on face with closest approach to origin y
        # IFF
        #   1. x . y = y . y
        # use polar reciprocal transformation to convert vertex x into plane x'
        #   2. x' = x / (x . x)
        #
        # in order for the transformation to be correct
        # vertex y' must lie on plane x'
        # IFF
        #   3. x' . y' = x' . x'
        # plugging Eq 2 into 3 for both x' and y' and simplifying
        #   x . y / (( y . y) ( x . x )) = x . x / (x . x) ^2
        #   x . y / (y . y) = 1
        #   x . y = y . y
        #
        # which proves that vertex x lies on plane y under the tranformation
        # so we've successfully shown that the polar reciprocal tranformation
        # can be used to tranform a polyhedral into its dual

        # normal can be computed using (v2 - v1) cross (v3 - v1)
        normal = np.cross(self.v2 - self.v1, self.v3 - self.v1)
        self.normal = normal/np.linalg.norm(normal)

        # plugging Eq 1. into Eq 2.
        # x' = x / (x . y)
        self.vertex = self.normal/(sum([ np.dot(v, self.normal) for v in vertices ])/len(vertices))
        self.vertex = self.vertex / np.linalg.norm(self.vertex)

    @classmethod
    def make(cls, vertices):
        id = len(cls.faces) + 1
        face = Face(vertices, id)
        if face.key in cls.faces:
            return cls.faces[face.key]
        else:
            cls.faces[face.key] = face
            return face

    @classmethod
    def clear(cls):
        cls.faces = {}

    @property
    def v1(self):
        return self.vertices[0].vertice

    @property
    def v2(self):
        return self.vertices[1].vertice

    @property
    def v3(self):
        return self.vertices[2].vertice

    @property
    def coordinates(self):
        coordinates = [vert.coordinates for vert in self.vertices]

        coordinates.append(coordinates[0])

        l = [c[0] for c in coordinates]
        dl = max(l) - min(l)
        if dl > 180:
            coordinates = [[c[0] + 360, c[1]] if c[0] < 0 else c for c in coordinates]

        return coordinates

    def __repr__(self):
        return 'Face(%d)' %(self.id)

class Icosahedron(object):

    def __init__(self, point=np.array([0, 1, PHI], dtype=np.float64)):
        """point = np.array([0, 1, PHI]"""
        self.point = point
        self.pairs = None # vertex and its 5 neighbours
        self.icosahedron_vertices()
        self.icosahedron_faces()
        self.vertices = Vertice.vertices
        self.faces = Face.faces

        #self.icosahedron_bisect(options.level)
        #self.icosahedron_dual()

    def getVertices(self):
        return self.vertices

    def getFaces(self):
        return self.faces

    def getPairs(self):
        return self.pairs

    def icosahedron_vertices(self):
        """
        generate icosahedron vertices
        """
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
        """
        generate icosahedron faces
        """
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
        """
        bisect the faces "level" times
        """

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
        """
        return geojson
        """
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
        """
        convert icosahedron into its dual
        """

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
    """
    run the program
    """
    icosahedron = Icosahedron( np.array([0, 1, PHI], dtype=np.float64) )
    #for key, value in icosahedron.getVertices().iteritems():
    #    print type(key), type(value)
    #for key, value in icosahedron.getFaces().iteritems():
    #    print key, value
    for key, value in icosahedron.getPairs().iteritems():
        print key, value

if __name__ == '__main__':
    main()
