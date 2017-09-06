##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

"""PyQtGraph representations of TubeTK data structures."""

# Avoid the local module of the same name.
from importlib import import_module
np = import_module('numpy')
pg = import_module('pyqtgraph')


def tubes_as_circles(tubes, point_colors=[0, 0, 1.0, 1.0]):
    """PyQtGraph MeshData representation of tube points.

    Each tube point is represented by a circle scaled by its radius.  It is
    centered at its position and orientated by its normals.

    Parameters
    ----------
    tubes : NumPy array
        NumPy representation of the tube points.
    point_colors : (len(tubes), 4) array or 4 element, 1D array_like
        RGBA colors for each point.
    """

    n_points = len(tubes)

    point_colors = np.array(point_colors)
    if point_colors.ndim == 1:
        point_colors = point_colors[np.newaxis, :]
        point_colors = np.repeat(point_colors, n_points, axis=0)

    if point_colors.shape[0] != n_points or point_colors.shape[1] != 4:
        raise ValueError('point_colors does not have the correct shape.')

    def xy_circle(resolution):
        "Create the vertexes and faces of a circle in the x-y plane."
        vertexes = np.empty((resolution + 1, 3), dtype=float)
        vertexes[:, 2] = 0.
        theta = np.arange(resolution) * 2 * np.pi / resolution
        vertexes[0, 0] = 0.
        vertexes[1:, 0] = np.cos(theta)
        vertexes[0, 1] = 0.
        vertexes[1:, 1] = np.sin(theta)

        return vertexes

    # should be function parameter?
    resolution = 16
    xy_vertexes = xy_circle(resolution)

    def rotate_vertexes(vertexes, point):
        """Rotate the vertexes in the xy plane by the normals and tangent defined
        in a tube point.  The vertexes will then lie in the plane define by the
        two normals."""

        rotation = np.empty((3, 3), dtype=float)
        rotation[0, :] = point['Normal1'] / np.linalg.norm(point['Normal1'])
        rotation[1, :] = point['Normal2'] / np.linalg.norm(point['Normal2'])
        rotation[2, :] = point['Tangent'] / np.linalg.norm(point['Tangent'])

        rotated = np.dot(vertexes, rotation)
        return rotated

    vertexes_per_point = resolution + 1
    vertexes = np.empty((vertexes_per_point * n_points, 3), dtype=float)
    faces_per_point = resolution
    faces = np.empty((faces_per_point * n_points, 3), dtype=int)
    face_colors = np.empty((faces_per_point * n_points, 4), dtype=float)
    for ii in range(n_points):
        vertex_start = ii * vertexes_per_point
        scaled_vertex = xy_vertexes * tubes[ii]['Radius']
        vertexes[vertex_start:vertex_start + vertexes_per_point, :] = \
            rotate_vertexes(scaled_vertex, tubes[ii])
        vertexes[vertex_start:vertex_start + vertexes_per_point, :] += \
            tubes[ii]['Position']

        # 'wagonwheel' set of faces.
        face_start = ii * faces_per_point
        faces[face_start:face_start + faces_per_point, 1] = vertex_start
        faces[face_start:face_start + faces_per_point, 0] = \
            np.arange(vertex_start + 1, vertex_start + resolution + 1, dtype=int)
        faces[face_start:face_start + faces_per_point, 2] = \
            np.arange(vertex_start + 2, vertex_start + resolution + 2, dtype=int)
        faces[face_start + faces_per_point - 1, 2] = vertex_start + 1

        face_colors[face_start:face_start + faces_per_point, :] = \
            point_colors[ii]

    gl = import_module('pyqtgraph.opengl')
    mesh_data = gl.MeshData(vertexes=vertexes, faces=faces,
                            faceColors=face_colors)

    return mesh_data
